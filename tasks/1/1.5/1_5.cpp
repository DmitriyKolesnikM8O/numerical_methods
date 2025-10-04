#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <utility>
#include <stdexcept>
#include "Vector.hpp"
#include "Matrix.hpp"

class Task5 {
public:
    static void CheckOrthogonality(const Matrix& Q) {
        int n = Q.GetLength();
        
        std::cout << "\n--- Проверка ортогональности Q (на последней итерации) ---\n";

        // 1. Проверка ортогональности Q (Ошибка Q = ||Q^T * Q - I||)
        Matrix Q_T = Q.Transpose(); 
        Matrix Q_T_Q = Q_T * Q;
        Matrix Identity = Matrix::GetSingularMatrix(n); 
        Matrix Residual_Q = Q_T_Q.Subtract(Identity);
        
        double error_Q = Residual_Q.GetNorm();

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Ошибка ортогональности ||Q^T*Q - I||_inf: " << error_Q << "\n";
        
        if (error_Q < 1e-8) {
            std::cout << "Ортогональность Q подтверждена. Ошибка мала.\n";
        } else {
            std::cerr << "Внимание: Ошибка ортогональности Q превышает допуск 1e-8. Возможно, проблемы с реализацией Householder.\n";
        }
        std::cout << "--------------------------------\n";
    }
    static void Do(Matrix A, double epsilon) {
        int k = 1;
        Matrix Q_last; 
        Matrix R_last; 
        do {
            auto pair = A.QR();
            Matrix Q = pair.first;
            Matrix R = pair.second;

            Q_last = Q; 
            R_last = R; 

            //получаем следующую матрицу в последовательности
            A = R * Q;
            
            // std::cout << "\nИтерация " << k << ":\nМатрица A:\n";
            // A.Show();
            // std::cout << "Press Enter to continue...\n";
            // std::cin.get();
            k++;
            
        } while (!Finish(A, epsilon));


        std::cout << "\nМатрица A:\n";
        A.Show();

        CheckOrthogonality(Q_last); 
    }

private:
    static bool Finish(const Matrix& matrix, double epsilon) {
        int m = matrix.GetLength();
        double sum1, sum2;
        //по каждому столбцу проверяем
        for (int j = 0; j < m; j++) {
            sum1 = matrix.SquareSumColumn(j, j + 1); // норма всех элементов под диагональю
            sum2 = matrix.SquareSumColumn(j, j + 2); // норма всех элементов ниже 2 элемента под диагональю

            if (sum2 > epsilon) {
                return false;
            } else if (sum1 <= epsilon) { //успех
                std::cout << "Лямбда" << j << ": " << matrix.Get(j, j) << "\n";
            } else if (sum1 > epsilon) { //пока не сошелся столбец
                if (j == 0) return false;

                /*
                проверка комплексно сопряженных корней
                когда матрица не сходится к диагональной, а оставляем блоки 2 на 2
                */
                double aii = matrix.Get(j, j);
                double ajj = matrix.Get(j + 1, j + 1);
                double aij = matrix.Get(j, j + 1);
                double aji = matrix.Get(j + 1, j);

                double x = (aii + ajj) / 2.0; //действительная часть корня
                double y = std::sqrt(-(aii + ajj) * (aii + ajj) + 4 * (aii * ajj - aij * aji)) / 2.0; //мнимая

                std::cout << "Лямбда  " << j << ": " << x << " + " << y << "i\n";
                std::cout << "Лямбда  " << (j + 1) << ": " << x << " - " << y << "i\n";
                j++;
            }
        }
        return true;
    }

    
};

int main(int argc, char* argv[]) {
    try
    {
        Matrix A;
        double epsilon;
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/1/1.5/1_5.txt";
            if (argc > 2) {
                filename = argv[2];
            }
            file_stream.open(filename);
            if (!file_stream.is_open()) {
                throw std::runtime_error("Не удалось открыть файл " + filename);
            }
            input_stream = &file_stream;
            std::cout << "Чтение данных из файла " << filename << std::endl;
            A = get_user_matrix_input_file(*input_stream);
            epsilon = get_user_epsilon_file(*input_stream);
        } else {
            std::cout << "Чтение данных с консоли" << std::endl;
            A = get_user_matrix_input_rows();
            epsilon = get_user_epsilon();
        }
        Task5::Do(A, epsilon);
        
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
    
}
