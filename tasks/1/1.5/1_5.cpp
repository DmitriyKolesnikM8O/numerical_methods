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
    static void Do(Matrix A, double epsilon) {
        int k = 1;
        do {
            auto pair = A.QR();
            Matrix Q = pair.first;
            Matrix R = pair.second;

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

                std::cout << "Лямбда " << j << ": " << x << " + " << y << "i\n";
                std::cout << "Лямбда " << (j + 1) << ": " << x << " - " << y << "i\n";
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