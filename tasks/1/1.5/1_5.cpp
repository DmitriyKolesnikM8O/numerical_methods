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
        for (int j = 0; j < m; j++) {
            sum1 = SquareSumColumn(matrix, j, j + 1);
            sum2 = SquareSumColumn(matrix, j, j + 2);

            if (sum2 > epsilon) {
                return false;
            } else if (sum1 <= epsilon) {
                std::cout << "Лямбда" << j << ": " << matrix.Get(j, j) << "\n";
            } else if (sum1 > epsilon) {
                if (j == 0) return false;

                double aii = matrix.Get(j, j);
                double ajj = matrix.Get(j + 1, j + 1);
                double aij = matrix.Get(j, j + 1);
                double aji = matrix.Get(j + 1, j);

                double x = (aii + ajj) / 2.0;
                double y = std::sqrt(-(aii + ajj) * (aii + ajj) + 4 * (aii * ajj - aij * aji)) / 2.0;

                std::cout << "Лямбда " << j << ": " << x << " + " << y << "i\n";
                std::cout << "Лямбда " << (j + 1) << ": " << x << " - " << y << "i\n";
                j++;
            }
        }
        return true;
    }

    static double SquareSumColumn(const Matrix& matrix, int column_number, int first_index) {
        int n = matrix.GetHeight();
        double sum = 0.0;
        for (int i = first_index; i < n; i++) {
            sum += matrix.Get(i, column_number) * matrix.Get(i, column_number);
        }
        return std::sqrt(sum);
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