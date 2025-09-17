#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task3 {
public:
    static void Do(const Matrix& A, const Vector& B, double epsilon) {
        int n = A.GetLength();
        
        // Проверка на диагональное преобладание (достаточное условие сходимости)
        if (!isDiagonallyDominant(A)) {
            std::cerr << "Внимание: Матрица не обладает свойством диагонального преобладания. Метод может не сойтись." << std::endl;
        }

        Vector X_current(n);
        Vector X_next(n);

        // Начальное приближение: X(0) = B / Aii
        for (int i = 0; i < n; i++) {
            if (std::abs(A.Get(i, i)) < 1e-9) {
                throw std::runtime_error("Диагональный элемент близок к нулю. Невозможно выполнить итерации.");
            }
            X_current.Set(i, B.Get(i) / A.Get(i, i));
        }

        int k = 0;
        double error;
        do {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        sum += A.Get(i, j) * X_current.Get(j);
                    }
                }
                X_next.Set(i, (B.Get(i) - sum) / A.Get(i, i));
            }

            error = 0.0;
            for (int i = 0; i < n; i++) {
                double diff = std::abs(X_next.Get(i) - X_current.Get(i));
                if (diff > error) {
                    error = diff;
                }
                X_current.Set(i, X_next.Get(i));
            }
            k++;

            // std::cout << "Error: " << error << "\n";
            // for (int i = 0; i < n; i++) {
            //     std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
            // }

            // std::cout << "Press Enter to continue...\n";
            // std::cin.get();

        } while (error > epsilon);

        std::cout << "\nРешение системы (за " << k << " итераций):\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i + 1 << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
        }
    }

private:
    static bool isDiagonallyDominant(const Matrix& matrix) {
        int n = matrix.GetLength();
        for (int i = 0; i < n; i++) {
            double diagonal_val = std::abs(matrix.Get(i, i));
            double sum_off_diagonal = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum_off_diagonal += std::abs(matrix.Get(i, j));
                }
            }
            if (diagonal_val < sum_off_diagonal) {
                return false;
            }
        }
        return true;
    }
};

int main(int argc, char* argv[]) {
    try {
        Matrix A;
        Vector B;
        double epsilon;
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/1/1.3/1_3.txt";
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
            B = get_free_members_vector_file(*input_stream);
            epsilon = get_user_epsilon_file(*input_stream);
        } else {
            std::cout << "Чтение данных с консоли" << std::endl;
            A = get_user_matrix_input_rows();
            B = get_free_members_vector_rows();
            epsilon = get_user_epsilon();
        }
        if (A.GetLength() != B.GetSize()) {
            throw std::invalid_argument("Размерность матрицы и вектора должны совпадать.");
        }

        Task3::Do(A, B, epsilon);
        
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}