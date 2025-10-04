#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task3 {
public:
    static void CheckSolution(const Matrix& A, const Vector& B, const Vector& X) {
        int n = A.GetLength();
        Vector Residual(n);

        //R - вектор невязки
        // 1. Вычисляем R = B - A * X
        for (int i = 0; i < n; i++) {
            double Ax_i = 0.0;
            for (int j = 0; j < n; j++) {
                Ax_i += A.Get(i, j) * X.Get(j);
            }
            Residual.Set(i, B.Get(i) - Ax_i);
        }

        // 2. Вычисляем норму невязки (например, бесконечную норму: max|R_i|)
        double residual_norm = 0.0;
        for (int i = 0; i < n; i++) {
            double abs_r = std::abs(Residual.Get(i));
            if (abs_r > residual_norm) {
                residual_norm = abs_r;
            }
        }

        std::cout << "\n--- Проверка решения ---\n";
        std::cout << "Максимальный элемент вектора невязки (||B - A*X||_inf): " 
                  << std::scientific << std::setprecision(2) << residual_norm << "\n";
        
        if (residual_norm > 1e-6) {
             std::cout << "Внимание: Норма невязки высока. Решение может быть неточным.\n";
        } else {
             std::cout << "Норма невязки мала. Решение считается надежным.\n";
        }
        std::cout << "----------------------\n";
    }

    static void Do(const Matrix& A, const Vector& B, double epsilon) {
        int n = A.GetLength();
        
        if (!A.isDiagonallyDominant()) {
            std::cerr << "Внимание: Матрица не обладает свойством диагонального преобладания. Метод может не сойтись." << std::endl;
        }

        Vector X_current(n);
        
        // Начальное приближение: X(0) = B / Aii
        for (int i = 0; i < n; i++) {
            if (std::abs(A.Get(i, i)) < 1e-9) {
                throw std::runtime_error("Диагональный элемент близок к нулю. Невозможно выполнить итерации.");
            }
            X_current.Set(i, B.Get(i) / A.Get(i, i));
        }

        int iterations_count = 0;
        double error;
        do {
            Vector X_prev = X_current;

            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        sum += A.Get(i, j) * X_current.Get(j);
                    }
                }
                X_current.Set(i, (B.Get(i) - sum) / A.Get(i, i));
            }

            error = 0.0;
            for (int i = 0; i < n; i++) {
                double diff = std::abs(X_current.Get(i) - X_prev.Get(i));
                if (diff > error) {
                    error = diff;
                }
            }
            iterations_count++;

            // std::cout << "Error: " << error << "\n";
            // for (int i = 0; i < n; i++) {
            //     std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
            // }

        } while (error > epsilon);

        std::cout << "\nРешение системы (за " << iterations_count << " итераций):\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i + 1 << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
        }

        Task3::CheckSolution(A, B, X_current);
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
