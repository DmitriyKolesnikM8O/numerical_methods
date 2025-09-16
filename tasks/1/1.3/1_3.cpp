#include <iostream>
#include <iomanip>
#include <cmath>
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

            std::cout << "Error: " << error << "\n";
            for (int i = 0; i < n; i++) {
                std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
            }

            std::cout << "Press Enter to continue...\n";
            std::cin.get();

        } while (error > epsilon);

        std::cout << "\nРешение системы (за " << k << " итераций):\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
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

int main() {
    try {
        Matrix A = get_user_matrix_input();
        Vector B = get_free_members_vector();
        double epsilon = get_user_epsilon();
        
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