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
        
        // Преобразование системы Ax=B к виду x = C*x + D
        Matrix C(n);
        Vector D(n);

        for (int i = 0; i < n; i++) {
            if (std::abs(A.Get(i, i)) < 1e-9) {
                throw std::runtime_error("Диагональный элемент близок к нулю. Невозможно выполнить итерации.");
            }
            D.Set(i, B.Get(i) / A.Get(i, i));
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    C.Set(i, j, 0.0);
                } else {
                    C.Set(i, j, -A.Get(i, j) / A.Get(i, i));
                }
            }
        }
        
        // Начальное приближение
        Vector X_current = D;
        Vector X_next(n);
        
        int k = 0;
        double error;
        do {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += C.Get(i, j) * X_current.Get(j);
                }
                X_next.Set(i, sum + D.Get(i));
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
            
            std::cout << "\nError: " << std::fixed << std::setprecision(6) << error << "\n";
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