#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task1 {
public:
    static void Do(Matrix A, Vector B) {
        int n = A.GetLength();
        
        // LU-разложение с частичным выбором главного элемента
        Matrix L = Matrix::GetSingularMatrix(n);
        Matrix U = A;
        Vector P(n);
        for (int i = 0; i < n; i++) {
            P.Set(i, i);
        }

        for (int k = 0; k < n; k++) {
            // Поиск главного элемента
            double max_val = std::abs(U.Get(k, k));
            int max_row = k;
            for (int i = k + 1; i < n; i++) {
                if (std::abs(U.Get(i, k)) > max_val) {
                    max_val = std::abs(U.Get(i, k));
                    max_row = i;
                }
            }

            // Перестановка строк в U и векторе перестановок P
            if (max_row != k) {
                for (int j = 0; j < n; j++) {
                    double temp_u = U.Get(k, j);
                    U.Set(k, j, U.Get(max_row, j));
                    U.Set(max_row, j, temp_u);
                }
                double temp_p = P.Get(k);
                P.Set(k, P.Get(max_row));
                P.Set(max_row, temp_p);

                // Перестановка строк в L
                for (int j = 0; j < k; j++) {
                    double temp_l = L.Get(k, j);
                    L.Set(k, j, L.Get(max_row, j));
                    L.Set(max_row, j, temp_l);
                }
            }

            // Создание L и U
            for (int i = k + 1; i < n; i++) {
                if (std::abs(U.Get(k, k)) < 1e-9) {
                    throw std::runtime_error("Матрица вырождена или плохо обусловлена.");
                }
                double factor = U.Get(i, k) / U.Get(k, k);
                L.Set(i, k, factor);
                for (int j = k; j < n; j++) {
                    U.Set(i, j, U.Get(i, j) - factor * U.Get(k, j));
                }
            }
        }
        
        std::cout << "\nМатрица U:\n";
        U.Show();
        std::cout << "\nМатрица L:\n";
        L.Show();

        // Перестановка вектора B
        Vector permuted_B(n);
        for (int i = 0; i < n; i++) {
            permuted_B.Set(i, B.Get(static_cast<int>(P.Get(i))));
        }

        // Решение Ly = B'
        Vector Y(n);
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L.Get(i, j) * Y.Get(j);
            }
            Y.Set(i, permuted_B.Get(i) - sum);
        }

        // Решение Ux = Y
        Vector X(n);
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += U.Get(i, j) * X.Get(j);
            }
            if (std::abs(U.Get(i, i)) < 1e-9) {
                throw std::runtime_error("Матрица вырождена или плохо обусловлена.");
            }
            X.Set(i, (Y.Get(i) - sum) / U.Get(i, i));
        }
        
        std::cout << "\nРешение системы:\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << X.Get(i) << "\n";
        }
    }
};

int main() {
    try {
        Matrix A = get_user_matrix_input();
        Vector B = get_free_members_vector();
        
        if (A.GetLength() != B.GetSize()) {
            throw std::invalid_argument("Размерность матрицы и вектора должны совпадать.");
        }

        Task1::Do(A, B);
        
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}