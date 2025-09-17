#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <fstream>
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
            // std::cout << "Press Enter to continue...\n";
            // std::cin.get();

        } while (error > epsilon);
        
        std::cout << "\nРешение системы (за " << k << " итераций):\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i + 1 << "] = " << std::fixed << std::setprecision(6) << X_current.Get(i) << "\n";
        }
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