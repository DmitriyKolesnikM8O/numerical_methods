#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <fstream>
#include <stdexcept>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task2 {
public:
    static void Do(const Matrix& A, const Vector& B) {
        int n = A.GetLength();
        
        // Проверка на трёхдиагональность
        if (!isTriDiagonal(A)) {
            throw std::invalid_argument("Матрица не является трехдиагональной.");
        }

        // Прямой ход (вычисление коэффициентов)
        Vector alpha(n);
        Vector beta(n);
        double z;

        alpha.Set(0, -A.Get(0, 1) / A.Get(0, 0));
        beta.Set(0, B.Get(0) / A.Get(0, 0));

        for (int i = 1; i < n - 1; i++) {
            z = A.Get(i, i) + A.Get(i, i - 1) * alpha.Get(i - 1);
            if (std::abs(z) < 1e-9) {
                 throw std::runtime_error("Матрица вырождена или плохо обусловлена.");
            }
            alpha.Set(i, -A.Get(i, i + 1) / z);
            beta.Set(i, (B.Get(i) - A.Get(i, i - 1) * beta.Get(i - 1)) / z);
        }

        // Обратный ход (вычисление решения)
        Vector X(n);
        z = A.Get(n - 1, n - 1) + A.Get(n - 1, n - 2) * alpha.Get(n - 2);
        if (std::abs(z) < 1e-9) {
            throw std::runtime_error("Матрица вырождена или плохо обусловлена.");
        }
        X.Set(n - 1, (B.Get(n - 1) - A.Get(n - 1, n - 2) * beta.Get(n - 2)) / z);

        for (int i = n - 2; i >= 0; i--) {
            X.Set(i, alpha.Get(i) * X.Get(i + 1) + beta.Get(i));
        }

        std::cout << "\nРешение системы:\n";
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i + 1 << "] = " << std::fixed << std::setprecision(6) << X.Get(i) << "\n";
        }
    }

private:
    static bool isTriDiagonal(const Matrix& matrix) {
        int n = matrix.GetLength();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (std::abs(i - j) > 1 && std::abs(matrix.Get(i, j)) > 1e-9) {
                    return false;
                }
            }
        }
        return true;
    }
};

int main(int argc, char* argv[]) {
    try {
        Matrix A;
        Vector B;
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/1/1.2/1_2.txt";
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
        } else {
            std::cout << "Чтение данных с консоли" << std::endl;
            A = get_user_matrix_input_rows();
            B = get_free_members_vector_rows();
        }        
        if (A.GetLength() != B.GetSize()) {
            throw std::invalid_argument("Размерность матрицы и вектора должны совпадать.");
        }

        Task2::Do(A, B);
        
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}