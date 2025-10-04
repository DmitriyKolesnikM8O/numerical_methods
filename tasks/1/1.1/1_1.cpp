#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <utility>
#include <stdexcept>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task1 {
public:
    static void Do(Matrix A, Vector B) {
        int n = A.GetLength();
        
        Matrix L, U;
        Vector P;
        int permutations;
        std::tie(L, U, P, permutations) = A.LU();

        std::cout << "\nМатрица U:\n";
        U.Show();
        std::cout << "\nМатрица L:\n";
        L.Show();
        


        double det = A.Determinant();
        std::cout << "\nОпределитель матрицы: " << std::fixed << std::setprecision(6) << det << "\n";

        Matrix inverse_A = A.Inverse();
        Matrix check = inverse_A * A;
        std::cout << "\nОбратная матрица:\n";
        inverse_A.Show();
        std::cout << "\nПроверка:\n";
        check.Show();
        
        //реализация решения Ax = b с помощью LU

        //переставляем исходный b
        Vector permuted_B(n);
        for (int i = 0; i < n; i++) {
            permuted_B.Set(i, B.Get(static_cast<int>(P.Get(i))));
        }

        //прямая подстановка: Ly = Pb
        Vector Y(n);
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L.Get(i, j) * Y.Get(j);
            }
            Y.Set(i, permuted_B.Get(i) - sum);
        }


        //обратная подстановка: Ux = y
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
            std::cout << "x[" << i + 1 << "] = " << std::fixed << std::setprecision(6) << X.Get(i) << "\n";
        }
    }
};

int main(int argc, char* argv[]) {
    try {
        Matrix A;
        Vector B;
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/1/1.1/1_1.txt";
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

        Task1::Do(A, B);
        
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
