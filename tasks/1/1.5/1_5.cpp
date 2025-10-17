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
    static void CheckOrthogonality(const Matrix& Q) {
        int n = Q.GetLength();
        
        std::cout << "\n--- Проверка ортогональности Q (на последней итерации) ---\n";

        // 1. Проверка ортогональности Q (Ошибка Q = ||Q^T * Q - I||)
        Matrix Q_T = Q.Transpose(); 
        Matrix Q_T_Q = Q_T * Q;
        Matrix Identity = Matrix::GetSingularMatrix(n); 
        Matrix Residual_Q = Q_T_Q - Identity;
        
        double error_Q = Residual_Q.GetNorm();

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Ошибка ортогональности ||Q^T*Q - I||_inf: " << error_Q << "\n";
        
        if (error_Q < 1e-8) {
            std::cout << "Ортогональность Q подтверждена. Ошибка мала.\n";
        } else {
            std::cerr << "Внимание: Ошибка ортогональности Q превышает допуск 1e-8. Возможно, проблемы с реализацией Householder.\n";
        }
        std::cout << "--------------------------------\n";
    }
    static void Do(Matrix A, double epsilon) {
        int k = 1;
        Matrix Q_last; 
        Matrix R_last; 
        do {
            auto pair = A.QR();
            Matrix Q = pair.first;
            Matrix R = pair.second;

            Q_last = Q; 
            R_last = R; 

            //получаем следующую матрицу в последовательности
            A = R * Q;
            
            // std::cout << "\nИтерация " << k << ":\nМатрица A:\n";
            // A.Show();
            // std::cout << "Press Enter to continue...\n";
            // std::cin.get();
            k++;
            
        } while (!Finish(A, epsilon));

        std::cout << "\nИтоговое число итераций: " << k << "\n";
        std::cout << "\nМатрица A:\n";
        A.Show();

        Matrix res = A.Subtract(Q_last * R_last);

        CheckOrthogonality(Q_last); 
    }

private:
    static bool Finish(const Matrix& matrix, double epsilon) {
        int m = matrix.GetLength();
        bool converged = true;
        
        int j = 0;
        while (j < m) {
            if (j == m - 1) {
                // Последний элемент - всегда вещественное собственное значение
                std::cout << "Лямбда " << j << ": " << matrix.Get(j, j) << "\n";
                j++;
            } else {
                double subdiag = std::abs(matrix.Get(j + 1, j));
                
                if (subdiag <= epsilon) {
                    // Поддиагональный элемент мал - вещественное собственное значение
                    std::cout << "Лямбда " << j << ": " << matrix.Get(j, j) << "\n";
                    j++;
                } else {
                    // Проверяем, сошелся ли блок 2x2
                    double next_subdiag = (j + 2 < m) ? std::abs(matrix.Get(j + 2, j + 1)) : 0.0;
                    
                    if (next_subdiag <= epsilon) {
                        // Блок 2x2 сошелся - вычисляем собственные значения
                        double a11 = matrix.Get(j, j);
                        double a12 = matrix.Get(j, j + 1);
                        double a21 = matrix.Get(j + 1, j);
                        double a22 = matrix.Get(j + 1, j + 1);
                        
                        // Характеристическое уравнение: λ² - (a11+a22)λ + (a11*a22 - a12*a21) = 0
                        double trace = a11 + a22;
                        double det = a11 * a22 - a12 * a21;
                        double discriminant = trace * trace - 4 * det;
                        
                        if (discriminant >= 0) {
                            // Два вещественных корня
                            double sqrt_d = std::sqrt(discriminant);
                            double lambda1 = (trace + sqrt_d) / 2.0;
                            double lambda2 = (trace - sqrt_d) / 2.0;
                            std::cout << "Лямбда " << j << ": " << lambda1 << "\n";
                            std::cout << "Лямбда " << (j + 1) << ": " << lambda2 << "\n";
                        } else {
                            // Комплексно-сопряженные корни
                            double real_part = trace / 2.0;
                            double imag_part = std::sqrt(-discriminant) / 2.0;
                            std::cout << "Лямбда " << j << ": " << real_part << " + " << imag_part << "i\n";
                            std::cout << "Лямбда " << (j + 1) << ": " << real_part << " - " << imag_part << "i\n";
                        }
                        j += 2; // Перескакиваем через обработанный блок
                    } else {
                        // Блок еще не сошелся - продолжаем итерации
                        converged = false;
                        break;
                    }
                }
            }
        }
        
        return converged;
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
