#include "Matrix.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <fstream>

using namespace std;

class Task2_2 {
private:
    
    // Увеличиваем TAU для ускорения сходимости (0.1 вместо 0.01)
    static constexpr double TAU = 0.1; 
    
    // Начальное приближение (определено графически)
    static constexpr double X1_INIT = 0.7; 
    static constexpr double X2_INIT = 1.0;

    // --- Функции невязки F(X) ---
    static double F1(double x1, double x2) {
        return x1 * x1 + 4.0 * x2 * x2 - 4.0;
    }

    static double F2(double x1, double x2) {
        return 2.0 * x2 - std::exp(x1) - x1;
    }
    
    // --- Функции для метода простой итерации G(X) (Метод релаксации) ---
    // G1 = x1 - TAU * F1 (x1, x2)
    static double G1(double x1, double x2) {
        return x1 - TAU * F1(x1, x2);
    }

    // G2 = x2 - TAU * F2 (x1, x2)
    static double G2(double x1, double x2) {
        return x2 - TAU * F2(x1, x2);
    }

    // --- Матрица Якоби J(X) и ее элементы ---
    static double J11(double x1) { return 2.0 * x1; }
    static double J12(double x2) { return 8.0 * x2; }
    static double J21(double x1) { return -std::exp(x1) - 1.0; }
    static double J22() { return 2.0; }

    // --- Норма вектора F (максимум модуля компоненты) ---
    static double norm_F(double x1, double x2) {
        return std::max(std::fabs(F1(x1, x2)), std::fabs(F2(x1, x2)));
    }

    // --- Норма приращения (максимум модуля компоненты) ---
    static double norm_delta(double delta1, double delta2) {
        return std::max(std::fabs(delta1), std::fabs(delta2));
    }


public:
    static void Do(double eps) {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "\n--- Решение системы для варианта 13 (a=2): x1^2+4x2^2-4=0, 2x2-e^x1-x1=0 ---\n";
        std::cout << "Начальное приближение: X(0) = (" << X1_INIT << ", " << X2_INIT << ")\n";

        // ---------------------------------------------------
        // Метод простой итерации (Метод релаксации)
        // ---------------------------------------------------
        
        
        // ПРОВЕРКА УСЛОВИЯ СХОДИМОСТИ (СЖАТИЯ) МПИ ***
        // Вычисляем норму Якоби итерационной функции G' = I - TAU * F' в начальной точке
        
        // 1. Элементы F'(X) в X0:
        double J11_init = J11(X1_INIT);
        double J12_init = J12(X2_INIT);
        double J21_init = J21(X1_INIT);
        double J22_init = J22();

        // 2. Элементы G'(X) = I - TAU * F'(X):
        double G_11 = 1.0 - TAU * J11_init;
        double G_12 = -TAU * J12_init;
        double G_21 = -TAU * J21_init;
        double G_22 = 1.0 - TAU * J22_init;
        
        // 3. Расчет нормы (Используем бесконечную норму: max сумма модулей по строкам)
        double norm_row1 = std::fabs(G_11) + std::fabs(G_12);
        double norm_row2 = std::fabs(G_21) + std::fabs(G_22);
        double norm_G_prime = std::max(norm_row1, norm_row2);
        
       if (norm_G_prime < 1.0) {
            std::cout << "Условие сходимости: ||G'(X0)||inf = " 
                      << norm_G_prime << " < 1. Сходимость ожидается.\n";
        } else {
             // Сохраняем вывод для информации, но убираем 'ВНИМАНИЕ' для чистоты
             std::cout << "Проверка сходимости: ||G'(X0)||inf = " 
                       << norm_G_prime << " >= 1 (достаточное условие не выполнено, но сходимость возможна).\n";
        }


        double x1_fp = X1_INIT;
        double x2_fp = X2_INIT;
        int iterations_fp = 0;
        vector<double> fp_errors;
        const int MAX_ITER_FP = 10000;

        try {
            double current_norm_F; // Норма невязки ||F(X)||
            
            // Первая итерация: просто расчет
            double x1_old = x1_fp;
            double x2_old = x2_fp;

            x1_fp = G1(x1_old, x2_old); 
            x2_fp = G2(x1_old, x2_old);
            
            double initial_err = std::max(std::fabs(x1_fp - x1_old), std::fabs(x2_fp - x2_old));
            fp_errors.push_back(initial_err);
            iterations_fp++;
            
            // Цикл: пока норма невязки больше точности (eps)
            do {
                x1_old = x1_fp;
                x2_old = x2_fp;

                // X(k+1) = G(X(k))
                x1_fp = G1(x1_old, x2_old); 
                x2_fp = G2(x1_old, x2_old);

                // Оценка погрешности (приращения)
                double err = std::max(std::fabs(x1_fp - x1_old), std::fabs(x2_fp - x2_old));
                fp_errors.push_back(err);
                
                // Главный критерий остановки: Норма невязки ||F(X)|| <= eps
                current_norm_F = std::max(std::fabs(F1(x1_fp, x2_fp)), std::fabs(F2(x1_fp, x2_fp)));

                iterations_fp++;
                if (iterations_fp > MAX_ITER_FP) { 
                    throw std::runtime_error("Метод простой итерации не сошелся за " + std::to_string(MAX_ITER_FP) + " шагов.");
                }
            } while (current_norm_F >= eps); // Использование невязки как критерия

        } catch (const std::exception& e) {
            std::cerr << "\nОшибка в методе простой итерации: " << e.what() << std::endl;
            if (iterations_fp > 0) {
                 std::cout << "Последнее приближение: X = (" << x1_fp << ", " << x2_fp << ")\n";
            }
            return;
        }

        std::cout << "\n--- Метод простой итерации (TAU=" << TAU << ") ---\n";
        std::cout << "Приближенный корень: X = (" << x1_fp << ", " << x2_fp << ")\n";
        std::cout << "Проверка F(X): (" << F1(x1_fp, x2_fp) << ", " << F2(x1_fp, x2_fp) << ")\n";
        std::cout << "Количество итераций: " << iterations_fp << endl;
        
        std::cout << "\nЗависимость погрешности от количества итераций:\n";
        size_t limit = std::min((size_t)10, fp_errors.size());
        for (size_t i = 0; i < limit; ++i) {
            std::cout << "Итерация " << (i + 1) << ": погрешность = " << fp_errors[i] << endl;
        }
        if (fp_errors.size() > 10) {
            std::cout << "...\n";
            std::cout << "Итерация " << fp_errors.size() << ": погрешность = " << fp_errors.back() << endl;
        }

        // ---------------------------------------------------
        // Метод Ньютона
        // ---------------------------------------------------
        double x1_newton = X1_INIT;
        double x2_newton = X2_INIT;
        int iterations_newton = 0;
        vector<double> newton_errors;
        
        try {
            do {
                double x1_old = x1_newton;
                double x2_old = x2_newton;

                double f1 = F1(x1_old, x2_old);
                double f2 = F2(x1_old, x2_old);

                double j11 = J11(x1_old);
                double j12 = J12(x2_old);
                double j21 = J21(x1_old); 
                double j22 = J22();

                //определитель матрицы якоби
                double det = j11 * j22 - j12 * j21;
                if (std::fabs(det) < 1e-10) {
                    throw std::runtime_error("Вырожденный якобиан в методе Ньютона.");
                }

                // Решение системы J * Delta = -F (по правилу Крамера)
                double delta1 = (-f1 * j22 + f2 * j12) / det;
                double delta2 = (f1 * j21 - f2 * j11) / det;


                //
                x1_newton = x1_old + delta1;
                x2_newton = x2_old + delta2;

                double err = std::max(std::fabs(delta1), std::fabs(delta2));
                newton_errors.push_back(err);

                iterations_newton++;
                if (iterations_newton > 100) {
                    throw std::runtime_error("Метод Ньютона не сошелся за 100 шагов.");
                }
            } while (newton_errors.back() >= eps);
        } catch (const std::exception& e) {
            std::cerr << "\nОшибка в методе Ньютона: " << e.what() << std::endl;
            return;
        }


        std::cout << "\n--- Метод Ньютона ---\n";
        std::cout << "Приближенный корень: X = (" << x1_newton << ", " << x2_newton << ")\n";
        std::cout << "Проверка F(X): (" << F1(x1_newton, x2_newton) << ", " << F2(x1_newton, x2_newton) << ")\n";
        std::cout << "Количество итераций: " << iterations_newton << endl;

        std::cout << "\nЗависимость погрешности от количества итераций:\n";
        for (size_t i = 0; i < newton_errors.size(); ++i) {
            std::cout << "Итерация " << (i + 1) << ": погрешность = " << newton_errors[i] << endl;
        }
    }
};

int main(int argc, char* argv[]) {
     try {
    
        std::cout << std::fixed << std::setprecision(10);
        double eps;
        
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/2/2.2/2_2.txt";
            if (argc > 2) {
                filename = argv[2];
            }
            file_stream.open(filename);
            if (!file_stream.is_open()) {
                throw std::runtime_error("Не удалось открыть файл " + filename);
            }
            input_stream = &file_stream;
            std::cout << "Чтение данных из файла " << filename << std::endl;
            eps = get_user_epsilon_file(*input_stream);
        } else {
            std::cout << "Чтение данных с консоли" << std::endl;
            eps = get_user_epsilon();
        }        
        
        Task2_2::Do(eps);
    } catch (const std::exception& e) {
        std::cerr << "Произошла фатальная ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
