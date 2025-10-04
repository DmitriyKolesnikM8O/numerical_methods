#include "Matrix.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <fstream>

using namespace std;


// Константы для варианта 13: ln(x+1) - 2x + 0.5 = 0


class Task2_1 {
private:
    // Интервал, содержащий положительный корень
    static constexpr double A = 0.0;
    static constexpr double B = 1.0;

    // Функция f(x)
    static double f(double x) {
        return std::log(x + 1) - 2.0 * x + 0.5;
    }

    // Функция g(x) для метода простой итерации: g(x) = 0.5 * (ln(x+1) + 0.5)
    static double g(double x) {
        return 0.5 * (std::log(x + 1) + 0.5);
    }

    // Производная f'(x) для метода Ньютона: f'(x) = 1/(x+1) - 2
    static double df(double x) {
        return 1.0 / (x + 1) - 2.0;
    }

public:
    static void Do(double eps) {
        std::cout << "\n--- Решение уравнения ln(x+1) - 2x + 0.5 = 0 на интервале [" 
                  << A << ", " << B << "] ---\n";
        
        // ---------------------------------------------------
        // Метод простой итерации
        // ---------------------------------------------------
        double x0_fp = (A + B) / 2.0; // Начальное приближение
        double x_fp = g(x0_fp);
        int iterations_fp = 1;
        vector<double> fp_errors;

        fp_errors.push_back(std::fabs(x_fp - x0_fp));

        // Проверка условия сходимости для x0
        if (std::fabs(g(x0_fp) - x0_fp) >= eps && std::fabs(1.0 / (2.0 * (x0_fp + 1))) >= 1.0) {
            std::cerr << "Внимание: Условие сходимости |g'(x)| < 1 может быть нарушено на границе. Продолжаем итерации." << std::endl;
        }

        while (std::fabs(x_fp - x0_fp) >= eps) {
            x0_fp = x_fp;
            x_fp = g(x0_fp);
            fp_errors.push_back(std::fabs(x_fp - x0_fp));
            iterations_fp++;
            if (iterations_fp > 100) { // Уменьшено ограничение для быстросходящегося метода
                std::cout << "Итерации не сошлись за 100 шагов (метод простой итерации)." << std::endl;
                return;
            }
        }

        std::cout << "\n--- Метод простой итерации ---\n";
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "Положительный корень примерно равен: " << x_fp << endl;
        std::cout << "Проверка: f(x) = " << f(x_fp) << endl;
        std::cout << "Количество итераций: " << iterations_fp << endl;

        std::cout << "\nЗависимость погрешности от количества итераций:" << endl;
        for (size_t i = 0; i < fp_errors.size(); ++i) {
            std::cout << "Итерация " << (i + 1) << ": погрешность = " << fp_errors[i] << endl;
        }

        // ---------------------------------------------------
        // Метод Ньютона
        // ---------------------------------------------------
        double x0_newton = (A + B) / 2.0; // Начальное приближение
        double x_newton;
        int iterations_newton = 0;
        vector<double> newton_errors;
        
        do {
            if (iterations_newton > 0) {
                x0_newton = x_newton;
            }
            
            // Проверка на деление на ноль
            if (std::fabs(df(x0_newton)) < 1e-9) {
                throw std::runtime_error("Производная f'(x) близка к нулю. Метод Ньютона не применим.");
            }

            x_newton = x0_newton - f(x0_newton) / df(x0_newton);
            
            if (iterations_newton > 0) {
                newton_errors.push_back(std::fabs(x_newton - x0_newton));
            } else {
                // Первая ошибка - это разность между первым итерационным значением и начальным приближением
                newton_errors.push_back(std::fabs(x_newton - x0_newton)); 
            }
            iterations_newton++;
            
            if (iterations_newton > 100) {
                std::cout << "Итерации не сошлись за 100 шагов (метод Ньютона)." << endl;
                return;
            }
        } while (std::fabs(x_newton - x0_newton) >= eps || iterations_newton == 1);


        std::cout << "\n--- Метод Ньютона ---\n";
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "Положительный корень примерно равен: " << x_newton << endl;
        std::cout << "Проверка: f(x) = " << f(x_newton) << endl;
        std::cout << "Количество итераций: " << iterations_newton << endl;

        std::cout << "\nЗависимость погрешности от количества итераций:" << endl;
        for (size_t i = 0; i < newton_errors.size(); ++i) {
            std::cout << "Итерация " << (i + 1) << ": погрешность = " << newton_errors[i] << endl;
        }
    }
};

int main(int argc, char* argv[]) {
     try {
        std::cout << std::ios::fixed << std::setprecision(10);
        double eps;
        std::istream* input_stream = &std::cin;
        std::ifstream file_stream;
        if(argc > 1 && std::string(argv[1]) == "-file") {
            std::string filename = "./tasks/2/2.1/2_1.txt";
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
        

        Task2_1::Do(eps);
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
