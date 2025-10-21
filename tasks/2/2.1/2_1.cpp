#include "Matrix.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <fstream>

using namespace std;


// Вариант 13: ln(x+1) - 2x + 0.5 = 0


class Task2_1 {
private:
    // Интервал, содержащий положительный корень
    // Выбран, так как функция на его концах имеет разные знаки
    // Соответственно в нем нечетное число корней (минимум 1)
    // static constexpr double A = 0.0;
    // static constexpr double B = 1.0;

    // Функция f(x)
    static double f(double x) {
        return std::log(x + 1) - 2.0 * x + 0.5;
    }

    // Функция g(x) для метода простой итерации: g(x) = 0.5 * (ln(x+1) + 0.5) (x выразили)
    static double g(double x) {
        return 0.5 * (std::log(x + 1) + 0.5);
    }

    // Производная f'(x) для метода Ньютона: f'(x) = 1/(x+1) - 2
    static double df(double x) {
        return 1.0 / (x + 1) - 2.0;
    }
    
    // Вторая производная f''(x) для проверки условия сходимости Ньютона: f''(x) = -1/(x+1)^2
    static double d2f(double x) {
        return -1.0 / ((x + 1) * (x + 1));
    }

public:

    /*
        start_x: начало поиска (например, 0.0)
        end_x: конец поиска (например, 10.0, или другое разумное число)
        h: шаг поиска (например, 0.1)
        теорема Больцано Коши:
            если функция f(x) непрерывна на отрезке [A,B] и на концах отрезка принимает значения разных знаков 
            (f(A)⋅f(B)<0), то внутри отрезка (A,B) существует хотя бы один корень ξ (то есть точка, где f(ξ)=0).
    */
    static std::pair<double, double> find_root_interval(double start_x, double end_x, double h) {
        
        double x_curr = start_x;
        double x_next;
        
        // Проверка, что логарифм определен: x+1 > 0
        if (start_x <= -1.0) {
            x_curr = -1.0 + h; // Сдвигаемся чуть правее -1
            if (x_curr >= end_x) throw std::runtime_error("Интервал поиска слишком мал.");
        }
        
        double f_curr = f(x_curr);

        while (x_curr < end_x) {
            x_next = x_curr + h;
            
            // Прерываем поиск, если вышли за end_x
            if (x_next > end_x) x_next = end_x;
            
            double f_next = f(x_next);

            // Условие изоляции корня: смена знака
            if (f_curr * f_next < 0) {
                // Корень найден в интервале [x_curr, x_next]
                return {x_curr, x_next};
            }

            // Если корень находится точно в узле (редко, но бывает)
            if (std::fabs(f_curr) < 1e-9) { 
                return {x_curr - h, x_curr}; 
            }

            // Сдвигаемся к следующему шагу
            x_curr = x_next;
            f_curr = f_next; 
        }

        // Если корень не найден в заданном диапазоне
        throw std::runtime_error("Корень не найден в заданном диапазоне [" + 
                                std::to_string(start_x) + ", " + std::to_string(end_x) + "].");
    }   


    static void Do(double eps) {

        double A = 0.0;
        double B = 0.0;
        try {
            auto interval = find_root_interval(0.0, 5.0, 0.1); // Ищем от 0 до 5 с шагом 0.1
            A = interval.first;
            B = interval.second;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка при поиске интервала: " << e.what() << std::endl;
            return;
        }
        // A = 0.4;
        // B = 0.5;


        std::cout << "\n--- Решение уравнения ln(x+1) - 2x + 0.5 = 0 на интервале [" 
                  << A << ", " << B << "] ---\n";
        
        // ---------------------------------------------------
        // Метод простой итерации
        // ---------------------------------------------------



        double x0_fp = (A + B) / 2.0; // Начальное приближение

        //проверка условия сходимости
        double q_1 = 1.0 / (2.0 * (B + 1.0));
        double q_2 = 1.0 / (2.0 * (A + 1.0)); // g'(x) = 1/(2(x+1)). Max на x=A=0

        double abs_q1 = std::fabs(q_1);
        double abs_q2 = std::fabs(q_2);
        double q_max = std::max(abs_q1, abs_q2);
        
        if (q_max >= 1.0) {
            throw std::runtime_error("Метод простой итерации: Условие сходимости |g'(x)| < 1 не выполнено (q_max=" + std::to_string(q_max) + ").");
        }
        std::cout << "Условие сжатия выполнено. Коэффициент сжатия q_max = " << q_max << endl;

        double x_fp = g(x0_fp); //первая итерация
        int iterations_fp = 1;
        vector<double> fp_errors_diff; 
        vector<double> fp_errors_bound;

        // Добавляем результаты первой итерации
        double diff_k1 = std::fabs(x_fp - x0_fp);
        fp_errors_diff.push_back(diff_k1);

        // Апостериорная оценка для k=1 (использует |x1 - x0|)
        double error_bound_k1 = (q_max / (1.0 - q_max)) * diff_k1;
        fp_errors_bound.push_back(error_bound_k1);

        
        while (std::fabs(error_bound_k1) >= eps) {
            x0_fp = x_fp;
            x_fp = g(x0_fp);
            
            
            if (x_fp < A || x_fp > B) {
                throw std::runtime_error("Метод простой итерации: Нарушено условие g(x) in [A, B]. Итерация вышла за границы.");
            }
            
            
            double diff_k = std::fabs(x_fp - x0_fp);
            error_bound_k1 = (q_max / (1.0 - q_max)) * diff_k;
            
            
            fp_errors_diff.push_back(diff_k);
            fp_errors_bound.push_back(error_bound_k1);
            
            
            iterations_fp++;
            
            if (iterations_fp > 1000) { 
                std::cout << "Итерации не сошлись за 1000 шагов (метод простой итерации)." << std::endl;
                return;
            }
        }

        std::cout << "\n--- Метод простой итерации ---\n";
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "Положительный корень примерно равен: " << x_fp << endl;
        std::cout << "Проверка: f(x) = " << f(x_fp) << endl;
        std::cout << "Количество итераций: " << iterations_fp << endl;

        std::cout << "\nЗависимость погрешности от количества итераций:" << endl;
        for (size_t i = 0; i < fp_errors_diff.size(); ++i) {
            std::cout << "Итерация " << (i + 1) << ": погрешность абсолютная = " << fp_errors_diff[i] 
            << "; погрешность из методички: " << fp_errors_bound[i] << endl;
        }

        

        

        // ---------------------------------------------------
        // Метод Ньютона
        // ---------------------------------------------------
        
        double x0_newton = (A + B) / 2.0; // Начальное приближение

        double f_at_x0 = f(x0_newton);
        double d2f_at_x0 = d2f(x0_newton);

        //проверка условия сходимости
        if (f_at_x0 * d2f_at_x0 <= 0.0) {
            std::cerr << "Предупреждение: Условие сходимости f(x0) * f''(x0) > 0 нарушено или близко к нулю.\n";
            std::cerr << "f(x0)=" << f_at_x0 << ", f''(x0)=" << d2f_at_x0 << ". Сходимость не гарантируется.\n";
        } else {
             std::cout << "Условие сходимости f(x0) * f''(x0) > 0 выполнено." << endl;
        }

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

            newton_errors.push_back(std::fabs(x_newton - x0_newton));

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
