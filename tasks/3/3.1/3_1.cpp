#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include <fstream>
#include "Vector.hpp"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Task3_1 {
public:
    Task3_1(std::function<double(double)> func, const Vector& nodes, double point)
        : f(func), X(nodes), x_star(point) {
        
        if (X.GetSize() == 0) {
            throw std::invalid_argument("Вектор узлов не может быть пустым.");
        }
        Y = Vector(X.GetSize());
        for (int i = 0; i < X.GetSize(); ++i) {
            Y.Set(i, f(X.Get(i)));
        }
    }

    void Do(const std::string& plot_filename) {
        std::cout << std::fixed << std::setprecision(8);
        
        PrintInitialData();

        double lagrange_result = SolveLagrange();
        
        // Вычисляем коэффициенты и значение Ньютона ОДИН РАЗ
        Vector newton_coeffs;
        double newton_result = SolveNewton(newton_coeffs);
        
        if (std::abs(lagrange_result - newton_result) > 1e-9) {
            std::cerr << "Внимание: Результаты методов Лагранжа и Ньютона не совпадают!" << std::endl;
        }

        // Передаем уже вычисленные коэффициенты в другие функции
        PrintNewtonPolynomial(newton_coeffs);
        GeneratePlotData(plot_filename, newton_coeffs);

        CalculateError(lagrange_result);
    }

private:
    std::function<double(double)> f;
    Vector X;
    Vector Y;
    double x_star;

    double EvaluateNewton(double x, const Vector& coeffs) {
        double result = coeffs.Get(0);
        double term = 1.0;
        for (int i = 1; i < coeffs.GetSize(); ++i) {
            term *= (x - X.Get(i - 1));
            result += coeffs.Get(i) * term;
        }
        return result;
    }

    void PrintInitialData() {
        std::cout << "--- Исходные данные ---\n";
        std::cout << "Функция: y = cos(x) + x\n";
        std::cout << "Точка для вычисления погрешности X* = " << x_star << "\n";
        std::cout << "Узлы интерполяции:\n";
        for (int i = 0; i < X.GetSize(); ++i) {
            std::cout << "  X[" << i << "] = " << std::setw(10) << X.Get(i) 
                      << "  |  Y[" << i << "] = " << std::setw(10) << Y.Get(i) << "\n";
        }
        std::cout << "------------------------\n\n";
    }

    double SolveLagrange() {
        std::cout << "--- 1. Метод Лагранжа ---\n";
        double sum = 0.0;
        for (int i = 0; i < X.GetSize(); ++i) {
            double basis_polynomial = 1.0;
            for (int j = 0; j < X.GetSize(); ++j) {
                if (i != j) {
                    basis_polynomial *= (x_star - X.Get(j)) / (X.Get(i) - X.Get(j));
                }
            }
            sum += Y.Get(i) * basis_polynomial;
        }
        std::cout << "Значение многочлена Лагранжа L(X*) в точке X* = " << x_star << " равно: " << sum << "\n";
        std::cout << "---------------------------\n\n";
        return sum;
    }
    
    double SolveNewton(Vector& coeffs) {
        std::cout << "--- 2. Метод Ньютона ---\n";
        int n = X.GetSize();
        
        std::vector<Vector> diffs(n);
        for(int i = 0; i < n; ++i) {
            //разделенная разность 0-го порядка
            diffs[i] = Vector(n - i);
            diffs[i].Set(0, Y.Get(i));
        }
        //j определяет порядок разделенной разности
        for (int j = 1; j < n; ++j) {
            for (int i = 0; i < n - j; ++i) {
                double val = (diffs[i+1].Get(j-1) - diffs[i].Get(j-1)) / (X.Get(i+j) - X.Get(i));
                diffs[i].Set(j, val);
            }
        }
        
        // ВОЗВРАЩАЕМ ВЫВОД ТАБЛИЦЫ РАЗДЕЛЕННЫХ РАЗНОСТЕЙ
        std::cout << "Таблица разделенных разностей:\n";
        for (int j = 0; j < n; ++j) {
            std::cout << "  f[..]_" << j << ": ";
            for (int i = 0; i < n - j; ++i) {
                std::cout << std::setw(12) << diffs[i].Get(j);
            }
            std::cout << std::endl;
        }
        
        //записываем разделенные разности для формулы
        coeffs = Vector(n);
        for (int i = 0; i < n; ++i) {
            coeffs.Set(i, diffs[0].Get(i));
        }
        
        double result = EvaluateNewton(x_star, coeffs);
        
        std::cout << "\nЗначение многочлена Ньютона P(X*) в точке X* = " << x_star << " равно: " << result << "\n";
        std::cout << "---------------------------\n\n";
        return result;
    }

    void PrintNewtonPolynomial(const Vector& coeffs) {
        std::cout << "--- Вид интерполяционного многочлена Ньютона ---\n";
        std::cout << "P(x) = " << coeffs.Get(0);
        for (int i = 1; i < coeffs.GetSize(); ++i) {
            double c = coeffs.Get(i);
            std::cout << (c >= 0 ? " + " : " - ") << std::abs(c);
            for (int j = 0; j < i; ++j) {
                double node_val = X.Get(j);
                std::cout << " * (x" << (node_val >= 0 ? " - " : " + ") << std::abs(node_val) << ")";
            }
        }
        std::cout << "\n-------------------------------------------------\n\n";
    }

    //точки графика для сравнения точной функции и интерполяционного многочлена
    void GeneratePlotData(const std::string& filename, const Vector& coeffs) {
        std::cout << "--- 4. Генерация данных для графика ---\n";
        
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Не удалось создать файл " << filename << std::endl;
            return;
        }

        outFile << "x,actual_f(x),polynomial_p(x)\n";
        
        double start = X.Get(0) - 0.1; 
        double end = X.Get(X.GetSize() - 1) + 0.1;
        int steps = 200;
        
        for (int i = 0; i <= steps; ++i) {
            double x = start + (end - start) * i / steps;
            double y_actual = f(x);
            double y_poly = EvaluateNewton(x, coeffs);
            outFile << x << "," << y_actual << "," << y_poly << "\n";
        }

        outFile.close();
        std::cout << "Данные для построения графика сохранены в файл: " << filename << "\n";
        std::cout << "-------------------------------------------\n\n";
    }

    double CalculateOmega(double x) {
        double omega = 1.0;
        for (int i = 0; i < X.GetSize(); ++i) {
            omega *= (x - X.Get(i));
        }
        return omega;
    }

    void CalculateError(double approx_value) {
        std::cout << "--- 5. Вычисление погрешности ---\n";
        
        // 1. Фактическая (апостериорная) погрешность
        double actual_value = f(x_star);
        double error_actual = std::abs(actual_value - approx_value);
        
        std::cout << "Точное значение функции f(X*):    " << actual_value << "\n";
        std::cout << "Приближенное значение P(X*):    " << approx_value << "\n";
        std::cout << "Абсолютная погрешность |f(X*) - P(X*)|: " << error_actual << "\n";

        // 2. Теоретическая (априорная) оценка погрешности

        // n - степень многочлена, n+1 - количество узлов
        int n_plus_1 = X.GetSize(); // 4
        
        // M_{n+1} = max |f^(n+1)(xi)|
        // Для f(x)=cos(x)+x, f^(4)(x) = cos(x). Max |cos(x)| на [0, pi/2] равен 1.0
        const double M_n_plus_1 = 1.0; 
        
        // (n+1)!
        double factorial_n_plus_1 = 1.0;
        for (int k = 1; k <= n_plus_1; ++k) {
            factorial_n_plus_1 *= k;
        }

        // |omega_{n+1}(X*)| = | product (X* - Xi) |
        double omega_x_star = CalculateOmega(x_star);
        double abs_omega_x_star = std::abs(omega_x_star);
        
        // Оценка: |E_n(X*)| <= M_{n+1} / (n+1)! * |omega_{n+1}(X*)|
        double error_theoretical_estimate = (M_n_plus_1 / factorial_n_plus_1) * abs_omega_x_star;

        std::cout << "\n--- Теоретическая (априорная) оценка погрешности ---\n";
        std::cout << "  n+1 (количество узлов) = " << n_plus_1 << "\n";
        std::cout << "  M_{n+1} = max |f^(" << n_plus_1 << ")(xi)| = " << M_n_plus_1 << "\n";
        std::cout << "  (n+1)! = " << factorial_n_plus_1 << "\n";
        std::cout << "  |omega_{n+1}(X*)| = " << abs_omega_x_star << "\n";
        std::cout << "  Оценка: E_n(X*) <= M_{n+1}/(n+1)! * |omega_{n+1}(X*)| \n";
        std::cout << "  Теоретическая оценка E_n(X*) <= " << error_theoretical_estimate << "\n";
        
        // Проверка: фактическая ошибка должна быть меньше или равна оценке
        if (error_actual > error_theoretical_estimate + 1e-10) {
            std::cerr << "Внимание: Фактическая ошибка (" << error_actual << ") ПРЕВЫШАЕТ теоретическую оценку (" << error_theoretical_estimate << ")!\n";
            std::cerr << "Проверьте вычисление M_{n+1} или узлы.\n";
        } else {
            std::cout << "  Проверка: Фактическая ошибка меньше теоретической оценки. (Корректно)\n";
        }
        
        std::cout << "------------------------------------\n";
    }
};

int main() {
    try {
        auto func = [](double x) { return cos(x) + x; };
        double x_star = 1.0;

        
        std::cout << "***************************************************\n";
        std::cout << "***     ВАРИАНТ 13, ПУНКТ а) РАВНОМЕРНАЯ СЕТКА    ***\n";
        std::cout << "***************************************************\n\n";
        
        Vector nodes_a(4);
        nodes_a.Set(0, 0.0);
        nodes_a.Set(1, M_PI / 6.0);
        nodes_a.Set(2, 2.0 * M_PI / 6.0); // то же, что и PI/3
        nodes_a.Set(3, 3.0 * M_PI / 6.0); // то же, что и PI/2

        Task3_1 task_a(func, nodes_a, x_star);
        task_a.Do("./tasks/3/3.1/plot_data_a.csv");


        
        std::cout << "\n\n***************************************************\n";
        std::cout << "***   ВАРИАНТ 13, ПУНКТ б) НЕРАВНОМЕРНАЯ СЕТКА   ***\n";
        std::cout << "***************************************************\n\n";

        Vector nodes_b(4);
        nodes_b.Set(0, 0.0);
        nodes_b.Set(1, M_PI / 6.0);
        nodes_b.Set(2, M_PI / 4.0);
        nodes_b.Set(3, M_PI / 2.0);

        Task3_1 task_b(func, nodes_b, x_star);
        task_b.Do("./tasks/3/3.1/plot_data_b.csv");


        std::cout << "\n\n***************************************************\n";
        std::cout << "***       СВОДНЫЕ РЕЗУЛЬТАТЫ И ВЫВОДЫ         ***\n";
        std::cout << "***************************************************\n\n";

        std::cout << "\n\nВывод:\n";
        std::cout << "1. Равномерная сетка (вариант а) дала значительно более точный результат.\n";
        std::cout << "   Фактическая ошибка для нее почти в 5 раз меньше, чем для неравномерной сетки\n";
        std::cout << "2. Теоретические оценки погрешности также подтверждают, что набор узлов 'а'\n";
        std::cout << "   является более удачным для интерполяции в точке X* = 1.0.\n";

    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}