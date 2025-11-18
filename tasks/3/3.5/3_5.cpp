#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <functional>
#include "Vector.hpp"

struct IntegralResults {
    double rect, trap, simp;
};

class Task3_5 {
public:
    Task3_5(std::function<double(double)> func, double start, double end, double h_1, double h_2)
        : f(func), x0(start), xk(end), h1(h_1), h2(h_2) {
    }

    void Do() {
        std::cout << std::fixed << std::setprecision(5);
        PrintInitialData();

        std::cout << "\n--- 1. Вычисление интеграла с шагом h1 = " << h1 << " ---\n";
        IntegralResults res_h1 = CalculateForStep(h1);

        std::cout << "\n--- 2. Вычисление интеграла с шагом h2 = " << h2 << " ---\n";
        IntegralResults res_h2 = CalculateForStep(h2);

        std::cout << "\n--- 3. Уточнение по методу Рунге-Ромберга ---\n";
        AnalyzeRungeRomberg(res_h1, res_h2);
    }

private:
    std::function<double(double)> f;
    double x0, xk, h1, h2;


    IntegralResults CalculateForStep(double h) {
        int n = static_cast<int>(round((xk - x0) / h));
        if (std::abs(n * h - (xk - x0)) > 1e-9) {
            throw std::runtime_error("Отрезок не делится нацело на шаг h = " + std::to_string(h));
        }

        Vector x(n + 1);
        Vector y(n + 1);
        for (int i = 0; i <= n; ++i) {
            double cx = x0 + i * h;
            x.Set(i, cx);
            y.Set(i, f(cx));
        }

        std::cout << std::setw(5) << "i" << std::setw(12) << "x_i" << std::setw(12) << "y_i"
                  << std::setw(20) << "           Прямоугольники" << std::setw(15) << "     Трапеции" << std::setw(15) << "    Симпсон" << "\n";
        std::cout << "-----------------------------------------------------------------------------\n";

        double rect_sum = 0, trap_sum = 0, simp_sum = 0;

        for (int i = 0; i <= n; ++i) {
            std::cout << std::setw(5) << i << std::setw(12) << x.Get(i) << std::setw(12) << y.Get(i);

            // ---= 1. МЕТОД ПРЯМОУГОЛЬНИКОВ =---
            // Формула: h * Σ f( (x_{i-1}+x_i)/2 )
            // Мы вычисляем значение функции в СЕРЕДИНЕ интервала, который НАЧИНАЕТСЯ с точки x_i.
            // Поэтому цикл накопления идет до предпоследнего узла (i < n).
            if (i < n) {
                rect_sum += f(x.Get(i) + h / 2.0);
                std::cout << std::setw(20) << h * rect_sum;
            } else {
                 std::cout << std::setw(20) << h * rect_sum;
            }
            
            // ---= 2. МЕТОД ТРАПЕЦИЙ =---
            // Формула: h * [ (y_0 + y_n)/2 + Σ y_i ]
            // Мы добавляем к сумме значения в узлах с разными весами.
            if (i == 0) trap_sum += y.Get(i) / 2.0;
            else if (i == n) trap_sum += y.Get(i) / 2.0;
            else trap_sum += y.Get(i);
            std::cout << std::setw(15) << h * trap_sum;

            // ---= 3. МЕТОД СИМПСОНА =---
            // Формула: (h/3) * [ y_0 + 4y_1 + 2y_2 + ... + 4y_{n-1} + y_n ]
            // Проверяем, что количество интервалов (n) - четное.
            if (n % 2 == 0) {
                if (i == 0 || i == n) simp_sum += y.Get(i);
                else if (i % 2 != 0) simp_sum += 4 * y.Get(i);
                else simp_sum += 2 * y.Get(i);
                std::cout << std::setw(15) << (h / 3.0) * simp_sum;
            } else {
                 if (i==0) std::cout << std::setw(15) << " (n нечетно)";
                 else std::cout << std::setw(15) << " ";
            }
            std::cout << "\n";
        }
        
        return {h * rect_sum, h * trap_sum, (n % 2 == 0) ? (h / 3.0) * simp_sum : NAN};
    }


    

    void AnalyzeRungeRomberg(const IntegralResults& coarse_res, const IntegralResults& fine_res) const {
        double k = h1 / h2;

        // --- Константы для теоретической оценки ---
        const double M2 = 0.00676; // max|f''(x)| на [-2, 2]
        const double M4 = 1.309;   // max|f''''(x)| на [-2, 2]
        const double b_minus_a = xk - x0;
        const double exact_value = -0.2040185; // Посчитано с высокой точностью

        
        std::cout << std::left << std::setw(20) << "Метод" 
                  << std::right << std::setw(18) << "Оценка по Р-Р"
                  << std::right << std::setw(22) << "        Теоретическая оценка"
                  << std::right << std::setw(20) << "         Истинная ошибка" << "\n";
        std::cout << "----------------------------------------------------------------------------------\n";
        
        // --- Прямоугольники (p=2) ---
        {
            double err_rr = std::abs((fine_res.rect - coarse_res.rect) / (pow(k, 2) - 1.0));
            double refined_value = fine_res.rect + err_rr;
            // Формула (3.24): R <= (b-a)/24 * h^2 * M2
            double err_theory = (b_minus_a / 24.0) * pow(h2, 2) * M2;
            double err_true = std::abs(exact_value - fine_res.rect);
            std::cout << std::left << std::setw(20) << "Прямоугольников"
                      << std::right << std::setw(18) << err_rr
                      << std::right << std::setw(22) << err_theory
                      << std::right << std::setw(20) << err_true
                      << std::right << std::setw(20) << refined_value << "\n";
        }

        
        {
               
            double err_rr = std::abs((fine_res.trap - coarse_res.trap) / (pow(k, 2) - 1.0));
            double refined_value = fine_res.trap + err_rr;
            // Формула (3.26): R <= (b-a)/12 * h^2 * M2 (в методичке опечатка, должно быть (b-a) * h^2 / 12)
            double err_theory = (b_minus_a / 12.0) * pow(h2, 2) * M2;
            double err_true = std::abs(exact_value - fine_res.trap);
            std::cout << std::left << std::setw(20) << "Трапеций"
                      << std::right << std::setw(18) << err_rr
                      << std::right << std::setw(22) << err_theory
                      << std::right << std::setw(20) << err_true
                      << std::right << std::setw(20) << refined_value << "\n";
        }

        // --- Симпсон (p=4) ---
        {
            
            double err_rr = std::abs((fine_res.simp - coarse_res.simp) / (pow(k, 4) - 1.0));
            double refined_value = fine_res.simp + err_rr;
            // Формула (3.29): R <= (b-a)/180 * h^4 * M4
            double err_theory = (b_minus_a / 180.0) * pow(h2, 4) * M4;
            double err_true = std::abs(exact_value - fine_res.simp);
            std::cout << std::left << std::setw(20) << "Симпсона"
                      << std::right << std::setw(18) << err_rr
                      << std::right << std::setw(22) << err_theory
                      << std::right << std::setw(20) << err_true
                      << std::right << std::setw(20) << refined_value << "\n";
            
        }
    }


    void PrintInitialData() const {
        std::cout << "--- Исходные данные (Вариант 13) ---\n";
        std::cout << "Функция: y = x^2 / (x^3 - 27)\n";
        std::cout << "Отрезок: [" << x0 << ", " << xk << "]\n";
        std::cout << "Шаги: h1 = " << h1 << ", h2 = " << h2 << "\n";
    }
};

int main() {
    try {
        auto func = [](double x) {
            return (x * x) / (x * x * x - 27.0);
        };
        
        double x_start = -2.0;
        double x_end = 2.0;
        double h1 = 1.0;
        double h2 = 0.5;

        Task3_5 task(func, x_start, x_end, h1, h2);
        task.Do();

    } catch (const std::exception& e) {
        std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}