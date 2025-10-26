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

        int num_intervals = static_cast<int>(round((xk - x0) / h));
        if (std::abs(num_intervals * h - (xk - x0)) > 1e-9) {
            throw std::runtime_error("Отрезок не делится нацело на шаг h = " + std::to_string(h));
        }


        Vector x(num_intervals + 1);
        Vector y(num_intervals + 1);
        for (int i = 0; i <= num_intervals; ++i) {
            double current_x = x0 + i * h;
            x.Set(i, current_x);
            y.Set(i, f(current_x));
        }


        double rect_res = MethodRectangles(h, num_intervals);
        double trap_res = MethodTrapezoids(h, num_intervals, y);
        double simp_res = MethodSimpson(h, num_intervals, y);

        std::cout << "  Метод прямоугольников: " << rect_res << "\n";
        std::cout << "  Метод трапеций:        " << trap_res << "\n";
        std::cout << "  Метод Симпсона:        " << simp_res << "\n";

        return {rect_res, trap_res, simp_res};
    }


    double MethodRectangles(double h, int n) const {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += f(x0 + i * h + h / 2.0);
        }
        return h * sum;
    }


    double MethodTrapezoids(double h, int n, const Vector& y) const {
        double sum = 0;
        for (int i = 1; i < n; ++i) {
            sum += y.Get(i);
        }
        sum += (y.Get(0) + y.Get(n)) / 2.0;
        return h * sum;
    }


    double MethodSimpson(double h, int n, const Vector& y) const {
        if (n % 2 != 0) {
            return NAN;
        }
        double sum = y.Get(0) + y.Get(n);
        for (int i = 1; i < n; i += 2) {
            sum += 4 * y.Get(i);
        }
        for (int i = 2; i < n; i += 2) {
            sum += 2 * y.Get(i);
        }
        return (h / 3.0) * sum;
    }
    

    void AnalyzeRungeRomberg(const IntegralResults& coarse_res, const IntegralResults& fine_res) const {
        double k = h1 / h2;

        std::cout << std::left << std::setw(25) << "Метод";
        std::cout << std::left << std::setw(40) << "Абс. погрешность"
                  << std::right << std::setw(30) << "Уточненное значение" << "\n";
        std::cout << "-------------------------------------------------------------------\n";
        
        // --- Прямоугольники (p=2) ---
        double p_rect = 2.0;
        double err_rect = (fine_res.rect - coarse_res.rect) / (pow(k, p_rect) - 1.0);
        std::cout << std::left << std::setw(20) << "Прямоугольников";
        std::cout << std::right << std::setw(19) << std::abs(err_rect) 
                  << std::right << std::setw(25) << (fine_res.rect + err_rect) << "\n";

        // --- Трапеции (p=2) ---
        double p_trap = 2.0;
        double err_trap = (fine_res.trap - coarse_res.trap) / (pow(k, p_trap) - 1.0);
        std::cout << std::left << std::setw(20) << "Трапеций";
        std::cout << std::right << std::setw(22) << std::abs(err_trap) 
                  << std::right << std::setw(25) << (fine_res.trap + err_trap) << "\n";

        // --- Симпсон (p=4) ---
        double p_simp = 4.0;
        double err_simp = (fine_res.simp - coarse_res.simp) / (pow(k, p_simp) - 1.0);
        std::cout << std::left << std::setw(20) << "Симпсона";
        std::cout << std::right << std::setw(22) << std::abs(err_simp) 
                  << std::right << std::setw(25) << (fine_res.simp + err_simp) << "\n";
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