#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include "Matrix.hpp"
#include "Vector.hpp"

class Task_4_2 {
private:
    static constexpr double A_INTERVAL = 0.0;
    static constexpr double B_INTERVAL = 1.0;
    static constexpr double Y_PRIME_A = 3.0 / 4.0; // 3/4
    static constexpr double Y_PRIME_B = (M_E * M_E * (M_E + 2.0)) / ((M_E + 1.0) * (M_E + 1.0)); // e^2(e + 2) / (e+1)^2

    // (eˣ+1)y'' - 2y' - eˣy = 0
    static Vector ode_system_rhs(double x, const Vector& y_vec) {
        Vector res(2);
        double y = y_vec.Get(0);
        double z = y_vec.Get(1); // z = y'
        res.Set(0, z);
        res.Set(1, (2.0 * z + exp(x) * y) / (exp(x) + 1.0));
        return res;
    }

    // (eˣ+1)y'' - 2y' - eˣy = 0
    static double exact_solution(double x) {
        return exp(x) - 1.0 + 1.0 / (exp(x) + 1.0);
    }

    // --- Вспомогательные методы ---

    // метод прогонки
    static Vector solve_tridiagonal_system(const Matrix& A, const Vector& B) {
        int n = A.GetLength();
        if (!A.isTriDiagonal()) throw std::invalid_argument("Матрица не является трехдиагональной.");

        Vector alpha(n), beta(n);
        double z = A.Get(0, 0);
        if (std::abs(z) < 1e-12) throw std::runtime_error("Деление на ноль в методе прогонки.");
        alpha.Set(0, -A.Get(0, 1) / z);
        beta.Set(0, B.Get(0) / z);

        for (int i = 1; i < n - 1; i++) {
            z = A.Get(i, i) + A.Get(i, i - 1) * alpha.Get(i - 1);
            if (std::abs(z) < 1e-12) throw std::runtime_error("Деление на ноль в методе прогонки.");
            alpha.Set(i, -A.Get(i, i + 1) / z);
            beta.Set(i, (B.Get(i) - A.Get(i, i - 1) * beta.Get(i - 1)) / z);
        }

        Vector X(n);
        z = A.Get(n - 1, n - 1) + A.Get(n - 1, n - 2) * alpha.Get(n - 2);
        if (std::abs(z) < 1e-12) throw std::runtime_error("Деление на ноль в методе прогонки.");
        X.Set(n - 1, (B.Get(n - 1) - A.Get(n - 1, n - 2) * beta.Get(n - 2)) / z);

        for (int i = n - 2; i >= 0; i--) X.Set(i, alpha.Get(i) * X.Get(i + 1) + beta.Get(i));
        return X;
    }

    // рунге-кутты 4 порядка
    static Vector solve_ode_system_rk4(double x0, double x_end, double h, const Vector& y0, const std::function<Vector(double, const Vector&)>& f) {

        // y0 - ВЕКТОР начальных условий [y(x0), z(x0), ...]
        Vector y_current = y0;
        double x_current = x0;
        int vec_size = y0.GetSize();
        int n_steps = static_cast<int>(round((x_end - x0) / h));

        for (int i = 0; i < n_steps; ++i) {

            // ---= РЕАЛИЗАЦИЯ ОДНОГО ШАГА МЕТОДА РУНГЕ-КУТТЫ 4-го ПОРЯДКА (Формулы 4.10, 4.16) =---
            // k1, k2, k3, k4 - это ВЕКТОРЫ коэффициентов.
            // k1.Get(0) соответствует K1 для 'y', k1.Get(1) - L1 для 'z' и т.д.

            // Шаг 1: Вычисляем ВЕКТОР k1 на основе наклона в НАЧАЛЬНОЙ точке.
            // Формула: k1_vec = f(x_k, y_k_vec)
            Vector k1 = f(x_current, y_current);

            // Шаг 2: Вычисляем ВЕКТОР k2 на основе наклона в СРЕДНЕЙ точке.
            // Сначала готовим временный вектор y_temp = y_k_vec + 0.5 * h * k1_vec
            Vector y_temp(vec_size);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + 0.5 * h * k1.Get(j));

            // Формула: k2_vec = f(x_k + h/2, y_temp)
            Vector k2 = f(x_current + 0.5 * h, y_temp);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + 0.5 * h * k2.Get(j));

            // Шаг 3: Вычисляем ВЕКТОР k3 на основе УТОЧНЕННОГО наклона в СРЕДНЕЙ точке.
            // Готовим y_temp = y_k_vec + 0.5 * h * k2_vec
            Vector k3 = f(x_current + 0.5 * h, y_temp);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + h * k3.Get(j));

            // Шаг 4: Вычисляем ВЕКТОР k4 на основе наклона в КОНЕЧНОЙ точке.
            // Готовим y_temp = y_k_vec + h * k3_vec
            Vector k4 = f(x_current + h, y_temp);

            // Финальный шаг: Обновляем ВЕКТОР y_current.
            // Это применение взвешенной суммы из формулы (4.10) к КАЖДОМУ компоненту системы.
            // y_{k+1}[j] = y_k[j] + (h/6) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
            for(int j=0; j < vec_size; ++j) y_current.Set(j, y_current.Get(j) + (h / 6.0) * (k1.Get(j) + 2.0*k2.Get(j) + 2.0*k3.Get(j) + k4.Get(j)));
            x_current += h;
        }
        return y_current;
    }
    
    static std::pair<std::vector<double>, Vector> solve_shooting_method(double h, bool silent = false) {
        if (!silent) std::cout << "\n\n--- 1. Метод стрельбы ---\n";


        // === ШАГ 1: Подготовка. Превращаем краевую задачу в задачу поиска корня ===
        //
        // ИДЕЯ: Мы не знаем начальное условие y(0), но знаем y'(0). Давайте угадаем y(0).
        // Обозначим наш "угадываемый" y(0) греческой буквой 'eta' (η).
        // Теперь у нас есть ЗАДАЧА КОШИ с начальными условиями:
        //   y(0) = η
        //   y'(0) = Y_PRIME_A (задано в условии)
        //
        // Решив эту задачу Коши (например, методом РК4), мы получим некоторое значение y'(1) на правом конце.
        // НАША ЦЕЛЬ: найти такое η*, чтобы полученное y'(1) совпало с заданным в условии Y_PRIME_B.
        //
        // Это эквивалентно поиску корня уравнения Φ(η) = 0, где:
        //   Φ(η) = y'(1, η) - Y_PRIME_B         (см. формулу 4.33 в методичке)
        //
        // Для поиска корня будем использовать МЕТОД СЕКУЩИХ (см. формулу 4.34),
        // так как производную dΦ/dη найти аналитически невозможно.

        double epsilon = 1e-7;


        // --- Делаем два "пристрелочных" выстрела с начальными догадками η₀ и η₁ ---
        
        // Первый "выстрел": берем произвольное η₀ = 1.0
        double eta0 = 1.0, eta1 = 0.0; 

        Vector y0_eta0(2);
        y0_eta0.Set(0, eta0); // y(0) = η₀
        y0_eta0.Set(1, Y_PRIME_A); // y'(0) - задано

        // Решаем задачу Коши от 0 до 1 и получаем вектор [y(1), y'(1)]
        Vector y_end_eta0 = solve_ode_system_rk4(A_INTERVAL, B_INTERVAL, h, y0_eta0, ode_system_rhs);

        // Вычисляем невязку Φ(η₀)
        double phi0 = y_end_eta0.Get(1) - Y_PRIME_B;

        // Второй "выстрел": берем произвольное η₁ = 0.0
        Vector y0_eta1(2);
        y0_eta1.Set(0, eta1); // y(0) = η₁
        y0_eta1.Set(1, Y_PRIME_A); // y'(0) - задано
        Vector y_end_eta1 = solve_ode_system_rk4(A_INTERVAL, B_INTERVAL, h, y0_eta1, ode_system_rhs);
        // Вычисляем невязку Φ(η₁)
        double phi1 = y_end_eta1.Get(1) - Y_PRIME_B;

        // === ШАГ 2: Итерационный процесс поиска корня (Метод секущих) ===
        
        int j = 1; // Счетчик итераций
        double eta_curr = eta1, eta_prev = eta0;
        double phi_curr = phi1, phi_prev = phi0;

        // Итерируем, пока невязка |Φ| не станет достаточно малой
        while (std::abs(phi_curr) > epsilon) {
            j++;

            // Защита от деления на ноль, если метод "застрял"
            if (std::abs(phi_curr - phi_prev) < 1e-12) throw std::runtime_error("Метод стрельбы не сходится.");

            // Формула метода секущих (аналог 4.34): η_{j+1} = η_j - Φ(η_j) * (η_j - η_{j-1}) / (Φ(η_j) - Φ(η_{j-1}))
            double eta_next = eta_curr - phi_curr * (eta_curr - eta_prev) / (phi_curr - phi_prev);
            
            // Обновляем значения для следующей итерации
            eta_prev = eta_curr;
            phi_prev = phi_curr;
            eta_curr = eta_next;

            Vector y0_eta_curr(2);
            y0_eta_curr.Set(0, eta_curr);
            y0_eta_curr.Set(1, Y_PRIME_A);
            Vector y_end_eta_curr = solve_ode_system_rk4(A_INTERVAL, B_INTERVAL, h, y0_eta_curr, ode_system_rhs);
            phi_curr = y_end_eta_curr.Get(1) - Y_PRIME_B;
            if (j > 100) throw std::runtime_error("Превышено число итераций в методе стрельбы.");
        }
        
        if (!silent) std::cout << "Найденное начальное значение y(0) = " << eta_curr << "\n";

        int n_steps = static_cast<int>(round((B_INTERVAL - A_INTERVAL) / h));
        std::vector<double> x_grid(n_steps + 1);
        Vector y_solution(n_steps + 1);
        Vector y_current(2);
        y_current.Set(0, eta_curr);
        y_current.Set(1, Y_PRIME_A);
        double x_current = A_INTERVAL;
        x_grid[0] = x_current;
        y_solution.Set(0, y_current.Get(0));
        
        for (int i = 0; i < n_steps; ++i) {
            y_current = solve_ode_system_rk4(x_current, x_current + h, h, y_current, ode_system_rhs);
            x_current += h;
            x_grid[i + 1] = x_current;
            y_solution.Set(i + 1, y_current.Get(0));
        }
        return {x_grid, y_solution};
    }



    static std::pair<std::vector<double>, Vector> solve_finite_difference_method(double h, bool silent = false) {
        if (!silent) std::cout << "\n\n--- 2. Конечно-разностный метод ---\n";

        // ---= ШАГ 1: Подготовка сетки и системы =---
        int N = static_cast<int>(round((B_INTERVAL - A_INTERVAL) / h));
        int n_system = N + 1; // Количество неизвестных y_0, y_1, ..., y_N
        
        // Приводим исходное уравнение к каноническому виду y'' + p(x)y' + q(x)y = f(x)
        // Исходное: (e^x+1)y'' - 2y' - e^x*y = 0
        // y'' - [2/(e^x+1)]*y' - [e^x/(e^x+1)]*y = 0
        auto p = [](double x) { return -2.0 / (exp(x) + 1.0); };
        auto q = [](double x) { return -exp(x) / (exp(x) + 1.0); };
        
        // Создаем матрицу A и вектор B для системы Ay=B
        Matrix A(n_system);
        Vector B(n_system);

        // --- Уравнение для k=0 (левая граница) ---
        // Используем y_{-1} = y_1 - 2*h*y'(0)
        // ---= ШАГ 2: Формирование уравнений системы =---
        // Мы будем использовать центральные разности (формула 4.37) для аппроксимации
        // производных y' и y'', которые имеют второй порядок точности O(h^2).
        // y'_k ≈ (y_{k+1} - y_{k-1}) / (2h)
        // y''_k ≈ (y_{k+1} - 2y_k + y_{k-1}) / h^2
        // Подстановка их в y''+p*y'+q*y=f дает уравнение (4.39) из методички:
        // (1-p*h/2)*y_{k-1} + (-2+h^2*q)*y_k + (1+p*h/2)*y_{k+1} = h^2*f

        // --- Уравнение для k=0 (левая граница) ---
        // Проблема: для центральной разности в точке k=0 нужна фиктивная точка y_{-1}.
        // Решение: выражаем y_{-1} из аппроксимации производной на границе.
        // y'(0) ≈ (y_1 - y_{-1}) / (2h)  =>  y_{-1} = y_1 - 2*h*y'(0)
        // Подставляем этот y_{-1} в общее уравнение для k=0 и группируем члены при y_0 и y_1.
        double x0 = A_INTERVAL;
        A.Set(0, 0, -2.0 + h*h*q(x0)); // Коэффициент при y_0
        A.Set(0, 1, 2.0); // Коэффициент при y_1 (получается из (1-p*h/2)*y_1 + (1+p*h/2)*y_1)
        // В правую часть переносятся все известные члены, включая y'(0) (Y_PRIME_A).
        // После упрощений правая часть B[0] для уравнения (1-p*h/2)*y_{-1} + ... = 0
        B.Set(0, h*h*p(x0)*Y_PRIME_A - p(x0)*h*( -2.0*h*Y_PRIME_A)); // Упрощается до (2*h - h*h*p(x0)) * p(x0)/(-2) ...
                                                                    // Проще: B[0] = 2*h*Y_PRIME_A*(1 - p(x0)*h/2) - h*h*...
        // Упрощенная и проверенная формула для правой части:
        B.Set(0, (2.0*h*Y_PRIME_A - h*h*p(x0)*Y_PRIME_A));

        // --- Внутренние узлы (k=1...N-1) ---
         // Здесь мы напрямую используем формулу (4.39) из методички.
        for (int k = 1; k < N; ++k) {
            double xk = A_INTERVAL + k * h;

            // Коэффициент при y_{k-1}
            A.Set(k, k - 1, 1.0 - p(xk) * h / 2.0);

            // Коэффициент при y_k
            A.Set(k, k, -2.0 + h * h * q(xk));

            // Коэффициент при y_{k+1}
            A.Set(k, k + 1, 1.0 + p(xk) * h / 2.0);

            // Правая часть (f(x) у нас равно 0)
            B.Set(k, 0.0);
        }

        // --- Уравнение для k=N (правая граница) ---
        // Используем y_{N+1} = y_{N-1} + 2*h*y'(N)
        double xN = B_INTERVAL;
        A.Set(N, N - 1, 2.0); // Коэффициент при y_0
        A.Set(N, N, -2.0 + h*h*q(xN)); // Коэффициент при y_1
        // Правая часть:
        B.Set(N, -(2.0*h*Y_PRIME_B + h*h*p(xN)*Y_PRIME_B));
        
        Vector y_solution = solve_tridiagonal_system(A, B);
        
        std::vector<double> x_grid(n_system);
        for(int i=0; i < n_system; ++i) x_grid[i] = A_INTERVAL + i * h;
        
        return {x_grid, y_solution};
    }
    
    static void print_error_report(const std::string& method_name, const std::vector<double>& x_grid, const Vector& y_numerical) {
        std::cout << "\nТаблица результатов и погрешности для '" << method_name << "':\n";
        std::cout << std::setw(15) << "x_k" << std::setw(20) << "y_numerical" << std::setw(20) << "y_exact" << std::setw(20) << "Error\n";
        std::cout << std::string(75, '-') << "\n";
        
        double max_error = 0.0;
        for (size_t i = 0; i < x_grid.size(); ++i) {
            double x = x_grid[i];
            double y_num = y_numerical.Get(i);
            double y_ex = exact_solution(x);
            double error = std::abs(y_num - y_ex);
            if (error > max_error) max_error = error;
            
            std::cout << std::fixed << std::setprecision(8) 
                      << std::setw(15) << x 
                      << std::setw(20) << y_num 
                      << std::setw(20) << y_ex 
                      << std::setw(20) << error << "\n";
        }
        std::cout << std::string(75, '-') << "\n";
        std::cout << "Максимальная погрешность (сравнение с точным решением): " << max_error << "\n";
    }

    static void runge_romberg_error(const std::string& method_name, double h, int p, 
        const std::function<std::pair<std::vector<double>, Vector>(double, bool)>& solver)
    {
        auto [x_h, y_h] = solver(h, true);
        auto [x_2h, y_2h] = solver(2*h, true);
        
        double max_rr_error = 0.0;
        for (int i = 0; i < y_2h.GetSize(); ++i) {
            double y_h_val = y_h.Get(i * 2);
            double y_2h_val = y_2h.Get(i);
            double error_estimate = std::abs(y_h_val - y_2h_val) / (pow(2.0, p) - 1.0);
            if (error_estimate > max_rr_error) {
                max_rr_error = error_estimate;
            }
        }
        std::cout << "Оценка погрешности по Рунге-Ромберга для '" << method_name << "' (p=" << p << "): " << max_rr_error << "\n";
    }

public:
    static void Do() {
        double h = 0.1;

        std::cout << "========== РЕШЕНИЕ КРАЕВОЙ ЗАДАЧИ ==========\n";
        std::cout << "Уравнение: (e^x + 1)y'' - 2y' - e^x*y = 0\n";
        std::cout << "Интервал: [" << A_INTERVAL << ", " << B_INTERVAL << "]\n";
        std::cout << "Граничные условия: y'(0)=" << Y_PRIME_A << ", y'(1)=" << Y_PRIME_B << "\n";
        std::cout << "Шаг h = " << h << "\n";

        auto [x_grid_shoot, y_shoot] = solve_shooting_method(h);
        print_error_report("Метод стрельбы", x_grid_shoot, y_shoot);
        runge_romberg_error("Метод стрельбы", h, 4, solve_shooting_method);

        auto [x_grid_fd, y_fd] = solve_finite_difference_method(h);
        print_error_report("Конечно-разностный метод", x_grid_fd, y_fd);
        runge_romberg_error("Конечно-разностный метод", h, 2, solve_finite_difference_method);
    }
};

int main() {
    try {
        Task_4_2::Do();
    } catch (const std::exception& e) {
        std::cerr << "\nПРОИЗОШЛА КРИТИЧЕСКАЯ ОШИБКА: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}