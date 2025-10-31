#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility> // for std::pair

// Подключаем ваши заголовочные файлы
#include "Matrix.hpp"
#include "Vector.hpp"

// --- Параметры задачи вынесены в пространство имен для удобства ---
namespace TaskConfig {
    const double A_INTERVAL = 0.0;
    const double B_INTERVAL = 1.0;
    const double Y_PRIME_A = 3.0 / 4.0;
    const double Y_PRIME_B = (exp(2.0) * (exp(1.0) + 2.0)) / pow(exp(1.0) + 1.0, 2.0);

    // Функция y'' = f(x, y, y') из (e^x+1)y'' - 2y' - e^x*y = 0
    Vector ode_system_rhs(double x, const Vector& y_vec) {
        Vector res(2);
        double y = y_vec.Get(0);
        double z = y_vec.Get(1); // z = y'
        res.Set(0, z);
        res.Set(1, (2.0 * z + exp(x) * y) / (exp(x) + 1.0));
        return res;
    }

    // Точное решение для сравнения
    double exact_solution(double x) {
        return exp(x) - 1.0 + 1.0 / (exp(x) + 1.0);
    }
}

class Task_4_2 {
private:
    // --- Вспомогательные методы, определены сразу ---

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

    static Vector solve_ode_system_rk4(double x0, double x_end, double h, const Vector& y0, const std::function<Vector(double, const Vector&)>& f) {
        Vector y_current = y0;
        double x_current = x0;
        int vec_size = y0.GetSize();
        int n_steps = static_cast<int>(round((x_end - x0) / h));

        for (int i = 0; i < n_steps; ++i) {
            Vector k1 = f(x_current, y_current);
            Vector y_temp(vec_size);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + 0.5 * h * k1.Get(j));
            Vector k2 = f(x_current + 0.5 * h, y_temp);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + 0.5 * h * k2.Get(j));
            Vector k3 = f(x_current + 0.5 * h, y_temp);
            for(int j=0; j < vec_size; ++j) y_temp.Set(j, y_current.Get(j) + h * k3.Get(j));
            Vector k4 = f(x_current + h, y_temp);
            for(int j=0; j < vec_size; ++j) y_current.Set(j, y_current.Get(j) + (h / 6.0) * (k1.Get(j) + 2.0*k2.Get(j) + 2.0*k3.Get(j) + k4.Get(j)));
            x_current += h;
        }
        return y_current;
    }
    
    static std::pair<std::vector<double>, Vector> solve_shooting_method(double h, bool silent = false) {
        if (!silent) std::cout << "\n\n--- 1. Метод стрельбы ---\n";

        double epsilon = 1e-7;
        double eta0 = 1.0, eta1 = 0.0; 

        Vector y0_eta0(2);
        y0_eta0.Set(0, eta0); 
        y0_eta0.Set(1, TaskConfig::Y_PRIME_A);
        Vector y_end_eta0 = solve_ode_system_rk4(TaskConfig::A_INTERVAL, TaskConfig::B_INTERVAL, h, y0_eta0, TaskConfig::ode_system_rhs);
        double phi0 = y_end_eta0.Get(1) - TaskConfig::Y_PRIME_B;

        Vector y0_eta1(2);
        y0_eta1.Set(0, eta1);
        y0_eta1.Set(1, TaskConfig::Y_PRIME_A);
        Vector y_end_eta1 = solve_ode_system_rk4(TaskConfig::A_INTERVAL, TaskConfig::B_INTERVAL, h, y0_eta1, TaskConfig::ode_system_rhs);
        double phi1 = y_end_eta1.Get(1) - TaskConfig::Y_PRIME_B;
        
        int j = 1;
        double eta_curr = eta1, eta_prev = eta0;
        double phi_curr = phi1, phi_prev = phi0;

        while (std::abs(phi_curr) > epsilon) {
            j++;
            if (std::abs(phi_curr - phi_prev) < 1e-12) throw std::runtime_error("Метод стрельбы не сходится.");
            double eta_next = eta_curr - phi_curr * (eta_curr - eta_prev) / (phi_curr - phi_prev);
            
            eta_prev = eta_curr;
            phi_prev = phi_curr;
            eta_curr = eta_next;

            Vector y0_eta_curr(2);
            y0_eta_curr.Set(0, eta_curr);
            y0_eta_curr.Set(1, TaskConfig::Y_PRIME_A);
            Vector y_end_eta_curr = solve_ode_system_rk4(TaskConfig::A_INTERVAL, TaskConfig::B_INTERVAL, h, y0_eta_curr, TaskConfig::ode_system_rhs);
            phi_curr = y_end_eta_curr.Get(1) - TaskConfig::Y_PRIME_B;
            if (j > 100) throw std::runtime_error("Превышено число итераций в методе стрельбы.");
        }
        
        if (!silent) std::cout << "Найденное начальное значение y(0) = " << eta_curr << "\n";

        int n_steps = static_cast<int>(round((TaskConfig::B_INTERVAL - TaskConfig::A_INTERVAL) / h));
        std::vector<double> x_grid(n_steps + 1);
        Vector y_solution(n_steps + 1);
        Vector y_current(2);
        y_current.Set(0, eta_curr);
        y_current.Set(1, TaskConfig::Y_PRIME_A);
        double x_current = TaskConfig::A_INTERVAL;
        x_grid[0] = x_current;
        y_solution.Set(0, y_current.Get(0));
        
        for (int i = 0; i < n_steps; ++i) {
            y_current = solve_ode_system_rk4(x_current, x_current + h, h, y_current, TaskConfig::ode_system_rhs);
            x_current += h;
            x_grid[i + 1] = x_current;
            y_solution.Set(i + 1, y_current.Get(0));
        }
        return {x_grid, y_solution};
    }

    static std::pair<std::vector<double>, Vector> solve_finite_difference_method(double h, bool silent = false) {
        if (!silent) std::cout << "\n\n--- 2. Конечно-разностный метод ---\n";

        int N = static_cast<int>(round((TaskConfig::B_INTERVAL - TaskConfig::A_INTERVAL) / h));
        int n_system = N + 1;
        
        auto p = [](double x) { return -2.0 / (exp(x) + 1.0); };
        auto q = [](double x) { return -exp(x) / (exp(x) + 1.0); };
        
        Matrix A(n_system);
        Vector B(n_system);

        double x0 = TaskConfig::A_INTERVAL;
        A.Set(0, 0, -2.0 + h * h * q(x0));
        A.Set(0, 1, 2.0);
        B.Set(0, -p(x0) * h * TaskConfig::Y_PRIME_A * 2.0);

        for (int k = 1; k < N; ++k) {
            double xk = TaskConfig::A_INTERVAL + k * h;
            A.Set(k, k - 1, 1.0 - p(xk) * h / 2.0);
            A.Set(k, k, -2.0 + h * h * q(xk));
            A.Set(k, k + 1, 1.0 + p(xk) * h / 2.0);
            B.Set(k, 0.0);
        }

        double xN = TaskConfig::B_INTERVAL;
        A.Set(N, N - 1, 2.0);
        A.Set(N, N, -2.0 + h * h * q(xN));
        B.Set(N, p(xN) * h * TaskConfig::Y_PRIME_B * 2.0);

        Vector y_solution = solve_tridiagonal_system(A, B);
        
        std::vector<double> x_grid(n_system);
        for(int i=0; i < n_system; ++i) x_grid[i] = TaskConfig::A_INTERVAL + i * h;
        
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
            double y_ex = TaskConfig::exact_solution(x);
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
        for (size_t i = 0; i < y_2h.GetSize(); ++i) {
            double y_h_val = y_h.Get(i * 2);
            double y_2h_val = y_2h.Get(i);
            double error_estimate = std::abs(y_h_val - y_2h_val) / (pow(2.0, p) - 1.0);
            if (error_estimate > max_rr_error) {
                max_rr_error = error_estimate;
            }
        }
        std::cout << "Оценка погрешности по методу Рунге-Ромберга (p=" << p << "): " << max_rr_error << "\n";
    }

public:
    static void Do() {
        double h = 0.1;

        std::cout << "========== РЕШЕНИЕ КРАЕВОЙ ЗАДАЧИ ==========\n";
        std::cout << "Уравнение: (e^x + 1)y'' - 2y' - e^x*y = 0\n";
        std::cout << "Интервал: [" << TaskConfig::A_INTERVAL << ", " << TaskConfig::B_INTERVAL << "]\n";
        std::cout << "Граничные условия: y'(0)=" << TaskConfig::Y_PRIME_A << ", y'(1)=" << TaskConfig::Y_PRIME_B << "\n";
        std::cout << "Шаг h = " << h << "\n";

        auto [x_grid_shoot, y_shoot] = solve_shooting_method(h);
        print_error_report("Метод стрельбы", x_grid_shoot, y_shoot);
        runge_romberg_error("Метод стрельбы", h, 4, solve_shooting_method);

        auto [x_grid_fd, y_fd] = solve_finite_difference_method(h);
        print_error_report("Конечно-разностный метод", x_grid_fd, y_fd);
        runge_romberg_error("Конечно-разностный метод", h, 2, solve_finite_difference_method);
    }
};


// --- Точка входа в программу ---
int main() {
    try {
        Task_4_2::Do();
    } catch (const std::exception& e) {
        std::cerr << "\nПРОИЗОШЛА КРИТИЧЕСКАЯ ОШИБКА: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}