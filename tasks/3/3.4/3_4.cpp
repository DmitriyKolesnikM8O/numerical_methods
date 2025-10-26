#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include "Vector.hpp"

class Task3_4 {
public:
    Task3_4(const Vector& nodes, const Vector& values, double point)
        : x(nodes), y(values), x_star(point) {
        if (x.GetSize() != y.GetSize() || x.GetSize() < 3) {
            throw std::invalid_argument("Для численного дифференцирования требуется минимум 3 узла.");
        }
    }

    void Do() {
        std::cout << std::fixed << std::setprecision(5);
        PrintInitialData();

        int idx = FindNodeIndex();
        if (idx == -1) {
            std::cout << "Точка X* = " << x_star << " не является узлом сетки.\n";
            return;
        }

        // --- 1. Первая производная (1-й порядок точности) ---
        std::cout << "\n--- 1. Первая производная (1-й порядок точности) ---\n";
        double y_left_deriv = NAN, y_right_deriv = NAN;

        if (idx > 0) {
            y_left_deriv = (y.Get(idx) - y.Get(idx-1)) / (x.Get(idx) - x.Get(idx-1));
            std::cout << "Левосторонняя производная y'(" << x_star << "): " << y_left_deriv << "\n";
        } else {
            std::cout << "Невозможно вычислить левостороннюю производную.\n";
        }

        if (idx < x.GetSize() - 1) {
            y_right_deriv = (y.Get(idx+1) - y.Get(idx)) / (x.Get(idx+1) - x.Get(idx));
            std::cout << "Правосторонняя производная y'(" << x_star << "): " << y_right_deriv << "\n";
        } else {
            std::cout << "Невозможно вычислить правостороннюю производную.\n";
        }

        // --- 2. Первая производная (2-й порядок точности) ---
        std::cout << "\n--- 2. Первая производная (2-й порядок точности, формула 3.20) ---\n";
        double y_deriv_2nd_order = NAN;
        if (idx > 0 && idx < x.GetSize() - 1) {
            double x_i = x.Get(idx-1), x_i1 = x.Get(idx), x_i2 = x.Get(idx+1);
            double y_i = y.Get(idx-1), y_i1 = y.Get(idx), y_i2 = y.Get(idx+1);

            double term1 = (y_i1 - y_i) / (x_i1 - x_i);
            double term2_num = (y_i2 - y_i1) / (x_i2 - x_i1) - term1;
            double term2_den = x_i2 - x_i;
            double term3 = (2 * x_star - x_i - x_i1);

            y_deriv_2nd_order = term1 + (term2_num / term2_den) * term3;
            std::cout << "Производная y'(" << x_star << ") по 3 точкам: " << y_deriv_2nd_order << "\n";
        } else {
            std::cout << "Невозможно вычислить (требуется 3 точки, X* не может быть крайней).\n";
        }

        // --- 3. Вторая производная ---
        std::cout << "\n--- 3. Вторая производная (формула 3.21) ---\n";
        double y_2nd_deriv = NAN;
         if (idx > 0 && idx < x.GetSize() - 1) {
            double x_i = x.Get(idx-1), x_i1 = x.Get(idx), x_i2 = x.Get(idx+1);
            double y_i = y.Get(idx-1), y_i1 = y.Get(idx), y_i2 = y.Get(idx+1);
            
            double term1 = (y_i1 - y_i) / (x_i1 - x_i);
            double term2 = (y_i2 - y_i1) / (x_i2 - x_i1);

            y_2nd_deriv = 2 * (term2 - term1) / (x_i2 - x_i);
            std::cout << "Вторая производная y''(" << x_star << ") по 3 точкам: " << y_2nd_deriv << "\n";
        } else {
            std::cout << "Невозможно вычислить (требуется 3 точки, X* не может быть крайней).\n";
        }

        
        VerifyFormulas(idx, y_left_deriv, y_right_deriv, y_deriv_2nd_order, y_2nd_deriv, 1e-7);
    }

private:
    Vector x, y;
    double x_star;

    // Проверка, является ли сетка равномерной
    bool IsGridUniform(double tolerance, double& step) const {
        if (x.GetSize() < 2) return true;
        step = x.Get(1) - x.Get(0);
        for (int i = 1; i < x.GetSize() - 1; ++i) {
            if (std::abs((x.Get(i+1) - x.Get(i)) - step) > tolerance) {
                return false;
            }
        }
        return true;
    }

    void VerifyFormulas(int idx, double left_d, double right_d, double d1, double d2, double tol) {
        std::cout << "\n--- 4. Проверка для равноотстоящих узлов ---\n";
        double step = 0;
        if (!IsGridUniform(tol, step)) {
            std::cout << "Сетка не является равномерной. Проверка упрощенных формул не выполняется.\n";
            return;
        }

        std::cout << "Сетка равномерная с шагом h = " << step << ". Выполняем проверку:\n";
        
        // Проверка 1: Первая производная
        double simplified_d1 = (left_d + right_d) / 2.0;
        bool d1_ok = std::abs(d1 - simplified_d1) < tol;
        std::cout << "  Проверка y' (полусумма): " << simplified_d1 << " vs " << d1 << " -> " << (d1_ok ? "OK" : "FAIL") << "\n";

        // Проверка 2: Вторая производная
        if (idx > 0 && idx < x.GetSize() - 1) {
            double simplified_d2 = (y.Get(idx+1) - 2*y.Get(idx) + y.Get(idx-1)) / (step * step);
            bool d2_ok = std::abs(d2 - simplified_d2) < tol;
            std::cout << "  Проверка y'' (центральная разность): " << simplified_d2 << " vs " << d2 << " -> " << (d2_ok ? "OK" : "FAIL") << "\n";
        }
    }
    
    int FindNodeIndex() const {
        for (int i = 0; i < x.GetSize(); ++i) {
            if (std::abs(x.Get(i) - x_star) < 1e-9) return i;
        }
        return -1;
    }

    void PrintInitialData() {
        std::cout << "--- Исходные данные (Вариант 13) ---\n";
        std::cout << "Точка для вычисления производных X* = " << x_star << "\n";
        std::cout << "Таблица значений:\n";
        for (int i = 0; i < x.GetSize(); ++i) {
            std::cout << "  i=" << i << " | x[" << i << "] = " << std::setw(8) << x.Get(i) 
                      << " | y[" << i << "] = " << std::setw(10) << y.Get(i) << "\n";
        }
    }
};

int main() {
    try {
        Vector nodes(5);
        nodes.Set(0, 0.2); nodes.Set(1, 0.5); nodes.Set(2, 0.8);
        nodes.Set(3, 1.1); nodes.Set(4, 1.4);

        Vector values(5);
        values.Set(0, 12.906); values.Set(1, 5.5273); values.Set(2, 3.8777);
        values.Set(3, 3.2692); values.Set(4, 3.0319);

        double x_star = 0.8;

        Task3_4 task(nodes, values, x_star);
        task.Do();

    } catch (const std::exception& e) {
        std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}