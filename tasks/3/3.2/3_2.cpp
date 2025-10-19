#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include "Vector.hpp"

class Task3_2 {
public:
    Task3_2(const Vector& nodes, const Vector& values, double point)
        : x(nodes), f(values), x_star(point) {
        if (x.GetSize() != f.GetSize() || x.GetSize() < 2) {
            throw std::invalid_argument("Некорректные входные данные для сплайна.");
        }
        
        int n = x.GetSize() - 1;
        a = Vector(n); b = Vector(n); c = Vector(n); d = Vector(n);
    }

    void Do() {
        std::cout << std::fixed << std::setprecision(8);
        PrintInitialData();
        int n = x.GetSize() - 1;
        Vector h(n);
        for (int i = 0; i < n; ++i) h.Set(i, x.Get(i+1) - x.Get(i));

        Vector c_full = SolveTridiagonalSystem(n, h);
        for(int i = 0; i < n; ++i) c.Set(i, c_full.Get(i));

        for (int i = 0; i < n; ++i) {
            a.Set(i, f.Get(i));
            if (i < n - 1) {
                d.Set(i, (c.Get(i+1) - c.Get(i)) / (3.0 * h.Get(i)));
                b.Set(i, (f.Get(i+1) - f.Get(i)) / h.Get(i) - (h.Get(i) / 3.0) * (c.Get(i+1) + 2.0 * c.Get(i)));
            } else {
                d.Set(i, -c.Get(i) / (3.0 * h.Get(i)));
                b.Set(i, (f.Get(i+1) - f.Get(i)) / h.Get(i) - (2.0 / 3.0) * h.Get(i) * c.Get(i));
            }
        }
        
        PrintSplineCoefficients();

        // <<<--- НОВЫЙ БЛОК: ВЫЗОВ ПРОВЕРКИ ---<<<
        VerifySpline(1e-7);

        EvaluateAtPoint();
    }

    // Проверка корректности сплайна
    void VerifySpline(double tolerance) {
        std::cout << "--- 4. Проверка корректности построенного сплайна ---\n";
        bool all_ok = true;

        // 1. Проверка прохождения через узлы
        std::cout << "1. Проверка прохождения через узлы S(x_i) = f_i:\n";
        for (int i = 0; i < x.GetSize(); ++i) {
            double val = Evaluate(x.Get(i));
            bool ok = std::abs(val - f.Get(i)) < tolerance;
            std::cout << "   Узел x[" << i << "]=" << x.Get(i) << ": S(x)=" << val << ", f=" << f.Get(i) << " -> " << (ok ? "OK" : "FAIL") << "\n";
            if (!ok) all_ok = false;
        }

        // 2. Проверка непрерывности производных во внутренних узлах
        std::cout << "\n2. Проверка непрерывности производных S'(x_i) и S''(x_i):\n";
        for (int i = 1; i < x.GetSize() - 1; ++i) {
            double node = x.Get(i);
            double d1_left = EvaluateDerivative(node, 1, i - 1);
            double d1_right = EvaluateDerivative(node, 1, i);
            double d2_left = EvaluateDerivative(node, 2, i - 1);
            double d2_right = EvaluateDerivative(node, 2, i);

            bool d1_ok = std::abs(d1_left - d1_right) < tolerance;
            bool d2_ok = std::abs(d2_left - d2_right) < tolerance;
            
            std::cout << "   Узел x[" << i << "]=" << node << ":\n";
            std::cout << "     S' слева=" << d1_left << ", S' справа=" << d1_right << " -> " << (d1_ok ? "OK" : "FAIL") << "\n";
            std::cout << "     S'' слева=" << d2_left << ", S'' справа=" << d2_right << " -> " << (d2_ok ? "OK" : "FAIL") << "\n";
            if (!d1_ok || !d2_ok) all_ok = false;
        }

        // 3. Проверка граничных условий
        std::cout << "\n3. Проверка граничных условий S''(x_0)=0, S''(x_n)=0:\n";
        double d2_start = EvaluateDerivative(x.Get(0), 2);
        double d2_end = EvaluateDerivative(x.Get(x.GetSize() - 1), 2);
        bool start_ok = std::abs(d2_start) < tolerance;
        bool end_ok = std::abs(d2_end) < tolerance;
        std::cout << "   S''(" << x.Get(0) << ") = " << d2_start << " -> " << (start_ok ? "OK" : "FAIL") << "\n";
        std::cout << "   S''(" << x.Get(x.GetSize()-1) << ") = " << d2_end << " -> " << (end_ok ? "OK" : "FAIL") << "\n";
        if (!start_ok || !end_ok) all_ok = false;
        
        std::cout << "\nИтог проверки: " << (all_ok ? "Все тесты пройдены успешно!" : "ОБНАРУЖЕНЫ НЕСООТВЕТСТВИЯ!") << "\n";
        std::cout << "---------------------------------------------------------\n\n";
    }


private:
    Vector x, f, a, b, c, d;
    double x_star;

    // Вспомогательные функции для вычислений и проверок
    int FindInterval(double point) const {
        for (int i = 0; i < x.GetSize() - 1; ++i) {
            if (point >= x.Get(i) && point <= x.Get(i+1)) {
                 // Особый случай для самой левой точки
                if (point == x.Get(0)) return 0;
                // Особый случай для самой правой точки
                if (point == x.Get(x.GetSize() - 1)) return x.GetSize() - 2;
                return i;
            }
        }
        return -1; // Точка вне диапазона
    }

    double Evaluate(double point) const {
        int i = FindInterval(point);
        if (i == -1) return NAN;
        double dx = point - x.Get(i);
        return a.Get(i) + b.Get(i) * dx + c.Get(i) * dx * dx + d.Get(i) * dx * dx * dx;
    }

    double EvaluateDerivative(double point, int order, int interval_idx = -1) const {
        int i = (interval_idx == -1) ? FindInterval(point) : interval_idx;
        if (i == -1) return NAN;
        double dx = point - x.Get(i);
        if (order == 1) { // Первая производная
            return b.Get(i) + 2.0 * c.Get(i) * dx + 3.0 * d.Get(i) * dx * dx;
        }
        if (order == 2) { // Вторая производная
            return 2.0 * c.Get(i) + 6.0 * d.Get(i) * dx;
        }
        return Evaluate(point); // По умолчанию (order=0) возвращаем значение
    }

    Vector SolveTridiagonalSystem(int n, const Vector& h) {
        int size = n - 1;
        if (size <= 0) return Vector(n);
        Vector A(size), B(size), C(size), F(size);
        for (int i = 0; i < size; ++i) {
            int j = i + 2;
            B.Set(i, 2.0 * (h.Get(j-2) + h.Get(j-1)));
            if (i > 0) A.Set(i, h.Get(j-2));
            if (i < size - 1) C.Set(i, h.Get(j-1));
            F.Set(i, 3.0 * (((f.Get(j) - f.Get(j-1)) / h.Get(j-1)) - ((f.Get(j-1) - f.Get(j-2)) / h.Get(j-2))));
        }
        Vector alpha(size), beta(size);
        alpha.Set(0, -C.Get(0) / B.Get(0)); beta.Set(0, F.Get(0) / B.Get(0));
        for (int i = 1; i < size; ++i) {
            double denom = B.Get(i) + A.Get(i) * alpha.Get(i-1);
            alpha.Set(i, -C.Get(i) / denom);
            beta.Set(i, (F.Get(i) - A.Get(i) * beta.Get(i-1)) / denom);
        }
        Vector result_c(size);
        result_c.Set(size-1, beta.Get(size-1));
        for (int i = size - 2; i >= 0; --i) {
            result_c.Set(i, alpha.Get(i) * result_c.Get(i+1) + beta.Get(i));
        }
        Vector full_c(n);
        for(int i = 0; i < size; ++i) full_c.Set(i+1, result_c.Get(i));
        std::cout << "--- Шаг 1: Найденные коэффициенты 'c' ---\n";
        for (int i = 0; i < n; ++i) std::cout << "c[" << i + 1 << "] = " << full_c.Get(i) << "\n";
        std::cout << "-------------------------------------------\n\n";
        return full_c;
    }
    void PrintInitialData() {
        std::cout << "--- Исходные данные для построения сплайна ---\n";
        std::cout << "Точка для вычисления X* = " << x_star << "\n";
        std::cout << "Узлы интерполяции:\n";
        for (int i = 0; i < x.GetSize(); ++i) std::cout << "  i=" << i << " | X[" << i << "] = " << std::setw(8) << x.Get(i) << " | f[" << i << "] = " << std::setw(10) << f.Get(i) << "\n";
        std::cout << "-----------------------------------------------\n\n";
    }
    void PrintSplineCoefficients() {
        std::cout << "--- Шаг 2: Коэффициенты кубических сплайнов ---\n";
        for (int i = 0; i < a.GetSize(); ++i) {
            std::cout << "Интервал [" << x.Get(i) << ", " << x.Get(i+1) << "] (i=" << i + 1 << "):\n";
            std::cout << "  a" << i+1 << " = " << a.Get(i) << "\n" << "  b" << i+1 << " = " << b.Get(i) << "\n" << "  c" << i+1 << " = " << c.Get(i) << "\n" << "  d" << i+1 << " = " << d.Get(i) << "\n";
        }
        std::cout << "--------------------------------------------------\n\n";
    }
    void EvaluateAtPoint() {
        std::cout << "--- Шаг 3 (ФИНАЛ): Вычисление значения в точке X* ---\n";
        int i = FindInterval(x_star);
        if (i == -1) { std::cout << "Точка X*=" << x_star << " находится вне диапазона.\n"; return; }
        double dx = x_star - x.Get(i);
        double result = a.Get(i) + b.Get(i) * dx + c.Get(i) * dx * dx + d.Get(i) * dx * dx * dx;
        std::cout << "Точка X* = " << x_star << " попадает в интервал [" << x.Get(i) << ", " << x.Get(i+1) << "].\n";
        std::cout << "Используем сплайн S" << i + 1 << "(x).\n";
        std::cout << "S(X*) = a + b(X*-x_i) + c(X*-x_i)^2 + d(X*-x_i)^3\n";
        std::cout << "S(" << x_star << ") = " << result << "\n";
        std::cout << "---------------------------------------------------------\n";
    }
};

int main() {
    try {
        Vector nodes(5);
        nodes.Set(0, 0.0); nodes.Set(1, 1.0); nodes.Set(2, 2.0); nodes.Set(3, 3.0); nodes.Set(4, 4.0);
        Vector values(5);
        values.Set(0, 1.0); values.Set(1, 1.5403); values.Set(2, 1.5839); values.Set(3, 2.01); values.Set(4, 3.3464);
        double x_star = 1.5;

        Task3_2 task(nodes, values, x_star);
        task.Do();

    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}