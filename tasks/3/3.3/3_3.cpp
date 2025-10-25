#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include "Vector.hpp"

class Task3_3 {
public:
    Task3_3(const Vector& nodes, const Vector& values)
        : x(nodes), y(values) {
        if (x.GetSize() != y.GetSize() || x.GetSize() == 0) {
            throw std::invalid_argument("Некорректные входные данные для МНК.");
        }
    }

    void Do() {
        std::cout << std::fixed << std::setprecision(5);
        PrintInitialData();
        SolveDegree1();
        SolveDegree2();
    }

private:
    Vector x, y;

    // Решение СЛАУ методом Гаусса
    Vector SolveLinearSystem(Vector& A_flat, Vector& B) {
        int n = B.GetSize();
        // Прямой ход
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A_flat.Get(j * n + i)) > std::abs(A_flat.Get(pivot * n + i))) {
                    pivot = j;
                }
            }
            for (int k = 0; k < n; ++k) {
                double temp = A_flat.Get(i * n + k);
                A_flat.Set(i * n + k, A_flat.Get(pivot * n + k));
                A_flat.Set(pivot * n + k, temp);
            }
            double temp_b = B.Get(i); B.Set(i, B.Get(pivot)); B.Set(pivot, temp_b);

            if (std::abs(A_flat.Get(i * n + i)) < 1e-10) throw std::runtime_error("Система вырождена.");

            for (int j = i + 1; j < n; ++j) {
                double factor = A_flat.Get(j * n + i) / A_flat.Get(i * n + i);
                B.Set(j, B.Get(j) - factor * B.Get(i));
                for (int k = i; k < n; ++k) {
                    double val = A_flat.Get(j * n + k) - factor * A_flat.Get(i * n + k);
                    A_flat.Set(j * n + k, val);
                }
            }
        }

        // Обратный ход
        Vector solution(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j) {
                sum += A_flat.Get(i * n + j) * solution.Get(j);
            }
            solution.Set(i, (B.Get(i) - sum) / A_flat.Get(i * n + i));
        }
        return solution;
    }
    
    // --- Часть a) Многочлен 1-й степени ---
    void SolveDegree1() {
        std::cout << "\n--- а) Построение многочлена 1-й степени F1(x) = a0 + a1*x ---\n";
        int n_points = x.GetSize();
        
        double sx = 0, sy = 0, sx2 = 0, sxy = 0;
        for (int i = 0; i < n_points; ++i) {
            double xi = x.Get(i), yi = y.Get(i);
            sx += xi; sy += yi; sx2 += xi * xi; sxy += xi * yi;
        }

        Vector A_flat(4);
        A_flat.Set(0, n_points); A_flat.Set(1, sx);
        A_flat.Set(2, sx);       A_flat.Set(3, sx2);

        Vector B(2); B.Set(0, sy); B.Set(1, sxy);
        
        std::cout << "Нормальная система МНК:\n";
        std::cout << "  " << A_flat.Get(0) << "*a0 + " << A_flat.Get(1) << "*a1 = " << B.Get(0) << "\n";
        std::cout << "  " << A_flat.Get(2) << "*a0 + " << A_flat.Get(3) << "*a1 = " << B.Get(1) << "\n\n";

        Vector coeffs = SolveLinearSystem(A_flat, B);
        double a0 = coeffs.Get(0);
        double a1 = coeffs.Get(1);
        
        std::cout << "Решение системы: a0 = " << a0 << ", a1 = " << a1 << "\n";
        std::cout << "Приближающий многочлен: F1(x) = " << a0 << (a1 >= 0 ? " + " : " - ") << std::abs(a1) << "*x\n\n";
        
        double error_sum = 0;
        for (int i = 0; i < n_points; ++i) {
            double fi = a0 + a1 * x.Get(i);
            error_sum += (fi - y.Get(i)) * (fi - y.Get(i));
        }
        std::cout << "Сумма квадратов ошибок Ф = " << error_sum << "\n";
    }

    // --- Часть б) Многочлен 2-й степени ---
    void SolveDegree2() {
        std::cout << "\n--- б) Построение многочлена 2-й степени F2(x) = a0 + a1*x + a2*x^2 ---\n";
        int n_points = x.GetSize();

        double sx=0, sy=0, sx2=0, sxy=0, sx3=0, sx4=0, sx2y=0;
        for (int i = 0; i < n_points; ++i) {
            double xi = x.Get(i), yi = y.Get(i), xi2 = xi * xi;
            sx += xi; sy += yi; sx2 += xi2;
            sxy += xi * yi; sx3 += xi * xi2; sx4 += xi2 * xi2;
            sx2y += xi2 * yi;
        }

        // Формируем матрицу 3x3 в виде 1D вектора размером 9
        Vector A_flat(9);
        A_flat.Set(0, n_points); A_flat.Set(1, sx);  A_flat.Set(2, sx2);
        A_flat.Set(3, sx);       A_flat.Set(4, sx2); A_flat.Set(5, sx3);
        A_flat.Set(6, sx2);      A_flat.Set(7, sx3); A_flat.Set(8, sx4);

        Vector B(3); B.Set(0, sy); B.Set(1, sxy); B.Set(2, sx2y);
        
        std::cout << "Нормальная система МНК:\n";
        std::cout << "  " << A_flat.Get(0) << "*a0 + " << A_flat.Get(1) << "*a1 + " << A_flat.Get(2) << "*a2 = " << B.Get(0) << "\n";
        std::cout << "  " << A_flat.Get(3) << "*a0 + " << A_flat.Get(4) << "*a1 + " << A_flat.Get(5) << "*a2 = " << B.Get(1) << "\n";
        std::cout << "  " << A_flat.Get(6) << "*a0 + " << A_flat.Get(7) << "*a1 + " << A_flat.Get(8) << "*a2 = " << B.Get(2) << "\n\n";

        Vector coeffs = SolveLinearSystem(A_flat, B);
        double a0 = coeffs.Get(0);
        double a1 = coeffs.Get(1);
        double a2 = coeffs.Get(2);

        std::cout << "Решение системы: a0 = " << a0 << ", a1 = " << a1 << ", a2 = " << a2 << "\n";
        std::cout << "Приближающий многочлен: F2(x) = " << a0 << (a1 >= 0 ? " + " : " - ") << std::abs(a1) 
                  << "*x" << (a2 >= 0 ? " + " : " - ") << std::abs(a2) << "*x^2\n\n";

        double error_sum = 0;
        for (int i = 0; i < n_points; ++i) {
            double xi = x.Get(i);
            double fi = a0 + a1 * xi + a2 * xi * xi;
            error_sum += (fi - y.Get(i)) * (fi - y.Get(i));
        }
        std::cout << "Сумма квадратов ошибок Ф = " << error_sum << "\n";
    }

    void PrintInitialData() {
        std::cout << "--- Исходные данные (Вариант 13) ---\n";
        for (int i = 0; i < x.GetSize(); ++i) {
            std::cout << "  i=" << i << " | x[" << i << "] = " << std::setw(8) << x.Get(i) 
                      << " | y[" << i << "] = " << std::setw(10) << y.Get(i) << "\n";
        }
    }
};

int main() {
    try {
        
        Vector nodes(6);
        nodes.Set(0, -1.0); nodes.Set(1, 0.0); nodes.Set(2, 1.0);
        nodes.Set(3, 2.0); nodes.Set(4, 3.0); nodes.Set(5, 4.0);

        Vector values(6);
        values.Set(0, -0.4597); values.Set(1, 1.0); values.Set(2, 1.5403);
        values.Set(3, 1.5839); values.Set(4, 2.010); values.Set(5, 3.3464);

        Task3_3 task(nodes, values);
        task.Do();

    } catch (const std::exception& e) {
        std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}