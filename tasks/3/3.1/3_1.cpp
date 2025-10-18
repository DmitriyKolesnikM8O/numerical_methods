#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Task3_1 {
public:
    Task3_1(std::function<double(double)> func, const std::vector<double>& nodes, double point)
        : f(func), X(nodes), x_star(point) {
        if (X.empty()) {
            throw std::invalid_argument("Вектор узлов не может быть пустым.");
        }
        for (double x_val : X) {
            Y.push_back(f(x_val));
        }
    }

    void Do() {
        std::cout << std::fixed << std::setprecision(8);
        
        PrintInitialData();

        double lagrange_result = SolveLagrange();
        
        // Вычисляем коэффициенты и значение Ньютона ОДИН РАЗ
        std::vector<double> newton_coeffs;
        double newton_result = SolveNewton(newton_coeffs);
        
        if (std::abs(lagrange_result - newton_result) > 1e-9) {
            std::cerr << "Внимание: Результаты методов Лагранжа и Ньютона не совпадают!" << std::endl;
        }

        // Передаем уже вычисленные коэффициенты в другие функции
        PrintNewtonPolynomial(newton_coeffs);
        GeneratePlotData("./tasks/3/3.1/plot_data.csv", newton_coeffs);

        CalculateError(lagrange_result);
    }

private:
    std::function<double(double)> f;
    std::vector<double> X;
    std::vector<double> Y;
    double x_star;

    double EvaluateNewton(double x, const std::vector<double>& coeffs) {
        double result = coeffs[0];
        double term = 1.0;
        for (size_t i = 1; i < coeffs.size(); ++i) {
            term *= (x - X[i-1]);
            result += coeffs[i] * term;
        }
        return result;
    }

    void PrintInitialData() {
        std::cout << "--- Исходные данные ---\n";
        std::cout << "Функция: y = cos(x) + x\n";
        std::cout << "Точка для вычисления погрешности X* = " << x_star << "\n";
        std::cout << "Узлы интерполяции:\n";
        for (size_t i = 0; i < X.size(); ++i) {
            std::cout << "  X[" << i << "] = " << std::setw(10) << X[i] 
                      << "  |  Y[" << i << "] = " << std::setw(10) << Y[i] << "\n";
        }
        std::cout << "------------------------\n\n";
    }

    double SolveLagrange() {
        std::cout << "--- 1. Метод Лагранжа ---\n";
        double sum = 0.0;
        for (size_t i = 0; i < X.size(); ++i) {
            double basis_polynomial = 1.0;
            for (size_t j = 0; j < X.size(); ++j) {
                if (i != j) {
                    basis_polynomial *= (x_star - X[j]) / (X[i] - X[j]);
                }
            }
            sum += Y[i] * basis_polynomial;
        }
        std::cout << "Значение многочлена Лагранжа L(X*) в точке X* = " << x_star << " равно: " << sum << "\n";
        std::cout << "---------------------------\n\n";
        return sum;
    }
    
    double SolveNewton(std::vector<double>& coeffs) {
        std::cout << "--- 2. Метод Ньютона ---\n";
        
        std::vector<std::vector<double>> diffs(X.size());
        for(size_t i = 0; i < X.size(); ++i) {
            diffs[i].resize(X.size() - i);
            diffs[i][0] = Y[i];
        }
        for (size_t j = 1; j < X.size(); ++j) {
            for (size_t i = 0; i < X.size() - j; ++i) {
                diffs[i][j] = (diffs[i+1][j-1] - diffs[i][j-1]) / (X[i+j] - X[i]);
            }
        }
        
        // ВОЗВРАЩАЕМ ВЫВОД ТАБЛИЦЫ РАЗДЕЛЕННЫХ РАЗНОСТЕЙ
        std::cout << "Таблица разделенных разностей:\n";
        for (size_t j = 0; j < X.size(); ++j) {
            std::cout << "  f[..]_" << j << ": ";
            for (size_t i = 0; i < X.size() - j; ++i) {
                std::cout << std::setw(12) << diffs[i][j];
            }
            std::cout << std::endl;
        }
        
        coeffs.clear();
        for (size_t i = 0; i < diffs.size(); ++i) {
            coeffs.push_back(diffs[0][i]);
        }
        
        double result = EvaluateNewton(x_star, coeffs);
        
        std::cout << "\nЗначение многочлена Ньютона P(X*) в точке X* = " << x_star << " равно: " << result << "\n";
        std::cout << "---------------------------\n\n";
        return result;
    }

    void PrintNewtonPolynomial(const std::vector<double>& coeffs) {
        std::cout << "--- Вид интерполяционного многочлена Ньютона ---\n";
        std::cout << "P(x) = " << coeffs[0];
        for (size_t i = 1; i < coeffs.size(); ++i) {
            std::cout << (coeffs[i] >= 0 ? " + " : " - ") << std::abs(coeffs[i]);
            for (size_t j = 0; j < i; ++j) {
                std::cout << " * (x" << (X[j] >= 0 ? " - " : " + ") << std::abs(X[j]) << ")";
            }
        }
        std::cout << "\n-------------------------------------------------\n\n";
    }

    // Теперь функция ПРИНИМАЕТ коэффициенты, а не вычисляет их
    void GeneratePlotData(const std::string& filename, const std::vector<double>& coeffs) {
        std::cout << "--- 4. Генерация данных для графика ---\n";
        
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Не удалось создать файл " << filename << std::endl;
            return;
        }

        outFile << "x,actual_f(x),polynomial_p(x)\n";
        
        double start = X.front() - 0.1; // Немного расширим диапазон для наглядности
        double end = X.back() + 0.1;
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

    void CalculateError(double approx_value) {
        std::cout << "--- 5. Вычисление погрешности ---\n"; // Изменил номер шага для порядка
        double actual_value = f(x_star);
        double error = std::abs(actual_value - approx_value);
        std::cout << "Точное значение функции f(X*):      " << actual_value << "\n";
        std::cout << "Приближенное значение P(X*):      " << approx_value << "\n";
        std::cout << "Абсолютная погрешность |f(X*) - P(X*)|: " << error << "\n";
        std::cout << "------------------------------------\n";
    }
};

int main() {
    try {
        auto func = [](double x) { return cos(x) + x; };
        std::vector<double> nodes = {0, M_PI / 6.0, M_PI / 4.0, M_PI / 2.0};
        double x_star = 1.0;

        Task3_1 task(func, nodes, x_star);
        task.Do();

    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}