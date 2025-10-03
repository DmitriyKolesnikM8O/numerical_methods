#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double f(double x) {
    return log(sqrt(1 - x * x) + 0.1) - x;
}

double g(double x) {
    return log(sqrt(1 - x * x) + 0.1);
}

double df(double x) {
    double denom = sqrt(1 - x * x) + 0.1;
    return (-x / (denom * sqrt(1 - x * x))) - 1;
}

int main() {
    double eps;
    cout << "Введите точность вычислений (eps): ";
    cin >> eps;
    double a = 0, b = 1;

    // Fixed-point iteration
    double x0_fp = (a + b) / 2.0;
    double x_fp = g(x0_fp);
    int iterations_fp = 1;
    vector<double> fp_errors;

    fp_errors.push_back(fabs(x_fp - x0_fp)); // Initial error
    while (fabs(x_fp - x0_fp) >= eps) {
        x0_fp = x_fp;
        x_fp = g(x0_fp);
        fp_errors.push_back(fabs(x_fp - x0_fp));
        iterations_fp++;
        if (iterations_fp > 10000) {
            cout << "Итерации не сошлись (метод простой итерации)." << endl;
            return 1;
        }
    }

    cout << "Положительный корень (метод простой итерации) примерно равен: " << x_fp << endl;
    cout << "Количество итераций (метод простой итерации): " << iterations_fp << endl;

    cout << "\nЗависимость погрешности от количества итераций (метод простой итерации):" << endl;
    for (int i = 0; i < iterations_fp; ++i) {
        cout << "Итерация " << (i + 1) << ": погрешность = " << fp_errors[i] << endl;
    }

    // Newton method
    double x0_newton = (a + b) / 2.0;
    double x_newton = x0_newton - f(x0_newton) / df(x0_newton);
    int iterations_newton = 1;
    vector<double> newton_errors;

    newton_errors.push_back(fabs(x_newton - x0_newton)); // Initial error
    while (fabs(x_newton - x0_newton) >= eps) {
        x0_newton = x_newton;
        x_newton = x0_newton - f(x0_newton) / df(x0_newton);
        newton_errors.push_back(fabs(x_newton - x0_newton));
        iterations_newton++;
        if (iterations_newton > 10000) {
            cout << "Итерации не сошлись (метод Ньютона)." << endl;
            return 1;
        }
    }

    cout << "\n----------------------------------------------------------------------\n" << endl;
    cout << "Положительный корень (метод Ньютона) примерно равен: " << x_newton << endl;
    cout << "Количество итераций (метод Ньютона): " << iterations_newton << endl;

    cout << "\nЗависимость погрешности от количества итераций (метод Ньютона):" << endl;
    for (int i = 0; i < iterations_newton; ++i) {
        cout << "Итерация " << (i + 1) << ": погрешность = " << newton_errors[i] << endl;
    }

    return 0;
}
