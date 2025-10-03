#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max

using namespace std;

int main() {
    double eps;
    cout << "Введите точность вычислений (eps): ";
    cin >> eps;

    double a = 4.0;
    double ah = a / 2.0;
    double a2 = a * a;
    double a3 = a2 * a;

    // Initial approximation determined graphically
    double x1_init = 4.0;
    double x2_init = 2.0;

    cout << "Метод простой итерации" << endl;

    vector<double> fp_errors;
    fp_errors.push_back(1.0); // Initial dummy error

    double x1_fp = x1_init;
    double x2_fp = x2_init;
    int iterations_fp = 0;
    int max_iter = 10000;

    do {
        double x1_old = x1_fp;
        double x2_old = x2_fp;

        // Check for valid sqrt argument
        double sq_arg = a2 - pow(x2_old - ah, 2);
        if (sq_arg < 0) {
            cout << "Ошибка: отрицательный аргумент под корнем в методе простой итерации." << endl;
            return 1;
        }

        x1_fp = ah + sqrt(sq_arg);
        x2_fp = a3 / (x1_old * x1_old + a2);

        double err = max(fabs(x1_fp - x1_old), fabs(x2_fp - x2_old));
        fp_errors.push_back(err);

        iterations_fp++;

        if (iterations_fp > max_iter) {
            cout << "Итерации не сошлись (метод простой итерации)." << endl;
            return 1;
        }
    } while (fp_errors.back() >= eps);

    cout << "Положительное решение (метод простой итерации): x1 = " << x1_fp << ", x2 = " << x2_fp << endl;
    cout << "Количество итераций (метод простой итерации): " << iterations_fp << endl;

    cout << "Зависимость погрешности от количества итераций:" << endl;
    for (int i = 1; i <= iterations_fp; ++i) {
        cout << "Итерация " << i << ": погрешность = " << fp_errors[i] << endl;
    }

    // Newton method
    cout << "\nМетод Ньютона" << endl;

    vector<double> newton_errors;
    newton_errors.push_back(1.0);

    double x1_newton = x1_init;
    double x2_newton = x2_init;
    int iterations_newton = 0;

    do {
        double f1 = (x1_newton * x1_newton + a2) * x2_newton - a3;
        double f2 = pow(x1_newton - ah, 2) + pow(x2_newton - ah, 2) - a2;

        double j11 = 2 * x1_newton * x2_newton;
        double j12 = x1_newton * x1_newton + a2;
        double j21 = 2 * (x1_newton - ah);
        double j22 = 2 * (x2_newton - ah);

        double det = j11 * j22 - j12 * j21;
        if (fabs(det) < 1e-10) {
            cout << "Ошибка: вырожденный якобиан в методе Ньютона." << endl;
            return 1;
        }

        double delta1 = (-f1 * j22 + f2 * j12) / det;
        double delta2 = (f1 * j21 - f2 * j11) / det;

        x1_newton += delta1;
        x2_newton += delta2;

        double err = max(fabs(delta1), fabs(delta2));
        newton_errors.push_back(err);

        iterations_newton++;

        if (iterations_newton > max_iter) {
            cout << "Итерации не сошлись (метод Ньютона)." << endl;
            return 1;
        }
    } while (newton_errors.back() >= eps);

    cout << "Положительное решение (метод Ньютона): x1 = " << x1_newton << ", x2 = " << x2_newton << endl;
    cout << "Количество итераций (метод Ньютона): " << iterations_newton << endl;

    cout << "Зависимость погрешности от количества итераций:" << endl;
    for (int i = 1; i <= iterations_newton; ++i) {
        cout << "Итерация " << i << ": погрешность = " << newton_errors[i] << endl;
    }

    return 0;
}
