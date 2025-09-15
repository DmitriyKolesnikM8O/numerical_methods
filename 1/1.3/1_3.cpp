#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Решение СЛАУ методом итераций
double* iter(double** a, double* y, int n) {
    double* res = new double[n];
    int i, j;

    // Инициализация начального приближения
    for (i = 0; i < n; i++) {
        res[i] = y[i] / a[i][i];
    }

    double eps = 0.0001;
    double* Xn = new double[n];
    bool flag;

    do {
        for (i = 0; i < n; i++) {
            Xn[i] = y[i] / a[i][i];
            for (j = 0; j < n; j++) {
                if (i == j)
                    continue;
                else {
                    Xn[i] -= a[i][j] / a[i][i] * res[j];
                }
            }
        }

        flag = true;
        for (i = 0; i < n; i++) {
            if (abs(Xn[i] - res[i]) > eps) {
                flag = false;
                break;
            }
        }

        for (i = 0; i < n; i++) {
            res[i] = Xn[i];
        }

    } while (!flag); // Продолжаем, пока не достигнута точность

    delete[] Xn; // Освобождение памяти
    return res;
}

int main() {
    setlocale(LC_ALL, "Russian");
    double** a;
    double* y;
    double* x;
    int n, i, j;

    cout << "Введите размер матрицы (n x n): ";
    cin >> n;

    y = new double[n];
    a = new double* [n];
    for (i = 0; i < n; i++) {
        a[i] = new double[n];
    }

    cout << "Введите элементы расширенной матрицы построчно (каждый ряд содержит " << n << " элемента матрицы A, затем элемент b):" << endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) {
            if (j != n) {
                cin >> a[i][j];
            }
            else {
                cin >> y[i];
            }
        }
    }

    cout << "Matrix A:" << endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cout << a[i][j] << "\t ";
        }
        cout << endl;
    }

    cout << "Matrix B:" << endl;
    for (i = 0; i < n; i++) {
        cout << y[i] << endl;
    }

    x = iter(a, y, n);

    cout << "Решение системы:" << endl;
    for (i = 0; i < n; i++) {
        cout << "x[" << setw(1) << i + 1 << "] = " << setw(6) << setprecision(3) << fixed << x[i] << endl;
    }

    // Освобождение памяти
    for (i = 0; i < n; i++) {
        delete[] a[i];
    }
    delete[] a;
    delete[] y;
    delete[] x;

    return 0;
}