#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    int i, n, k, n1;
    double z;
    double eps[50];
    double et[50];
    cout << "Введите размер матрицы (n x n): ";
    cin >> n1;
    double A[n1][n1];
    double B[n1];
    double X[n1];
    cout << "Введите элементы расширенной матрицы построчно (каждый ряд содержит " << n1 << " элемента матрицы A, затем элемент b):" << endl;
    for (i = 0; i < n1; i++) {
        for (k = 0; k < n1; k++) {
            cin >> A[i][k];
        }
        cin >> B[i];
    }

    // (Опционально) Проверка на трехдиагональность
    for (i = 0; i < n1; i++) {
        for (k = 0; k < n1; k++) {
            if ((k > i + 1 || k < i - 1) && i != k && abs(A[i][k]) > 1e-10) {
                cout << "Ошибка: матрица не является трехдиагональной!" << endl;
                return 1;
            }
        }
    }

    cout << "Матрица A:" << endl;
    for (i = 0; i < n1; i++) {
        for (k = 0; k < n1; k++) {
            cout << A[i][k] << "\t ";
        }
        cout << endl;
    }

    cout << "Матрица B:" << endl;
    for (i = 0; i < n1; i++) {
        cout << B[i] << endl;
    }

    n = n1 - 1;
    eps[0] = -A[0][1] / A[0][0];
    et[0] = B[0] / A[0][0];

    for (i = 1; i < n; i++) {
        z = A[i][i] + A[i][i - 1] * eps[i - 1];
        eps[i] = -A[i][i + 1] / z;
        et[i] = (B[i] - A[i][i - 1] * et[i - 1]) / z;
    }

    X[n] = (B[n] - A[n][n - 1] * et[n - 1]) / (A[n][n] + A[n][n - 1] * eps[n - 1]);

    for (i = n - 1; i >= 0; i--) {
        X[i] = eps[i] * X[i + 1] + et[i];
    }

    cout << "Решение системы:" << endl;
    for (i = 0; i < n1; i++) {
        cout << "x[" << setw(1) << i + 1 << "] = " << setw(6) << setprecision(3) << fixed << X[i] << endl;
    }

    return 0;
}