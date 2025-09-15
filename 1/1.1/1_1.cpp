#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <cmath>

using namespace std;

// Функция для вычисления строки верхней треугольной матрицы U с выбором главного элемента
void urow(vector<vector<float>>& u, vector<vector<float>>& l, vector<vector<float>>& a, vector<int>& p, int i, int n) {
    int max_row = i;
    float max_val = abs(a[i][i]);

    // Поиск максимального элемента в столбце i ниже текущей строки
    for (int k = i + 1; k < n; k++) {
        if (abs(a[k][i]) > max_val) {
            max_val = abs(a[k][i]);
            max_row = k;
        }
    }

    // Перестановка строк, если найден больший элемент
    if (max_row != i) {
        swap(a[i], a[max_row]);
        swap(p[i], p[max_row]);
        if (i > 0) {
            for (int k = 0; k < i; k++) {
                swap(l[i][k], l[max_row][k]);
            }
        }
    }

    float s;
    for (int j = i; j < n; j++) {
        s = 0;
        for (int k = 0; k < i; k++) {
            s += u[k][j] * l[i][k];
        }
        u[i][j] = a[i][j] - s;
    }
}

// Функция для вычисления столбца нижней треугольной матрицы L
void lcol(vector<vector<float>>& u, vector<vector<float>>& l, vector<vector<float>>& a, int j, int n) {
    float s;
    int i, k;
    for (i = j + 1; i < n; i++) {
        s = 0;
        for (k = 0; k < i; k++) {
            s += u[k][j] * l[i][k];
        }
        l[i][j] = (a[i][j] - s) / u[j][j];
    }
}

// Функция для вывода матрицы
void printmat(const vector<vector<float>>& x, int n) {
    cout << fixed; // Установка фиксированного формата
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << setprecision(3) << x[i][j];
        }
        cout << endl;
    }
    cout.unsetf(ios::fixed); // Сброс fixed для последующих выводов, если нужно
}

// Функция для перестановки вектора b в соответствии с P
void permute_b(vector<float>& b, const vector<int>& p, int n) {
    vector<float> temp = b;
    for (int i = 0; i < n; i++) {
        b[i] = temp[p[i]];
    }
}

int main() {
    int n;
    cout << "Введите размер матрицы (n x n): ";
    cin >> n;

    if (n <= 0) {
        cout << "Ошибка: размер матрицы должен быть положительным." << endl;
        return 1;
    }

    // Динамическое выделение векторов
    vector<vector<float>> a(n, vector<float>(n, 0.0));
    vector<vector<float>> l(n, vector<float>(n, 0.0));
    vector<vector<float>> u(n, vector<float>(n, 0.0));
    vector<float> b(n, 0.0);
    vector<float> x(n, 0.0);
    vector<float> v(n, 0.0);
    vector<int> p(n); // Вектор перестановок

    // Инициализация диагонали L и вектора перестановок
    for (int i = 0; i < n; i++) {
        l[i][i] = 1.0;
        p[i] = i; // Изначально тождественная перестановка
    }

    cout << "Введите элементы расширенной матрицы построчно (каждый ряд содержит " << n << " элемента матрицы A, затем элемент b):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    // LU разложение с выбором главного элемента
    for (int m = 0; m < n; m++) {
        urow(u, l, a, p, m, n);
        if (abs(u[m][m]) < 1e-6) {
            cout << "Ошибка: нулевой элемент на диагонали U[" << m << "][" << m << "]. Матрица может быть вырожденной." << endl;
            return 1;
        }
        if (m < n - 1) {
            lcol(u, l, a, m, n);
        }
    }

    // Перестановка вектора b в соответствии с перестановками
    permute_b(b, p, n);

    cout << setw(14) << "Матрица U:" << endl;
    printmat(u, n);
    cout << endl;
    cout << setw(14) << "Матрица L:" << endl;
    printmat(l, n);

    // Решение системы: Ly = b
    for (int i = 0; i < n; i++) {
        float s = 0;
        for (int j = 0; j < i; j++) {
            s += l[i][j] * v[j];
        }
        v[i] = b[i] - s;
    }

    // Решение системы: Ux = y
    for (int i = n - 1; i >= 0; i--) {
        float s = 0;
        for (int j = i + 1; j < n; j++) {
            s += u[i][j] * x[j];
        }
        if (abs(u[i][i]) < 1e-6) {
            cout << "Ошибка: нулевой элемент на диагонали U[" << i << "][" << i << "]" << endl;
            return 1;
        }
        x[i] = (v[i] - s) / u[i][i];
    }

    cout << fixed; // Установка фиксированного формата для вывода решения
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << setw(1) << i + 1 << "] = " << setw(6) << setprecision(3) << x[i] << endl;
    }
    cout.unsetf(ios::fixed); // Сброс fixed

    return 0;
}