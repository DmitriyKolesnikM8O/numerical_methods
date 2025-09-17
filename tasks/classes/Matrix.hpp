#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include <vector>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <iomanip>


class Matrix {
private:
    int n;
    std::vector<std::vector<double>> data;
    
public:
    Matrix() : n(0), data(0, std::vector<double>(0, 0.0)) {}
    Matrix(int size) : n(size), data(size, std::vector<double>(size, 0.0)) {}

    Matrix(const Matrix& other) : n(other.n), data(other.data) {}

    Matrix& operator=(const Matrix& other) {
        if (this != &other) { // Проверка на самоприсваивание
            n = other.n;
            data = other.data;
        }
        return *this;
    }

    double Get(int i, int j) const {
        if (i < 0 || i >= n || j < 0 || j >= n) throw std::out_of_range("Индекс не в границах матрицы");
        return data[i][j];
    }

    void Set(int i, int j, double value) {
        if (i < 0 || i >= n || j < 0 || j >= n) throw std::out_of_range("Индекс не в границах матрицы");
        data[i][j] = value;
    }

    int GetLength() const { return n; }

    int GetHeight() const { return n; }

    Vector GetColumn(int k) const {
        if (k < 0 || k >= n) throw std::out_of_range("Столбец не в границах матрицы");
        Vector column(n);
        for (int i = 0; i < n; i++) {
            column.Set(i, data[i][k]);
        }
        return column;
    }

    void Show() const {
        std::cout << std::fixed << std::setprecision(4); 
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << std::setw(12) << data[i][j];
            }
            std::cout << "\n";
        }
        std::cout.unsetf(std::ios::fixed);
    }

    static Matrix GetSingularMatrix(int size) {
        Matrix result(size);
        for (int i = 0; i < size; i++) {
            result.Set(i, i, 1.0);
        }
        return result;
    }

    Matrix operator+(const Matrix& other) const {
        if (n != other.n) throw std::invalid_argument("Недопустимые размерности");
        Matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.Set(i, j, data[i][j] + other.data[i][j]);
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.Set(i, j, data[i][j] * scalar);
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (n != other.n) throw std::invalid_argument("Недопустимые размерности");
        Matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += data[i][k] * other.data[k][j];
                }
                result.Set(i, j, sum);
            }
        }
        return result;
    }

    std::pair<Matrix, Matrix> QR() {
        Matrix R = *this;
        Matrix Q = GetSingularMatrix(n);
        Vector v(n);

        for (int k = 0; k < n - 1; k++) {
            for (int i = 0; i < n; i++) {
                if (i < k) {
                    v.Set(i, 0.0);
                } else if (i == k) {
                    double norma = 0.0;
                    Vector column = R.GetColumn(k);
                    for (int j = k; j < n; j++) {
                        norma += column.Get(j) * column.Get(j);
                    }
                    norma = std::sqrt(norma);
                    v.Set(i, R.Get(i, i) + std::copysign(norma, R.Get(i, i)));
                } else if (i > k) {
                    v.Set(i, R.Get(i, k));
                }
            }

            Matrix temp(n);
            double p = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    temp.Set(i, j, v.Get(i) * v.Get(j));
                }
                p += v.Get(i) * v.Get(i);
            }

            p = -2.0 / p;
            Matrix H = GetSingularMatrix(n) + (temp * p);

            R = H * R;
            Q = Q * H;
        }

        return std::make_pair(Q, R);
    }
};

Matrix get_user_matrix_input() { 
    size_t rows;
    std::cout << "Введите размерность матрицы: "; 
    while (!(std::cin >> rows) || rows <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите целое число больше 0: "; 
        std::cin.clear(); 
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    } 
    Matrix matrix(rows);
    std::cout << "Введите элементы матрицы (" << rows << "x" << rows << "):\n";

    for (size_t i = 0; i < rows; ++i) { 
        for (size_t j = 0; j < rows; ++j) { 
            std::cout << "Элемент [" << i << "][" << j << "]: "; 
            double n;
            while (!(std::cin >> n)) {
                
                std::cout << "Некорректный ввод. Пожалуйста, введите число: "; 
                std::cin.clear(); 
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
            }
            matrix.Set(i, j, n); 
        } 
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    matrix.Show();
    return matrix;
}

Matrix get_user_matrix_input_rows() {
    size_t rows;
    std::cout << "Введите размерность матрицы: "; 
    while (!(std::cin >> rows) || rows <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите целое число больше 0: "; 
        std::cin.clear(); 
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    } 
    Matrix matrix(rows);
    std::cout << "Введите элементы матрицы (" << rows << "x" << rows << "):\n";
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "Введите элементы " << i << " строки матрицы (формат: 1 2 3): ";
        for (size_t j = 0; j < rows; ++j) {
            double n;
            std::cin >> n;
            matrix.Set(i, j, n);
        }
    }

    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    matrix.Show();
    return matrix;
} 

Matrix get_user_matrix_input_file(std::istream& is) {
    size_t n;
    is >> n;
    if (is.fail() || n <= 0) {
        throw std::invalid_argument("Некорректная размерность матрицы в потоке ввода.");
    }
    Matrix matrix(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double value;
            is >> value;
            if (is.fail()) {
                throw std::invalid_argument("Некорректные данные матрицы в потоке ввода.");
            }
            matrix.Set(i, j, value);
        }
    }
    matrix.Show();
    return matrix;
}

Vector get_free_members_vector() {
    size_t rows;
    std::cout << "Введите размерность вектора: ";
    while (!(std::cin >> rows) || rows <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите целое число больше 0: "; 
        std::cin.clear(); 
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    }
    Vector vector(rows);
    std::cout << "Введите элементы вектора (" << rows << "):\n";
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "Элемент [" << i << "]: ";
        double n;
        while (!(std::cin >> n)) {
                std::cout << "Некорректный ввод. Пожалуйста, введите число: "; 
                std::cin.clear(); 
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
        }
        vector.Set(i, n);
    }

    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return vector;
}

Vector get_free_members_vector_rows() {
        size_t rows;
    std::cout << "Введите размерность вектора: ";
    while (!(std::cin >> rows) || rows <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите целое число больше 0: "; 
        std::cin.clear(); 
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    }
    Vector vector(rows);
    std::cout << "Введите элементы вектора (" << rows << "): ";
    for (size_t i = 0; i < rows; ++i) {
        double n;
        std::cin >> n;
        vector.Set(i, n);

    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return vector;
}

Vector get_free_members_vector_file(std::istream& is) {
    size_t n;
    is >> n;
    if (is.fail() || n <= 0) {
        throw std::invalid_argument("Некорректная размерность вектора в потоке ввода.");
    }
    Vector vector(n);
    for (size_t i = 0; i < n; i++) {
        double value;
        is >> value;
        if (is.fail()) {
            throw std::invalid_argument("Некорректные данные вектора в потоке ввода.");
        }
        vector.Set(i, value);
    }
    return vector;
}

double get_user_epsilon() { 
    double epsilon; std::cout << "\nВведите точность вычислений (например, 1e-6): "; 
    while (!(std::cin >> epsilon) || epsilon <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите положительное число: "; 
        std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    } 
    return epsilon; 
}

double get_user_epsilon_file(std::istream& is) {
    double epsilon;
    is >> epsilon;
    if (is.fail() || epsilon <= 0) {
        throw std::invalid_argument("Некорректное значение точности в потоке ввода.");
    }
    return epsilon;
}


#endif //MATRIX_H