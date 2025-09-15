#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>

class Vector {
private:
    std::vector<double> data;
    int size;

public:
    Vector(int n) : size(n), data(n, 0.0) {}

    double Get(int i) const {
        if (i < 0 || i >= size) throw std::out_of_range("Vector index out of range");
        return data[i];
    }

    void Set(int i, double value) {
        if (i < 0 || i >= size) throw std::out_of_range("Vector index out of range");
        data[i] = value;
    }

    int GetSize() const { return size; }
};

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int n; // Assuming square matrix for simplicity

public:
    Matrix(int size) : n(size), data(size, std::vector<double>(size, 0.0)) {}

    Matrix(const Matrix& other) : n(other.n), data(other.data) {}

    double Get(int i, int j) const {
        if (i < 0 || i >= n || j < 0 || j >= n) throw std::out_of_range("Matrix index out of range");
        return data[i][j];
    }

    void Set(int i, int j, double value) {
        if (i < 0 || i >= n || j < 0 || j >= n) throw std::out_of_range("Matrix index out of range");
        data[i][j] = value;
    }

    int GetLength() const { return n; }

    int GetHeight() const { return n; }

    Vector GetColumn(int k) const {
        if (k < 0 || k >= n) throw std::out_of_range("Column index out of range");
        Vector column(n);
        for (int i = 0; i < n; i++) {
            column.Set(i, data[i][k]);
        }
        return column;
    }

    void Show() const {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << data[i][j] << "\t";
            }
            std::cout << "\n";
        }
    }

    static Matrix GetSingularMatrix(int size) {
        Matrix result(size);
        for (int i = 0; i < size; i++) {
            result.Set(i, i, 1.0);
        }
        return result;
    }

    Matrix operator+(const Matrix& other) const {
        if (n != other.n) throw std::invalid_argument("Matrix dimensions must match");
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
        if (n != other.n) throw std::invalid_argument("Matrix dimensions must match");
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

class Task5 {
public:
    static void Do(Matrix A, double epsilon) {
        int k = 1;
        do {
            auto pair = A.QR();
            Matrix Q = pair.first;
            Matrix R = pair.second;

            A = R * Q;
            
            std::cout << "\nIteration " << k << ":\nMatrix A:\n";
            A.Show();
            // std::cout << "Press Enter to continue...\n";
            // std::cin.get();
            k++;

            
        } while (!Finish(A, epsilon));

        std::cout << "\nMatrix A:\n";
        A.Show();
    }

private:
    static bool Finish(const Matrix& matrix, double epsilon) {
        int m = matrix.GetLength();
        double sum1, sum2;
        for (int j = 0; j < m; j++) {
            sum1 = SquareSumColumn(matrix, j, j + 1);
            sum2 = SquareSumColumn(matrix, j, j + 2);

            if (sum2 > epsilon) {
                return false;
            } else if (sum1 <= epsilon) {
                std::cout << "lambda" << j << ": " << matrix.Get(j, j) << "\n";
            } else if (sum1 > epsilon) {
                if (j == 0) return false;

                double aii = matrix.Get(j, j);
                double ajj = matrix.Get(j + 1, j + 1);
                double aij = matrix.Get(j, j + 1);
                double aji = matrix.Get(j + 1, j);

                double x = (aii + ajj) / 2.0;
                double y = std::sqrt(-(aii + ajj) * (aii + ajj) + 4 * (aii * ajj - aij * aji)) / 2.0;

                std::cout << "lambda" << j << ": " << x << " + " << y << "i\n";
                std::cout << "lambda" << (j + 1) << ": " << x << " - " << y << "i\n";
                j++;
            }
        }
        return true;
    }

    static double SquareSumColumn(const Matrix& matrix, int column_number, int first_index) {
        int n = matrix.GetHeight();
        double sum = 0.0;
        for (int i = first_index; i < n; i++) {
            sum += matrix.Get(i, column_number) * matrix.Get(i, column_number);
        }
        return std::sqrt(sum);
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

    return matrix;
}

double get_user_epsilon() { 
    double epsilon; std::cout << "\nВведите точность вычислений (например, 1e-6): "; 
    while (!(std::cin >> epsilon) || epsilon <= 0) { 
        std::cout << "Некорректный ввод. Пожалуйста, введите положительное число: "; 
        std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    } 
    return epsilon; 
}


int main() {
    Matrix A = get_user_matrix_input();

    // A.Set(0, 0, 1.0); A.Set(0, 1, 3.0); A.Set(0, 2, 1.0);
    // A.Set(1, 0, 1.0); A.Set(1, 1, 1.0); A.Set(1, 2, 4.0);
    // A.Set(2, 0, 4.0); A.Set(2, 1, 3.0); A.Set(2, 2, 1.0);

    double epsilon = get_user_epsilon();
    Task5::Do(A, epsilon);
    return 0;
}