#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include <vector>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <iomanip>

//всегда квадратная
class Matrix {
private:
    int n;
    std::vector<std::vector<double>> data;
    
public:
    Matrix() : n(0), data(0, std::vector<double>(0, 0.0)) {}
    Matrix(int size) : n(size), data(size, std::vector<double>(size, 0.0)) {}

    Matrix(const Matrix& other) : n(other.n), data(other.data) {}

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
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

    //вычсляем норму части столбца
    double SquareSumColumn(int column_number, int first_index) const {
        int n = this->GetHeight();
        double sum = 0.0;
        for (int i = first_index; i < n; i++) {
            sum += this->Get(i, column_number) * this->Get(i, column_number);
        }
        return std::sqrt(sum);
    }

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

    double Determinant() const {
        Matrix L, U;
        Vector P;
        int permutations;
        std::tie(L, U, P, permutations) = this->LU();
        
        double det = 1.0;
        for (int i = 0; i < n; i++) {
            det *= U.Get(i, i);
        }
        
        if (permutations % 2 != 0) {
            det = -det;
        }
        
        return det;
    }

    //вычисление обратной матрицы
    Matrix Inverse() const {
        Matrix L, U;
        Vector P;
        int permutations;
        std::tie(L, U, P, permutations) = this->LU();
        
        Matrix inverse(n);

        //на каждой итерации вычисляем j столбец обратной матрицы
        for (int j = 0; j < n; j++) {
            Vector e(n);
            e.Set(j, 1.0); //столбец единичной матрицы

            //переставляем e согласно вектору перестановок
            Vector permuted_e(n);
            for (int i = 0; i < n; i++) {
                permuted_e.Set(i, e.Get(static_cast<int>(P.Get(i))));
            }

            //прямая подстановка: Ly = b (b = permutated_e)
            Vector Y(n);
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int k = 0; k < i; k++) {
                    sum += L.Get(i, k) * Y.Get(k);
                }
                Y.Set(i, permuted_e.Get(i) - sum);
            }

            //обратная подстановка: Ux = y
            Vector X(n);
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0.0;
                for (int k = i + 1; k < n; k++) {
                    sum += U.Get(i, k) * X.Get(k);
                }
                if (std::abs(U.Get(i, i)) < 1e-9) {
                    throw std::runtime_error("Матрица вырождена или плохо обусловлена. Обратная матрица не существует.");
                }
                X.Set(i, (Y.Get(i) - sum) / U.Get(i, i));
            }
            
            for (int i = 0; i < n; i++) {
                inverse.Set(i, j, X.Get(i));
            }
        }
        
        return inverse;
    }

    //возвращает единичную матрицу
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

    /*
    Проверка на трехдиагональность
    - true: есть
    - false: нет
    */
    bool isTriDiagonal() const {
        int n = this->GetLength();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (std::abs(i - j) > 1 && std::abs(this->Get(i, j)) > 1e-9) {
                    return false;
                }
            }
        }
        return true;
    }

    /*
    Проверка на диагональное преобладение
    - true: есть
    - false: нет
    */
    bool isDiagonallyDominant() const {
        int n = this->GetLength();
        for (int i = 0; i < n; i++) {
            double diagonal_val = std::abs(this->Get(i, i));
            double sum_off_diagonal = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum_off_diagonal += std::abs(this->Get(i, j));
                }
            }
            if (diagonal_val <= sum_off_diagonal) {
                return false;
            }
        }
        return true;
    }

    /*
    возврат: (Q, R)
    Q - унитарная или ортогональная матрица
    R - верхнетреугольная матрица
    используется Хаусхолдер:
        - строим вектор v
        - делаем внешнее произведение
        - считаем матрицу Хаусхолдера H
        - делаем преобразования
    */
    std::pair<Matrix, Matrix> QR() {
        Matrix R = *this;
        Matrix Q = GetSingularMatrix(n);
        Vector v(n); //вектор Хаусхолдера

        //цикл по столбцам
        //каждая итерация k - обнуление элементов под главной диагональю в k-ом столбце
        for (int k = 0; k < n - 1; k++) {
            //строим вектор для текущего столбца k
            for (int i = 0; i < n; i++) {
                //выше главной диагонали в 0
                if (i < k) {
                    v.Set(i, 0.0);
                } else if (i == k) {
                    double norma = 0.0;
                    Vector column = R.GetColumn(k);  
                    //проходим по элементам столбца с диагонали   
                    for (int j = k; j < n; j++) {
                        norma += column.Get(j) * column.Get(j);
                    }
                    norma = std::sqrt(norma);
                    v.Set(i, R.Get(i, i) + std::copysign(norma, R.Get(i, i)));
                } else if (i > k) {
                    v.Set(i, R.Get(i, k));
                }
            }

            //это v * v^t (внешнее произведение)
            Matrix temp(n);
            double p = 0.0; //скалярный кэф из формулы Хаусхолдера
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    temp.Set(i, j, v.Get(i) * v.Get(j));
                }
                p += v.Get(i) * v.Get(i);
            }

            p = -2.0 / p;
            Matrix H = GetSingularMatrix(n) + (temp * p); //сама матрица

            R = H * R; //элементы обнуляются
            Q = Q * H; //накапливаем преобразование
        }

        return std::make_pair(Q, R);
    }

    /*
    возврат: (L, U, P, permutations)
    P - вектор индексов. Например: [2, 1, 0] -> 0 строка в конечной матрице изначально была второй и т.д.
    permutations - для определителя. det(A) = (-1)^k * det(U)
    принцип работы: 
            - проходим по всем столбцам
            - ищем главный элемент (если не диагональ, то перестановка)
            - для каждой строки под диагональю вычисляем кэф и обнуляем
    */
    std::tuple<Matrix, Matrix, Vector, int> LU() const {
        Matrix temp = *this;
        Matrix L = GetSingularMatrix(n);
        Vector P(n);
        int permutations = 0;
        for (int i = 0; i < n; i++) {
            P.Set(i, i);
        }

        //цикл по столбцам
        for (int k = 0; k < n; k++) {
            double max_val = std::abs(temp.Get(k, k));
            int max_row = k;
            //ищем настоящий главный элемент в столбце со следующей строки
            for (int i = k + 1; i < n; i++) {
                if (std::abs(temp.Get(i, k)) > max_val) {
                    max_val = std::abs(temp.Get(i, k));
                    max_row = i;
                }   
            }

            //если нашли главный элемент не на диагонали
            if (max_row != k) {
                //цикл для перестановки строк матриц поэлементно
                for (int j = 0; j < n; j++) {
                    double temp_u = temp.Get(k, j);
                    temp.Set(k, j, temp.Get(max_row, j));
                    temp.Set(max_row, j, temp_u);
                }
                //переставляем инфу в векторе перестановок
                double temp_p = P.Get(k);
                P.Set(k, P.Get(max_row));
                P.Set(max_row, temp_p);

                //переставляем элементы в L
                for (int j = 0; j < k; j++) {
                    double temp_l = L.Get(k, j);
                    L.Set(k, j, L.Get(max_row, j));
                    L.Set(max_row, j, temp_l);
                }
                permutations++;
            }

            //прямой ход метода Гаусса: строки под главной диагональю
            for (int i = k + 1; i < n; i++) {
                if (std::abs(temp.Get(k, k)) < 1e-9) {
                    throw std::runtime_error("Матрица вырождена или плохо обусловлена.");
                }
                //вычисляем кэф для домножения
                double factor = temp.Get(i, k) / temp.Get(k, k);
                L.Set(i, k, factor);
                //обнуляем
                for (int j = k; j < n; j++) {
                    temp.Set(i, j, temp.Get(i, j) - factor * temp.Get(k, j));
                }
            }
        }
        
        return std::make_tuple(L, temp, P, permutations);
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
    std::cout << "Исходная матрица:\n"; 
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

    std::cout << "Исходная матрица:\n"; 
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
    std::cout << "Исходная матрица:\n"; 
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