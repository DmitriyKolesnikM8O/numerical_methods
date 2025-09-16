#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>
#include "Vector.hpp"
#include "Matrix.hpp"

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

int main() {
    Matrix A = get_user_matrix_input();
    double epsilon = get_user_epsilon();

    // A.Set(0, 0, 1.0); A.Set(0, 1, 3.0); A.Set(0, 2, 1.0);
    // A.Set(1, 0, 1.0); A.Set(1, 1, 1.0); A.Set(1, 2, 4.0);
    // A.Set(2, 0, 4.0); A.Set(2, 1, 3.0); A.Set(2, 2, 1.0);

    Task5::Do(A, epsilon);
    return 0;
}