#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Matr
{
private:
    int size;
    double **mas;
    double *mas1;

public:
    Matr()
    {
        size = 0;
        mas = nullptr;
        mas1 = nullptr;
    }

    Matr(int l)
    {
        size = l;
        mas = new double*[l];
        for (int i = 0; i < l; i++)
            mas[i] = new double[l];
        mas1 = new double[l];
    }

    void Add()
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                cin >> mas[i][j];
            cin >> mas1[i];
        }
    }

    void Print()
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                cout << setw(4) << mas[i][j] << " ";
            }
            cout << " " << mas1[i] << endl;
        }
    }

    void Preob()
    {
        double temp = 0;
        for (int k = 0; k < size; k++)
        {
            for (int i = 0; i < size; i++)
            {
                temp = mas[i][i] * (-1);
                mas1[i] /= temp;
                for (int j = 0; j < size; j++)
                {
                    mas[i][j] /= temp;
                }
            }
        }
        for (int i = 0; i < size; i++)
        {
            mas1[i] *= -1;
            for (int j = 0; j < size; j++)
                mas[i][i] = 0;
        }
    }

    double Pogr(double **mas, double epsilon)
    {
        double eps = 0, sum = 0, max = 0;
        double norm1 = 0, norm2 = 0;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < i; j++)
            {
                sum += fabs(mas[i][j]);
                if (sum > norm1) norm1 = sum;
            }
            sum = 0;
            for (int j = i + 1; j < size; j++)
            {
                sum += fabs(mas[i][j]);
                if (sum > norm2) norm2 = sum;
            }
            sum = 0;
        }
        if (norm1 >= 1 || norm2 >= 1)
        {
            cerr << "Норма матрицы больше или равна 1." << endl;
            exit(1);
        }
        eps = ((1 - norm1) / norm2) * epsilon;
        return eps;
    }

    void Itera()
    {
        auto *x = new double[size];
        double *p = new double[size];
        double *a = new double[size];
        double *E = new double[size];
        double per = Pogr(mas, 0.0001), max = 0; // Фиксированное значение epsilon = 0.0001
        for (int i = 0; i < size; i++)
        {
            x[i] = mas1[i];
            p[i] = 0;
        }
        double var = 0;
        for (int i = 0; i < size; i++)
        {
            var = 0;
            for (int k = 0; k < size; k++)
                var += mas[i][k] * mas1[k];
            x[i] = var;
        }
        for (int i = 0; i < size; i++)
            p[i] = x[i] + mas1[i];
        int counter = 0;
        do
        {
            counter++;
            for (int i = 0; i < size; i++)
            {
                var = 0;
                for (int j = 0; j < i; j++)
                    var += (mas[i][j] * p[j]);
                for (int j = i + 1; j < size; j++)
                    var += (mas[i][j] * x[j]);
                a[i] = var;
                x[i] = mas1[i] + a[i];
            }
            max = 0;
            for (int i = 0; i < size; i++)
            {
                E[i] = fabs(x[i] - p[i]);
                if (max < E[i]) max = E[i];
                p[i] = x[i];
            }
        } while (max > per);

        cout << "Решение системы:" << endl;
        for (int i = 0; i < size; i++)
        {
            cout << "x[" << setw(1) << i + 1 << "] = " << setw(6) << setprecision(3) << fixed << x[i] << endl;
        }

        delete[] x;
        delete[] p;
        delete[] E;
        delete[] a;
    }

    ~Matr()
    {
        for (int i = 0; i < size; i++)
            delete[] mas[i];
        delete[] mas;
        delete[] mas1;
    }
};

int main()
{
    int n;
    cout << "Введите размер матрицы (n x n): ";
    cin >> n;
    Matr a(n);
    cout << "Введите элементы расширенной матрицы построчно (каждый ряд содержит " << n << " элемента матрицы A, затем элемент b):" << endl;
    a.Add();
    cout << endl << "Расширенная матрица:" << endl;
    a.Print();
    a.Preob();
    cout << endl << "Преображенная матрица:" << endl;
    a.Print();
    cout << endl;
    a.Itera();
    return 0;
}