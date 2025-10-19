#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <stdexcept>
#include <initializer_list>


class Vector {
private:
    int size;
    std::vector<double> data;

public:
    Vector(std::initializer_list<double> il) : size(il.size()), data(il) {}
    Vector() : size(0), data(0, 0.0) {}
    Vector(int n) : size(n), data(n, 0.0) {}

    double Get(int i) const {
        if (i < 0 || i >= size) throw std::out_of_range("Индекс не в границах вектора");
        return data[i];
    }

    void Set(int i, double value) {
        if (i < 0 || i >= size) throw std::out_of_range("Индекс не в границах вектора");
        data[i] = value;
    }

    int GetSize() const { return size; }
};

#endif // VECTOR_H