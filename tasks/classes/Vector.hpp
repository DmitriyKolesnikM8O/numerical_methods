#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <stdexcept>


class Vector {
private:
    int size;
    std::vector<double> data;

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

#endif // VECTOR_H