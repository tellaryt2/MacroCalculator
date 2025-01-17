#pragma once

#ifndef COMPLEXARRAY_H
#define COMPLEXARRAY_H

#include "complex.h"
#include <vector>

class ComplexArray {
private:
    std::vector<complex> array;

public:
    // ����������� �� ���������
    ComplexArray() {}

    // ���������� ��������
    void add(const complex& c) {
        array.push_back(c);
    }

    // ��������� �������� �� �������
    complex& operator[](size_t index) {
        return array[index];
    }

    // ��������� ������� �������
    size_t size() const {
        return array.size();
    }

    // �������� ���������� ��������
    void removeLast() {
        if (!array.empty()) {
            array.pop_back();
        }
    }

    //ComplexArray arr;
    //arr.add(complex(1.0, 2.0));
    //arr.add(complex(3.0, 4.0));
    //complex sum = arr[0] + arr[1];

};

#endif 