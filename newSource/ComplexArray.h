#pragma once

#ifndef COMPLEXARRAY_H
#define COMPLEXARRAY_H

#include "complex.h"
#include <vector>

class ComplexArray {
private:
    std::vector<complex> array;

public:
    // Конструктор по умолчанию
    ComplexArray() {}

    // Добавление элемента
    void add(const complex& c) {
        array.push_back(c);
    }

    // Получение элемента по индексу
    complex& operator[](size_t index) {
        return array[index];
    }

    // Получение размера массива
    size_t size() const {
        return array.size();
    }

    // Удаление последнего элемента
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