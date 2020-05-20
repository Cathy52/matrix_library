#ifndef PROJECT_COMPLEX_HPP
#define PROJECT_COMPLEX_HPP

#include "iostream"

using namespace std;

class complex {
public:
    complex() : real(0), img(0) {};

    complex(double r) : real(r), img(0) {};

    complex(double r, double i) : real(r), img(i) {};

    friend ostream &operator<<(ostream &os, complex &c) {
        if (c.real == 0 && c.img == 0) {
            os << 0;
        } else if (c.real == 0) {
            os << c.img << "i";
        } else if (c.img == 0) {
            os << c.real;
        } else {
            if(c.img > 0) os << c.real << " + " << c.img << "i";
            else  os << c.real << c.img << "i";
        }
        return os;
    };

    friend complex operator+(complex &c1, complex &c2) {
        return {c1.real + c2.real, c1.img + c2.img};
    }

    friend complex operator+(double i, complex &c2) {
        return {i + c2.real, c2.img};
    }

    friend complex operator+(complex &c2, double i) {
        return {i + c2.real, c2.img};
    }

    friend complex operator-(complex &c1, complex &c2) {
        return {c1.real - c2.real, c1.img - c2.img};
    }

    friend complex operator-(double i, complex &c2) {
        return {i - c2.real, -c2.img};
    }

    friend complex operator-(complex &c2, double i) {
        return {-i + c2.real, c2.img};
    }

    friend complex operator*(complex &c1, complex &c2) {
        return {c1.real * c2.real - c1.img * c2.img, c1.img * c2.real + c1.real * c2.img};
    }

    friend complex operator*(double i, complex &c2) {
        return {i * c2.real, i * c2.img};
    }

    friend complex operator*(complex &c2, double i) {
        return {i * c2.real, i * c2.img};
    }

    friend complex operator/(complex &c2, double i) {
        return {c2.real/i, i * c2.img/i};
    }

    friend complex operator~(complex & c){
        return {c.real,c.img*-1};
    }

private:
    double real;
    double img;
};


#endif //PROJECT_COMPLEX_HPP
