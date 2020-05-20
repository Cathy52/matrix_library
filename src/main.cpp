#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "complex.hpp"

using namespace std;

void test_int();

void test_double();

void test_complex();

int main() {
    int op = 0;
    cout << "Please enter an integer to test different types(1-int,2-double,3-complex): " << endl;
    cin >> op;
    if (op == 1) {
        test_int();
    } else if (op == 2) {
        test_double();
    } else {
        test_complex();
    }
}

void test_int() {
    int arr1[2] = {1, 2};
    int arr2[2] = {1, 2};
    int a1[1] = {1};
    int a2[1] = {2};
    int *arr3[2] = {arr1, arr2};
    int *arr4[2] = {a1, a2};
    matrix<int> m1(2, 2, arr3);
    matrix<int> m2(2, 2, arr3);
    matrix<int> m3(2, 1, arr4);
    matrix<int> m4 = m1 + m2;
    matrix<int> m5 = m4 - m2;
    matrix<int> m6 = m1 * m3;
    matrix<int> m7 = 2 * m3;
    matrix<int> m8 = m3 * 3;
    matrix<int> m0 = m8 / 3;
    cout << "m1 = \n" << m1;
    cout << "m2 = \n" << m2;
    cout << "m3 = \n" << m3;
    cout << "m4 = m1+m2 = \n" << m4;
    cout << "m5 = m4-m2 = \n" << m5;
    cout << "m6 = m1*m3 = \n" << m6;
    cout << "m7 = 2*m3 = \n" << m7;
    cout << "m8 = m3*3 = \n" << m8;
    cout << "m0 = m8/3 = \n" << m0;
    matrix<int> m9 = m5.tran_matrix();
    cout << "m9 = m5's trans = \n" << m9;
    matrix<int> m10 = m5.conjugation_matrix();
    cout << "m10 = m5's conju = \n" << m10;
    matrix<int> ma1 = m9.elem_wise_mul(m10);
    cout << "ma1 = m9 elem_wise_mul m10 = \n" << ma1;
    vector<int> v;
    v.push_back(1);
    v.push_back(2);
    matrix<int> ma2 = ma1.vector_mul(v);
    cout << "ma2 = vector(1,2) * ma1 = \n" << ma2 << endl;
    matrix<int> ma3 = m0.mul_vector(v);
    cout << "ma3 = m0 * vector(1,2) = \n" << ma3 << endl;
    int ar1[3] = {1,2,3};
    int ar2[3] = {1,-1,2};
    matrix<int> mat1(3,1,ar1);
    matrix<int> mat2(3,1,ar2);
    int result = mat1.dot_product(mat2);
    matrix<int> mat3 = mat1.cross_product(mat2);
    cout << "mat1 = \n" << mat1 << endl;
    cout << "mat2 = \n" << mat2 << endl;
    cout << "result = mat1 dot product mat2 = " << result << endl;
    cout << "mat3 = mat1 cross product mat2 = \n" << mat3 << endl;
};

void test_double() {
    double arr1[2] = {1.1, 2.3};
    double arr2[2] = {1.2, 2.2};
    double a1[1] = {1.2};
    double a2[1] = {2.4};
    double *arr3[2] = {arr1, arr2};
    double *arr4[2] = {a1, a2};
    matrix<double> m1(2, 2, arr3);
    matrix<double> m2(2, 2, arr3);
    matrix<double> m3(2, 1, arr4);
    matrix<double> m4 = m1 + m2;
    matrix<double> m5 = m4 - m2;
    matrix<double> m6 = m1 * m3;
    matrix<double> m7 = 2 * m3;
    matrix<double> m8 = m3 * 3;
    cout << "m1 = \n" << m1;
    cout << "m2 = \n" << m2;
    cout << "m3 = \n" << m3;
    cout << "m4 = m1+m2 = \n" << m4;
    cout << "m5 = m4-m2 = \n" << m5;
    cout << "m6 = m1*m3 = \n" << m6;
    cout << "m7 = 2*m3 = \n" << m7;
    cout << "m8 = m3*3 = \n" << m8;
};

void test_complex() {
    complex c1(1, 1);
    complex c2(2, 1);
    complex c3 = c1 + c2;
    complex c4 = c2 - c1;
    complex c5 = c1 * c2;
    cout << "c1 = " << c1 << endl;
    cout << "c2 = " << c2 << endl;
    cout << "c3 = c1+c2 = " << c3 << endl;
    cout << "c4 = c2-c1 = " << c4 << endl;
    cout << "c5 = c1*c2 = " << c5 << endl;
    complex arr1[2] = {c1, c2};
    complex arr2[2] = {c3, c4};
    complex *p[2] = {arr1, arr2};
    matrix<complex> m1(2, 2, p);
    matrix<complex> m2(2, 2, p);
    matrix<complex> m3 = m1 + m2;
    matrix<complex> m4 = m1 * m2;
    matrix<complex> m5 = 2 * m2;
    matrix<complex> m6 = m2 * 4;
    matrix<complex> m0 = m6 / 2;

    cout << "m1 = \n" << m1 << endl;
    cout << "m2 = \n" << m2 << endl;
    cout << "m3 = m1+m2 = \n" << m3 << endl;
    cout << "m4 = m1*m2 = \n" << m4 << endl;
    cout << "m5 = 2 *m2 = \n" << m5 << endl;
    cout << "m6 = m2* 4 = \n" << m6 << endl;
    cout << "m0 = m6/ 2 = \n" << m0 << endl;
    matrix<complex> m7 = m6.tran_matrix();
    cout << "m7 = m6'trans = \n" << m7 << endl;
    matrix<complex> m8 = m6.conjugation_matrix();
    cout << "m8 = m6'conju = \n" << m8 << endl;
    matrix<complex> m9 = m8.elem_wise_mul(m6);
    cout << "m9 = m8 elem_wise_mul m6 = \n" << m9 << endl;
    vector<complex> v;
    v.push_back(1);
    v.push_back(2);
    matrix<complex> ma1 = m8.vector_mul(v);
    cout << "ma1 = vector(1,2) * ma1 = \n" << ma1 << endl;
    complex a1[1]={c1};
    complex a2[1]={c2};
    complex *p1[2]={a1,a2};
    matrix<complex> ma2(2,1,p1);
    matrix<complex> ma3 = ma2.mul_vector(v);
    cout << "ma2 = \n" << ma2 << endl;
    cout << "ma3 = ma1 * vector(1,2) = \n" << ma3 << endl;
    complex ar1[3] = {complex(1,2),complex(1,3),complex(1,4)};
    complex ar2[3] = {complex(1,-1),complex(1,0),complex(1,1)};
    matrix<complex> mat1(3,1,ar1);
    matrix<complex> mat2(3,1,ar2);
    cout << "mat1 = \n" << mat1 << endl;
    cout << "mat2 = \n" << mat2 << endl;
    complex result = mat1.dot_product(mat2);
    cout << "result = mat1 dot product mat2 = " << result << endl;
    matrix<complex> mat3 = mat1.cross_product(mat2);
    cout << "mat3 = mat1 cross product mat2 = \n" << mat3 << endl;
}
