#ifndef PROJECT_MATRIX_HPP
#define PROJECT_MATRIX_HPP

#include "iostream"
#include "complex.hpp"

using namespace std;

template<typename T>
struct is_complex {
    operator bool() {
        return false;
    }
};

template<>
struct is_complex<complex> {
    operator bool() {
        return true;
    }
};

template<typename T>
class matrix {
public:
    matrix(int r, int c, T **d) : rows(r), cols(c), data(d) {};

    matrix(int r, int c, T *d) : rows(r), cols(c) {
        T **rs = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *line = new T[cols];
            for (int j = 0; j < cols; j++) {
                line[j] = d[i*cols+j];
            }
            rs[i] = line;
        }
        data = rs;
    };

    friend ostream &operator<<(ostream &os, matrix<T> &m) {
        for (int i = 0; i < m.rows; i++) {
            if (i == 0) cout << "[";
            else cout << " ";
            for (int j = 0; j < m.cols; j++) {
                os << m.data[i][j];
                if (j != m.cols - 1) cout << ", ";
            }
            if (i == m.rows - 1) os << "]\n";
            else os << ";\n";
        }
        return os;
    };

    matrix<T> operator+(matrix<T> &m) {
        if (m.rows == rows && m.cols == cols) {
            T **rs = new T *[rows];
            for (int i = 0; i < rows; i++) {
                T *line = new T[cols];
                for (int j = 0; j < cols; j++) {
                    line[j] = data[i][j] + m.data[i][j];
                }
                rs[i] = line;
            }
            matrix<T> re(rows, cols, rs);
            return re;
        } else {
            cout << "ERROR: Size is different\n" << endl;
            exit(1);
        }
    };

    matrix<T> operator-(matrix<T> &m) {
        if (m.rows == rows && m.cols == cols) {
            T **rs = new T *[rows];
            for (int i = 0; i < rows; i++) {
                T *line = new T[cols];
                for (int j = 0; j < cols; j++) {
                    line[j] = data[i][j] - m.data[i][j];
                }
                rs[i] = line;
            }
            matrix<T> re(rows, cols, rs);
            return re;
        } else {
            cout << "ERROR: Size is different\n" << endl;
            exit(1);
        }
    };

    matrix<T> operator*(matrix<T> &m) {
        if (m.rows == cols) {
            T **rs = new T *[rows];
            for (int i = 0; i < rows; i++) {
                T *line = new T[m.cols];
                for (int j = 0; j < m.cols; j++) {
                    // line[j] = line i * col j
                    T tmp = 0;
                    T temp = 0;
                    for (int z = 0; z < cols; z++) {
                        temp = data[i][z] * m.data[z][j];
                        tmp = tmp + temp;
                    }
                    line[j] = tmp;
                }
                rs[i] = line;
            }
            matrix<T> re(rows, m.cols, rs);
            return re;
        } else {
            cout << "ERROR: Size is not match\n" << endl;
            exit(1);
        }
    };

    friend matrix<T> operator*(double num, matrix<T> &m) {
        T **rs = new T *[m.rows];
        for (int i = 0; i < m.rows; i++) {
            T *line = new T[m.cols];
            for (int j = 0; j < m.cols; j++) {
                T tmp = 0;
                line[j] = m.data[i][j] * num;
            }
            rs[i] = line;
        }
        matrix<T> re(m.rows, m.cols, rs);
        return re;
    };

    friend matrix<T> operator*(matrix<T> &m, double num) {
        T **rs = new T *[m.rows];
        for (int i = 0; i < m.rows; i++) {
            T *line = new T[m.cols];
            for (int j = 0; j < m.cols; j++) {
                T tmp = 0;
                line[j] = m.data[i][j] * num;
            }
            rs[i] = line;
        }
        matrix<T> re(m.rows, m.cols, rs);
        return re;
    };

    friend matrix<T> operator/(matrix<T> &m, double num) {
        T **rs = new T *[m.rows];
        for (int i = 0; i < m.rows; i++) {
            T *line = new T[m.cols];
            for (int j = 0; j < m.cols; j++) {
                T tmp = 0;
                line[j] = m.data[i][j] / num;
            }
            rs[i] = line;
        }
        matrix<T> re(m.rows, m.cols, rs);
        return re;
    };

    void free_matrix(matrix<T> &m) {
        delete m.data;
    }

    matrix<T> copy_matrix() {
        T **t = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *t2 = new T[cols];
            for (int j = 0; j < cols; j++) {
                t2[j] = data[i][j];
            }
            t[i] = t2;
        }
        return {rows, cols, t};
    }

    matrix<T> tran_matrix() {
        T **t = new T *[cols];
        for (int i = 0; i < cols; i++) {
            T *t2 = new T[rows];
            for (int j = 0; j < rows; j++) {
                t2[j] = data[j][i];
            }
            t[i] = t2;
        }
        return {cols, rows, t};
    }

    matrix<T> conjugation_matrix() {
        if (is_complex<T>() == 1) {
            auto **t = new T *[rows];
            for (int i = 0; i < rows; i++) {
                auto *t2 = new T[cols];
                for (int j = 0; j < cols; j++) {
                    t2[j] = ~data[i][j];
                }
                t[i] = t2;
            }
            return {rows, cols, t};
        } else {
            return copy_matrix();
        }
    }

    matrix<T> elem_wise_mul(matrix<T> &m) {
        if (cols == m.cols && rows == m.rows) {
            T **t = new T *[rows];
            for (int i = 0; i < rows; i++) {
                T *t2 = new T[rows];
                for (int j = 0; j < cols; j++) {
                    t2[j] = m.data[i][j] * data[i][j];
                }
                t[i] = t2;
            }
            return {rows, cols, t};
        } else {
            cout << "ERROR: Size is different\n" << endl;
            exit(1);
        }
    }

    // default: assuming vector is a row vector
    matrix<T> vector_mul(vector<T> v) {
        if (rows != v.size()) {
            cout << "ERROR: The size is not match! The rows must equal to size of vector" << endl;
            exit(1);
        }
        T **t = new T *[1];
        T *t2 = new T[cols];
        for (int i = 0; i < cols; i++) {
            T tmp = 0;
            T temp = 0;
            for (int j = 0; j < rows; j++) {
                temp = v[j] * data[j][i];
                tmp = tmp + temp;
            }
            t2[i] = tmp;
        }
        t[0] = t2;
        return {1, cols, t};
    }

    matrix<T> mul_vector(vector<T> v) {
        if (cols != 1) {
            cout << "ERROR: The size is not match! col size must be 1" << endl;
            exit(1);
        }
        T ** t = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *t2 = new T[v.size()];
            for (int j = 0; j < v.size(); j++) {
                t2[j] = data[i][0] * v[j];
            }
            t[i] = t2;
        }
        return {rows, v.size(), t};
    }

    // 只适用于：维数相同的列矩阵，且结果是标量
    T dot_product(matrix<T> &m) {
        if (rows == m.rows && cols == 1 && m.cols == 1) {
            T tmp = 0;
            T temp = 0;
            for (int i = 0; i < rows; i++) {
                temp = data[i][0] * m.data[i][0];
                tmp = tmp + temp;
            }
            return tmp;
        } else {
            cout << "ERROR: dot product only applicable for two column matrix with the same dimension" << endl;
            exit(1);
        }
    }

    // 只适用于：两个三维列矩阵，且结果为三维列矩阵
    matrix<T> cross_product(matrix<T> &m) {
        if (m.rows == 3 && rows == 3 && m.cols == 1 && cols == 1) {
            T **t = new T *[3];

            T *t1 = new T[1];
            T *t2 = new T[1];
            T *t3 = new T[1];
            T tmp1 = data[1][0] * m.data[2][0]; // y * z'
            T tmp2 = m.data[1][0] * data[2][0]; // y' * z
            t1[0] = tmp1 - tmp2;
            tmp1 = data[2][0] * m.data[0][0]; // z * x'
            tmp2 = data[0][0] * m.data[2][0]; // x * z'
            t2[0] = tmp1 - tmp2;
            tmp1 = data[0][0] * m.data[1][0]; // x * y'
            tmp2 = data[1][0] * m.data[0][0]; // y * x'
            t3[0] = tmp1 - tmp2;
            t[0] = t1;
            t[1] = t2;
            t[2] = t3;
            return {3, 1, t};
        } else {
            cout << "ERROR: cross product only applicable for two 3-dimensional column matrix" << endl;
            exit(1);
        }
    }

private:
    T **data;
    int cols;
    int rows;
};

#endif //PROJECT_MATRIX_HPP
