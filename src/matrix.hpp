#ifndef PROJECT_MATRIX_HPP
#define PROJECT_MATRIX_HPP

#include <cstring>
#include "iostream"
#include "vector"
#include <opencv.hpp>
#include "my_exception.cpp"

using namespace std;
using namespace cv;

template<typename T>
class matrix {
private:
    T **data;
    int cols;
    int rows;
    int channels=3;

public:
    /** generate matrix by size and fill it with all 0 */
    matrix(int r, int c) {
        T **d = new T *[r];
        for (int i = 0; i < r; i++) {
            T *line = new T[c];
            for (int j = 0; j < c; j++) {
                line[j] = 0;
            }
            d[i] = line;
        }
        data = d;
        rows = r;
        cols = c;
    }

    /** generate matrix by size and a 2-dimensional array */
    matrix(int r, int c, T **d, int ch = 3) : rows(r), cols(c), data(d), channels(ch) {};

    /** generate matrix by size and a 1-dimensional array */
    matrix(int r, int c, T *d) : rows(r), cols(c) {
        T **rs = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *line = new T[cols];
            for (int j = 0; j < cols; j++) {
                line[j] = d[i * cols + j];
            }
            rs[i] = line;
        }
        data = rs;
    };

    /** generate matrix by Mat in opencv
     * assume picture is BGR color space */
    matrix(Mat mat) {
        channels = mat.channels();
        rows = mat.rows;
        cols = channels * mat.cols;
        int **arr = new int *[rows];
        for (int i = 0; i < rows; i++) {
            int *row = new int[cols];
            for (int j = 0; j <cols; j++) {
                row[j] = mat.at<uchar>(i, j);
            }
            arr[i] = row;
        }
        data = arr;
    };

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

    void free_matrix(matrix<T> &m) {
        delete m.data;
    }

    /** convert to Mat*/
    Mat to_Mat() {
        Mat ans(rows, cols / channels, CV_8UC(channels),
                Scalar::all(0));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols ; j++) {
                ans.at<uchar>(i, j) = data[i][j];
            }
        }
        return ans;
    }

    /** print the content of matrix */
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

    /** Addition */
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
            throw my_exception("Addition: two matrixes should have same size!\n");
        }
    };

    /** Subtraction */
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
            throw my_exception("Substraction: two matrixes should have same size!\n");
        }
    };

    /** Scalar multiplication */
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

    /** Scalar multiplication */
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

    /** Scalar division */
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

    /** Transposition */
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

    // these four functions are used for function below - conjugation_matrix
    void conjugation(int i, int *arr, int j) { arr[j] = i; }

    void conjugation(double i, double *arr, int j) { arr[j] = i; }

    void conjugation(complex<int> i, complex<int> *arr, int j) {
        arr[j] = conj(i);
    }

    void conjugation(complex<double> i, complex<double> *arr, int j) {
        arr[j] = conj(i);
    }

    /** Conjugation: for complex number */
    matrix<T> conjugation_matrix() {
        // if the type is complex, then deal with it, or just copy
        if (strcmp(typeid(T).name(), "St7complexIiE") == 0 || strcmp(typeid(T).name(), "St7complexIdE") == 0) {
            auto **t = new T *[rows];
            for (int i = 0; i < rows; i++) {
                auto *t2 = new T[cols];
                for (int j = 0; j < cols; j++) {
                    conjugation(data[i][j], t2, j);
                }
                t[i] = t2;
            }
            return {rows, cols, t};
        } else {
            return copy_matrix();
        }
    }

    /** Element-wise multiplication: each element * element in the same position */
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
            throw my_exception("Element-wise multiplication: two matrixes should have same size!\n");
        }
    }

    /** Matrix-matrix multiplication */
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
            throw my_exception("Matrix-matrix multiplication: the size of two matrixes are not match!\n");
        }
    };

    /** Matrix * vector: assuming vector is a row vector */
    matrix<T> vector_mul(vector<T> v) {
        if (rows != v.size()) {
            throw my_exception("Matrix * vector: the size of two matrixes are not match!\n"
                               "The rows must equal to size of vector");
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

    /** Vector * matrix: assuming vector is a row vector */
    matrix<T> mul_vector(vector<T> v) {
        if (cols != 1) {
            throw my_exception("Vector * matrix : The size is not match! col size must be 1!\n");
        }
        int temp = v.size();
        T **t = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *t2 = new T[temp];
            for (int j = 0; j < v.size(); j++) {
                t2[j] = data[i][0] * v[j];
            }
            t[i] = t2;
        }
        return {rows, temp, t};
    }

    /** Dot product: same size matrix and col is 1, result is a num */
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
            throw my_exception("Dot product: dot product only applicable for two column matrix with the same dimension!\n");
        }
    }

    /** Cross product: appliable for two 3*1 matrix, a 3*1 matrix*/
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
            throw my_exception("Cross product only applicable for two 3-dimensional column matrix!\n");
        }
    }

    /** Convolution */
    matrix<T> convolution(matrix<T> &original_core,int div = 1) {
        if (original_core.rows != original_core.cols || original_core.rows % 2 != 1) {
            throw my_exception("Convolution: convolution core must be n*n and n must be an odd.\n");
        }

        // convert 180
        matrix<T> core(original_core.rows, original_core.cols);
        for (int i = 0; i < core.rows * core.cols; i++) {
            core.data[i / core.cols][i % core.cols] = original_core.data[core.rows - 1 - i / core.cols]
            [core.cols - 1 -i % core.cols];
        }

        // center
        int cen_row = (core.rows + 1) / 2 - 1;
        int cen_col = (core.cols + 1) / 2 - 1;

        // calculate the result
        T **d = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *line = new T[cols];
            for (int j = 0; j < cols; j++) {
                T result = 0;
                // calculate (i,j), traverse every element to calculate the result
                for (int m = 0; m < core.rows; m++) {
                    for (int n = 0; n < core.cols; n++) {
                        if (i + (m - cen_row) >= 0 && i + m - cen_row <= rows - 1 && j + 3*(n - cen_col) >= 0 &&
                            j + 3*(n - cen_col) <= cols - 1) {
                            result += core.data[m][n] * data[i + m - cen_row][j + 3*(n - cen_col)];
                        }
                    }
                }
                // set value
                result = result/div;
                if (result > 255) result = 255;
                if (result < 0) result = 0;
                line[j] = result;
            }
            d[i] = line;
        }
        matrix<T> m(rows, cols, d, channels);
        free_matrix(core);
        return m;
    }

    /** Find the minimum value */
    T max_value(int &loc_r, int &loc_c) const {
        T max = data[0][0];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if(max < data[i][j]){
                    max = data[i][j];
                    loc_r = i;
                    loc_c = j;
                }
            }
        }
        return max;
    }

    T max_value(char c, int n, int &loc) const {
        if(c=='r'){
            if(n<0 || n>=rows){
                throw my_exception("max_value: value range of n should be 0<=n<rows.\n");
            }
            T max = data[n][0];
            for (int i = 0; i < cols; i++) {
                if(max < data[n][i]){
                    max = data[n][i];
                    loc = i;
                }
            }
        } else if(c=='c'){
            if(n<0 || n>=cols){
                throw my_exception("min_value: value range of n should be 0<=n<cols.\n");
            }
            T max = data[0][n];
            for (int i = 0; i < rows; i++) {
                if(max < data[i][n]){
                    max = data[i][n];
                    loc = i;
                }
            }
        } else{
            throw my_exception("max_value: char c should be 'r' or 'c' to indicate row or column.\n");
        }
        return max;
    }

    T min_value(int &loc_r, int &loc_c) const {
        T min = data[0][0];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                for (int j = 0; j < cols; j++) {
                    if(min > data[i][j]){
                        min = data[i][j];
                        loc_r = i;
                        loc_c = j;
                    }
                }
            }
        }
        return min;
    }

    T min_value(char c, int n, int &loc) const {
        if(c=='r'){
            if(n<0 || n>=rows){
                throw my_exception("min_value: value range of n should be 0<=n<rows.\n");
            }
            T min = data[n][0];
            for (int i = 0; i < cols; i++) {
                if(min > data[n][i]){
                    min = data[n][i];
                    loc = i;
                }
            }
        } else if(c=='c'){
            if(n<0 || n>=cols){
                throw my_exception("min_value: value range of n should be 0<=n<cols.\n");
            }
            T min = data[0][n];
            for (int i = 0; i < rows; i++) {
                if(min > data[i][n]){
                    min = data[i][n];
                    loc = i;
                }
            }
        } else{
            throw my_exception("max_value: char c should be 'r' or 'c' to indicate row or column.\n");
        }
        return min;
    }

/** Sum all items*/
    T sum() const {
        T sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum += data[i][j];
            }
        }
        return sum;
    }

    T sum(char c, int n) const {
        T sum = 0;
        if(c=='r'){
            if(n<0 || n>=rows){
                throw my_exception("sum: value range of n should be 0<=n<rows.\n");
            }
            for (int i = 0; i < cols; i++) {
                sum += data[n][i];
            }
        } else if(c=='c'){
            if(n<0 || n>=cols){
                throw my_exception("sum: value range of n should be 0<=n<cols.\n");
            }
            for (int i = 0; i < rows; i++) {
                sum += data[i][n];
            }
        } else{
            throw my_exception("sum: char c should be 'r' or 'c' to indicate row or column.\n");
        }
        return sum;
    }

/** Average value*/
    T average() const {
        T sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum += data[i][j];
            }
        }
        return sum/(rows*cols);
    }

    T average(char c, int n) const {
        T sum = 0;
        if(c=='r'){
            if(n<0 || n>=rows){
                throw my_exception("average: value range of n should be 0<=n<rows.\n");
            }
            for (int i = 0; i < cols; i++) {
                sum += data[n][i];
            }
            return sum/cols;
        } else if(c=='c'){
            if(n<0 || n>=cols){
                throw my_exception("average: value range of n should be 0<=n<cols.\n");
            }
            for (int i = 0; i < rows; i++) {
                sum += data[i][n];
            }
            return sum/rows;
        } else{
            throw my_exception("average: char c should be 'r' or 'c' to indicate row or column.\n");
        }
    }

    /** Determinant*/
    T determinant(){
        try {
            if (rows != cols) {
                throw my_exception("determinant: the matrix is not a square matrix.\n");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }

        T m[rows][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                m[i][j] = 0;
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                m[i][j] = data[i][j];
            }
        }

        int t;
        int n_exch = 0;
        T max;
        T a=0,b=0;
        for(int i=0; i<rows-1; i++){
            t=i;
            max=m[i][i];
            for(int j=i+1; j<rows; j++){
                if(abs(max)<abs(m[j][i])){
                    max=m[j][i];
                    t=j;
                }
            }
            if(abs(max) == 0){
                return 0;
            }

            if(t!=i){
                n_exch++;
                for(int j=i; j<rows; j++){
                    a = m[i][j];
                    m[i][j] = m[t][j];
                    m[t][j] = a;
                }
            }

            for(int j=i+1; j<rows; j++){
                b = m[j][i]/m[i][i];
                for(int k=i; k<rows; k++){
                    m[j][k]=m[j][k]-m[i][k]*b;
                }
            }
        }
        b=1;
        for(int i=0; i<rows; i++)
            b*=m[i][i];
        if(n_exch%2==0)
            return b;
        else
            return -b;
    }

    /** Trace */
    T trace(){
        try {
            if (rows != cols) {
                throw my_exception("trace: the matrix is not a square matrix.\n");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }
        T tra = 0;
        for (int i = 0; i < rows; i++) {
            tra+=data[i][i];
        }
        return tra;
    }

    /** Inverse matrix*/
    matrix<T> inverse() {
        try {
            if (rows != cols || this->determinant() == 0) {
                throw my_exception("inverse: the matrix has no inverse matrix.\n");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }
        T m[rows][2*rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < 2*rows; j++) {
                m[i][j] = 0;
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                m[i][j] = data[i][j];
            }
            for (int j = rows; j < 2*rows; j++) {
                m[i][j] = 0;
            }
        }
        for (int i = 0; i < rows; i++) {
            m[i][i + rows] = 1;
        }
        int t;
        T max;
        T a, b;
        for (int i = 0; i < rows - 1; i++) {
            t = i;
            max = m[i][i];
            for (int j = i + 1; j < rows; j++) {
                if (abs(max) < abs(m[j][i])) {
                    max = m[j][i];
                    t = j;
                }
            }
            if (t != i) {
                for (int j = i; j < 2*rows; j++) {
                    a = m[i][j];
                    m[i][j] = m[t][j];
                    m[t][j] = a;
                }
            }
            for (int j = i + 1; j < rows; j++) {
                b = m[j][i]/m[i][i];
                for (int k = i; k < 2*rows; k++){
                    m[j][k] = m[j][k] - b*m[i][k];
                }
            }
        }

        for (int i = rows-1; i > 0; i--) {
            for (int k = i; k > 0; k--){
                b = m[k-1][i]/m[i][i];
                for (int j = i; j < 2*rows; j++) {
                    m[k-1][j] = m[k-1][j] - b*m[i][j];
                }
            }

        }
        for (int i=0; i<rows; i++){
            b = m[i][i];
            for (int j = i; j < 2*rows; j++) {
                m[i][j] = m[i][j]/b;
            }
        }
        auto **m2 = new T *[rows];
        for (int i = 0; i < rows; i++) {
            T *t2 = new T[rows];
            for (int j = 0; j < rows; j++){
                t2[j] = m[i][j+rows];
            }
            m2[i] = t2;
        }
        return {rows ,rows, m2};
    }

    /** Eigenvalues and eigenvectors*/
    void eigen(matrix<T> &eValue, matrix<T> &eVector){
        try {
            if(rows != cols){
                throw my_exception("eigen: the matrix is not a square matrix.\n");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }
        int N = 100;
        matrix<T> temp = this->copy_matrix();
        matrix<T> temp_Q(rows,rows);
        matrix<T> temp_R(rows,rows);
        for (int k = 0; k < N; ++k){
            temp.QR(temp_Q, temp_R);
            temp = temp_R*temp_Q;
        }
        for (int k = 0; k < rows; k++){
            eValue.data[k][0] = temp.data[k][k];
        }
        double e_Value;
        for (int count = 0; count < cols; count++){
            e_Value = eValue.data[count][0];
            temp = this->copy_matrix();
            for (int i = 0; i < rows; i++){
                temp.data[i][i] -= e_Value;
            }
            for (int i = 0; i < rows-1; i++){
                double coe = temp.data[i][i];
                for (int j = i; j < rows; j++){
                    temp.data[i][j] /= coe;
                }
                for (int i1 = i + 1; i1 < rows; i1++){
                    coe = temp.data[i1][i];
                    for (int j1 = i; j1 < rows; j1++){
                        temp.data[i1][j1] -= coe * temp.data[i][j1];
                    }
                }
            }
            double sum1 = eVector.data[rows - 1][count] = 1;
            for (int i2 = rows-2; i2 >= 0; i2--){
                double sum2 = 0;
                for (int j2 = i2 + 1; j2 < rows; j2++){
                    sum2 += temp.data[i2][j2] * eVector.data[j2][count];
                }
                sum2 = -sum2/temp.data[i2][i2];
                sum1 += sum2*sum2;
                eVector.data[i2][count] = sum2;
            }
            sum1 = sqrt(sum1);
            for (int i = 0; i < rows; i++){
                eVector.data[i][count] /= sum1;
            }
        }
        free_matrix(temp);
        free_matrix(temp_Q);
        free_matrix(temp_R);
    }

    void QR(matrix<T> &Q, matrix<T> &R)
    {
        T col_A[rows], col_Q[rows];

        if (rows != cols)
            printf("A is not a square matrix!");

        for (int j = 0; j < cols; j++){
            for (int i = 0; i < cols; i++){
                col_A[i] = data[i][j];
                col_Q[i] = data[i][j];
            }
            for (int k = 0; k < j; k++){
                R.data[k][j] = 0;
                for (int i = 0; i < rows; i++){
                    R.data[k][j] += col_A[i] * Q.data[i][k];
                }
                for (int i = 0; i < rows; i++){
                    col_Q[i] -= R.data[k][j] * Q.data[i][k];
                }
            }

            double temp = 0;
            for (int i = 0; i < rows; i++)
                temp+=(col_Q[i])*(col_Q[i]);
            temp = sqrt(temp);
            R.data[j][j] = temp;
            for (int i = 0; i < rows; i++){
                Q.data[i][j] = col_Q[i]/temp;
            }
        }
    }


    /** displaying detailed information  */
    void findInfo() {
        int op = -1;
        cout << "You have choose to check the information of the martix.\n"
             << "1--check for rows\n"
             << "2--check for columns\n"
             << "anything else--check for matrix\n";
        cin >> op;
        if (op == -1) {
            cout << "The specific row you want to check in row (integer number):";
            cin >> op;
            while (op <= 0 || op > rows) {
                cout << "Please resume load the correct number(the row serial number): ";
                cin >> op;
            }
            double max = data[op - 1][0], min = data[op - 1][0], sum = data[op - 1][0], aver;
            int site[2] = {1, 1};
            for (int i = 1; i < cols; ++i) {
                if (max < data[op - 1][i]) {
                    max = data[op - 1][i];
                    site[0] = i + 1;
                }
                if (min > data[op - 1][i]) {
                    min = data[op - 1][i];
                    site[1] = i + 1;
                }
                sum += data[op - 1][i];
            }
            aver = sum / cols;
            cout << "For all items of the row " << op << " :\n"
                 << "The sum is " << sum << endl
                 << "The average is " << aver << endl
                 << "The maximum item is " << max
                 << " lies in column " << site[0] << endl
                 << "The minimum item is " << min
                 << " lies in column " << site[1] << endl;
        } else if (op == 2) {
            cout << "The specific row you want to check in column (integer number):";
            cin >> op;
            while (op <= 0 || op > cols) {
                cout << "Please resume load the correct number";
                cin >> op;
            }
            double max = data[0][op - 1], min = data[0][op - 1], sum = data[0][op - 1], aver;
            int site[2] = {1, 1};
            for (int i = 1; i < rows; ++i) {
                if (max < data[i][op - 1]) {
                    max = data[i][op - 1];
                    site[0] = i + 1;
                }
                if (min > data[i][op - 1]) {
                    min = data[i][op - 1];
                    site[1] = i + 1;
                }
                sum += data[i][op - 1];
            }
            aver = sum / rows;
            cout << "For all items of the row " << op << " :\n"
                 << "The sum is " << sum << endl
                 << "The average is " << aver << endl
                 << "The maximum item is " << max
                 << " lies in row " << site[0] << endl
                 << "The minimum item is " << min
                 << " lies in row " << site[1] << endl;
        } else {
            double max = data[0][0], min = data[0][0], sum, aver;
            int site[2][2] = {(0, 0), (0, 0)};
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (max < data[i][j]) {
                        max = data[i][j];
                        site[0][0] = i;
                        site[0][1] = j;
                    }
                    if (min > data[i][j]) {
                        min = data[i][j];
                        site[1][0] = i;
                        site[1][1] = j;
                    }
                    sum += data[i][j];
                }
            }
            aver = sum / rows / cols;

            cout << "For all items of the whole matrix:\n"
                 << "The sum is " << sum << endl
                 << "The average is " << aver << endl
                 << "The maximum item is " << max
                 << " lies in (" << site[0][0] << "," << site[0][1] << ")" << endl
                 << "The minimum item is " << min
                 << " lies in (" << site[1][0] << "," << site[1][1] << ")" << endl;
        }
    }

    /** Slice matrix by choosing small block */
    matrix<T> BlockMatrix(int new_rows, int new_cols, int ser){
        try {
            if (rows % new_rows != 0 || cols % new_cols != 0) {
                throw my_exception("Slice the whole matrix fail for unreasonable block size.");
            }
        }
        catch (const my_exception &e1) {
            cout<<e1.what();
        }
        int r = rows / new_rows, c = cols / new_cols;
        int total = r * c;
        int r_ser = ser / r, c_ser = ser % c;
        try {
            if (ser >= total || ser < 0) {
                throw my_exception("Unreasonable serial of the block matrix.");
            }
        }
        catch (const my_exception &e2) {
            cout<<e2.what();
        }
        T **newData = new T*[new_rows];
        int a, b = 0;
        for (int i = r_ser * new_rows; i < (r_ser+1) * new_rows; ++i) {
            T *lin = new T[new_cols];
            for (int j = c_ser * new_cols; j < (c_ser+1) * new_cols; ++j) {
                lin[a] = data[i][j];
                a++;
            }
            a = 0;
            newData[b] = lin;
            b++;
        }
        return matrix<T>(new_rows, new_cols, newData);
    }

    /** Slice matrix */
    matrix<T> SilceMatrix
    (int rs, int cs, int re, int ce) {
        if(re<rs){
            int temp = rs;
            rs = re;
            re = temp;
        }
        if(ce<cs){
            int temp = cs;
            cs = ce;
            ce = temp;
        }
        try {
            if (rs < 0 || cs < 0 || re > rows || ce > cols) {
                throw my_exception("Out of boundary while choosing column or row.");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }
        T **newData = new T*[re - rs + 1];
        int lin_label = 0, nd_label = 0;
        for (int i = rs - 1; i < re; ++i) {
            T *lin = new T [ce - cs + 1];
            for (int j = cs - 1; j < ce; ++j) {
                lin[lin_label]= data[i][j];
                lin_label++;
            }
            lin_label = 0;
            newData[nd_label] = lin;
            nd_label++;
        }
        return matrix<T>(re - rs + 1, ce - cs + 1, newData);
    }

    /** Reshape enlarging matrix */
    matrix<T> ReshapeMatrix
    (int new_rows, int new_cols, int row_start, int col_start) {
        try {
            if (row_start + rows - 1 > new_rows || col_start + cols - 1 > new_cols) {
                throw my_exception("Reshaping matrix fail for unreasonable size for new matrix\n");
            }
        }
        catch (const my_exception &e) {
            cout<<e.what();
        }
        T **newData = new T *[new_rows];
        int r_data = 0;
        for (int i = 0; i < new_rows; ++i) {
            T *lin = new T [new_cols];
            int c_data = 0;
            for (int j = 0; j < new_cols; ++j) {
                if (i >= row_start - 1 && j >= col_start - 1
                    && i <= row_start + rows - 2 && j <= col_start + cols - 2) {
                    lin[j] = data[r_data][c_data];
                    if(c_data == 7){
                        r_data++;
                    }
                    c_data++;
                }else{
                    lin[j] = 0;
                }
            }
            newData[i] = lin;
        }
        return matrix<T>(new_rows, new_cols, newData);
    }

    matrix<T> ReshapeMatrix(int new_rows, int new_cols) {
        try {
            if (new_rows * new_cols != cols * rows) {
                throw my_exception("Reshaping matrix fail for unreasonable size for new matrix.");
            }
        }
        catch(const my_exception &e){
            cout<<e.what();
        }

        T *arr = new T [rows*cols];
        int arr_label = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                arr[arr_label] = data[i][j];
                arr_label++;
            }
        }
        arr_label = 0;
        T **d = new T *[new_rows];
        for (int i = 0; i < new_rows; ++i) {
            T *lin = new T [new_cols];
            for (int j = 0; j < new_cols; ++j) {
                lin[j] = arr[arr_label];
                arr_label++;
            }
            d[i] = lin;
        }
        return matrix<T> (new_rows,new_cols,d);
    }

    /** get data */
    int getCols(){return cols;}
    int getRows(){return rows;}
    T** getData(){return data;}
};




#endif //PROJECT_MATRIX_HPP
