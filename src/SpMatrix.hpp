//
// Created by Trian on 15/6/2020.
//

#ifndef PROJ_MATRIX_SPMATRIX_HPP
#define PROJ_MATRIX_SPMATRIX_HPP
#include <iostream>
#include <cstring>
#include "matrix.hpp"
using namespace std;
template<typename T>
class spmatrix {
private:
    /** CSR */
    int zero_num, ro_len, data_len;
    int *row_offset;
    int *col_indice;
    T *data;
public:
    spmatrix();
    spmatrix(int *r, int *c, T *d, int zn, int ro_l, int da_l)
    : row_offset(r), col_indice(c), data(d), zero_num(zn), ro_len(ro_l), data_len(da_l){
        cout<<"Sparse matrix create successfully"<<endl;
    };
    ~spmatrix(){
        delete row_offset;
        delete col_indice;
        delete data;
    }

    /** print the content of the spmatrix */
    friend ostream &operator<<(ostream &os, spmatrix<T> &spm) {
        int r = spm.ro_len;
        int ci_len = spm.data_len, da_len = spm.data_len;
        cout << "row_offset:[ ";
        for (int i = 0; i < r; i++) {
            cout<<spm.row_offset[i];
            if(i == r-1){
                cout<<"]\n";
                break;
            }
            cout<<", ";
        }
        cout << "col_indice:[ ";
        for (int j = 0; j < ci_len; ++j) {
            cout<<spm.col_indice[j];
            if(j == ci_len-1){
                cout<<"]\n";
                break;
            }
            cout<<", ";
        }
        cout<<"value:[ ";
        for (int k = 0; k < da_len; ++k) {
            cout<<spm.data[k];
            if(k == da_len-1){
                cout<<"]\n";
                break;
            }
            cout<<", ";
        }
        cout<<"Number of zero elements: "<<spm.zero_num<<endl;
        return os;
    };

    /** vector multiply -> one column vector */
    matrix<T> sp_vector_mul (matrix<T> &m){
        if((zero_num+data_len)%(ro_len-1)!=0){
            throw my_exception("Data in Sparse matrix is wrong, please check.");
        }
        int c = (zero_num+data_len)/(ro_len-1);
        if(m.getRows() != c || m.getCols() !=1){
            throw my_exception("Unreasonable size for vector");
        }
        T **d = new T*[ro_len-1];
        int value_label = 0;
        for (int i = 0; i < ro_len-1; ++i) {
            T *lin = new T [1];
            int count = 0;
            T sum = 0;
            for (int j = 0; j < c; ++j) {
                if(col_indice[value_label] == j && count < row_offset[i+1]-row_offset[i]){
                    sum += data[value_label]*m.getData()[j][0];
                    value_label++;
                    count ++;
                }
            }
            lin[0] = sum;
            d[i] = lin;
        }
        return matrix<T>(ro_len-1,1,d);
    }

    /** recover the sparse matrix to matrix */
    matrix<T> recover_matrix()  {
        if((zero_num+data_len)%(ro_len-1)!=0){
            throw my_exception("Data in Sparse matrix is wrong, please check.");
        }
        int c = (zero_num+data_len)/(ro_len-1);
        T **d = new T *[ro_len-1];
        int value_label = 0;
        for (int i = 0; i < ro_len-1; ++i) {
            T *lin = new T [c];
            int count = 0;
            for (int j = 0; j < c; ++j) {
                if(j==col_indice[value_label] && row_offset[i]!= -1
                   && count<row_offset[i+1]-row_offset[i] ){
                    lin[j] = data[value_label];
                    value_label++;
                    count++;
                    continue;
                }
                lin[j] = 0;
            }
            d[i] = lin;
        }
        return matrix<T>(ro_len-1, c, d);
    }
};
/** compress matrix into sparse matrix */
template<typename T>
spmatrix<T> compress_matrix(matrix<T> &m) {
    int non_zero_num = 0;
    for (int i = 0; i < m.getRows(); ++i) {
        for (int j = 0; j < m.getCols(); ++j) {
            if(m.getData()[i][j]!=0){
                non_zero_num++;
            }
        }
    }
    if(non_zero_num*100/m.getCols()/m.getRows()>10){
        cout<<"This matrix may not be suitable to compress.\n";
    }
    T *d = new T [non_zero_num];
    int *ro = new int[m.getRows()+1];
    int *ci = new int[non_zero_num];
    int zero_num = m.getRows()*m.getCols()-non_zero_num;
    int d_label = 0, ro_label = 0;
    int ro_value = 0;
    for (int i = 0; i < m.getRows(); ++i) {
        int temp = 0;
        for (int j = 0; j < m.getCols(); ++j) {
            if(m.getData()[i][j]!=0){
                if(temp == 0){
                    ro[ro_label] = ro_value;
                    ro_label++;
                }
                ci[d_label] = j;
                d[d_label] = m.getData()[i][j];
                d_label++;
                temp++;
            }
        }
        if(temp == 0){
            ro[ro_label] = ro_value;
            ro_label++;
        }
        ro_value +=temp;
    }
    ro[m.getRows()] = non_zero_num;
    return spmatrix<T>(ro, ci, d, zero_num, m.getRows()+1, non_zero_num);
}
#endif //PROJ_MATRIX_SPMATRIX_HPP
