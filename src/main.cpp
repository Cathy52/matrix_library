#include <iostream>
#include <vector>
#include <complex>
#include "opencv.hpp"
#include "SpMatrix.hpp"
#include <opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

void test_img(Mat, Mat);

void test_cal();

void test_matrix_computation_int();

void test_matrix_computation_double();

void test_matrix_computation_complex();

void pic_process_gray();

void pic_process();

void pic_blend();

void test_sparse_matrix();

void test_SliceAndReshape();

void test_detailed_matrix();

void test_exception_handling();

int main() {
    /** 支持多种数据类型：基本数据类型，以及自定义数据类型
     *  实现了矩阵的多种运算，包括：
     *  Addition, subtraction, scalar multiplication, scalar division, transposition, conjugation,
     *  element-wise multiplication, matrix-matrix multiplication, matrix-vector multiplication, dot product
     *  and cross product.
     */
//    test_matrix_computation_int();
//    test_matrix_computation_double();
//    test_matrix_computation_complex();
//    pic_blend();

    /** 支持稀疏矩阵和普通矩阵互换*/
//    test_sparse_matrix();

    /** 支持一些基本归约操作：找到最大值，找到最小值，求和，算平均值等 */
//    test_detailed_matrix();

    /** 支持计算特征值和特征向量，求逆和行列式 */
//    test_cal();

    /** 支持改变形状和切割 */
//    test_SliceAndReshape();

    /** 支持矩阵卷积以及矩阵卷积在图像处理中的运用*/
   pic_process();
//    pic_process_gray();

    /** 支持 Opencv Mat 与 Matrix 类的转换 */
//    Mat img1 = imread("../img1.jpg"); // 3 channel BGR color image
//    Mat img2 = imread("../img1.jpg",COLOR_BGR2GRAY);
//    if (img1.empty() || img2.empty())                     // Check for invalid input
//    {
//        cout << "Could not open or find the image" << std::endl;
//        return -1;
//    }
//    test_img(img1, img2);

    /** 代码中包含异常处理 */
//    test_exception_handling();
}

void test_cal() {
    double arr1[9] = {1, 2, 3, 4, 5, 6, 7, 4, 5};
    matrix<double> m1(3, 3, arr1);
    matrix<double> eValue(3, 1);
    matrix<double> eVector(3, 3);
    m1.eigen(eValue, eVector);
    matrix<double> inv(3, 3);
    inv = m1.inverse();
    cout << m1;
    cout << "eigenvalues:\n" << eValue << endl;
    cout << "eigenvectors:\n" << eVector << endl;
    cout << "inverse matrix:\n" << inv << endl;
    cout << "determinant:\n" << m1.determinant() << endl;
    cout << "trace:\n" << m1.trace() << endl;
    cout << "average value:\n" << m1.average() << endl;

    int r, c;
    cout << "max value:\n" << m1.max_value(r, c) << "\nat (" << r << "," << c << ")" << endl;
}

void pic_process_gray() {
    // FYI: https://zhuanlan.zhihu.com/p/55013828
    // get picture
    Mat M = imread("../img1.jpg");
    cvtColor(M, M, COLOR_BGR2GRAY);
    matrix<int> mat(M);
    imshow("original", mat.to_Mat());

    // 模糊
    int arr1[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    matrix<int> core1(5, 5, arr1);
    matrix<int> pic1 = mat.convolution(core1, 25);
    pic1 = pic1.convolution(core1, 25);
    pic1 = pic1.convolution(core1, 25);
    imshow("box filter", pic1.to_Mat());

    // 边缘检测
    int arr2[9] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};
    matrix<int> core2(3, 3, arr2);
    matrix<int> pic2 = mat.convolution(core2);
    imshow("edge Detection", pic2.to_Mat());

    // 锐化
    int arr3[25] = {1, 1, 1, 1, -7, 1, 1, 1, 1};
    matrix<int> core3(5, 5, arr3);
    matrix<int> pic3 = mat.convolution(core3);
    imshow("Sharpness Filter", pic3.to_Mat());

    // 浮雕滤波器
    int arr4[25] = {-1, -1, -1, -1, 0, -1, -1, -1, 0, 1, -1, -1, 0, 1, 1, -1, 0, 1, 1, 1, 0, 1, 1, 1, 1};
    matrix<int> core4(5, 5, arr4);
    matrix<int> pic4 = mat.convolution(core4);
    imshow("Embossing Filter", pic4.to_Mat());

    waitKey();
}

void pic_process() {
    // FYI: https://zhuanlan.zhihu.com/p/55013828
    // get picture
    Mat M = imread("../img1.jpg");
    matrix<int> mat(M);
    imshow("original", mat.to_Mat());

    // 模糊
    int arr1[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    matrix<int> core1(5, 5, arr1);
    matrix<int> pic1 = mat.convolution(core1, 25);
    pic1 = pic1.convolution(core1, 25);
    pic1 = pic1.convolution(core1, 25);
    imshow("box filter", pic1.to_Mat());

    // 边缘检测
    int arr2[9] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};
    matrix<int> core2(3, 3, arr2);
    matrix<int> pic2 = mat.convolution(core2);
    imshow("edge Detection", pic2.to_Mat());

    // 锐化
    int arr3[25] = {1, 1, 1, 1, -7, 1, 1, 1, 1};
    matrix<int> core3(5, 5, arr3);
    matrix<int> pic3 = mat.convolution(core3);
    imshow("Sharpness Filter", pic3.to_Mat());

    // 浮雕滤波器
    int arr4[25] = {-1, -1, -1, -1, 0, -1, -1, -1, 0, 1, -1, -1, 0, 1, 1, -1, 0, 1, 1, 1, 0, 1, 1, 1, 1};
    matrix<int> core4(5, 5, arr4);
    matrix<int> pic4 = mat.convolution(core4);
    imshow("Embossing Filter", pic4.to_Mat());

    waitKey();
}

void test_img(Mat img1, Mat img2) {
    /** use opencv to get the image and convert it to our matrix to do computation
     *  after computation, convert matrix to Mat and use opencv to show it
     * */
    matrix<int> mat1(img1);
    matrix<int> mat2(img2);
    matrix<int> mat = mat1 + mat2;
    imshow("mat1+mat2", mat.to_Mat());
    mat = mat1;
    imshow("mat1", mat.to_Mat());
    mat = mat1 - mat2;
    imshow("mat1-mat2", mat.to_Mat());
    mat = 2 * mat1;
    imshow("2*mat1", mat.to_Mat());
    mat = mat1 / 2;
    imshow("mat1/2", mat.to_Mat());
    waitKey();
}

void test_matrix_computation_int() {
    int arr1[4] = {1, 2, 1, 2};
    int arr2[2] = {1, 2};
    matrix<int> m1(2, 2, arr1);
    matrix<int> m2(2, 2, arr1);
    matrix<int> m3(2, 1, arr2);
    cout << "m1 = \n" << m1;
    cout << "m2 = \n" << m2;
    cout << "m3 = \n" << m3;
    matrix<int> m4 = m1 + m2;
    matrix<int> m5 = m4 - m2;
    matrix<int> m6 = m1 * m3;
    cout << "m4 = m1+m2 = \n" << m4;
    cout << "m5 = m4-m2 = \n" << m5;
    cout << "m6 = m1*m3 = \n" << m6;
    matrix<int> m7 = 2 * m3;
    matrix<int> m8 = m7 / 2;
    cout << "m7 = 2*m3 = \n" << m7;
    cout << "m8 = m7/2 = \n" << m8;
    matrix<int> m9 = m7.tran_matrix();
    cout << "m9 = m7's trans = \n" << m9;
    matrix<int> m10 = m7.conjugation_matrix();
    cout << "m10 = m7's conju = \n" << m10;
    matrix<int> m11 = m9.elem_wise_mul(m10);
    cout << "m11 = m9 elem_wise_mul m10 = \n" << m11;
    vector<int> v;
    v.push_back(1);
    v.push_back(2);
    matrix<int> m12 = m11.vector_mul(v);
    cout << "m12 = vector(1,2) * m11 = \n" << m12 << endl;
    matrix<int> m13 = m8.mul_vector(v);
    cout << "m13 = m0 * vector(1,2) = \n" << m13 << endl;
    int ar1[3] = {1, 2, 3};
    int ar2[3] = {1, -1, 2};
    matrix<int> mat1(3, 1, ar1);
    matrix<int> mat2(3, 1, ar2);
    int result = mat1.dot_product(mat2);
    matrix<int> mat3 = mat1.cross_product(mat2);
    cout << "mat1 = \n" << mat1 << endl;
    cout << "mat2 = \n" << mat2 << endl;
    cout << "result = mat1 dot product mat2 = " << result;
    cout << "mat3 = mat1 cross product mat2 = \n" << mat3;
};

void test_matrix_computation_double() {
    double arr1[4] = {1.1, 2.3, 1.2, 2.2};
    double arr2[2] = {1.2, 2.4};
    matrix<double> m1(2, 2, arr1);
    matrix<double> m2(2, 2, arr1);
    matrix<double> m3(2, 1, arr2);
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

void test_matrix_computation_complex() {
    complex<int> ci1(1, 1);
    complex<int> ci2 = conj(ci1);
    complex<int> t[4] = {ci1, ci2, ci1, ci2};
    matrix<complex<int>> mm1(2, 2, t);
    matrix<complex<int>> mm2 = mm1.conjugation_matrix();
    cout << mm1 << endl;
    cout << mm2 << endl;
    complex<double> c1(1, 1);
    complex<double> c2(2, 1);
    complex<double> c3 = c1 + c2;
    complex<double> c4 = c2 - c1;
    complex<double> c5 = c1 * c2;
    cout << "c1 = " << c1 << endl;
    cout << "c2 = " << c2 << endl;
    cout << "c3 = c1+c2 = " << c3 << endl;
    cout << "c4 = c2-c1 = " << c4 << endl;
    cout << "c5 = c1*c2 = " << c5 << endl;
    complex<double> arr1[4] = {c1, c2, c3, c4};
    matrix<complex<double>> m1(2, 2, arr1);
    matrix<complex<double>> m2(2, 2, arr1);
    matrix<complex<double>> m3 = m1 + m2;
    matrix<complex<double>> m4 = m1 * m2;
    matrix<complex<double>> m5 = 2 * m2;
    matrix<complex<double>> m6 = m2 * 4;
    matrix<complex<double>> m0 = m6 / 2;

    cout << "m1 = \n" << m1 << endl;
    cout << "m2 = \n" << m2 << endl;
    cout << "m3 = m1+m2 = \n" << m3 << endl;
    cout << "m4 = m1*m2 = \n" << m4 << endl;
    cout << "m5 = 2 *m2 = \n" << m5 << endl;
    cout << "m6 = m2* 4 = \n" << m6 << endl;
    cout << "m0 = m6/ 2 = \n" << m0 << endl;
    matrix<complex<double>> m7 = m6.tran_matrix();
    cout << "m7 = m6'trans = \n" << m7 << endl;
    matrix<complex<double>> m8 = m6.conjugation_matrix();
    cout << "m8 = m6'conju = \n" << m8 << endl;
    matrix<complex<double>> m9 = m8.elem_wise_mul(m6);
    cout << "m9 = m8 elem_wise_mul m6 = \n" << m9 << endl;
    vector<complex<double>> v;
    v.push_back(1);
    v.push_back(2);
    matrix<complex<double>> ma1 = m8.vector_mul(v);
    cout << "ma1 = vector(1,2) * ma1 = \n" << ma1 << endl;
    complex<double> a1[2] = {c1, c2};
    matrix<complex<double>> ma2(2, 1, a1);
    matrix<complex<double>> ma3 = ma2.mul_vector(v);
    cout << "ma2 = \n" << ma2 << endl;
    cout << "ma3 = ma1 * vector(1,2) = \n" << ma3 << endl;
    complex<double> ar1[3] = {complex<double>(1, 2), complex<double>(1, 3), complex<double>(1, 4)};
    complex<double> ar2[3] = {complex<double>(1, -1), complex<double>(1, 0), complex<double>(1, 1)};
    matrix<complex<double>> mat1(3, 1, ar1);
    matrix<complex<double>> mat2(3, 1, ar2);
    cout << "mat1 = \n" << mat1 << endl;
    cout << "mat2 = \n" << mat2 << endl;
    complex<double> result = mat1.dot_product(mat2);
    cout << "result = mat1 dot product mat2 = " << result << endl;
    matrix<complex<double>> mat3 = mat1.cross_product(mat2);
    cout << "mat3 = mat1 cross product mat2 = \n" << mat3 << endl;
}

void pic_blend() {
    Mat img1 = imread("../img1.jpg");
    Mat img2 = imread("../img2.jpg");
    matrix<int> mat1(img1);
    matrix<int> mat2(img2);
    matrix<int> mat_1 = mat1;
    matrix<int> mat_2 = mat2;
    matrix<int> mat = mat_1 + mat_2;
    imshow("picture blend", mat.to_Mat());
    waitKey(0);
}

void test_sparse_matrix() {
    int arr1[8] = {2, 0, 0, 0, 5, 4, 1, 0};
    int arr2[8] = {0, 1, 0, 0, 0, 0, 0, 0};
    int arr3[8]{0, 0, 0, 14, 0, 0, 0, 0};
    int arr4[8]{0, 0, 5, 0, 0, 0, 0, 0};
    int arr5[8]{0, 0, 0, 0, 5, 0, 0, 0};
    int arr6[8]{0, 12, 0, 0, 0, 0, 0, 0};
    int arr7[8]{0, 0, 2, 1, 0, 1, 0, 0};
    int arr8[8] = {0, 0, 0, 0, 0, 0, 4, 0};
    int *arr0[8] = {arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8};
    int arra[1] = {3};
    int arrb[1] = {4};
    int arrc[1] = {1};
    int arrd[1] = {7};
    int *arrz[4] = {arra, arrb, arrc, arrd};
    matrix<int> m1(8, 8, arr0);
    matrix<int> ma(4, 1, arrz);
    matrix<int> m4 = m1.BlockMatrix(4, 4, 2);

    cout << m4 << endl;
    cout << ma << endl;
    spmatrix<int> m5 = compress_matrix(m4);

    matrix<int> m6 = m5.sp_vector_mul(ma);
    matrix<int> m7 = m5.recover_matrix();

    cout << "Sparse matrix information:" << endl;
    cout << m5 << endl;
    cout << "The result of vector multiply in spmatrix: " << endl;
    cout << m6 << endl;
    cout << "Recover spmatrix into standard matrix:" << endl;
    cout << m7 << endl;
}

void test_SliceAndReshape() {
    int arr1[8] = {2, 0, 0, 0, 5, 4, 1, 0};
    int arr2[8] = {0, 1, 0, 0, 0, 0, 0, 0};
    int arr3[8]{0, 0, 0, 14, 0, 0, 0, 0};
    int arr4[8]{0, 0, 5, 0, 0, 0, 0, 0};
    int arr5[8]{0, 0, 0, 0, 5, 0, 0, 0};
    int arr6[8]{0, 12, 0, 0, 0, 0, 0, 0};
    int arr7[8]{0, 0, 2, 1, 0, 1, 0, 0};
    int arr8[8] = {0, 0, 0, 0, 0, 0, 4, 0};
    int *arr0[8] = {arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8};
    matrix<int> m1(8, 8, arr0);
    matrix<int> m2 = m1.ReshapeMatrix(12, 12, 2, 3);
    matrix<int> m5 = m1.ReshapeMatrix(4, 16);
    matrix<int> m3 = m1.SilceMatrix(3, 4, 7, 7);
    matrix<int> m4 = m1.BlockMatrix(4, 4, 2);

    cout << "The origin matrix:" << endl;
    cout << m1 << endl;
    cout << "Reshaping the origin matrix:" << endl;
    cout << m2 << endl;
    cout << m5 << endl;
    cout << "Slice the origin matrix:" << endl;
    cout << m3 << endl;
    cout << "Choosing the thrid block of the matrix:(left bottom):" << endl;
    cout << m4 << endl;
}

void test_exception_handling() {
    int arr1[8] = {2, 0, 0, 0, 5, 4, 1, 0};
    int arr2[8] = {0, 1, 0, 0, 0, 0, 0, 0};
    int arr3[8]{0, 0, 0, 14, 0, 0, 0, 0};
    int arr4[8]{0, 0, 5, 0, 0, 0, 0, 0};
    int arr5[8]{0, 0, 0, 0, 5, 0, 0, 0};
    int arr6[8]{0, 12, 0, 0, 0, 0, 0, 0};
    int arr7[8]{0, 0, 2, 1, 0, 1, 0, 0};
    int arr8[8] = {0, 0, 0, 0, 0, 0, 4, 0};
    int *arr0[8] = {arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8};
    matrix<int> m1(8, 8, arr0);
    cout << m1 << endl;
    cout << "Reshape matrix into 4*15 exception out.";
    matrix<int> m2 = m1.ReshapeMatrix(4, 15);
}
