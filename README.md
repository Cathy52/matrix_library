# matrix_library


done:

- [x] 矩阵创建：支持用一维数组，二维数组，或者 Mat 对象（默认为：8UC3）创建矩阵
- [x] 矩阵转换：支持 matrix 转换为 opencv Mat对象（8UC3）。（颜色空间：RGB，数据类型：unsigned char，3 channel BGR color image）
- [x] matrix 的基本计算。
- [x] matrix 支持 int, double, complex。

unsolved:

- [ ] 要求4，5，6，7，9
- [ ] 稀疏矩阵
- [ ] matrix 计算方式与 opencv 不一致。比如：opencv 中加法超过255就直接为255，两像素值差为负，则直接为0

