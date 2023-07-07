#include "matrix.hpp"

int main(){
    Matrix<double> A = {
        {1, 1, 1},
        {-1, 1, 1},
        {1, 1, -1}
    };
    
    // Identity matrix
    std::cout <<"Identity matrix of size 3" << std::endl;
    Matrix<double> identity = Matrix<double>::IdentityMatrix(3);
    identity.print();
    
    // Rotation matrix !! Parametr phi select in radians !!
    std::cout <<"Rotation matrix of size 2" << std::endl;
    Matrix<double> rotation = Matrix<double>::RotationMatrix(3, 0, 1, 0.52);
    rotation.print();
    
    // rowsCount and colsCount
    std::size_t size_m = A.rowsCount();
    std::cout << "Number of rows: " << size_m << std::endl;
    std::size_t size_n = A.colsCount();
    std::cout << "Number of cols: " << size_n << std::endl;
    
    // Print Matrix
    std::cout << "Matrix A: \n";
    A.print();
    
    // getRow
    std::size_t row_index = 1; 
    std::vector<double> row = A.getRow(row_index);
    std::cout << "Row " << row_index << ": ";
    for (const auto& element : row) {
        std::cout << element << " ";
    }
    std::cout << std::endl;

    // getCol
    std::size_t col_index = 2; 
    std::vector<double> col = A.getCol(col_index);
    std::cout << "Column " << col_index << ": ";
    for (const auto& element : col) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    
    /*
    // swapRows
    std::cout << "Matrix A after rows swapping: \n";
    A.swapRows(0, 2);
    A.print();

    // swapCols
    std::cout << "Matrix A after cols swapping: \n";
    A.swapCols(0, 1);
    A.print();
    
    // appendCol
    std::cout << "Matrix A after append new col: \n";
    std::vector<double> col1 = {10, 11, 12};
    A.appendCol(col1);
    A.print();

    // appendRow
    std::cout << "Matrix A after append new row: \n";
    std::vector<double> row1 = {13, 14, 15, 16};
    A.appendRow(row1);
    A.print();
    
    //getSubMatrix
    std::cout << "SubMatrix A: \n";
    Matrix<double> SubMatrixA = A.getSubMatrix(0, 2, 0, 1);
    SubMatrixA.print();
    
    // getTrace
    double trace = A.getTrace();
    std::cout << "Matrix A trace: " << trace << std::endl;
    
    // getDet
    double determinant = A.getDet();
    std::cout << "Matrix A determinant: " << determinant << std::endl;

    // getNorm
    double norm = A.getNorm();
    std::cout << "Matrix A norm: " << norm << std::endl;
    
    std::cout << "Matrix B: \n";
    Matrix<double> B = {
        {17, 18},
        {19, 20},
        {21, 22}
    };
    
    B.print();
    std::cout << "Matrix A after appending new matrix B: \n";
    A.appendMatrixToRight(B);
    A.print();
    */
    // getUpTrapezoidal
    std::cout << "Up Trapezoidal matrix A: \n";
    Matrix<double> UpTrapezoidal = A.getUpTrapezoidal();
    UpTrapezoidal.print();

    // getDownTrapezoidal
    std::cout << "Down Trapezoidal matrix A: \n";
    Matrix<double> DownTrapezoidal = A.getDownTrapezoidal();
    DownTrapezoidal.print();
    
    // isSquare
    std::cout << "Matrix A is square? \n";
    bool check_square = A.isSquare();
    std::cout << check_square << std::endl;
    
    // isSymmetric
    std::cout << "Matrix A is symmetric? \n";
    bool check_symmetric = A.isSymmetric();
    std::cout << check_symmetric << std::endl;
    
    // isUpperTriangular
    std::cout << "Matrix A is Upper Triangular? \n";
    bool check_upper = A.isUpperTriangular();
    std::cout << check_upper << std::endl;

    // isLowerTriangular
    std::cout << "Matrix A is Lower Triangular? \n";
    bool check_lower = A.isLowerTriangular();
    std::cout << check_lower << std::endl;
    
    // getInversedMatrix
    std::cout << "Inverse Matrix A: \n";
    Matrix<double> InversedMatrix = A.getInversedMatrix();
    InversedMatrix.print();
    
    // isOrthogonal
    std::cout << "Matrix A is Orthogonal? \n";
    bool check_orthogonal = A.isOrthogonal();
    std::cout << check_orthogonal << std::endl;
    
    // getLUDecomposition
    std::pair<Matrix<double>, Matrix<double>> LU_decomposition = A.getLUDecomposition();
    Matrix<double> L = LU_decomposition.first;
    Matrix<double> U = LU_decomposition.second;

    std::cout << "Lower Matrix (L):" << std::endl;
    L.print();

    std::cout << "Upper Matrix (U):" << std::endl;
    U.print();
    
    // getCholeskyDecomposition
    std::cout << "Cholesky decomposition A: \n";
    Matrix<double> cholesky = A.getCholeskyDecomposition();
    cholesky.print();
    
    // getQRDecomposition using Gram_Schmidt method
    std::cout << "QR decomposition of Matrix A using Gram-Schmidt orthogonalization process: \n";
    std::pair<Matrix<double>, Matrix<double>> QR_decomposition_Gram_Schmidt = A.getQRDecomposition_Gram_Schmidt();
    Matrix<double> Q = QR_decomposition_Gram_Schmidt.first;
    Matrix<double> R = QR_decomposition_Gram_Schmidt.second;

    std::cout << "Orthogonal Matrix (Q):" << std::endl;
    Q.print();
    std::cout << "Upper Triangular Matrix (R):" << std::endl;
    R.print();

    // getQRDecomposition using Givens method (rotation method)
    std::cout << "QR decomposition of Matrix A using Givens method (rotation method): \n";
    std::pair<Matrix<double>, Matrix<double>> QR_decomposition_rotation = A.getQRDecomposition_rotation();
    Matrix<double> Q1 = QR_decomposition_rotation.first;
    Matrix<double> R1 = QR_decomposition_rotation.second;

    std::cout << "Orthogonal Matrix (Q):" << std::endl;
    Q1.print();
    std::cout << "Upper Triangular Matrix (R):" << std::endl;
    R1.print();

    // getQRDecomposition using Housholder method (reflection method)
    std::cout << "QR decomposition of Matrix A using Housholder method (reflection method): \n";
    std::pair<Matrix<double>, Matrix<double>> QR_decomposition_reflection = A.getQRDecomposition_reflection();
    Matrix<double> Q2 = QR_decomposition_reflection.first;
    Matrix<double> R2 = QR_decomposition_reflection.second;

    std::cout << "Orthogonal Matrix (Q):" << std::endl;
    Q2.print();
    std::cout << "Upper Triangular Matrix (R):" << std::endl;
    R2.print();
    
    // Complex matrix
    Matrix<std::complex<double>> Complex_Matrix(2, 2);
    Complex_Matrix[0][0] = std::complex<double>(1.0, 0.0);
    Complex_Matrix[0][1] = std::complex<double>(3, -2);
    Complex_Matrix[1][0] = std::complex<double>(3, 2);
    Complex_Matrix[1][1] = std::complex<double>(2.0, 0);
    
    std::cout << "Complex_Matrix: \n";
    Complex_Matrix.print();

    // Conjugate Matrix
    std::cout << "Conjugate Complex_Matrix: \n";
    Matrix<std::complex<double>> conjugate = Complex_Matrix.getConjugate();
    conjugate.print();
    
    std::cout << "Matrix A is Conjugate? \n";
    bool check_conjugate = Complex_Matrix.isConjugate();
    std::cout << check_conjugate << std::endl;
    
    // Hermitian Matrix
    std::cout << "Matrix A is Hermitian? \n";
    bool check_hermitian = Complex_Matrix.isHermitian();
    std::cout << check_hermitian << std::endl;
    
    // Normal Matrix
    std::cout << "Matrix A is Normal? \n";
    bool check_normal = Complex_Matrix.isNormal();
    std::cout << check_normal << std::endl;

    return 0;
}