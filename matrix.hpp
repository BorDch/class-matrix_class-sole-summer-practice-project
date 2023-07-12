#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <ctime>
#include <random>

template<typename T = double>
class Matrix {
public:
	// class Matrix Constructors 
	Matrix();
	Matrix(std::size_t, std::size_t, const T& = T{});
	Matrix(const std::initializer_list<std::initializer_list<T>>&);
	Matrix(const Matrix<T>&);

	// STATIC methods for matrix transformation
	static Matrix<T> IdentityMatrix(std::size_t);
	static Matrix<T> RotationMatrix(std::size_t, std::size_t, std::size_t, const T&);

	// Main methods for working with matrices
	void print() const;
	std::size_t rowsCount() const;
	std::size_t colsCount() const;
	std::vector<T> getRow(std::size_t) const;
	std::vector<T> getCol(std::size_t) const;

	void swapRows(std::size_t, std::size_t);
	void swapCols(std::size_t, std::size_t);

	void appendCol(const std::vector<T>&);
	void appendRow(const std::vector<T>&);

	void appendMatrixToRight(const Matrix<T>&);

	Matrix<T> getSubMatrix(std::size_t, std::size_t, std::size_t, std::size_t) const;

	T getDet() const;
	T getTrace() const;
	T getNorm() const;

	Matrix<T> getTransposed() const;
	Matrix<T> getUpTrapezoidal() const;
	Matrix<T> getDownTrapezoidal() const;
	Matrix<T> getConjugate() const;
	std::pair<Matrix<T>, Matrix<T>> getLUDecomposition() const;
	Matrix<T> getCholeskyDecomposition() const;
	Matrix<T> getInversedMatrix() const;
	std::pair<Matrix<T>, Matrix<T>> getQRDecomposition_Gram_Schmidt() const;
	std::pair<Matrix<T>, Matrix<T>> getQRDecomposition_rotation() const;
	std::pair<Matrix<T>, Matrix<T>> getQRDecomposition_reflection() const;

	// Eigen values (Non-const method)
	std::vector<T> findEigenvalues(const Matrix<T>&, std::size_t);
	
	// Eigen values (Const method) Extra modificaton
	std::vector<T> QR_findEigenvalues(std::size_t numIterations) const {
    	Matrix<T> Ak(*this);
        std::size_t n = Ak.rowsCount();

        for (std::size_t iteration = 0; iteration < numIterations; ++iteration) {
            std::pair<Matrix<T>, Matrix<T>> QR = Ak.getQRDecomposition_reflection();
            Matrix<T> Q = QR.first;
            Matrix<T> R = QR.second;
            Ak = R * Q;
        }

        std::vector<T> eigenvalues(n);

        for (std::size_t i = 0; i < n; ++i) {
            eigenvalues[i] = Ak[i][i];
        }

        return eigenvalues;
    } 

	// Random Matrix Generation
	static Matrix<T> FLOAT_RandomMatrix(std::size_t, std::size_t);
	static Matrix<T> INT_RandomMatrix(std::size_t, std::size_t);

	// Matrix check methods
	bool isSquare() const;
	bool isSymmetric() const;
	bool isUpperTriangular() const;
	bool isLowerTriangular() const;
	bool isDiagonal() const;
	bool isOrthogonal() const;
	bool isHermitian() const;
	bool isConjugate() const;
	bool isNormal() const;

	// Overloading operators and several friendly functions to working with matrices
	Matrix<T>& operator=(const Matrix<T>&);
	std::vector<T>& operator[](std::size_t);
	const std::vector<T>& operator[](std::size_t) const;

	// Operator *
	template<typename V>
	friend std::vector<V> operator*(const Matrix<V>&, const std::vector<V>&);

	template<typename V>
	friend std::vector<V> operator*(const std::vector<V>&, const Matrix<V>&);

	template<typename V>
	friend Matrix<V> operator*(const Matrix<V>&, const Matrix<V>&);

	template<typename V>
	friend Matrix<V> operator*(const Matrix<V>&, const V&);
	
	template<typename V>
	friend Matrix<V> operator*(const V&, const Matrix<V>&);

	// Operator /
	template<typename V>
	friend Matrix<V> operator/(const Matrix<V>&, const V&);

	// Operator +
	template<typename V>
	friend Matrix<V> operator+(const Matrix<V>&, const Matrix<V>&);

	// Operator - 
	template<typename V>
	friend Matrix<V> operator-(const Matrix<V>&, const Matrix<V>&);

	// Other arithmetic operators
	Matrix<T>& operator+=(const Matrix<T>&);
	Matrix<T>& operator-=(const Matrix<T>&);
	Matrix<T>& operator*=(const Matrix<T>&);
	Matrix<T>& operator*=(const T&);
	Matrix<T>& operator/=(const T&);

	// Boolean Comparison Operators
	template<typename V>
	friend bool operator==(const Matrix<V>&, const Matrix<V>&);
	
	template<typename V>
	friend bool operator!=(const Matrix<V>&, const Matrix<V>&);

	// Friendly functions to Input/Output
	template<typename V>
	friend std::istream& operator>>(std::istream&, Matrix<V>&);
	
	template<typename V>
	friend std::ostream& operator<<(std::ostream&, const Matrix<V>&);

private:
	std::size_t m, n;
	std::vector<std::vector<T>> data;
};

#include "matrix.cpp"
