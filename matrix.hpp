#pragma once

#include <iostream>
#include <vector>

template<typename T = double>
class Matrix {
public:
	static Matrix<T> IdentityMatrix(std::size_t);
	static Matrix<T> RotationMatrix(std::size_t, std::size_t, std::size_t, const T&);
public:
	Matrix();
	Matrix(std::size_t, std::size_t, const T& = T{});
	Matrix(const std::initializer_list<std::initializer_list<T>>&);
	Matrix(const Matrix<T>&);
	Matrix<T>& operator=(const Matrix<T>&);
public:
	std::size_t rowsCount() const;
	std::size_t colsCount() const;
public:
	std::vector<T>& operator[](std::size_t);
	const std::vector<T>& operator[](std::size_t) const;
public:
	std::vector<T> getRow(std::size_t) const;
	std::vector<T> getCol(std::size_t) const;

	void swapRows(std::size_t, std::size_t);
	void swapCols(std::size_t, std::size_t);

	void appendCol(const std::vector<T>&);
	void appendRow(const std::vector<T>&);
public:
	void appendMatrixToRight(const Matrix<T>&);

	Matrix<T> getSubMatrix(std::size_t, std::size_t, std::size_t, std::size_t) const;
public:
	T getDet() const;
	T getTrace() const;
	T getNorm() const;
public:
	Matrix<T> getUpTrapezoidal() const;
	Matrix<T> getDownTrapezoidal() const;
	std::pair<Matrix<T>, Matrix<T>> getLUDecomposition() const;
	Matrix<T> getCholeskyDecomposition() const;
	Matrix<T> getInversedMatrix() const;
	std::pair<Matrix<T>, Matrix<T>> getQRDecomposition() const;
public:
	bool isSquare() const;
	bool isSymmetric() const;
	bool isUpperTriangular() const;
	bool isLowerTriangular() const;
	bool isDiagonal() const;
	bool isOrthognal() const;
	bool isHermitian() const;
	bool isConjugate() const;
	bool isNormal() const;
public:
	Matrix<T>& operator+=(const Matrix<T>&);
	Matrix<T>& operator-=(const Matrix<T>&);
	Matrix<T>& operator*=(const Matrix<T>&);
	Matrix<T>& operator*=(const T&);
	Matrix<T>& operator/=(const T&);
public:
	template<typename V>
	friend bool operator==(const Matrix<V>&, const Matrix<V>&);
	template<typename V>
	friend bool operator!=(const Matrix<V>&, const Matrix<V>&);
public:
	template<typename V>
	friend Matrix<V> operator+(const Matrix<V>&, const Matrix<V>&);

	template<typename V>
	friend Matrix<V> operator-(const Matrix<V>&, const Matrix<V>&);

	template<typename V>
	friend Matrix<V> operator*(const Matrix<V>&, const Matrix<V>&);

	template<typename V>
	friend std::vector<V> operator*(const Matrix<V>&, const std::vector<V>&);

	template<typename V>
	friend std::vector<V> operator*(const std::vector<V>&, const Matrix<V>&);

	template<typename V>
	friend Matrix<V> operator*(const Matrix<V>&, const V&);
	template<typename V>
	friend Matrix<V> operator*(const V&, const Matrix<V>&);
	template<typename V>
	friend Matrix<V> operator/(const Matrix<V>&, const V&);
public:
	template<typename V>
	friend std::istream& operator>>(std::istream&, Matrix<V>&);
	template<typename V>
	friend std::ostream& operator<<(std::ostream&, const Matrix<V>&);
private:
	std::size_t m, n;
	std::vector<std::vector<T>> data;
};


















#include "matrix.cpp"
