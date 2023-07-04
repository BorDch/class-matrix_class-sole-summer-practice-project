#include <cstdlib>
#include <vector>
#include <cmath>

// Constructors
template<typename T>
Matrix<T>::Matrix() : m(0), n(0) {}

template<typename T>
Matrix<T>::Matrix(std::size_t m, std::size_t n, const T& value) : m(m), n(n) {
	data = std::vector<std::vector<T>>(m);
	for (std::size_t i{}; i < m; i++) {
		data[i] = std::vector<T>(n, value);
	}
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& other) {
	m = other.size();
	data = std::vector<std::vector<T>>();
	for (const auto& row : other) {
		data.push_back(row);
	}
	n = data[0].size();
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& other) : m(other.m), n(other.n), data(other.data) {}

// IdentityMatrix
template<typename T>
Matrix<T> Matrix<T>::IdentityMatrix(std::size_t size) {
	Matrix<T> identity(size, size, 0);
	for (std::size_t i = 0; i < size; i++) {
		identity[i][i] = 1;
	}
	return identity;
}

// RotationMatrix
template<typename T>
Matrix<T> Matrix<T>::RotationMatrix(std::size_t n, std::size_t i, std::size_t j, const T& phi) {
	Matrix<T> res = IdentityMatrix(n);
	res.data[i][i] = std::cos(phi);
	res.data[i][j] = -std::sin(phi);
	res.data[j][i] = std::sin(phi);
	res.data[j][j] = std::cos(phi);
	return res;
}

// MAIN METHODS

// Print Matrix
template<typename T>
void Matrix<T>::print() const {
    for (const auto& row : data) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Get size
template<typename T>
std::size_t Matrix<T>::rowsCount() const {
	return m;
}

template<typename T>
std::size_t Matrix<T>::colsCount() const {
	return n;
}

// Get Row
template<typename T>
std::vector<T> Matrix<T>::getRow(std::size_t row) const {
	if (row >= m) {
		throw std::runtime_error("Invalid row index");
	}
	return data[row];
}

// Get Col
template<typename T>
std::vector<T> Matrix<T>::getCol(std::size_t col) const {
	if (col >= n) {
		throw std::runtime_error("Invalid column index");
	}
	std::vector<T> column(m);
	
	for (std::size_t i = 0; i < m; i++) {
		column[i] = data[i][col];
	}
	return column;
}

// swapRows
template<typename T>
void Matrix<T>::swapRows(std::size_t row1, std::size_t row2) {
	if (row1 >= m || row2 >= m) {
		throw std::runtime_error("Invalid row index");
	}

	for (std::size_t j = 0; j < n; j++) {
		T temp = data[row1][j];
		data[row1][j] = data[row2][j];
		data[row2][j] = temp;
	}
}

// swapCols
template<typename T>
void Matrix<T>::swapCols(std::size_t col1, std::size_t col2) {
	if (col1 >= n || col2 >= n) {
		throw std::runtime_error("Invalid column index");
	}

	for (std::size_t i = 0; i < m; i++) {
		T temp = data[i][col1];
		data[i][col1] = data[i][col2];
		data[i][col2] = temp;
	}
}

// appendCol
template<typename T>
void Matrix<T>::appendCol(const std::vector<T>& col) {
	if (col.size() != m) {
		throw std::runtime_error("Invalid size of cols");
	}

	for (std::size_t i = 0; i < m; i++) {
		data[i].push_back(col[i]);
	}

	n++;
}

// appendRow
template<typename T>
void Matrix<T>::appendRow(const std::vector<T>& row) {
	if (row.size() != n) {
		throw std::runtime_error("Invalid size of rows");
	}

	data.push_back(row);
	m++;
}

// appendMatrixToRight
template<typename T>
void Matrix<T>::appendMatrixToRight(const Matrix<T>& matrix) {
	if (matrix.m != m) {
		throw std::runtime_error("Invalid size of matrices");
	}

	for (std::size_t i = 0; i < m; i++) {
		data[i].insert(data[i].end(), matrix.data[i].begin(), matrix.data[i].end());
	}

	n += matrix.n;
}

// getSubMatrix
template<typename T>
Matrix<T> Matrix<T>::getSubMatrix(std::size_t startRow, std::size_t endRow, std::size_t startCol, std::size_t endCol) const {
    std::size_t numRows = endRow - startRow + 1;
    std::size_t numCols = endCol - startCol + 1;

    if (endRow >= m || endCol >= n) {
        throw std::runtime_error("Invalid submatrix indices");
    }

    Matrix<T> subMatrix(numRows, numCols);

    for (std::size_t i = startRow; i <= endRow; i++) {
        for (std::size_t j = startCol; j <= endCol; j++) {
            subMatrix[i - startRow][j - startCol] = data[i][j];
        }
    }

    return subMatrix;
}

// getDet using the Gauss Method
template<typename T>
T Matrix<T>::getDet() const {
    if (m != n) {
        throw std::runtime_error("Determinant is only defined for square matrices");
    }

    Matrix<T> Upper_Triangular = *this;
    T determinant = 1;

    for (std::size_t i = 0; i < n; i++) {
        // Find max row
        std::size_t max_Row = i;
        T max_Val = std::abs(Upper_Triangular[i][i]);

        for (std::size_t j = i + 1; j < n; j++) {
            T abs_Val = std::abs(Upper_Triangular[j][i]);
            if (abs_Val > max_Val) {
                max_Row = j;
                max_Val = abs_Val;
            }
        }

        // Swap rows
        if (max_Row != i) {
            Upper_Triangular.swapRows(i, max_Row);
            determinant *= -1; 
        }

        T pivot = Upper_Triangular[i][i];
        if (std::abs(pivot) < 1e-10) {
            return 0;
        }

        determinant *= pivot;

        for (std::size_t j = i + 1; j < n; j++) {
            T temp = Upper_Triangular[j][i] / pivot;
            for (std::size_t k = i; k < n; k++) {
                Upper_Triangular[j][k] -= temp * Upper_Triangular[i][k];
            }
        }
    }

    return determinant;
}

// getTrace
template<typename T>
T Matrix<T>::getTrace() const {
	if (m != n) {
		throw std::runtime_error("Trace is only defined for square matrices");
	}

	T trace = 0;
	for (std::size_t i = 0; i < m; i++) {
		trace += data[i][i];
	}

	return trace;
}

// getNorm
template<typename T>
T Matrix<T>::getNorm() const {
	T norm = 0;

	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < n; j++) {
			norm += std::pow(std::abs(data[i][j]), 2);
		}
	}

	return std::sqrt(norm);
}

// getTransposed
template<typename T>
Matrix<T> Matrix<T>::getTransposed() const {
	Matrix res(n, m);
	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < n; j++) {
			res.data[j][i] = data[i][j];
		}
	}
	return res;
}

// getUpTrapezoidal
template<typename T>
Matrix<T> Matrix<T>::getUpTrapezoidal() const {
    Matrix<T> result(*this);

	for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < i; j++) {
            result.data[i][j] = 0;
        }
    }
    
    return result;
}

// getDownTrapezoidal
template<typename T>
Matrix<T> Matrix<T>::getDownTrapezoidal() const {
    Matrix<T> result(*this);

    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = i + 1; j < n; j++) {
            result.data[i][j] = 0;
        }
    }

    return result;
}

// getConjugate
template<typename T>
Matrix<T> Matrix<T>::getConjugate() const {
    Matrix<T> conjugate(m, n);
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            if (std::imag(data[i][j]) == 0) {
                conjugate[i][j] = std::conj(std::real(data[i][j]));
            } else {
                conjugate[i][j] = std::conj(data[i][j]);
            }
        }
    }
    return conjugate;
}

// getLUDecomposition
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::getLUDecomposition() const {
    if (!isSquare()) {
        throw std::runtime_error("Invalid size of matrix");
    }

    Matrix<T> L = IdentityMatrix(m);
    Matrix<T> U = *this;

    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = i + 1; j < m; j++) {
            T temp = U[j][i] / U[i][i];
            L[j][i] = temp;
            for (std::size_t k = i; k < m; k++) {
                U[j][k] -= temp * U[i][k];
            }
        }
    }

    return std::make_pair(L, U);
}

// getCholeskyDecomposition
template<typename T>
Matrix<T> Matrix<T>::getCholeskyDecomposition() const {
    if (!isSquare()) {
        throw std::runtime_error("Matrix is not square");
    }

    Matrix<T> L(m, n);

    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j <= i; j++) {
            T sum = 0;

            if (i == j) {
                for (std::size_t k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                L[j][j] = std::sqrt(data[i][j] - sum);
            } else {
                for (std::size_t k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (data[i][j] - sum) / L[j][j];
            }
        }
    }

    return L;
}

// getInversedMatrix
template<typename T>
Matrix<T> Matrix<T>::getInversedMatrix() const {
	if (!isSquare()) {
		throw std::runtime_error("Matrix is not square");
	}

	Matrix<T> augmentedMatrix(*this);
	Matrix<T> identityMatrix = IdentityMatrix(m);
	augmentedMatrix.appendMatrixToRight(identityMatrix);

	std::size_t num_Rows = augmentedMatrix.rowsCount();
	std::size_t num_Cols = augmentedMatrix.colsCount();

	for (std::size_t i = 0; i < num_Rows; i++) {
		std::size_t pivot_Row = i;
		T pivot = augmentedMatrix[i][i];

		for (std::size_t j = i + 1; j < num_Rows; j++) {
			if (std::abs(augmentedMatrix[j][i]) > std::abs(pivot)) {
				pivot_Row = j;
				pivot = augmentedMatrix[j][i];
			}
		}

		if (pivot_Row != i) {
			augmentedMatrix.swapRows(i, pivot_Row);
		}

		for (std::size_t j = i + 1; j < num_Cols; j++) {
			augmentedMatrix[i][j] /= pivot;
		}

		for (std::size_t j = 0; j < num_Rows; j++) {
			if (j != i) {
				T temp = augmentedMatrix[j][i];
				for (std::size_t k = i + 1; k < num_Cols; k++) {
					augmentedMatrix[j][k] -= temp * augmentedMatrix[i][k];
				}
			}
		}
	}

	Matrix<T> inversedMatrix(m, m);

	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < m; j++) {
			inversedMatrix[i][j] = augmentedMatrix[i][j + m];
		}
	}

	return inversedMatrix;
}

// getQRDecomposition
/* First realization using rotation method
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::getQRDecomposition() const {
	Matrix R(*this);
	Matrix Q = IdentityMatrix(n);
	for (std::size_t j = 0; j < n - 1; j++) {
		for (std::size_t i = j + 1; i < n; i++) {
			Matrix V = RotationMatrix(n, i, j, std::atan2(R.data[i][j], R.data[j][j]));
			R = V * R;
			Q = Q * V;
			// Extra modification: Check the signs of elements in R and Q
            for (std::size_t k = 0; k < n; k++) {
                if (R.data[k][k] < 0) {
                    R.data[k][k] = -R.data[k][k];
                    Q.data[k][0] = -Q.data[k][0];
                    Q.data[k][1] = -Q.data[k][1];
                }
            }
		}
	}
	
	// Checking values close to zero for the matrix R
	for (std::size_t i = 0; i < n; i++) {
    	for (std::size_t j = 0; j < n; j++) {
        	if (std::abs(R.data[i][j]) <= 1e-10) {
            	R.data[i][j] = 0.0;
        	}
    	}
	}

	Q = Q.getTransposed();
	// Extra modification: Check if R and Q are 2x2 matrices
    if (R.n == 2 && Q.n == 2) {
        // If the second row of R is all zeros, set the second column of Q to zeros
        if (R.data[1][0] == 0.0 && R.data[1][1] == 0.0) {
            Q.data[0][1] = 0.0;
            Q.data[1][1] = 0.0;
        }
    }
	
    std::cout << Q << "\n" << R << "\n";   
	std::cout << Q * R << "\n";  
	return std::make_pair(Q, R); 
} */

// Second realization: getQRDecomposition using Gram-Schmidt orthogonalization process 
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::getQRDecomposition() const {
    Matrix Q(*this);
    Matrix R(n, n, T(0));

    for (std::size_t j = 0; j < n; j++) {
        for (std::size_t i = 0; i < j; i++) {
            T dotProduct = T(0);
            T normSquared = T(0);

            for (std::size_t k = 0; k < n; k++) {
                dotProduct += Q.data[k][i] * Q.data[k][j];
                normSquared += Q.data[k][i] * Q.data[k][i];
            }

            R.data[i][j] = dotProduct;
            for (std::size_t k = 0; k < n; k++) {
                Q.data[k][j] -= dotProduct * Q.data[k][i] / normSquared;
            }
        }

        T norm = T(0);
        for (std::size_t k = 0; k < n; k++) {
            norm += Q.data[k][j] * Q.data[k][j];
        }
        norm = std::sqrt(norm);

        R.data[j][j] = norm;
        for (std::size_t k = 0; k < n; k++) {
            Q.data[k][j] /= norm;
        }
    }

	Q = Q.getTransposed();
//	std::cout << Q * R << "\n"; 
    return std::make_pair(Q, R);
}

// isSquare
template<typename T>
bool Matrix<T>::isSquare() const {
	return (m == n); 
}

// isSymmetric
template<typename T>
bool Matrix<T>::isSymmetric() const {
	if (m != n) 
		return false; // Symmetric matrix should be a square
	for (std::size_t i = 0; i < m; i++) {
		for(std::size_t j = 0; j < n; j++) {
			if (data[i][j] != data[j][i]) 
				return false;
		}
	}
	
	return true;
}

// isUpperTriangular
template<typename T>
bool Matrix<T>::isUpperTriangular() const {
	if (m != n) 
		return false; // Matrix should be a square
	for (std::size_t i = 1; i < m; i++) {
		for(std::size_t j = 0; j < i; j++) {
			if (data[i][j] != 0)
				return false; 
		}
	}
	
	return true;
}		

// isLowerTriangular
template<typename T>
bool Matrix<T>::isLowerTriangular() const {
	if (m != n) 
		return false; // Matrix should be a square
	for (std::size_t i = 0; i < m - 1; i++) {
		for(std::size_t j = i + 1; j < n; j++) {
			if (data[i][j] != 0)
				return false; 
		}
	}
	
	return true;
}

// isDiagonal
template<typename T>
bool Matrix<T>::isDiagonal() const {
	if (m != n) 
		return false;
	for (std::size_t i = 0; i < m ; i++) {
		for(std::size_t j = 0; j < n; j++){
			if (i != j && data[i][j] != 0)
				return false; 
		}
	}
	
	return true;
}

// isOrthogonal ---> (A * A^T = I && A^T == A^-1)
template<typename T>
bool Matrix<T>::isOrthogonal() const {
    if (!isSquare())
        return false;
	
	Matrix<T> transposed = getTransposed();
    Matrix<T> identity = IdentityMatrix(n);

    Matrix<T> product = (*this) * transposed;
    Matrix<T> inverse = getInversedMatrix();

	// Rounding cycle of matrix product elements
	for (std::size_t i = 0; i < n; i++) {
        if (std::abs(product[i][i] - identity[i][i]) < 1.0)
        // Diff between an element approximated to 1 from a final matrix with 1
            product[i][i] = T(1);
    }

	// Rounding cycle of inverse matrix elements
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            T round_Inverse = static_cast<T>(static_cast<int>(inverse[i][j] * 10)) / 10.0;
            inverse[i][j] = round_Inverse;
        }
    }
	
    if (inverse == transposed && product == identity)
        return true;
        
    return false;
}

// isHermitian
template<typename T>
bool Matrix<T>::isHermitian() const {
    if (m != n)
    	return false; // Hermitian matrix should be a square
    Matrix<T> transposedMatrix = getTransposed();
	Matrix<T> conjugateMatrix = transposedMatrix.getConjugate();

    return (*this == conjugateMatrix);
}

// isConjugate
template<typename T>
bool Matrix<T>::isConjugate() const {
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            if (std::imag(data[i][j]) != 0) 
            	return true;
        }
    }
    return false; 
}

// isNormal
template<typename T>
bool Matrix<T>::isNormal() const {
    if (!isSquare()) {
        return false;
    }

    Matrix<T> transposedMatrix = getTransposed();
	Matrix<T> conjugateMatrix = transposedMatrix.getConjugate();

    return ((*this) * conjugateMatrix == conjugateMatrix * (*this));
}

// Assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
	m = other.m;
	n = other.n;
	data = other.data;
	return *this;
}

// Element Access
template<typename T>
std::vector<T>& Matrix<T>::operator[](std::size_t pos) {
	return data[pos];
}

template<typename T>
const std::vector<T>& Matrix<T>::operator[](std::size_t pos) const {
	return data[pos];
}

// Matrix operators

// Operator * Operands: matrices
template<typename V>
Matrix<V> operator*(const Matrix<V>& left, const Matrix<V>& right) {
	if (left.n != right.m) {
		throw std::runtime_error("Invlid size of operands");
	}

	Matrix<V> res(left.n, right.m);
	for (std::size_t i = 0; i < left.m; i++) {
		for (std::size_t j = 0; j < right.n; j++) {
			for (std::size_t k = 0; k < left.n; k++) {
				res.data[i][j] += left.data[i][k] * right.data[k][j];
			}
		}
	}
	return res;
}

// Operator * Operands: matrix - vector
template<typename V>
std::vector<V> operator*(const Matrix<V>& mat, const std::vector<V>& vec) {
	if (mat.n != vec.size()) {
		throw std::runtime_error("Invalid size of operands");
	}

	std::vector<V> res(vec.size());

	for (std::size_t i = 0; i < mat.m; i++) {
		for (std::size_t j = 0; j < mat.n; j++) {
			res[i] += mat.data[i][j] * vec[j];
		}
	}
	return res;
}

// Operator * Operands: vector - matrix
template<typename V>
std::vector<V> operator*(const std::vector<V>& mat, const Matrix<V>& vec) {
	if (mat.n != vec.size()) {
		throw std::runtime_error("Invalid size of operands");
	}

	std::vector<V> res(vec.size());

	for (std::size_t i = 0; i < mat.m; i++) {
		for (std::size_t j = 0; j < mat.n; j++) {
			res[i] += vec[j] * mat.data[i][j];
		}
	}
	return res;
}

// Operator * Operands: matrix - scalar
template<typename V>
Matrix<V> operator*(const Matrix<V>& matrix, const V& scalar) {
    Matrix<V> result(matrix.m, matrix.n);
    for (std::size_t i = 0; i < matrix.m; i++) {
        for (std::size_t j = 0; j < matrix.n; j++) {
            result.data[i][j] = matrix.data[i][j] * scalar;
        }
    }
    return result;
}

// Operator * Operands: scalar - matrix
template<typename V>
Matrix<V> operator*(const V& scalar, const Matrix<V>& matrix) {
	 Matrix<V> result(matrix.m, matrix.n);
    for (std::size_t i = 0; i < matrix.m; i++) {
        for (std::size_t j = 0; j < matrix.n; j++) {
            result.data[i][j] = scalar * matrix.data[i][j];
        }
    }
    return result;
}

// Operator /
template<typename V>
Matrix<V> operator/(const Matrix<V>& matrix, const V& scalar) {
    Matrix<V> result(matrix.m, matrix.n);
    for (std::size_t i = 0; i < matrix.m; i++) {
        for (std::size_t j = 0; j < matrix.n; j++) {
            result.data[i][j] = matrix.data[i][j] / scalar;
        }
    }
    return result;
}

// Operator +
template<typename V>
Matrix<V> operator+(const Matrix<V>& matrix1, const Matrix<V>& matrix2) {
    if (matrix1.m != matrix2.m || matrix1.n != matrix2.n)
		throw std::runtime_error("Invalid size of operands");
    Matrix<V> result(matrix1.m, matrix1.n);
    for (std::size_t i = 0; i < matrix1.m; i++) {
        for (std::size_t j = 0; j < matrix1.n; j++) {
            result.data[i][j] = matrix1.data[i][j] + matrix2.data[i][j];
        }
    }
    return result;
}

// Operator -
template<typename V>
Matrix<V> operator-(const Matrix<V>& matrix1, const Matrix<V>& matrix2) {
    if (matrix1.m != matrix2.m || matrix1.n != matrix2.n)
		throw std::runtime_error("Invalid size of operands");
    Matrix<V> result(matrix1.m, matrix1.n);
    for (std::size_t i = 0; i < matrix1.m; i++) {
        for (std::size_t j = 0; j < matrix1.n; j++) {
            result.data[i][j] = matrix1.data[i][j] - matrix2.data[i][j];
        }
    }
    return result;
}

// Operator +=
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
    if (m != other.m || n != other.n)
		throw std::runtime_error("Invalid size of operands");
    for (std::size_t i = 0; i < m; i++) {
    	for (std::size_t j = 0; j < n; j++) {
        	data[i][j] += other.data[i][j];
        }
    }
    return *this;
}

// Operator -=
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
    if (m != other.m || n != other.n)
		throw std::runtime_error("Invalid size of operands");
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            data[i][j] -= other.data[i][j];
        }
    }
    return *this;
}

// Operator *= with other matrix
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
    if (n == other.n)
		throw std::runtime_error("Invalid size of operands");
    Matrix<T> result(m, other.n);
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < other.n; j++) {
            result.data[i][j] = 0;
            for (std::size_t k = 0; k < n; k++) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    *this = result;
    return *this;
}

// Operator *= with scalar
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& scalar) {
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            data[i][j] *= scalar;
        }
    }
    return *this;
}

// Operator /=
template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& scalar) {
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            data[i][j] /= scalar;
        }
    }
    return *this;
}

// Boolean Comparison Operators

// Operator ==
template<typename V>
bool operator==(const Matrix<V>& left, const Matrix<V>& right) {
    if (left.m != right.m || left.n != right.n) {
        return false;
    }
    
    for (std::size_t i = 0; i < left.m; i++) {
        for (std::size_t j = 0; j < left.n; j++) {
            if (left.data[i][j] != right.data[i][j]) {
                return false;
            }
        }
    }
    
    return true;
}

// Operator !=
template<typename V>
bool operator!=(const Matrix<V>& left, const Matrix<V>& right) {
	if (left.m == right.m && left.n == right.n) {
        return true;
    }
    
    for (std::size_t i = 0; i < left.m; i++) {
        for (std::size_t j = 0; j < left.n; j++) {
            if (left.data[i][j] == right.data[i][j]) {
                return true;
            }
        }
    }
    
    return false;
}

// In/Out operators

template<typename T>
std::istream& operator>>(std::istream& in, Matrix<T>& self) {
	for (auto& row : self.data) {
		for (auto& elem : row) {
			in >> elem;
		}
	}
	return in;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& self) {
	for (const auto& row : self.data) {
		for (const auto& elem : row) {
			out << elem << " ";
		}
		out << std::endl;
	}
	return out;
}
