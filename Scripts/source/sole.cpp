#include <cmath>

// Constructors
template<typename T>
SOLE<T>::SOLE(std::size_t n, const Matrix<T>& A, const std::vector<T>& b) : n(n), A(A), b(b) {}

// Several additional methods
template<typename T>
std::vector<T> operator+(const std::vector<T>& left, const std::vector<T>& right) {
	if (left.size() != right.size()) {
		throw std::runtime_error("Invalid size of operands");
	}
	
    std::vector<T> res(left);
	for (std::size_t i = 0; i < res.size(); i++) {
		res[i] += right[i];
	}
	
    return res;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& left, const std::vector<T>& right) {
	if (left.size() != right.size()) {
		throw std::runtime_error("Invalid size of operands");
	}

	std::vector<T> res(left);
	for (std::size_t i = 0; i < res.size(); i++) {
		res[i] -= right[i];
	}

	return res;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T& scalar) {
	std::vector<T> res(vec);
	for (auto& elem : res) {
		elem *= scalar;
	}

	return res;
}

template<typename T>
double norm(const std::vector<T>& vec) {
	double sum = 0;
	for (const auto& elem : vec) {
		sum += elem * elem;
	}

	return std::sqrt(sum);
}

// SOLE solution methods

// MAIN METHODS

// GaussSolve
template<typename T>
std::vector<T> SOLE<T>::GaussSolve() const {
    std::size_t n = b.size();

    // Create copies of the matrix and vector
    Matrix<T> A_copy = A;
    std::vector<T> b_copy(b);
    
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            A_copy[i][j] = A[i][j];
        }
    }
    // Forward move
    for (std::size_t k = 0; k < n; k++) {
        // Find max element in col k
        std::size_t maxRow = k;
        T maxElement = A_copy[k][k];
        for (std::size_t i = k + 1; i < n; i++) {
            if (A_copy[i][k] > maxElement) {
                maxElement = A_copy[i][k];
                maxRow = i;
            }
        }

        // swap rows
        std::swap(A_copy[k], A_copy[maxRow]);
        std::swap(b_copy[k], b_copy[maxRow]);

        // Transform to triangular form
        for (std::size_t i = k + 1; i < n; i++) {
            T temp = A_copy[i][k] / A_copy[k][k];
            A_copy[i][k] = 0.0;
            for (std::size_t j = k + 1; j < n; j++) {
                A_copy[i][j] -= temp * A_copy[k][j];
            }
            b_copy[i] -= temp * b_copy[k];
        }
    }

    // Reverse move
    std::vector<T> x(n);
    for (long k = n - 1; k >= 0; --k) {
        x[k] = b_copy[k];
        for (long i = k + 1; i < n; ++i) {
            x[k] -= A_copy[k][i] * x[i];
        }
        x[k] /= A_copy[k][k];
    }

    return x;
}   

// LUSolve 
template<typename T>
std::vector<T> SOLE<T>::LUSolve() const {
    Matrix<T> A_copy = A;
    std::vector<T> b_copy = b;

    // Sizes check
    if (A_copy.colsCount() != A_copy.rowsCount()) {
        throw std::runtime_error("Matrix A must be square");
    }

    if (A_copy.colsCount() != b_copy.size()) {
        throw std::runtime_error("Invalid size of vector b");
    }

    std::size_t n = A_copy.colsCount();

    // Initialize the matrices L and Y as a vector of vectors
    std::vector<std::vector<T>> L(n, std::vector<T>(n, 0.0));
    std::vector<std::vector<T>> U(n, std::vector<T>(n, 0.0));

    // LU decomposition
    for (std::size_t i = 0; i < n; i++) {
        // U elements
        for (std::size_t k = i; k < n; k++) {
            T sum = 0.0;
            for (std::size_t j = 0; j < i; j++) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A_copy[i][k] - sum;
        }

        // L elements
        for (std::size_t k = i; k < n; k++) {
            if (i == k) {
                L[i][i] = 1.0;
            } else {
                T sum = 0.0;
                for (std::size_t j = 0; j < i; j++) {
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (A_copy[k][i] - sum) / U[i][i];
            }
        }
    }

    // Solve system Ly = b_copy (forward move)
    std::vector<T> y(n);
    for (std::size_t i = 0; i < n; i++) {
        T sum = 0.0;
        for (std::size_t j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b_copy[i] - sum) / L[i][i];
    }

    // Solve system Ux = y (reverse move)
    std::vector<T> x(n);
    for (long i = n - 1; i >= 0; --i) {
        T sum = 0.0;
        for (long j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
        if (x[i] == -0.0)
            x[i] = 0.0;
    }

    return x;
}

// CholeskySolve
template<typename T>
std::vector<T> SOLE<T>::CholeskySolve() const {
    Matrix<T> A_copy = A;
    std::vector<T> b_copy = b;

    // Sizes check
    if (A_copy.colsCount() != A_copy.rowsCount())
        throw std::runtime_error("Matrix A must be square");

    if (A_copy.colsCount() != b_copy.size())
        throw std::runtime_error("Invalid size of operands");

    if (A_copy.getDet() <= 0)
        throw std::runtime_error("Matrix is not positive definite"); 
    
    std::size_t n = A_copy.colsCount();

    // Perform Cholesky decomposition
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < i; j++) {
            for (std::size_t k = 0; k < j; k++) {
                A_copy[i][j] -= A_copy[i][k] * A_copy[j][k];
            }
            A_copy[i][j] /= A_copy[j][j];
        }

        for (std::size_t k = 0; k < i; k++) {
            A_copy[i][i] -= A_copy[i][k] * A_copy[i][k];
        }

        A_copy[i][i] = std::sqrt(A_copy[i][i]);
        if (std::isnan(A_copy[i][i]))
            throw std::runtime_error("Cholesky decomposition failed: NaN encountered");
    }

    // Solve system L*L^T*x = b_copy
    std::vector<T> y(n);
    for (std::size_t i = 0; i < n; i++) {
        T sum = 0.0;
        for (std::size_t j = 0; j < i; j++) {
            sum += A_copy[i][j] * y[j];
        }
        y[i] = (b_copy[i] - sum) / A_copy[i][i];
    }

    std::vector<T> x(n);
    for (std::size_t i = n - 1; i >= 0; i--) {
        T sum = 0.0;
        for (std::size_t j = i + 1; j < n; j++) {
            sum += A_copy[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / A_copy[i][i];
    }

    return x;
}

// InversedSolve
template<typename T>
std::vector<T> SOLE<T>::InversedSolve() const {
    Matrix<T> A_copy = A;
    std::vector<T> b_copy = b;

    // Sizes check
    if (A_copy.colsCount() != A_copy.rowsCount()) {
        throw std::runtime_error("Matrix A must be square");
    }

    if (A_copy.colsCount() != b_copy.size()) {
        throw std::runtime_error("Invalid size of operands");
    }

    std::size_t n = A_copy.colsCount();

    // Calculate the inverse matrix
    Matrix<T> A_inverse = A_copy.getInversedMatrix();

    // Solve the system Ax = b -> x = A^-1 * b
    std::vector<T> x = A_inverse * b_copy;

    return x;
}

// CramerSolve
template<typename T>
std::vector<T> SOLE<T>::CramerSolve() const {
    Matrix<T> A_copy = A;
    std::vector<T> b_copy = b;

    // Sizes check
    if (A_copy.colsCount() != A_copy.rowsCount()) {
        throw std::runtime_error("Matrix A must be square");
    }

    if (A_copy.colsCount() != b_copy.size()) {
        throw std::runtime_error("Invalid size of operands");
    }

    std::size_t n = A_copy.colsCount();

    // Calculate the determinant of A
    T detA = A_copy.getDet();

    // Check if the determinant is zero
    if (detA == 0) {
        throw std::runtime_error("The system has no unique solution");
    }

    std::vector<T> x(n);

    for (std::size_t i = 0; i < n; i++) {
        Matrix<T> A_i_copy = A_copy;

        // Replace the i-th column of A_i_copy with b_copy
        for (std::size_t j = 0; j < n; j++) {
            A_i_copy[j][i] = b_copy[j];
        }

        // Calculate the determinant of A_i_copy
        T detA_i = A_i_copy.getDet();

        // Calculate the i-th component of the solution x
        x[i] = detA_i / detA;
    }

    return x;
}

// ThomasAlgorithm
template <typename T>
std::vector<T> SOLE<T>::ThomasAlgorithm() const {
     for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            if ((i == j || i == j + 1 || i == j - 1) && A[i][j] == 0)
            	throw std::logic_error("Matrix is not tridiagonal");
            else if (!((i == j || i == j + 1 || i == j - 1) || A[i][j] == 0))
            	throw std::logic_error("Matrix is not tridiagonal");
    	}
    }

    std::vector<T> x(n);
    std::vector<T> alpha(n + 1);
    std::vector<T> beta(n + 1);
    alpha[1] = -A[0][1] / A[0][0];
    beta[1] = b[0] / A[0][0];

    for (std::size_t i = 1; i < n; i++) {
        beta[i + 1] = (b[i] - A[i][i - 1] * beta[i]) / (A[i][i - 1] * alpha[i] + A[i][i]);
        alpha[i + 1] = -A[i][i + 1] / (A[i][i - 1] * alpha[i] + A[i][i]);
    }

    x[n - 1] = beta[n];
    for (std::size_t i = n - 2; i >= 0; i--) {
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
    }
    
	return x;
}

// Iterative Methods

// JacobiIteration --- https://en.wikipedia.org/wiki/Jacobi_method
template<typename T>
std::vector<T> SOLE<T>::JacobiIteration(const std::vector<T>& x0, double eps) const {
    std::vector<T> x(x0);
    std::vector<T> xk(x0);
    std::vector<T> remain = b - A * x;

    while (norm(remain) > eps) {
        for (std::size_t i = 0; i < n; i++) {
            T sum = 0.0;
            for (std::size_t j = 0; j < n; j++) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            xk[i] = (b[i] - sum) / A[i][i];
        }

        remain = b - A * xk;
        x = xk;
    }

    return x;
}

// SeidelIteration --- https://en.wikipedia.org/wiki/Gauss–Seidel_method
template<typename T>
std::vector<T> SOLE<T>::SeidelIteration(const std::vector<T>& x0, double epsilon) const {
    std::size_t n = b.size();
    std::vector<T> x(x0);
    std::vector<T> xPrev(x);

    T error = epsilon + 1.0;
    T iteration = 0;

    while (error > epsilon) {
        xPrev = x;
        error = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            T sum1 = 0.0;
            T sum2 = 0.0;

            for (std::size_t j = 0; j < i; ++j)
                sum1 += A[i][j] * x[j];

            for (std::size_t j = i + 1; j < n; ++j)
                sum2 += A[i][j] * xPrev[j];

            x[i] = (b[i] - sum1 - sum2) / A[i][i];

            T diff = std::abs(x[i] - xPrev[i]);
            if (diff > error)
                error = diff;
        }

        iteration++;
    }
    
    std::cout << "Number of steps: " << iteration << std::endl;
    return x;
}

// UpperRelaxationIteration --- https://en.wikipedia.org/wiki/Successive_over-relaxation
template<typename T>
std::vector<T> SOLE<T>::UpperRelaxationIteration(const std::vector<T>& x0, double omega, const T& epsilon) const {
    std::size_t n = b.size();
    std::vector<T> x(x0);
    std::vector<T> xPrev(x);

    T error = epsilon + 1.0;

    while (error > epsilon) {
        for (std::size_t i = 0; i < n; i++) {
            // sum1 of the right side of the equation starting from element i + 1
            T sum1 = 0.0;
            for (std::size_t j = i + 1; j < n; j++) {
                sum1 += A[i][j] * x[j];
            }
            // sum2 of the right side of the equation starting from 0 to i - 1
            T sum2 = 0.0;
            for (std::size_t j = 0; j < i; j++) {
                sum2 += A[i][j] * x[j];
            }

            x[i] = (1.0 - omega) * xPrev[i] + (omega / A[i][i]) * (b[i] - sum1 - sum2);
        }

        // Checking the convergence condition
        error = 0.0;
        for (std::size_t i = 0; i < n; i++) {
            error += std::abs(x[i] - xPrev[i]);
        }

        if (error < epsilon)
            return x;
        xPrev = x; 
    }

    return x;
}

// SimpleIteration --- https://en.wikipedia.org/wiki/Fixed-point_iteration
template<typename T>
std::vector<T> SOLE<T>::SimpleIteration(const std::vector<T>& x0, double eps, const T& tau) const {
    std::vector<T> x(x0);
	std::vector<T> xk(x0);
	
    while(true) {
		xk = (b - A * x) * tau + x;
		// std::cout << norm(xk - x) << "\n";
		if (norm(xk - x) < eps) {
			break;
		}
		
        x = xk;
	}

	return x;
}