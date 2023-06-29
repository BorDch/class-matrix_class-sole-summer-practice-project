#include <cmath>

template<typename T>
SOLE<T>::SOLE(std::size_t n, const Matrix<T>& A, const std::vector<T>& b) : n(n), A(A), b(b) {}

template<typename T>
double norm(const std::vector<T>& vec) {
	double sum = 0;
	for (const auto& elem : vec) {
		sum += elem * elem;
	}
	return std::sqrt(sum);
}

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
std::vector<T> SOLE<T>::SimpleIteration(const std::vector<T>& x0, double eps, const T& tau) const {

	std::vector<T> x(x0);
	std::vector<T> xk(x0);
	while(true) {
		xk = (b - A * x) * tau + x;
		std::cout << norm(xk - x) << "\n";
		if (norm(xk - x) < eps) {
			break;
		}
		x = xk;
	}

	return x;
}
