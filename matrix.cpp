#include <cstdlib>

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

// Assignment

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
	m = other.m;
	n = other.n;
	data = other.data;
	return *this;
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

// Element Access

template<typename T>
std::vector<T>& Matrix<T>::operator[](std::size_t pos) {
	return data[pos];
}

template<typename T>
const std::vector<T>& Matrix<T>::operator[](std::size_t pos) const {
	return data[pos];
}

//Matrix operations

template<typename T>
std::vector<T> operator*(const Matrix<T>& mat, const std::vector<T>& vec) {
	if (mat.n != vec.size()) {
		throw std::runtime_error("Invalid size of operands");
	}

	std::vector<T> res(vec.size());

	for (std::size_t i = 0; i < mat.m; i++) {
		for (std::size_t j = 0; j < mat.n; j++) {
			res[i] += mat.data[i][j] * vec[j];
		}
	}
	return res;
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
