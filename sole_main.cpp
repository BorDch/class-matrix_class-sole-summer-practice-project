#include "sole.hpp"

int main() {
    Matrix<double> A = {
        {4, 1, 0}, 
        {2, 2, 1}, 
        {0, 1, 2}
    };

    std::cout << "Matrix A: \n";
    A.print();

    std::vector<double> b = {1, 2, -3};
 
    SOLE<double> sole(3, A, b);
    
    // GaussSolve
    std::cout << "System solve using Gauss method:" << std::endl;
    std::vector<double> Gauss_solution = sole.GaussSolve();
    for (std::size_t i = 0; i < Gauss_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Gauss_solution[i] << std::endl;
    }
    
    /*
    // LUSolve
    std::cout << "System solve using LU decomposition:" << std::endl;
    std::vector<double> LU_solution = sole.LUSolve();
    for (std::size_t i = 0; i < LU_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << LU_solution[i] << std::endl;
    } 

    //CholeskySolve
    std::cout << "System solve using Cholesky method:" << std::endl;
    std::vector<double> Cholesky_solution = sole.CholeskySolve();
    for (std::size_t i = 0; i < Cholesky_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Cholesky_solution[i] << std::endl;
    } 
    
    // InversedSolve
    std::cout << "System solve using Inverse Matrix:" << std::endl;
    std::vector<double> Inverse_solution = sole.InversedSolve();
    for (std::size_t i = 0; i < Inverse_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Inverse_solution[i] << std::endl;
    } 

    // CramerSolve
    std::cout << "System solve using Cramer method:" << std::endl;
    std::vector<double> Cramer_solution = sole.CramerSolve();
    for (std::size_t i = 0; i < Cramer_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Cramer_solution[i] << std::endl;
    } 

    // ThomasAlgorithm
    std::cout << "System solve using Thomas method:" << std::endl;
    std::vector<double> solution = sole.ThomasAlgorithm();
    for (std::size_t i = 0; i < solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << solution[i] << std::endl;
    } 
    */

    // Iterative methods
    
    // JacobiIteration
    std::cout << "System solve using Jacobi method:" << std::endl;
    std::vector<double> x0(3, 0.0);
    double epsilon = 0.000001;
    std::vector<double> Jacobi_solution = sole.JacobiIteration(x0, epsilon);
    for (std::size_t i = 0; i < Jacobi_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Jacobi_solution[i] << std::endl;
    } 

    // SeidelIteration
    std::cout << "System solve using Seidel method:" << std::endl;
    std::vector<double> x1(3, 0.0);
    double epsilon1 = 0.000001;
    std::vector<double> Seidel_solution = sole.SeidelIteration(x1, epsilon1);
    for (std::size_t i = 0; i < Seidel_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << Seidel_solution[i] << std::endl;
    }

    // UpperRelaxationIteration
    std::cout << "System solve using Upper Relaxation method:" << std::endl;
    std::vector<double> x2(3, 0.0);
    double epsilon2 = 0.000001;
    double omega = 1.25;  // omega -> range(0, 2)
    std::vector<double> UpperRelaxation_solution = sole.UpperRelaxationIteration(x2, omega, epsilon2);
    for (std::size_t i = 0; i < UpperRelaxation_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << UpperRelaxation_solution[i] << std::endl;
    }  

    // SimpleIteration
    std::cout << "System solve using Simple Iteration method:" << std::endl;
    std::vector<double> x3(3, 0.0);
    double epsilon3 = 0.0000001;
    double tau = 0.1;
    std::vector<double> SimpleIteration_solution = sole.SimpleIteration(x3, epsilon3, tau);
    for (std::size_t i = 0; i < SimpleIteration_solution.size(); i++) {
        std::cout << "x" << i + 1 << " = " << SimpleIteration_solution[i] << std::endl;
    }
    
    return 0;
}
