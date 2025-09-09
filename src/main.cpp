#include <iostream>
#include <Eigen/Dense>

int main(int argc, char** argv) {
    // Create two 3x3 matrices with some example values
    Eigen::Matrix3d matrixA;
    matrixA << 1, 2, 3,
               4, 5, 6,
               7, 8, 9;
    
    Eigen::Matrix3d matrixB;
    matrixB << 9, 8, 7,
               6, 5, 4,
               3, 2, 1;
    
    // Display the input matrices
    std::cout << "Matrix A:\n" << matrixA << "\n\n";
    std::cout << "Matrix B:\n" << matrixB << "\n\n";
    
    // Perform matrix multiplication
    Eigen::Matrix3d result = matrixA * matrixB;
    
    // Display the result
    std::cout << "Result (A * B):\n" << result << "\n\n";
    
    // Demonstrate with different sized matrices
    Eigen::Matrix<double, 2, 3> matrixC;
    matrixC << 1, 2, 3,
               4, 5, 6;
    
    Eigen::Matrix<double, 3, 2> matrixD;
    matrixD << 1, 2,
               3, 4,
               5, 6;
    
    std::cout << "Matrix C (2x3):\n" << matrixC << "\n\n";
    std::cout << "Matrix D (3x2):\n" << matrixD << "\n\n";
    
    // Multiply 2x3 * 3x2 = 2x2 result
    Eigen::Matrix2d result2 = matrixC * matrixD;
    std::cout << "Result (C * D):\n" << result2 << "\n\n";
    
    // Demonstrate dynamic matrices
    Eigen::MatrixXd dynamicA = Eigen::MatrixXd::Random(4, 3);
    Eigen::MatrixXd dynamicB = Eigen::MatrixXd::Random(3, 4);
    
    std::cout << "Random Matrix A (4x3):\n" << dynamicA << "\n\n";
    std::cout << "Random Matrix B (3x4):\n" << dynamicB << "\n\n";
    
    Eigen::MatrixXd dynamicResult = dynamicA * dynamicB;
    std::cout << "Result (A * B) - 4x4:\n" << dynamicResult << "\n";
    
    return 0;
}