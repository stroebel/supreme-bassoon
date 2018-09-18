#include <iostream>

#include <eigen3/Eigen/Core>
#include "include/typedefs.h"
#include "include/jacobi.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{
  MatrixXd A(4, 4);
  VectorXd b(4), m(4);
  VectorXd x0 = VectorXd::Zero(4);

  A(0, 0) = 10; A(0, 1) = -1; A(0, 2) = 2; A(0, 3) = 0;
  A(1, 0) = -1; A(1, 1) = 11; A(1, 2) = -1; A(1, 3) = 3;
  A(2, 0) = 2; A(2, 1) = -1; A(2, 2) = 10; A(2, 3) = -1;
  A(3, 0) = 0; A(3, 1) = 3; A(3, 2) = -1; A(3, 3) = 8;

  b(0) = 6; b(1) = 25; b(2) = -11; b(3) = 15;

  m = jacobi(A, b, x0, 20, 0.005);

  std::cout << m << std::endl;

  return 0;
}
