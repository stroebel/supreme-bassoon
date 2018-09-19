#include <iostream>
#include <armadillo>

#include "include/jacobi.h"

using namespace arma;

int main()
{
  mat A = {
    {10, -1, 2, 0},
    {-1, 11, -1, 3},
    {2, -1, 10, -1},
    {0, 3, -1, 8}
  };

  vec b = { 6, 25, -11, 15};
  std::cout << "Armadillo solver solution" << std::endl;
  std::cout << solve(A, b) << std::endl;

  mat s(4,4);
  vec x(4, fill::zeros);

  s = jacobi(A, b, x, 20, 0.0001);

  std:: cout << "Jacobi iterative solver" << std:: endl;
  std:: cout << s << std:: endl;

  return 0;
}
