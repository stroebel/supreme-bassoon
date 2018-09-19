#include <iostream>
#include <armadillo>
#include <tuple>

#include "include/jacobi.h"
#include "include/Gauss-Seidel.h"

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
  vec x(4, fill:: zeros);

  std:: string output;

  std:: cout << "Armadillo solver solution" << std:: endl;
  std:: cout << solve(A, b) << std:: endl;

  vec s1(4);

  std:: tie(s1, output) = jacobi(A, b, x, 50, 0.005);

  std:: cout << "Jacobi iterative solver" << std:: endl;
  std:: cout << output << std:: endl;
  std:: cout << s1 << std:: endl;

  vec s2(4);

  std:: tie(s2, output) = gauss_seidel(A, b, x, 50, 0.005);

  std:: cout << "Gauss_seidel iterative solver" << std:: endl;
  std:: cout << output << std:: endl;
  std:: cout << s2 << std:: endl;

  return 0;
}
