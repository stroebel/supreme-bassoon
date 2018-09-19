#ifndef GAUSS_H
#define GAUSS_H

#include <armadillo>
#include <tuple>

using namespace arma;

/*
  Linear equation solution by Gauss-Seidel iterative technique.
  Input: Number of equations; a matrix A; a vector b; initial guess;
  maximum number of iterations; error tolerance
  Output: Approximate solution of the system or an error message

  Ref: Burden, R., Faires, D., Burden, A. 2016. Numerical Analysis, 10th ed.
  Cengage: Somewhere awful.
*/

std:: tuple<vec, std:: string>
  gauss_seidel(const mat &A, const mat &b, const vec &XO, const int N, const double TOL)
{
  /* First draft: Too much copying of vectors */
  int k = 0;

  int n = b.size() - 1;
  vec x(size(XO));
  /* This is not a great copy to do, but it is there fore now cf. previous comment */
  vec x0 = XO;
  double sum1, sum2;
  vec temp;

  while (k <= N)
    {
      for (int i = 0; i <= n; i++)
        {
          sum1 = 0; sum2 = 0;
          /* Using access with range checking here cf. Armadillo documentation */
          for (int j = 0; j <= i -1 ; j++)
            sum1 = sum1 + A(i, j)*x0(j);
          for (int j = i + 1; j <= n; j++)
            sum2 = sum2 + A(i, j)*x0(j);
          x(i) = (b(i) - sum1 - sum2)/A(i, i);
        }
      if (norm(x - x0, "inf") < TOL)
        return std:: make_tuple(x, "Tolerance reached");
      k++;
      x0 = x;
    }
  return std:: make_tuple(x, "Number of iterations exceeded");
}
#endif
