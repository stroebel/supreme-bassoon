#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <tuple>

using namespace arma;

/*
 Linear equation solution by Jacobi iterative method.
 Input: Number of equations; a matrix A; a vector b; initial guess;
        maximum number of iterations; error tolerance
 Output: Approximate solution of the system or an error message

 Ref: Burden, R., Faires, D., Burden, A. 2016. Numerical Analysis, 10th ed.
      Cengage: Somewhere awful.
 */
std:: tuple<vec, std:: string>
  jacobi(const mat &A, const mat &b, const vec &XO, const int N, const double TOL)
{
  /* First draft: Too much copying of vectors */
  int k = 0;

  int n = b.size() - 1;
  vec x(4);
  vec x0 = XO;
  double sum;
  vec temp;

  while (k <= N)
    {
      for (int i = 0; i <= n; i++)
        {
          sum = 0;
          /* Using access with range checking here cf. armadillo documentation */
          for (int j = 0; j <=n; j++)
            if (j != i)
              sum = sum + A(i, j)*x0(j);
          x(i) = (b(i) - sum)/A(i, i);
        }
      if (norm(x - x0, "inf") < TOL)
        return std:: make_tuple(x, "Tolerance reached");
      k++;
      x0 = x;
    }
  return std:: make_tuple(x, "Number of iterations exceeded");
}

#endif
