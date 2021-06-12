#ifndef UTILS_NUMERICAL_H
#define UTILS_NUMERICAL_H

#include "Utils.h"
#include <cmath>
#include <vector>
using namespace std;

// someuseful definitions
const double MAX_DOUBLE_VAL = 1.0e100;
const double MAX_NEG_DOUBLE_VAL = -1.0 * MAX_DOUBLE_VAL;
const double MIN_POS_VAL = 1.0e-40;

// Some matrix utilities

// template<class T>
// T MatrixPermanent(const vector<T>& A, int n); // expects n by n matrix
// encoded as vector
inline int *dec2binarr(long n, int dim) {
  // note: res[dim] will save the sum res[0]+...+res[dim-1]
  int *res = (int *)calloc(dim + 1, sizeof(int));
  int pos = dim - 1;

  // note: this will crash if dim < log_2(n)...
  while (n > 0) {
    res[pos] = n % 2;
    res[dim] += res[pos];
    n = n / 2; // integer division
    pos--;
  }

  return res;
}

template <class T> T MatrixPermanent(const vector<T> &A, int n) {
  // cout << "MatrixPermanent: n = " << n << endl;
  // expects n by n matrix encoded as vector
  T sum = 0;
  T rowsumprod, rowsum;
  // int* chi = new int[n + 1];
  int *chi;
  double C = (double)pow((double)2, n);

  // loop all 2^n submatrices of A
  for (int k = 1; k < C; k++) {
    // cout << "k = " << k << endl;
    rowsumprod = 1;
    chi = dec2binarr(k, n); // characteristic vector

    // loop columns of submatrix #k
    for (int m = 0; m < n; m++) {
      // cout << "m = " << m << endl;
      rowsum = 0;

      // loop rows and compute rowsum
      for (int p = 0; p < n; p++) {
        // cout << "p = " << p << endl;
        YW_ASSERT_INFO(m * n + p < (int)A.size(), "array out of bound");
        rowsum += chi[p] * A[m * n + p];
      }
      // update product of rowsums
      rowsumprod *= rowsum;

      // (optional -- use for sparse matrices)
      // if (rowsumprod == 0) break;
    }

    sum += (T)pow((double)-1, n - chi[n]) * rowsumprod;
    free(chi);
  }

  // delete [] chi;

  return sum;
}

// compute the product
template <class T> T CalcProductOfVec(const vector<T> &A) {
  YW_ASSERT_INFO(A.size() > 0, "Must have at least one item");
  T res = A[0];
  for (int i = 1; i < (int)A.size(); ++i) {
    res *= A[i];
  }
  return res;
}

// useful algorithms like Brent's method
class NumericalAlgoUtils {
public:
  virtual double EvaluateAt(double pt, void *pParam) = 0;
  double Func1DMinBrent(double ax, double bx, double cx, double tol,
                        double *xmin);
  virtual bool IsSignificantlyLarge(double v1, double v2) const;
  static bool IsLikeliSignificantlyLargeThresNum(double valLikeli1,
                                                 double valLikeli2,
                                                 int numItems, double thres);
};

// statistics related
double CalcApproxCDFStdNormal(double val);
double CalcBinomialProb(double p, int n, int k);

// other numerical related methods
double RoundDoubleValTo(double val, int numFractionDigits);
int GetCeilingPowerOf(int val, int base);
int RoundToInt(double val);
bool IsConvergedWithin(double valCurr, double valPre, double maxDiffFrac);
void NormalizeVec(vector<double> &vecDoubles);
double CalcSumOfSquareError(const vector<double> &vecDoubles1,
                            const vector<double> &vecDoubles2);
double CalcFactorial(int n);

#endif
