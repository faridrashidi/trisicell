#include "UtilsNumerical.h"
#include "Utils3.h"
#include <cmath>
#include <cstdlib>

// Some matrix utilities
// YW: seem to be some risk of memory issue: not freeing???

#if 0 // move implementaiton to header file
inline int* dec2binarr(long n, int dim)
{
    // note: res[dim] will save the sum res[0]+...+res[dim-1]
    int* res = (int*)calloc(dim + 1, sizeof(int));
    int pos = dim - 1;

    // note: this will crash if dim < log_2(n)...
    while (n > 0)
    {
        res[pos] = n % 2;
        res[dim] += res[pos];
        n = n / 2; // integer division
        pos--;
    }

    return res;
}

template<class T>
T MatrixPermanent(const vector<T>& A, int n)
{
	// expects n by n matrix encoded as vector
    T sum = 0;
    T rowsumprod, rowsum;
    int* chi = new int[n + 1];
    double C = (double)pow((double)2, n);

    // loop all 2^n submatrices of A
    for (int k = 1; k < C; k++)
    {
        rowsumprod = 1;
        chi = dec2binarr(k, n); // characteristic vector

        // loop columns of submatrix #k
        for (int m = 0; m < n; m++)
        {
            rowsum = 0;

            // loop rows and compute rowsum
            for (int p = 0; p < n; p++)
                rowsum += chi[p] * A[m * n + p];

            // update product of rowsums
            rowsumprod *= rowsum;

            // (optional -- use for sparse matrices)
            // if (rowsumprod == 0) break;
        }

        sum += (T)pow((double)-1, n - chi[n]) * rowsumprod;
    }

	//delete [] chi;

    return sum;
}

#endif

///////////////////////////////////////////////////////////////////////////////////////

double NumericalAlgoUtils ::Func1DMinBrent(double ax, double bx, double cx,
                                           double tol, double *xmin) {
  // cout << "Func1DMinBrent: " << ", [" << ax << ", " << bx << ", " << cx  <<
  // ", tol  " << tol << "], \n";
  // YW: this function is based Numerical Receipe in C book.
  // search for best 1 D function (in this case, the likelihood) using Brent's
  // method
  // Given a function f, and given a bracketing triplet of abscissas ax, bx, cx
  // (such that bx is between ax and cx, and f(bx) is less than both f(ax) and
  // f(cx)), this routine isolates the minimum to a fractional precision of
  // about tol using Brent�s method. The abscissa of the minimum is returned as
  // xmin, and the minimum function value is returned as brent, the returned
  // function value.
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
// Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden
// ratio; ZEPS is a small number that protects against trying to achieve
// fractional accuracy for a minimum that happens to be exactly zero.
#define SHFT(a, b, c, d)                                                       \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

  // cout << "Func1DMinBrent: ax=" << ax << ", bx=" << bx << ", cx=" << cx << ",
  // tol=" << tol << endl;
  int iter;
  double a, b, d = 0.0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x,
               xm;
  double e = 0.0; // This will be the distance moved on the step before last.
  a = (ax < cx ? ax : cx); // a and b must be in ascending order,
  b = (ax > cx ? ax : cx); // but input abscissas need not be.
  x = w = v = bx;          // Initializations...
  fw = fv = fx = EvaluateAt(x, NULL);

#if 0
	// in case f(a) < f(b) < f(c), stop
	double fa1 = EvaluateAt( a, NULL );
	double fb1 = EvaluateAt( b, NULL );
cout << "fa1 = " << fa1 << " for a = " << a << ", fb1= " << fb1 << " for b = " << b << ", fx = " << fx << ", x = " << x << endl;
	if( IsSignificantlyLarge(fa1, fx) == false || IsSignificantlyLarge( fb1, fx) == false )
	{
		//  take the minimum
		*xmin = a;
		double res1 = fa1;
		if( fa1 > fb1)
		{
			*xmin = b;
			res1 = fb1;
		}
		return res1;
	}
#endif

  for (iter = 1; iter <= ITMAX; iter++) { // Main program loop.
    // cout << "iteration " << iter << endl;
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) { // Test for done here.
      *xmin = x;
      // cout << "x = " << x << ", xm = " << xm << ", tol2 = " << tol2 << ", b =
      // " << b << ", a = " << a << endl; cout << "Here: STOP EARLY\n";
      return fx;
    }
    if (fabs(e) > tol1) { // Construct a trial parabolic fit.
                          // cout << "here...\n";
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) ||
          p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));
      // The above conditions determine the acceptability of the parabolic fit.
      // Here we take the golden section step into the larger of the two
      // segments.
      else {
        d = p / q; // Take the parabolic step.
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      }
    } else {
      // cout << "here2\n";
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }
    u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
    // cout << "u=" << u << endl;
    fu = EvaluateAt(u, NULL);
    // This is the one function evaluation per iteration.
    if (fu <= fx) { // Now decide what to do with our func
      if (u >= x)
        a = x;
      else
        b = x;         // tion evaluation.
      SHFT(v, w, x, u) // Housekeeping follows:
      SHFT(fv, fw, fx, fu)
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    } // Done with housekeeping. Back for
    // cout << "** -fx = " << -1.0*fx << endl;
  } // another iteration.
  // YW_ASSERT_INFO(false, "Too many iterations in brent");
  cout << "WARNING: Too many iterations in brent.\n";
  *xmin = x; // Never get here.
  return fx;
}

bool NumericalAlgoUtils ::IsSignificantlyLarge(double v1, double v2) const {
  // is v1 significantly larger than v2 (i.e. larger by some threshold)?
  // by default, the computed values are in log-space, and thus we ask then to
  // differ by at least 5%
  const double thresDef = log(1.05);
  return v1 >= v2 + thresDef;
}

bool NumericalAlgoUtils ::IsLikeliSignificantlyLargeThresNum(double valLikeli1,
                                                             double valLikeli2,
                                                             int numItems,
                                                             double thres) {
  // assume both are log-likelihood; thres: log(1.05) say
  // is likeli1 (per item) is significantly larger than likeli2 (per item)?
  double valLikeli1Ave = valLikeli1 / numItems;
  double valLikeli2Ave = valLikeli2 / numItems;
  return valLikeli1Ave >= valLikeli2Ave + thres;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
double RoundDoubleValTo(double val, int numFractionDigits) {
  // numFractiondigits: how many digits after . we want to keep
  YW_ASSERT_INFO(numFractionDigits >= 0, "numFracDigits:; must be positive");
  double ratioInc = pow(10.0, numFractionDigits);
  return round(val * ratioInc) / ratioInc;
}

int GetCeilingPowerOf(int val, int base) {
  // given a value, find the smallest (positive) power of base that is at least
  // this value
  int res = 1;
  while (val > res) {
    res *= base;
  }
  return res;
}

// statistics related
double CalcApproxCDFStdNormal(double val) {
  //
  const double pi = 3.1415926535897;
  double sign = 1.0;
  if (val < 0) {
    sign = -1.0;
  }
  return 0.5 * (1.0 + sign * (sqrt(1.0 - exp(-2.0 * val * val / pi))));
}

double CalcBinomialProb(double p, int n, int k) {
  YW_ASSERT_INFO(k <= n, "CalcBinomialProb: k must be smaller than n");
  double res = pow(p, k) * pow(1.0 - p, n - k) * CalcNumNChooseK(n, k);
  return res;
}

int RoundToInt(double val) { return (int)(val + 0.5); }

bool IsConvergedWithin(double valCurr, double valPre, double maxDiffFrac) {
  double valDiff = std::abs(valCurr - valPre);
  double valBase1 = std::abs(valCurr);
  double valBase2 = std::abs(valPre);
  double valBase = std::max(valBase1, valBase2);
  return valDiff <= maxDiffFrac * valBase;
}

void NormalizeVec(vector<double> &vecDoubles) {
  double sum = GetSumOfElements(vecDoubles);
  YW_ASSERT_INFO(sum > 0.0, "Cannot normalize a zero vector");
  for (int i = 0; i < (int)vecDoubles.size(); ++i) {
    vecDoubles[i] = vecDoubles[i] / sum;
  }
}

double CalcSumOfSquareError(const vector<double> &vecDoubles1,
                            const vector<double> &vecDoubles2) {
  //
  double res = 0.0;
  YW_ASSERT_INFO(vecDoubles1.size() == vecDoubles2.size(), "Sizes don't match");
  for (int i = 0; i < (int)vecDoubles1.size(); ++i) {
    double diff = vecDoubles1[i] - vecDoubles2[i];
    res += diff * diff;
  }
  return res;
}

double CalcFactorial(int n) {
  double res = 1.0;
  for (int i = 2; i <= n; ++i) {
    res *= i;
  }
  return res;
}
