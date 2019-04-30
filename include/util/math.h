#ifndef MATH_H_
#define MATH_HPP_

// STL
#include <cmath>
#include <complex>
#include <limits>

// Boost
#include <boost/math/special_functions/sign.hpp>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

#ifndef PI
  #define PI 3.1415926535897932384
#endif

template<typename Scalar>
inline bool is_negligible(Scalar x) {
    return std::abs(x) <= std::numeric_limits<Scalar>::epsilon();
}

Eigen::Vector4d hesseNormalForm(const Eigen::Vector4d& p) {
    Eigen::Vector4d h = p;

    h.head<3>() = h.head<3>().normalized();

    if (!is_negligible(p(3)))
        h *= -boost::math::sign(p(3)); // negative d

    return h;
}

Eigen::Vector4d transform_plane(const Eigen::Vector4d& p, const Eigen::Isometry3d& T) {
    Eigen::Vector4d h = hesseNormalForm(p);

    h.head<3>() = T.linear() * h.head<3>(); // Rotate normal vector

    const Eigen::Vector3d t = T.translation();
    const double d = h.dot(Eigen::Vector4d(-t(0), -t(1), -t(2), 1.0));
    h(3) = d;

    if (!is_negligible(d))
        h *= -boost::math::sign(d); // negative d

    return h;
}

/* Solve for real roots of the standard quadratic equation,
 * returning the number of real roots found.
 */
/* solve_quadratic.c - finds the real roots of a x^2 + b x + c = 0 */
// Adapted from gsl_poly.h (https://github.com/ampl/gsl/blob/master/poly/gsl_poly.h)
int solveQuadratic(double a, double b, double c, double &x0, double &x1) {
    if (is_negligible(a)) { /* Handle linear case */
        if (is_negligible(b))
            return 0;
        else {
          x0 = -c / b;
          return 1;
        }
    }

    double disc = b * b - 4.0 * a * c;
    if (disc > 0.0) {
        if (is_negligible(b)) {
            double r = std::sqrt(-c / a);
            x0 = -r;
            x1 =  r;
        } else {
            // What Every Computer Scientist Should Know About Floating-Point Arithmetic, Goldberg, D. in Computing Surveys, 1991 (Sec. 1.4 - Cancellation)
            // www.itu.dk/~sestoft/bachelor/IEEE754_article.pdf
            // Stable quadratic roots according to BKP Horn.
            // http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
            double sgnb = (b > 0 ? 1 : -1);
            double temp = -0.5 * (b + sgnb * std::sqrt(disc));
            double r1 = temp / a ;
            double r2 = c / temp ;

            if (r1 < r2) {
                x0 = r1;
                x1 = r2;
            } else {
                x0 = r2;
                x1 = r1;
            }
        }

        return 2;
    } else if (is_negligible(disc)) {
        x0 = -0.5 * b / a;
        x1 = -0.5 * b / a;

        return 2;
    } else
        return 0;
}


// Solves the right nullspace from QR decomposition,
// returning the size of the kernel
template<typename Derived, typename Scalar>
int solveNullspace(const Eigen::MatrixBase<Derived> &A, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &k) {
    Eigen::ColPivHouseholderQR<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr(A.transpose());
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ();

    int n = qr.dimensionOfKernel();
    k.resize(Q.rows(), n);

    k = Q.block(0, Q.cols() - n, Q.rows(), n);
    return qr.dimensionOfKernel();
}

template<typename T>
Eigen::Matrix<T,4,4> QuaternionMultMatLeft(const Eigen::Quaternion<T>& q)
{
    return (Eigen::Matrix<T,4,4>() << q.w(), -q.z(), q.y(), q.x(),
                                      q.z(), q.w(), -q.x(), q.y(),
                                      -q.y(), q.x(), q.w(), q.z(),
                                      -q.x(), -q.y(), -q.z(), q.w()).finished();
}

template<typename T>
Eigen::Matrix<T,4,4> QuaternionMultMatRight(const Eigen::Quaternion<T>& q)
{
    return (Eigen::Matrix<T,4,4>() << q.w(), q.z(), -q.y(), q.x(),
                                      -q.z(), q.w(), q.x(), q.y(),
                                      q.y(), -q.x(), q.w(), q.z(),
                                      -q.x(), -q.y(), -q.z(), q.w()).finished();
}

template<typename T>
void mat2RPY(const Eigen::Matrix<T, 3, 3>& m, T& roll, T& pitch, T& yaw)
{
    roll = std::atan2(m(2,1), m(2,2));
    pitch = std::atan2(-m(2,0), std::sqrt(m(2,1) * m(2,1) + m(2,2) * m(2,2)));
    yaw = std::atan2(m(1,0), m(0,0));
}

bool solveQuadraticEquation(double a, double b, double c, double& x1, double& x2)
{
    if (std::abs(a) < 1e-12)
    {
        x1 = x2 = -c / b;
        return true;
    }
    double delta2 = b * b - 4.0 * a * c;

    if (delta2 < 0.0)
    {
        return false;
    }

    double delta = std::sqrt(delta2);

    x1 = (-b + delta) / (2.0 * a);
    x2 = (-b - delta) / (2.0 * a);

    return true;
}

namespace {

// Solve depressed cubic using Cardano's method.
int SolveDepressedCubic(const double p, const double q,
                        std::complex<double>* roots) {
  if (p == 0.0) {
    roots[0] = std::pow(-1.0 * q, 1.0 / 3.0);
    return 1;
  }

  std::complex<double> cubic_root_of_unity(-0.5, 0.5 * std::sqrt(3.0));
  std::complex<double> temp = q * q / 4.0 + p * p * p / 27.0;
  std::complex<double> sqrt_t = std::sqrt(temp);
  std::complex<double> u = std::pow(-0.5 * q + sqrt_t, 1.0 / 3.0);
  std::complex<double> v = std::pow(-0.5 * q - sqrt_t, 1.0 / 3.0);
  roots[0] = u + v;
  roots[1] =
      u * cubic_root_of_unity + v * cubic_root_of_unity * cubic_root_of_unity;
  roots[2] =
      u * cubic_root_of_unity * cubic_root_of_unity + v * cubic_root_of_unity;
  return 3;
}

} // namespace

// Provides solutions to the equation a*x^2 + b*x + c = 0.
int SolveQuadratic(const double a, const double b, const double c,
                   std::complex<double>* roots) {
  // If the equation is actually linear.
  if (a == 0.0) {
    roots[0] = -1.0 * c / b;
    return 1;
  }

  const double D = b * b - 4 * a * c;
  const double sqrt_D = std::sqrt(std::abs(D));

  // Real roots.
  if (D >= 0) {
    // Stable quadratic roots according to BKP Horn.
    // http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
    if (b >= 0) {
      roots[0] = (-b - sqrt_D) / (2.0 * a);
      roots[1] = (2.0 * c) / (-b - sqrt_D);
    } else {
      roots[0] = (2.0 * c) / (-b + sqrt_D);
      roots[1] = (-b + sqrt_D) / (2.0 * a);
    }
    return 2;
  }

  // Use the normal quadratic formula for the complex case.
  roots[0].real(-b / (2.0 * a));
  roots[1].real(-b / (2.0 * a));
  roots[0].imag(sqrt_D / (2.0 * a));
  roots[1].imag(-sqrt_D / (2.0 * a));
  return 2;
}

int SolveQuadraticReals(const double a, const double b, const double c,
                        double* roots) {
  std::complex<double> complex_roots[2];
  int num_complex_solutions = SolveQuadratic(a, b, c, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    roots[num_real_solutions++] = complex_roots[i].real();
  }
  return num_real_solutions;
}

int SolveQuadraticReals(const double a, const double b, const double c,
                        const double tolerance, double* roots) {
  std::complex<double> complex_roots[2];
  int num_complex_solutions = SolveQuadratic(a, b, c, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    if (std::abs(complex_roots[i].imag()) < tolerance) {
      roots[num_real_solutions++] = complex_roots[i].real();
    }
  }
  return num_real_solutions;
}

// Provides solutions to the equation a*x^3 + b*x^2 + c*x + d = 0 using Cardan's
// method.
int SolveCubic(const double a, const double b, const double c, const double d,
               std::complex<double>* roots) {
  if (is_negligible(a)) {
      throw std::runtime_error("Quadratic!");
    //return SolveQuadratic(b, c, d, roots);
  }

  // Solve by first reducing the problem to a depressed cubic.
  double p = (3.0 * a * c - b * b) / (3.0 * a * a);
  double q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) /
             (27.0 * a * a * a);
  int num_solutions = SolveDepressedCubic(p, q, roots);
  // Transform solution back to normal params.
  roots[0] -= b / (3.0 * a);
  roots[1] -= b / (3.0 * a);
  roots[2] -= b / (3.0 * a);
  return num_solutions;
}

int SolveCubicReals(const double a, const double b, const double c,
                    const double d, double* roots) {
  std::complex<double> complex_roots[3];
  int num_complex_solutions = SolveCubic(a, b, c, d, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    roots[num_real_solutions++] = complex_roots[i].real();
  }
  return num_real_solutions;
}

int SolveCubicReals(const double a, const double b, const double c,
                    const double d, const double tolerance, double* roots) {
  std::complex<double> complex_roots[3];
  int num_complex_solutions = SolveCubic(a, b, c, d, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    if (std::abs(complex_roots[i].imag()) < tolerance) {
      roots[num_real_solutions++] = complex_roots[i].real();
    }
  }
  return num_real_solutions;
}

#endif // MATH_HPP_