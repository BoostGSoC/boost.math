#include <cstdint>
#include <ctime>
#include <vector>
#include <iostream>
#include <iomanip>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>


// Investigate high precision calculations of the tgamma function with Boost.Multiprecision.

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000>, boost::multiprecision::et_off>      float_cpp_dec_et_off_type;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> >                                    float_cpp_dec_et_on_type;

typedef float_cpp_dec_et_off_type float_mp_t;
//typedef float_cpp_dec_et_on_type  float_mp_t;

namespace
{
  const float_mp_t& zero() { static const float_mp_t value_zero(0U);    return value_zero; }
  const float_mp_t& half() { static const float_mp_t value_half("0.5"); return value_half; }
  const float_mp_t& one () { static const float_mp_t value_one (1U);    return value_one; }

  const std::size_t highest_bernoulli_index()
  {
    // Find the high index n for Bn to produce the desired precision in Stirling's calculation.
    // TBD: This max-index is easy to template when considering a generic tgamma.
    return static_cast<std::size_t>(18.0 + (0.6 * static_cast<double>(std::numeric_limits<float_mp_t>::digits10)));
  }
}

namespace my
{
  float_mp_t tgamma(const float_mp_t& x);

  void generate_bernoulli_table(std::vector<float_mp_t>& bn, const std::uint32_t n)
  {
    // See reference "Computing Bernoulli and Tangent Numbers", Richard P. Brent.
    // See also the book Richard P. Brent and Paul Zimmermann, "Modern Computer Arithmetic",
    // Cambridge University Press, 2010, 237 pp.

    const std::uint32_t nn = (((n % static_cast<std::uint32_t>(2)) != static_cast<std::uint32_t>(0U)) ? static_cast<std::uint32_t>(n + 1U) : n);

    const std::int32_t m = nn / 2;

    std::vector<float_mp_t> tangent_numbers(static_cast<std::vector<float_mp_t>::size_type>(m + 1));

    tangent_numbers[0U] = zero();
    tangent_numbers[1U] = one();

    for(std::int32_t k = 2; k <= m; k++)
    {
      tangent_numbers[k] = (k - 1) * tangent_numbers[k - 1];
    }

    for(std::int32_t k = 2; k <= m; k++)
    {
      for(std::int32_t j = k; j <= m; j++)
      {
        tangent_numbers[j] = (tangent_numbers[j - 1] * (j - k)) + (tangent_numbers[j] * (j - k + 2));
      }
    }

    float_mp_t two_pow_two_m(4);

    bn.clear();
    bn.resize(static_cast<std::vector<float_mp_t>::size_type>(m + 1));

    for(std::int32_t i = 1; i < static_cast<std::int32_t>(tangent_numbers.size()); i++)
    {
      const std::int32_t two_i = static_cast<std::int32_t>(static_cast<std::int32_t>(2) * i);

      const float_mp_t b = (tangent_numbers[i] * two_i) / (two_pow_two_m * (two_pow_two_m - 1));

      const bool  b_neg = (static_cast<std::int32_t>(two_i % static_cast<std::int32_t>(4)) == static_cast<std::int32_t>(0));

      bn[i] = ((!b_neg) ? b : -b);

      two_pow_two_m *= 4;
    }

    bn[0U] =  one();
  }

  const float_mp_t& bernoulli_bn(const std::uint32_t n)
  {
    // TBD: Needs a better mechanism for warm-caching.
    static std::vector<float_mp_t> bernoulli_b2n;

    if(bernoulli_b2n.empty())
    {
      generate_bernoulli_table(bernoulli_b2n, static_cast<std::uint32_t>(highest_bernoulli_index()));
    }

    if(n == static_cast<std::uint32_t>(0U))
    {
      static const float_mp_t minus_half = -half();
      return minus_half;
    }
    else if((n % static_cast<std::uint32_t>(2)) != static_cast<std::uint32_t>(0U))
    {
      return zero();
    }
    else
    {
      const std::uint32_t n_half = static_cast<std::uint32_t>(n / static_cast<std::uint32_t>(2U));

      return ((n_half < bernoulli_b2n.size()) ? bernoulli_b2n[n_half] : zero());
    }
  }
}

float_mp_t my::tgamma(const float_mp_t& x)
{
  // TBD: small argument approximation?
  // TBD: NaN for zero or very near zero.
  // TBD: +infinity for negative integers or very close thereto.

  if(!boost::math::isfinite(x))
  {
    return x;
  }

  const bool b_neg = (x < zero());

  // Make a local, unsigned copy of the input argument.

  float_mp_t xx((!b_neg) ? x : -x);

  float_mp_t G;

  // Check if the argument should be scaled up for the Bernoulli series expansion.
  static const double digit_scale_of_argument = static_cast<double>(std::numeric_limits<float_mp_t>::digits10 * 1.7);

  static const std::int32_t min_arg_n = static_cast<std::int32_t>(digit_scale_of_argument);
  static const float_mp_t   min_arg_x = float_mp_t(min_arg_n);

  const std::int32_t n_recur = ((xx < min_arg_x) ? static_cast<std::int32_t>((min_arg_n - xx.convert_to<std::int32_t>()) + 1)
                                                 : static_cast<std::int32_t>(0));

  // Scale the argument up and use downward recursion later for the final result.
  if(n_recur != static_cast<std::int32_t>(0))
  {
    xx += n_recur;
  }

        float_mp_t one_over_x_pow_two_n_minus_one = one() / xx;
  const float_mp_t one_over_x2                    = one_over_x_pow_two_n_minus_one * one_over_x_pow_two_n_minus_one;
        float_mp_t sum                            = (bernoulli_bn(static_cast<std::uint32_t>(2U)) / static_cast<std::int32_t>(2)) * one_over_x_pow_two_n_minus_one;

  // Perform the Bernoulli series expansion.
  for(std::int32_t n2 = static_cast<std::int32_t>(4); n2 < static_cast<std::int32_t>(highest_bernoulli_index()); n2 += static_cast<std::int32_t>(2))
  {
    one_over_x_pow_two_n_minus_one *= one_over_x2;

    const float_mp_t term = (bernoulli_bn(static_cast<std::uint32_t>(n2)) * one_over_x_pow_two_n_minus_one) / static_cast<std::int32_t>(n2 * (n2 - static_cast<std::int32_t>(1)));

    // TBD: Use a break condition based on the first increasing term after
    // the end of the asymptotic decrease, but with more than 10 terms or so.

    sum += term;
  }

  static const float_mp_t half_ln_two_pi = log(boost::math::constants::two_pi<float_mp_t>()) / static_cast<std::int32_t>(2);

  G = exp(((((xx - half()) * log(xx)) - xx) + half_ln_two_pi) + sum);

  // Rescale the result using downward recursion if necessary.
  for(std::int32_t k = static_cast<std::int32_t>(0); k < n_recur; k++)
  {
    xx -= 1;
    G /= xx;
  }

  // Return the result, accounting for possible negative arguments.
  return ((!b_neg) ? G : -boost::math::constants::pi<float_mp_t>() / (xx * G * sin(boost::math::constants::pi<float_mp_t>() * xx)));
}

int main(int, char**)
{
  ::clock_t t0_first = ::clock();
  const float_mp_t g_first = my::tgamma(half());
  ::clock_t t1_first = ::clock();

  std::cout << std::setprecision(std::numeric_limits<float_mp_t>::digits10)
            << g_first
            << std::endl;

  ::clock_t t0_second = ::clock();
  const float_mp_t g_second = my::tgamma(half());
  ::clock_t t1_second = ::clock();

  std::cout << std::setprecision(std::numeric_limits<float_mp_t>::digits10)
            << g_second
            << std::endl;

  // Test with tgamma(1/2) == sqrt(pi).
  std::cout << std::setprecision(std::numeric_limits<float_mp_t>::digits10)
            << sqrt(boost::math::constants::pi<float_mp_t>())
            << std::endl;

  std::cout << "Time: Gamma with Bernoulli generation:  " << std::setprecision(7) << (static_cast<double>(t1_first)  - t0_first ) / CLOCKS_PER_SEC << std::endl;
  std::cout << "Time: Gamma with warm-cached constants: " << std::setprecision(7) << (static_cast<double>(t1_second) - t0_second) / CLOCKS_PER_SEC << std::endl;
}
