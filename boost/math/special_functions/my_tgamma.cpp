
#include <cstdint>
#include <ctime>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>

// Investigate high precision calculations of the tgamma function with
// Bernoulli numbers (Stirling's approximation) and Boost.Multiprecision.

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000>,
                                      boost::multiprecision::et_off> mp_type;

namespace
{
  const std::size_t highest_bernoulli_index()
  {
    // Find the high index n for Bn to produce the desired precision in Stirling's calculation.
    // TBD: This max-index is easy to template when considering a generic tgamma.
    return static_cast<std::size_t>(18.0 + (0.6 * static_cast<double>(std::numeric_limits<mp_type>::digits10)));
  }
}

namespace my
{
  const mp_type& zero()
  {
    static mp_type z0(0);

    return z0;
  }

  mp_type tgamma(const mp_type& x);

  const mp_type& bernoulli_table(const std::uint32_t n)
  {
    static std::vector<mp_type> bernoulli_b2n_data;

    if(bernoulli_b2n_data.empty())
    {
      bernoulli_b2n_data.resize(highest_bernoulli_index());

      boost::math::bernoulli_b2n<mp_type>(0U,
                                          unsigned(bernoulli_b2n_data.size()),
                                          bernoulli_b2n_data.begin());
    }

    return ((n < bernoulli_b2n_data.size()) ? bernoulli_b2n_data[n] : my::zero());
  }
}

mp_type my::tgamma(const mp_type& x)
{
  // Make a local copy of the input argument.

  mp_type xx(x);

  // Check if the argument should be scaled up for the Bernoulli series expansion.

  static const mp_type min_arg_for_recursion(float(std::numeric_limits<mp_type>::digits10 * 1.7F));

  const mp_type n_recur = ((xx < min_arg_for_recursion) ? (floor(min_arg_for_recursion - xx) + 1)
                                                        : mp_type(0));

  // Scale the argument up and use downward recursion later for the final result.
  if(n_recur != 0)
  {
    xx += n_recur;
  }

        mp_type one_over_x_pow_two_n_minus_one = 1 / xx;
  const mp_type one_over_x2                    = one_over_x_pow_two_n_minus_one * one_over_x_pow_two_n_minus_one;
        mp_type sum                            = (bernoulli_table(1) / static_cast<std::int32_t>(2)) * one_over_x_pow_two_n_minus_one;

  // Perform the Bernoulli series expansion of Stirling's approximation.
  for(std::int32_t n2 = static_cast<std::int32_t>(4); n2 < static_cast<std::int32_t>(highest_bernoulli_index()); n2 += static_cast<std::int32_t>(2))
  {
    one_over_x_pow_two_n_minus_one *= one_over_x2;

    const mp_type term = (bernoulli_table(static_cast<std::uint32_t>(n2 / 2)) * one_over_x_pow_two_n_minus_one) / static_cast<std::int32_t>(n2 * (n2 - static_cast<std::int32_t>(1)));

    sum += term;
  }

  static const mp_type half_ln_two_pi = log(boost::math::constants::two_pi<mp_type>()) / static_cast<std::int32_t>(2);

  const mp_type log_gamma_value = ((((xx - boost::math::constants::half<mp_type>()) * log(xx)) - xx) + half_ln_two_pi) + sum;

  mp_type gamma_value = exp(log_gamma_value);

  // Rescale the result using downward recursion if necessary.
  for(mp_type kx = mp_type(0); kx < n_recur; ++kx)
  {
    xx -= 1;
    gamma_value /= xx;
  }

  // Return the result, accounting for possible negative arguments.
  return gamma_value;
}

int main(int, char**)
{
  ::clock_t t0_first = ::clock();
  const mp_type g_first = my::tgamma(boost::math::constants::half<mp_type>());
  ::clock_t t1_first = ::clock();

  std::cout << std::setprecision(std::numeric_limits<mp_type>::digits10)
            << g_first
            << std::endl;

  ::clock_t t0_second = ::clock();
  const mp_type g_second = my::tgamma(boost::math::constants::half<mp_type>());
  ::clock_t t1_second = ::clock();

  std::cout << std::setprecision(std::numeric_limits<mp_type>::digits10)
            << g_second
            << std::endl;

  // Test with tgamma(1/2) == sqrt(pi).
  std::cout << std::setprecision(std::numeric_limits<mp_type>::digits10)
            << sqrt(boost::math::constants::pi<mp_type>())
            << std::endl;

  std::cout << "Time: Gamma with Bernoulli generation:  " << std::setprecision(7) << (static_cast<double>(t1_first)  - t0_first ) / CLOCKS_PER_SEC << std::endl;
  std::cout << "Time: Gamma with warm-cached constants: " << std::setprecision(7) << (static_cast<double>(t1_second) - t0_second) / CLOCKS_PER_SEC << std::endl;
}
