#include <cstdint>
#include <ctime>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>

/*
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
  // TBD: small argument approximation?
  // TBD: NaN for zero or very near zero.
  // TBD: +infinity for negative integers or very close thereto.

  if(!boost::math::isfinite(x))
  {
    return x;
  }

  const bool b_neg = (x < 0);

  // Make a local, unsigned copy of the input argument.

  mp_type xx((!b_neg) ? x : -x);

  // Check if the argument should be scaled up for the Bernoulli series expansion.
  static const double digit_scale_of_argument = static_cast<double>(std::numeric_limits<mp_type>::digits10 * 1.7);

  static const std::int32_t min_arg_n = static_cast<std::int32_t>(digit_scale_of_argument);
  static const mp_type   min_arg_x = mp_type(min_arg_n);

  const std::int32_t n_recur = ((xx < min_arg_x) ? static_cast<std::int32_t>((min_arg_n - xx.convert_to<std::int32_t>()) + 1)
                                                 : static_cast<std::int32_t>(0));

  // Scale the argument up and use downward recursion later for the final result.
  if(n_recur != static_cast<std::int32_t>(0))
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
  for(std::int32_t k = static_cast<std::int32_t>(0); k < n_recur; k++)
  {
    xx -= 1;
    gamma_value /= xx;
  }

  // Return the result, accounting for possible negative arguments.
  return ((!b_neg) ? gamma_value : -boost::math::constants::pi<mp_type>() / (xx * gamma_value * sin(boost::math::constants::pi<mp_type>() * xx)));
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
*/

#define BOOST_TEST_MODULE MyTest

#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "bernoulli.hpp"

using namespace boost::test_tools;
typedef boost::multiprecision::number<cpp_dec_float<50>, boost::multiprecision::et_off> cpp_dec_float_50_et_off;

template <class RealType>
void generic_test_built_in()
{
  RealType tolerance = (std::max)
    (boost::math::tools::epsilon<RealType>(),
    static_cast<RealType>(std::numeric_limits<double>::epsilon()));
  tolerance *= 100;

  std::cout << "Tolerance for type " << typeid(RealType).name()  << " is "
    << std::setprecision(3) << tolerance  << " (or " << tolerance * 100 << "%)." << std::endl;

  BOOST_CHECK_CLOSE_FRACTION(boost::math::bernoulli_b2n<RealType>(12),
                             static_cast<RealType>(-86580.253113553113553113553113553113553114L),
                             tolerance);

  typename std::vector<RealType> bn(5);
  typename std::vector<RealType>::iterator it=bn.begin();

  boost::math::bernoulli_b2n<RealType>(20,5,it);

  BOOST_CHECK_CLOSE_FRACTION(bn[0],static_cast<RealType>(-1.9296579341940068148632668144863266814486e16L),tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[1],static_cast<RealType>(8.4169304757368261500055370985603543743079e17L),tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[2],static_cast<RealType>(-4.0338071854059455413076811594202898550725e19L),tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[3],static_cast<RealType>(2.1150748638081991605601453900709219858156e21L),tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[4],static_cast<RealType>(-1.2086626522296525934602731193708252531782e23L),tolerance);
}

template <>
void generic_test_built_in<cpp_dec_float_50_et_off>()
{
  std::cout << "Test cpp_dec_float<50>" << std::endl;

  cpp_dec_float_50_et_off tolerance=boost::math::tools::epsilon<cpp_dec_float_50_et_off>() * 100;

   std::cout << "Tolerance is "
    << std::setprecision(3) << tolerance  << " (or " << tolerance * 100 << "%)." << std::endl;

  output_test_stream output;

  // TBD: Use BOOST_CHECK_CLOSE_FRACTION in all of these tests.
  // TBD: Weird compiler error. Chris should ask John.

  //should call unchecked bernoulli


//*** the following test fails
  BOOST_CHECK_CLOSE_FRACTION(boost::math::bernoulli_b2n<cpp_dec_float_50_et_off>(12),
                             static_cast<cpp_dec_float_50_et_off>(-86580.2531135531135531135531135531135531135531136),
                             tolerance);

  output << std::setprecision(std::numeric_limits<cpp_dec_float_50_et_off>::digits10 - 2)
         << boost::math::bernoulli_b2n<cpp_dec_float_50_et_off>(12);
  BOOST_CHECK(output.is_equal("-86580.2531135531135531135531135531135531135531136"));

  //should be calculated by tangent numbers algorithm


//*** compiler returns a warning warning: floating constant exceeds range of ‘double’
//*** therefore not converting any more tests to close fraction for the moment
  cpp_dec_float_50_et_off x(-5.31870446941552203648291374376708554520936600522e+1769);
  BOOST_CHECK_CLOSE_FRACTION(boost::math::bernoulli_b2n<cpp_dec_float_50_et_off>(500),
                             x,
                             tolerance);

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << boost::math::bernoulli_b2n<boost::multiprecision::cpp_dec_float_50>(500);
  BOOST_CHECK(output.is_equal("-5.31870446941552203648291374376708554520936600522e+1769"));

  std::vector<boost::multiprecision::cpp_dec_float_50> bn(5);
  std::vector<boost::multiprecision::cpp_dec_float_50>::iterator it = bn.begin();

  //series retrieved from unchecked bernoulli
  boost::math::bernoulli_b2n<boost::multiprecision::cpp_dec_float_50>(10, 5, it);

  output<<std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)<<bn[0];
  BOOST_CHECK(output.is_equal("-529.124242424242424242424242424242424242424242424"));

  output<<std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)<<bn[1];
  BOOST_CHECK(output.is_equal("6192.12318840579710144927536231884057971014492754"));

  output<<std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)<<bn[2];
  BOOST_CHECK(output.is_equal("-86580.2531135531135531135531135531135531135531136"));

  output<<std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)<<bn[3];
  BOOST_CHECK(output.is_equal("1425517.16666666666666666666666666666666666666667"));

  output<<std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)<<bn[4];
  BOOST_CHECK(output.is_equal("-27298231.0678160919540229885057471264367816091954"));

  //series retrieved from tangent numbers algorithm
  boost::math::bernoulli_b2n<boost::multiprecision::cpp_dec_float_50>(50, 5, it);

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << bn[0];
  BOOST_CHECK(output.is_equal("-2.8382249570693706959264156336481764738284680928e+78"));

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << bn[1];
  BOOST_CHECK(output.is_equal("7.40642489796788506297508271409209841768797317881e+80"));

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << bn[2];
  BOOST_CHECK(output.is_equal("-2.00964548027566044834656196727153631868672708225e+83"));

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << bn[3];
  BOOST_CHECK(output.is_equal("5.66571700508059414457193460305193569614194682875e+85"));

  output << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10 - 2)
         << bn[4];
  BOOST_CHECK(output.is_equal("-1.65845111541362169158237133743199123014949626147e+88"));
}

/*template <>
void generic_test_built_in<boost::multiprecision::number<cpp_dec_float<9> > >()
{
  std::cout << "Test cpp_dec_float<9>" << std::endl;

  output_test_stream output;

  // TBD: Use BOOST_CHECK_CLOSE_FRACTION in all of these tests.
  // TBD: Weird compiler error. Chris should ask John.

  //should call unchecked bernoulli
  output << std::setprecision(std::numeric_limits<boost::multiprecision::number<cpp_dec_float<9> > >::digits10 - 2)
         << boost::math::bernoulli_b2n<boost::multiprecision::number<cpp_dec_float<9> > >(12);
  BOOST_CHECK(output.is_equal("-86580.25"));

  //should be calculated by tangent numbers algorithm
  output << std::setprecision(std::numeric_limits<boost::multiprecision::number<cpp_dec_float<9> > >::digits10 - 2)
         << boost::math::bernoulli_b2n<boost::multiprecision::number<cpp_dec_float<9> > >(500);
  BOOST_CHECK(output.is_equal("-5.318704e+1769"));

}*/

BOOST_AUTO_TEST_CASE(generic_built_in_test)
{
  generic_test_built_in<float>();
  generic_test_built_in<double>();
  generic_test_built_in<long double>();
  generic_test_built_in<cpp_dec_float_50_et_off>();
  generic_test_built_in<boost::multiprecision::number<cpp_dec_float<9>, boost::multiprecision::et_off> >();
}
