
#include <iostream>
#include <iomanip>
#include <limits>

#include <boost/cstdint.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>

// Investigate high precision calculations.

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>,
                                      boost::multiprecision::et_off> mp_type;

                                      template<typename float_type>
void test_order_17_argument_pi_over_11()
{
  const float_type pg = boost::math::polygamma(17, boost::math::constants::pi<float_type>() / 11);

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << pg
            << std::endl;
}

template<typename float_type>
void test_order_17_argument_pi_plus_12()
{
  const float_type pg = boost::math::polygamma(17, 12 + boost::math::constants::pi<float_type>());

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << pg
            << std::endl;
}

template<typename float_type>
void test_order_17_argument_pi_plus_401()
{
  const float_type pg = boost::math::polygamma(17, 401 + boost::math::constants::pi<float_type>());

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << pg
            << std::endl;
}

int main(int, char**)
{
  test_order_17_argument_pi_over_11<float>();
  test_order_17_argument_pi_over_11<double>();
  test_order_17_argument_pi_over_11<mp_type>();

  test_order_17_argument_pi_plus_12<float>();
  test_order_17_argument_pi_plus_12<double>();
  test_order_17_argument_pi_plus_12<mp_type>();

  test_order_17_argument_pi_plus_401<float>();
  test_order_17_argument_pi_plus_401<double>();
  test_order_17_argument_pi_plus_401<mp_type>();
}
