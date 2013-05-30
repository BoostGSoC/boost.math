
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_BERNOULLI_B2N_2013_05_30_HPP_
  #define _BOOST_BERNOULLI_B2N_2013_05_30_HPP_

  #include <limits>
  #include <vector>
  #include <cmath>
  #include <boost/cstdint.hpp>
  #include <boost/math/special_functions/pow.hpp>
  #include <boost/math/policies/policy.hpp>
  #include <boost/static_assert.hpp>

  struct pol{};

  namespace boost { namespace math { namespace detail {

  template <class T>
  T tangent_numbers(const unsigned nn)
  {
    if(nn % 2)
    {
      T x(0U);
      return x;
    }

    const std::int32_t m = nn / 2;

    std::vector<T> tangent_numbers(m + 1);

    tangent_numbers[0U] = T(0U);
    tangent_numbers[1U] = T(1U);

    T power_four(4);
    T power_two(2);

    for(boost::int32_t k = 2; k <= m; k++)
    {
      power_four*=4;
      power_two*=2;
      tangent_numbers[k] = (k - 1) * tangent_numbers[k - 1];
    }

    power_four *= power_four;
    power_two *= power_two;

    for(boost::int32_t k = 2; k <= m; k++)
    {
      for(boost::int32_t j = k; j <= m; j++)
      {
        tangent_numbers[j] = (tangent_numbers[j - 1] * (j - k)) + (tangent_numbers[j] * (j - k + 2));
      }
    }

    T x(nn);
    x = x / (power_four - power_two);
    x *= tangent_numbers[nn / 2];

    return (nn % 4 == 0) ? -x : x;
  }

  template <class T>
  void tangent_numbers_series(std::vector<T>& bn, const unsigned m)
  {
    std::vector<T> tangent_numbers(m + 1);

    tangent_numbers[0U] = T(0U);
    tangent_numbers[1U] = T(1u);

    for(boost::int32_t k = 2; k <= m; k++)
    {
      tangent_numbers[k] = (k - 1) * tangent_numbers[k - 1];
    }

    for(boost::int32_t k = 2; k <= m; k++)
    {
      for(boost::int32_t j = k; j <= m; j++)
      {
        tangent_numbers[j] = (tangent_numbers[j - 1] * (j - k)) + (tangent_numbers[j] * (j - k + 2));
      }
    }

    T two_pow_two_m(4);

    bn.clear();
    bn.resize(m + 1);

    for(std::int32_t i = 1; i < static_cast<std::int32_t>(tangent_numbers.size()); i++)
    {
      const std::int32_t two_i = static_cast<std::int32_t>(static_cast<std::int32_t>(2) * i);

      const T b = (tangent_numbers[i] * two_i) / (two_pow_two_m * (two_pow_two_m - 1));

      const bool  b_neg = (static_cast<std::int32_t>(two_i % static_cast<std::int32_t>(4)) == static_cast<std::int32_t>(0));

      bn[i] = ((!b_neg) ? b : -b);

      two_pow_two_m *= 4;
    }

    bn[0U] =  T(1U);
  }

  template <class T, class Policy>
  T bernoulli_number_imp(const unsigned n, const Policy &pol)
  {
    if((n/2)<=max_bernoulli_index<T>::value)
    {
      return unchecked_bernoulli_b2n<T>(n/2);
    }
    else
    {
      return tangent_numbers<T>(n);
    }
  }

  template <class T, class OutputIterator, class Policy>
  inline OutputIterator bernoulli_series_imp(unsigned start_index,
                                      unsigned number_of_bernoulis_bn,
                                      OutputIterator out_it,
                                      const Policy& pol)
  {
    if((start_index + number_of_bernoulis_bn) < max_bernoulli_index<T>::value)
    {
      OutputIterator last= out_it + number_of_bernoulis_bn;

      while(out_it!=last)
      {
        *out_it=unchecked_bernoulli_b2n<T>(start_index);
        ++out_it;
        ++start_index;
      }

      //return one past the last element
      return out_it;
    }

    std::vector<T> bn;

    tangent_numbers_series(bn, start_index + number_of_bernoulis_bn);

    OutputIterator last = out_it + number_of_bernoulis_bn;

    while(out_it != last)
    {
      *out_it=bn[start_index];
      ++out_it;
      ++start_index;
    }

    return out_it;
  }

} } } // namespace boost::math::detail

#endif // _BOOST_BERNOULLI_B2N_2013_05_30_HPP_
