
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
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/type_traits/is_convertible.hpp>
  #include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::multiprecision;
using std::size_t;


namespace boost { namespace math { namespace detail {

  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
        return x; //just a test return value
  }


} } } // namespace boost::math::detail

#endif // _BOOST_BERNOULLI_B2N_2013_05_30_HPP_
