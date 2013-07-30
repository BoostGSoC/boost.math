
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_POLYGAMMA_2013_07_30_HPP_
  #define _BOOST_POLYGAMMA_2013_07_30_HPP_

  #include <boost/array.hpp>
  #include <boost/cstdint.hpp>
  #include "detail/polygamma.hpp"

  namespace boost { namespace math {

  template<class T, class Policy>
  inline T polygamma(const int n, T x, const Policy &pol)
  {
        return boost::math::detail::polygamma_imp(n,x,pol);
  }

  template<class T>
  inline T polygamma(const int n, T x)
  {
      return boost::math::polygamma(n,x,policies::policy<>());
  }

  template<class T, class Policy>
  inline T digamma(T x, const Policy &pol)
  {
      return boost::math::polygamma(2,x,pol);
  }

  template<class T>
  inline digamma(T x)
  {
      return boost::math::digamma(x,policies::policy<>());
  }

  template<class T, class Policy>
  inline T trigamma(T x, const Policy &pol)
  {
      return boost::math::polygamma(3,x,pol);
  }

  template<class T>
  inline trigamma(T x)
  {
      return boost::math::trigamma(x,policies::policy<>());
  }


} } // namespace boost::math

#endif // _BOOST_BERNOULLI_2013_05_30_HPP_
