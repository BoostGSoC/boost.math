
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_BERNOULLI_2013_05_30_HPP_
  #define _BOOST_BERNOULLI_2013_05_30_HPP_

  #include <boost/array.hpp>
  #include <boost/cstdint.hpp>
  #include "detail/bernoulli_b2n.hpp"

  namespace boost { namespace math {

  template <class T, class Policy>
  inline T bernoulli_b2n(const unsigned i, const Policy &pol)
  {
    const unsigned i_2 = 2 * i;

    return boost::math::detail::bernoulli_number_imp<T,Policy>(i_2, pol);
  }

  template <class T>
  inline T bernoulli_b2n(const unsigned i)
  {
    return boost::math::bernoulli_b2n<T>(i, policies::policy<>());
  }

  template <class T, class OutputIterator, class Policy>
  inline OutputIterator bernoulli_b2n(unsigned start_index,
                                      unsigned number_of_bernoulis_b2n,
                                      OutputIterator out_it,
                                      const Policy& pol)
  {
    return boost::math::detail::bernoulli_series_imp<T, OutputIterator, Policy>(start_index,
                                                                                number_of_bernoulis_b2n,
                                                                                out_it,
                                                                                policies::policy<>());
  }

  template <class T, class OutputIterator>
  inline OutputIterator bernoulli_b2n(unsigned start_index,
                                      unsigned number_of_bernoulis_b2n,
                                      OutputIterator out_it)
  {
    return boost::math::bernoulli_b2n<T, OutputIterator>(start_index,
                                                         number_of_bernoulis_b2n,
                                                         out_it,
                                                         policies::policy<>());
  }



} } // namespace boost::math

#endif // _BOOST_BERNOULLI_2013_05_30_HPP_
