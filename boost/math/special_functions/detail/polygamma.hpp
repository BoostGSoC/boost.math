
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_POLYGAMMA_2013_07_30_CPP_
  #define _BOOST_POLYGAMMA_2013_07_30_CPP_

  #include <limits>
  #include <vector>
  #include <cmath>
  #include <boost/cstdint.hpp>
  #include <boost/math/special_functions/trunc.hpp>
  #include <boost/math/special_functions/zeta.hpp>
  #include <boost/math/policies/policy.hpp>
  #include <boost/static_assert.hpp>
  #include "../bernoulli.hpp"
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/type_traits/is_convertible.hpp>
  #include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::multiprecision;
using std::size_t;


namespace boost { namespace math { namespace detail {

  template<class T>
  struct max_iteration
  {
     //TODO Derive a suitable formula based on precision of T
     static const int value=1000;
  };

  template<class T, class Policy>
  T digamma_atinfinityplus(const int n, const T &x, const Policy &pol)
  {

	  BOOST_MATH_STD_USING

	  // calculate a high bernoulli number upfront to make use of cache
	    unsigned int bernoulli_index = 100;
	    boost::math::bernoulli_b2n<T>(bernoulli_index);
        T z=x;
        T log_z=T(log(z));
        T one_over_2z= T(1)/(2*z);
        T sum=T(0);

        for(int two_k=2; two_k < max_iteration<T>::value; two_k+=2)
        {
                if(two_k/2 > bernoulli_index)
                {
                    try
                    {
                        int temp = bernoulli_index * 1.5;
                        boost::math::bernoulli_b2n<T>(temp);
                        bernoulli_index = temp;
                    }
                    catch(...)
                    {
                        break;
                    }
                }
                T term=T(1);
                T one_over_two_k=T(1)/two_k;
                T z_pow_two_k=pow(z,static_cast<boost::int32_t>(two_k));
                T one_over_z_pow_two_k = T(1)/z_pow_two_k;
                T bernoulli_term= boost::math::bernoulli_b2n<T>(two_k/2);

                term = bernoulli_term * one_over_two_k * one_over_z_pow_two_k;

                if(term == 0 ) continue;

                sum += term;

                T term_base_10_exp = term < 0 ? -term: term;
                T sum_base_10_exp  = sum < 0 ? -sum: sum;

                int ll;

                T significand=frexp(term_base_10_exp,&ll);
                term_base_10_exp=ll*0.303;

                significand=frexp(sum_base_10_exp,&ll);
                sum_base_10_exp=ll*0.303;

                long int order_check =  boost::math::ltrunc(term_base_10_exp)-boost::math::ltrunc(sum_base_10_exp);
                long int tol         =  std::numeric_limits<T>::digits10;


                if((two_k > static_cast<boost::int32_t>(50)) && (order_check < -tol))
                {
                    break;
                }
        }

        T answer = log_z - one_over_2z -sum;

        return answer;
  }

  template<class T, class Policy>
  T polygamma_atinfinityplus(const int n, const T &x, const Policy &pol) // for large values of x such as for x> 400
  {
     BOOST_MATH_STD_USING

     if(n==0)
        return digamma_atinfinityplus(n,x,pol);

     unsigned int bernoulli_index = 100;
     boost::math::bernoulli_b2n<T>(bernoulli_index);

     const bool b_negate = (n % 2 == 0);

     const T n_minus_one_fact            = boost::math::factorial<T>(n - 1);
     const T nn                          = T(n);
     const T n_fact                      = n_minus_one_fact * nn;
     const T one_over_z                  = T(1)/ x;
     const T one_over_z2                 = one_over_z * one_over_z;
     const T one_over_z_pow_n            = T(1) / pow(x, static_cast<boost::int64_t>(n));
           T one_over_x_pow_two_k_plus_n = one_over_z_pow_n * one_over_z2;
           T two_k_plus_n_minus_one      = nn + T(1);
           T two_k_plus_n_minus_one_fact = n_fact * T(static_cast<boost::uint64_t>(n + 1)); //(n+3)! ?
           T one_over_two_k_fact         = T(1)/2;
           T sum                         = (  (boost::math::bernoulli_b2n<T>(static_cast<boost::uint64_t>(1u)) * two_k_plus_n_minus_one_fact)
                                                  * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

     // Perform the Bernoulli series expansion.
     for(boost::int32_t two_k = 4; two_k < max_iteration<T>::value; two_k += 2)
     {
        if(two_k/2 > bernoulli_index)
        {
           try
           {
              int temp = bernoulli_index * 1.5;
              boost::math::bernoulli_b2n<T>(temp);
              bernoulli_index = temp;
            }
            catch(...)
            {
               break;
            }
        }
        one_over_x_pow_two_k_plus_n *= one_over_z2;
        two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
        two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
        one_over_two_k_fact         /= static_cast<boost::int32_t>(two_k * static_cast<boost::int32_t>(two_k - static_cast<boost::int32_t>(1)));

        const T term = (  (boost::math::bernoulli_b2n<T>(two_k/2) * two_k_plus_n_minus_one_fact)
                              * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

        if(term == 0 ) continue;

	    sum += term;

	    T term_base_10_exp = term < 0 ? -term: term;
	    T sum_base_10_exp  = sum < 0 ? -sum: sum;

	    int ll;

	    T significand=frexp(term_base_10_exp,&ll);
	    term_base_10_exp=ll*0.303;

	    significand=frexp(sum_base_10_exp,&ll);
	    sum_base_10_exp=ll*0.303;

	    long int order_check =  boost::math::ltrunc(term_base_10_exp)-boost::math::ltrunc(sum_base_10_exp);
	    long int tol         =  std::numeric_limits<T>::digits10;


        if((two_k > static_cast<boost::int32_t>(50)) && (order_check < -tol))
        {
          break;
        }

     }

     sum += ((((n_minus_one_fact * (nn + (x * static_cast<boost::int32_t>(2)))) * one_over_z_pow_n) * one_over_z) / 2);

     return !b_negate ? sum : -sum;
  }

  template<class T, class Policy>
  T polygamma_attransitionplus(const int n, const T &x, const Policy &pol)
  {
  // this doesn't work for digamma either

      // Use Euler-Maclaurin summation.

    // Use N = (0.4 * digits) + (4 * n)
    BOOST_MATH_STD_USING
    static const boost::int64_t prec = static_cast<boost::int64_t>(std::numeric_limits<T>::digits10);
    static const boost::int64_t d4d  = static_cast<boost::int64_t>(double(0.4) * prec);
           const boost::int64_t N4dn = static_cast<boost::int64_t>(d4d + 4 * n);
    static const boost::int64_t n32m = static_cast<boost::int64_t>(std::numeric_limits<boost::int32_t>::max());
           const boost::int32_t N    = static_cast<boost::int32_t>(std::min(N4dn, n32m));
           const boost::int32_t m    = n;

    const boost::int64_t minus_m_minus_one = static_cast<boost::int64_t>(static_cast<boost::int32_t>(-m)
                                                                        - static_cast<boost::int32_t>(1));

    T z                              = x;
    T sum0                           = T(0);
    T z_plus_k_pow_minus_m_minus_one = T(0);

    for(boost::int32_t k = 1; k <= N; k++)
    {
      z_plus_k_pow_minus_m_minus_one = pow(z, minus_m_minus_one);

      sum0 += z_plus_k_pow_minus_m_minus_one;

      ++z;
    }

    const T one_over_z_plus_N_pow_minus_m           = pow(z, static_cast<boost::int64_t>(-m));
    const T one_over_z_plus_N_pow_minus_m_minus_one = one_over_z_plus_N_pow_minus_m / z;

    const T term0 = one_over_z_plus_N_pow_minus_m_minus_one / static_cast<boost::int32_t>(2);
    const T term1 = one_over_z_plus_N_pow_minus_m           / m;

          T sum1                                      = T(0);
          T one_over_two_k_fact                       = T(1)/2;
          boost::int32_t   mk                         = m + 1;
          T am                                        = T(mk);
    const T one_over_z_plus_N_squared                 = T(1) / (z * z);
          T one_over_z_plus_N_pow_minus_m_minus_two_k = one_over_z_plus_N_pow_minus_m * one_over_z_plus_N_squared;

    for(boost::int32_t k = 1; k < max_iteration<T>::value; k++)
    {
      const boost::int32_t two_k = 2 * k; // k << 1

      const T term = ((boost::math::bernoulli_b2n<T>(two_k/2) * am) * one_over_two_k_fact) * one_over_z_plus_N_pow_minus_m_minus_two_k;

      /*const boost::int64_t order_check = term.order() - sum1.order();*/

      //TODO Devise a good breaking condition

      if(k > 80 /*&& (order_check < -ef::tol())*/)
      {
        break;
      }

      sum1 += term;

      one_over_two_k_fact /= static_cast<boost::int32_t>(two_k + static_cast<boost::int32_t>(1));
      one_over_two_k_fact /= static_cast<boost::int32_t>(two_k + static_cast<boost::int32_t>(2));

      am *= static_cast<boost::int32_t>(++mk);
      am *= static_cast<boost::int32_t>(++mk);

      one_over_z_plus_N_pow_minus_m_minus_two_k *= one_over_z_plus_N_squared;
    }

    const T pg = (((sum0 + term0) + term1) + sum1) * factorial<T>(static_cast<boost::uint32_t>(m));

    const bool b_negate = (m % 2 == 0);

    return (!b_negate ? pg : -pg);
  }

  template<class T, class Policy>
  T polygamma_nearzero(const int n, const T &x, const Policy &pol)
  {
    BOOST_MATH_STD_USING
    // not defined for digamma

    // Use a series expansion for x near zero which uses poly_gamma(m, 1) which,
    // in turn, uses the Riemann zeta function for integer arguments.
    // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/01/03/01/02/
    const bool b_negate = (( n % 2 ) == 0 ) ;

    const T n_fact               =  boost::math::factorial<T>(n);
    const T z_pow_n_plus_one     =  pow(x, static_cast<boost::int64_t>(n + 1));
    const T n_fact_over_pow_term =  n_fact / z_pow_n_plus_one;
    const T term0                =  !b_negate ? n_fact_over_pow_term : -n_fact_over_pow_term;

          T one_over_k_fact      =  T(1);
          T z_pow_k              =  T(1);
          T k_plus_n_fact        =  boost::math::factorial<T>(n);
          T k_plus_n_plus_one    =  T(n + 1);
    const T pg_kn                =  k_plus_n_fact * boost::math::zeta<T>(k_plus_n_plus_one);
          bool    b_neg_term     =  ((n % 2) == 0);
          T sum                  =  !b_neg_term ? pg_kn : -pg_kn;

    for(int k = 1; k < max_iteration<T>::value; k++)
    {
      k_plus_n_fact   *= k_plus_n_plus_one++;
      one_over_k_fact /= k;
      z_pow_k         *= x;

      const T pg = k_plus_n_fact * boost::math::zeta<T>(k_plus_n_plus_one);

      const T term = (pg * z_pow_k) * one_over_k_fact;

      //const boost::int64_t order_check = static_cast<boost::int64_t>(term.order() - sum.order());

      //TODO devise a good breaking condition
      if(k > 500 /*&& (order_check < -ef::tol())*/)
      {
        break;
      }

      b_neg_term = !b_neg_term;

      !b_neg_term ? sum += term : sum -= term;
    }

    return term0 + sum;

  }

  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
//	  std::cout<<typeid(T).name()<<std::endl;
        if(x < 0.5)
            return polygamma_nearzero(n,x,pol);
        else if (x > 400)
            return polygamma_atinfinityplus(n,x,pol); //just a test return value
        else
            return polygamma_attransitionplus(n,x,pol);
  }


} } } // namespace boost::math::detail

#endif // __BOOST_POLYGAMMA_2013_07_30_CPP_
