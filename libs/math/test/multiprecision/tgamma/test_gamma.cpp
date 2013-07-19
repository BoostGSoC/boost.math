//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "test_gamma.hpp"
#include <boost/multiprecision/mpfr.hpp>


void expected_results()
{
   //
   // Define the max and mean errors expected for
   // various compilers and platforms.
   //
   const char* largest_type;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   if(boost::math::policies::digits<double, boost::math::policies::policy<> >() == boost::math::policies::digits<long double, boost::math::policies::policy<> >())
   {
      largest_type = "(long\\s+)?double";
   }
   else
   {
      largest_type = "long double";
   }
#else
   largest_type = "(long\\s+)?double";
#endif
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                // test type(s)
      ".*",                // test data group
      ".*", 20, 5);        // test function

   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

BOOST_AUTO_TEST_CASE( test_main )
{
   expected_results();
   BOOST_MATH_CONTROL_FP;
#ifndef DIGITS10
#error "No precision specified."
#endif

   typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<DIGITS10> > mp_type;

   std::cout << "Lanczos type in use is: " << typeid(boost::math::lanczos::lanczos<mp_type, boost::math::policies::policy<> >::type).name() << std::endl;
   std::cout << "Precision is: " << DIGITS10 << " decimal digits" << std::endl;

   test_gamma(mp_type(0), "number<mpfr_float_backend<" BOOST_STRINGIZE(DIGITS10) "> >");
}



