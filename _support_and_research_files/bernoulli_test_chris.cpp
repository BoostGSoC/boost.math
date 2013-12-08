#define BOOST_TEST_MODULE my_bernoulli_test

#include <vector>
#include <iostream>
#include <typeinfo>

#include <boost/math/concepts/real_concept.hpp> // for real_concept

#include <boost/math/special_functions.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "bernoulli.hpp"
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>, boost::multiprecision::et_off> cpp_dec_float_50_noet;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 9>, boost::multiprecision::et_off> cpp_dec_float_9_noet;

template <class RealType>
void generic_bernoulli_test()
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

  BOOST_CHECK_CLOSE_FRACTION(bn[0], static_cast<RealType>(-1.9296579341940068148632668144863266814486e16L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[1], static_cast<RealType>(8.4169304757368261500055370985603543743079e17L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[2], static_cast<RealType>(-4.0338071854059455413076811594202898550725e19L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[3], static_cast<RealType>(2.1150748638081991605601453900709219858156e21L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[4], static_cast<RealType>(-1.2086626522296525934602731193708252531782e23L), tolerance);
}

template <>
void generic_bernoulli_test<cpp_dec_float_50_noet>()
{
  static const cpp_dec_float_50_noet tolerance = std::numeric_limits<cpp_dec_float_50_noet>::epsilon() * 100;

  std::cout << "Tolerance for type " << typeid(cpp_dec_float_50_noet).name()  << " is "
    << std::setprecision(3) << tolerance  << " (or " << tolerance * 100 << "%)." << std::endl;

  // test for negative argument
  BOOST_CHECK_THROW(boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(-1), std::domain_error);

  //should call unchecked bernoulli
  BOOST_CHECK_CLOSE_FRACTION(boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(12),
                             cpp_dec_float_50_noet("-86580.253113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553"),
                             tolerance);

  //should be calculated by tangent numbers algorithm
  BOOST_CHECK_CLOSE_FRACTION(boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(500),
                             cpp_dec_float_50_noet("-5.3187044694155220364829137437670855452093660052166637474739092287348819200049898438096886914613452320729737902e+1769"),
                             tolerance);

  std::vector<cpp_dec_float_50_noet> bn(5);
  std::vector<cpp_dec_float_50_noet>::iterator it = bn.begin();

  // test for negative argument
  BOOST_CHECK_THROW(boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(-1,5,it), std::domain_error);

  //series retrieved from unchecked bernoulli
  boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(10, 5, it);

  BOOST_CHECK_CLOSE_FRACTION(bn[0U],
                             cpp_dec_float_50_noet("-529.12424242424242424242424242424242424242424242424242424242424242424242424242424242424242424242424242424242424"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[1U],
                             cpp_dec_float_50_noet("6192.1231884057971014492753623188405797101449275362318840579710144927536231884057971014492753623188405797101449"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[2U],
                             cpp_dec_float_50_noet("-86580.253113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553113553"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[3U],
                             cpp_dec_float_50_noet("1425517.1666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[4U],
                             cpp_dec_float_50_noet("-27298231.067816091954022988505747126436781609195402298850574712643678160919540229885057471264367816091954022989"),
                             tolerance);

  //series retrieved from tangent numbers algorithm
  boost::math::bernoulli_b2n<cpp_dec_float_50_noet>(50, 5, it);

  BOOST_CHECK_CLOSE_FRACTION(bn[0U],
                             cpp_dec_float_50_noet("-2.8382249570693706959264156336481764738284680928012882128228531714464865111070281341434143414341434143414341434e+78"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[1U],
                             cpp_dec_float_50_noet("7.4064248979678850629750827140920984176879731788088706673116100348748532844121085501410078594544613962089690245e+80"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[2U],
                             cpp_dec_float_50_noet("-2.0096454802756604483465619672715363186867270822532876624346130198921356500977969888305220125786163522012578616e+83"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[3U],
                             cpp_dec_float_50_noet("5.6657170050805941445719346030519356961419468287510420621387564452152460861972277798400157320872274143302180685e+85"),
                             tolerance);

  BOOST_CHECK_CLOSE_FRACTION(bn[4U],
                             cpp_dec_float_50_noet("-1.6584511154136216915823713374319912301494962614725464727402466815589878137712650743149939341946471014554066220e+88"),
                             tolerance);
}

BOOST_AUTO_TEST_CASE(bernoulli_test)
{
  generic_bernoulli_test<float>();
  generic_bernoulli_test<double>();
  generic_bernoulli_test<long double>();
  generic_bernoulli_test<cpp_dec_float_9_noet>();
  generic_bernoulli_test<cpp_dec_float_50_noet>();
}
