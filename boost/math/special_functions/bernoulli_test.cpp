#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace boost::test_tools;

template <class RealType>
void generic_test_built_in(RealType)
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

BOOST_AUTO_TEST_CASE(generic_built_in_test)
{
    generic_test_built_in(0.0F);
    generic_test_built_in(0.0);
    generic_test_built_in(0.0L);
}

BOOST_AUTO_TEST_CASE( cpp_dec_float_50_test )
{
  //should call unchecked bernoulli
  output_test_stream output;
  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<boost::math::bernoulli_b2n<cpp_dec_float_50>(12);
  BOOST_CHECK(output.is_equal("-86580.253113553113553113553113553113553113553113553"));

  //should be calculated by tangent numbers algorithm

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<boost::math::bernoulli_b2n<cpp_dec_float_50>(100);
  BOOST_CHECK(output.is_equal("-3.6470772645191354362138308865549944904868234686191e+215"));

  std::vector<cpp_dec_float_50> bn(5);
  std::vector<cpp_dec_float_50>::iterator it=bn.begin();

  //series retrieved from unchecked bernoulli
  boost::math::bernoulli_b2n<cpp_dec_float_50>(10,5,it);

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[0];
  BOOST_CHECK(output.is_equal("-529.12424242424242424242424242424242424242424242424"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[1];
  BOOST_CHECK(output.is_equal("6192.1231884057971014492753623188405797101449275362"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[2];
  BOOST_CHECK(output.is_equal("-86580.253113553113553113553113553113553113553113553"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[3];
  BOOST_CHECK(output.is_equal("1425517.1666666666666666666666666666666666666666667"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[4];
  BOOST_CHECK(output.is_equal("-27298231.067816091954022988505747126436781609195402"));

  //series retrieved fro tangent numbers algorithm
  boost::math::bernoulli_b2n<cpp_dec_float_50>(50,5,it);

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[0];
  BOOST_CHECK(output.is_equal("-2.8382249570693706959264156336481764738284680928013e+78"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[1];
  BOOST_CHECK(output.is_equal("7.4064248979678850629750827140920984176879731788089e+80"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[2];
  BOOST_CHECK(output.is_equal("-2.0096454802756604483465619672715363186867270822533e+83"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[3];
  BOOST_CHECK(output.is_equal("5.665717005080594144571934603051935696141946828751e+85"));

  output<<std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)<<bn[4];
  BOOST_CHECK(output.is_equal("-1.6584511154136216915823713374319912301494962614725e+88"));
}
