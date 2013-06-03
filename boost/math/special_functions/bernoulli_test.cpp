#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include "bernoulli.hpp"

using namespace boost::test_tools;

BOOST_AUTO_TEST_CASE( float_test )
{

	output_test_stream output;
    output<<std::setprecision(std::numeric_limits<float>::digits10)<<boost::math::bernoulli_b2n<float>(12);
	BOOST_CHECK(output.is_equal("-86580.3"));

    std::vector<float> bn(5);
    std::vector<float>::iterator it=bn.begin();

    //calculate series
    boost::math::bernoulli_b2n<float>(20,5,it);

    output<<std::setprecision(std::numeric_limits<float>::digits10)<<bn[0];
	BOOST_CHECK(output.is_equal("-1.92966e+16"));

    output<<std::setprecision(std::numeric_limits<float>::digits10)<<bn[1];
	BOOST_CHECK(output.is_equal("8.41693e+17"));

    output<<std::setprecision(std::numeric_limits<float>::digits10)<<bn[2];
	BOOST_CHECK(output.is_equal("-4.03381e+19"));

    output<<std::setprecision(std::numeric_limits<float>::digits10)<<bn[3];
	BOOST_CHECK(output.is_equal("2.11507e+21"));

    output<<std::setprecision(std::numeric_limits<float>::digits10)<<bn[4];
	BOOST_CHECK(output.is_equal("-1.20866e+23"));




}

BOOST_AUTO_TEST_CASE( double_test )
{

	output_test_stream output;
    output<<std::setprecision(std::numeric_limits<double>::digits10)<<boost::math::bernoulli_b2n<double>(12);
	BOOST_CHECK(output.is_equal("-86580.2531135531"));

    //no rounding error here as opposed to float
    output<<std::setprecision(std::numeric_limits<double>::digits10)<<boost::math::bernoulli_b2n<double>(5);
	BOOST_CHECK(output.is_equal("0.0757575757575758"));

    std::vector<double> bn(5);
    std::vector<double>::iterator it=bn.begin();

    //calculate series
    boost::math::bernoulli_b2n<double>(100,5,it);

    output<<std::setprecision(std::numeric_limits<double>::digits10)<<bn[0];
	BOOST_CHECK(output.is_equal("-3.64707726451914e+215"));

    output<<std::setprecision(std::numeric_limits<double>::digits10)<<bn[1];
	BOOST_CHECK(output.is_equal("3.75087554364544e+218"));

    output<<std::setprecision(std::numeric_limits<double>::digits10)<<bn[2];
	BOOST_CHECK(output.is_equal("-3.9345867296439e+221"));

    output<<std::setprecision(std::numeric_limits<double>::digits10)<<bn[3];
	BOOST_CHECK(output.is_equal("4.20882111481901e+224"));

    output<<std::setprecision(std::numeric_limits<double>::digits10)<<bn[4];
	BOOST_CHECK(output.is_equal("-4.59022962206179e+227"));



}

BOOST_AUTO_TEST_CASE( long_double_test )
{

	output_test_stream output;
    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<boost::math::bernoulli_b2n<long double>(12);
	BOOST_CHECK(output.is_equal("-86580.2531135531136"));

    //no error rounding here either. Problem only in float
    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<boost::math::bernoulli_b2n<long double>(6);
	BOOST_CHECK(output.is_equal("-0.253113553113553114"));

	std::vector<long double> bn(5);
    std::vector<long double>::iterator it=bn.begin();

    //calculate series
    boost::math::bernoulli_b2n<long double>(100,5,it);

    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<bn[0];
	BOOST_CHECK(output.is_equal("-3.64707726451913544e+215"));

    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<bn[1];
	BOOST_CHECK(output.is_equal("3.75087554364544091e+218"));

    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<bn[2];
	BOOST_CHECK(output.is_equal("-3.93458672964390283e+221"));

    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<bn[3];
	BOOST_CHECK(output.is_equal("4.2088211148190082e+224"));

    output<<std::setprecision(std::numeric_limits<long double>::digits10)<<bn[4];
	BOOST_CHECK(output.is_equal("-4.59022962206179187e+227"));

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

