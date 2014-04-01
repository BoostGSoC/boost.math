#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <boost/math/fixed_point/fixed_point.hpp>

typedef math::fixed_point::negatable<32, -24> fixed_point_type;

int main()
{
  bool   success = true;
  fixed_point_type    fxp_a, fxp_b;
  float  flo_a, flo_b;
  double dou_a, dou_b;

  fxp_a = fixed_point_type( 64)   + fixed_point_type(1)   / fixed_point_type(3);
  flo_a =      64.0F +     1.0F /     3.0F;
  dou_a =      64.0  +     1.0  /     3.0;

  fxp_b = fixed_point_type(127)   + fixed_point_type(1)   / fixed_point_type(3);
  flo_b =     127.0F +     1.0F /     3.0F;
  dou_b =     127.0  +     1.0  /     3.0;

  std::cout << "Testing range for math::fixed_point::negatable<32, -24> is: [64, 128)" << std::endl;
  std::cout << "Precision Comparison:" << std::setprecision(std::numeric_limits<long double>::digits10) << std::endl;
  std::cout << "Fixed: a=" << fxp_a << "; b=" << fxp_b << std::endl;
  std::cout << "Float: a=" << flo_a << "; b=" << flo_b << std::endl;
  std::cout << "Doubl: a=" << dou_a << "; b=" << dou_b << std::endl;


  //Test cases bzgl ctors und Präzision:
  if( !( fixed_point_type(71) == 71.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) != 71.0 ) ) { success = false; }
  if( !( fixed_point_type(71) <= 71.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) <  71.0 ) ) { success = false; }
  if( !( fixed_point_type(71) >= 71.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) >  71.0 ) ) { success = false; }

  if(  ( fixed_point_type(71) == 84.0 ) ) { success = false; }
  if( !( fixed_point_type(71) != 84.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) <= 69.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) <  69.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) >= 84.0 ) ) { success = false; }
  if(  ( fixed_point_type(71) >  84.0 ) ) { success = false; }

  if( !( fixed_point_type(71) <= 84.0 ) ) { success = false; }
  if( !( fixed_point_type(71) <  84.0 ) ) { success = false; }
  if( !( fixed_point_type(71) >= 69.0 ) ) { success = false; }
  if( !( fixed_point_type(71) >  69.0 ) ) { success = false; }

  if( !( fixed_point_type(71.17) == 71.17 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) != 71.17 ) ) { success = false; }
  if( !( fixed_point_type(71.17) <= 71.17 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) <  71.17 ) ) { success = false; }
  if( !( fixed_point_type(71.17) >= 71.17 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) >  71.17 ) ) { success = false; }

  if(  ( fixed_point_type(71.17) == 84.0 ) ) { success = false; }
  if( !( fixed_point_type(71.17) != 84.0 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) <= 69.0 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) <  69.0 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) >= 84.0 ) ) { success = false; }
  if(  ( fixed_point_type(71.17) >  84.0 ) ) { success = false; }

  if( !( fixed_point_type(71.17) <= 84.0 ) ) { success = false; }
  if( !( fixed_point_type(71.17) <  84.0 ) ) { success = false; }
  if( !( fixed_point_type(71.17) >= 69.0 ) ) { success = false; }
  if( !( fixed_point_type(71.17) >  69.0 ) ) { success = false; }


  //Test cases for all operators with data types: fixed_point_type, int, double:
  //LH: negatable
  if( !( (fixed_point_type(5) +  fixed_point_type(7)   ) == ( 5.0 +  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) +      7    ) == ( 5.0 +  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) +      7.0  ) == ( 5.0 +  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) -  fixed_point_type(7)   ) == ( 5.0 -  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) -      7    ) == ( 5.0 -  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) -      7.0  ) == ( 5.0 -  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) *  fixed_point_type(7)   ) == ( 5.0 *  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) *      7    ) == ( 5.0 *  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) *      7.0  ) == ( 5.0 *  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) /  fixed_point_type(7)   ) == ( 5.0 /  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) /      7    ) == ( 5.0 /  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) /      7.0  ) == ( 5.0 /  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <  fixed_point_type(7)   ) == ( 5.0 <  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <      7    ) == ( 5.0 <  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <      7.0  ) == ( 5.0 <  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >  fixed_point_type(7)   ) == ( 5.0 >  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >      7    ) == ( 5.0 >  7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >      7.0  ) == ( 5.0 >  7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <= fixed_point_type(7)   ) == ( 5.0 <= 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <=     7    ) == ( 5.0 <= 7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) <=     7.0  ) == ( 5.0 <= 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >= fixed_point_type(7)   ) == ( 5.0 >= 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >=     7    ) == ( 5.0 >= 7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) >=     7.0  ) == ( 5.0 >= 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) == fixed_point_type(7)   ) == ( 5.0 == 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) ==     7    ) == ( 5.0 == 7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) ==     7.0  ) == ( 5.0 == 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) != fixed_point_type(7)   ) == ( 5.0 != 7.0 ) ) ) { success = false; }
  if( !( (fixed_point_type(5) !=     7    ) == ( 5.0 != 7   ) ) ) { success = false; }
  if( !( (fixed_point_type(5) !=     7.0  ) == ( 5.0 != 7.0 ) ) ) { success = false; }

  //LH: arithmetic
  if( !( ( fixed_point_type(7)   +  fixed_point_type(5) ) == ( 7.0 +  5.0 ) ) ) { success = false; }
  if( !( (     7    +  fixed_point_type(5) ) == ( 7   +  5.0 ) ) ) { success = false; }
  if( !( (     7.0  +  fixed_point_type(5) ) == ( 7.0 +  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   -  fixed_point_type(5) ) == ( 7.0 -  5.0 ) ) ) { success = false; }
  if( !( (     7    -  fixed_point_type(5) ) == ( 7   -  5.0 ) ) ) { success = false; }
  if( !( (     7.0  -  fixed_point_type(5) ) == ( 7.0 -  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   *  fixed_point_type(5) ) == ( 7.0 *  5.0 ) ) ) { success = false; }
  if( !( (     7    *  fixed_point_type(5) ) == ( 7   *  5.0 ) ) ) { success = false; }
  if( !( (     7.0  *  fixed_point_type(5) ) == ( 7.0 *  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   /  fixed_point_type(5) ) == ( 7.0 /  5.0 ) ) ) { success = false; }
  if( !( (     7    /  fixed_point_type(5) ) == ( 7   /  5.0 ) ) ) { success = false; }
  if( !( (     7.0  /  fixed_point_type(5) ) == ( 7.0 /  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   <  fixed_point_type(5) ) == ( 7.0 <  5.0 ) ) ) { success = false; }
  if( !( (     7    <  fixed_point_type(5) ) == ( 7   <  5.0 ) ) ) { success = false; }
  if( !( (     7.0  <  fixed_point_type(5) ) == ( 7.0 <  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   >  fixed_point_type(5) ) == ( 7.0 >  5.0 ) ) ) { success = false; }
  if( !( (     7    >  fixed_point_type(5) ) == ( 7   >  5.0 ) ) ) { success = false; }
  if( !( (     7.0  >  fixed_point_type(5) ) == ( 7.0 >  5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   <= fixed_point_type(5) ) == ( 7.0 <= 5.0 ) ) ) { success = false; }
  if( !( (     7    <= fixed_point_type(5) ) == ( 7   <= 5.0 ) ) ) { success = false; }
  if( !( (     7.0  <= fixed_point_type(5) ) == ( 7.0 <= 5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   >= fixed_point_type(5) ) == ( 7.0 >= 5.0 ) ) ) { success = false; }
  if( !( (     7    >= fixed_point_type(5) ) == ( 7   >= 5.0 ) ) ) { success = false; }
  if( !( (     7.0  >= fixed_point_type(5) ) == ( 7.0 >= 5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   == fixed_point_type(5) ) == ( 7.0 == 5.0 ) ) ) { success = false; }
  if( !( (     7    == fixed_point_type(5) ) == ( 7   == 5.0 ) ) ) { success = false; }
  if( !( (     7.0  == fixed_point_type(5) ) == ( 7.0 == 5.0 ) ) ) { success = false; }
  if( !( ( fixed_point_type(7)   != fixed_point_type(5) ) == ( 7.0 != 5.0 ) ) ) { success = false; }
  if( !( (     7    != fixed_point_type(5) ) == ( 7   != 5.0 ) ) ) { success = false; }
  if( !( (     7.0  != fixed_point_type(5) ) == ( 7.0 != 5.0 ) ) ) { success = false; }

  std::cout << '\n'
            << std::boolalpha
            << "All results are correct: "
            << success
            << '\n'
            << std::endl;

  std::cout << "Sqrt: "  << '\n' << sqrt (fixed_point_type(0.5))      << '\n' << std::sqrt(0.5)           << std::endl;
  std::cout << "Exp: "   << '\n' << exp  (fixed_point_type(3.14159))  << '\n' << std::exp(3.14159)        << std::endl;
  std::cout << "Log: "   << '\n' << log  (fixed_point_type(3.14159))  << '\n' << std::log(3.14159)        << std::endl;
  std::cout << "Log: "   << '\n' << log  (fixed_point_type(0.5))      << '\n' << std::log(0.5)            << std::endl;
  std::cout << "Sin: "   << '\n' << sin  (fixed_point_type(2))        << '\n' << std::sin(2.0)            << std::endl;
  std::cout << "Cos: "   << '\n' << cos  (fixed_point_type(2))        << '\n' << std::cos(2.0)            << std::endl;
  std::cout << "Tan: "   << '\n' << tan  (fixed_point_type(2))        << '\n' << std::tan(2.0)            << std::endl;
  std::cout << "Asin: "  << '\n' << asin (fixed_point_type(0.5))      << '\n' << std::asin(0.5)           << std::endl;
  std::cout << "Acos: "  << '\n' << acos (fixed_point_type(0.5))      << '\n' << std::acos(0.5)           << std::endl;
  std::cout << "Atan: "  << '\n' << atan (fixed_point_type(3.14159))  << '\n' << std::atan(3.14159)       << std::endl;
  std::cout << "Sinh: "  << '\n' << sinh (fixed_point_type(3.14159))  << '\n' << std::sinh(3.14159)       << std::endl;
  std::cout << "Cosh: "  << '\n' << cosh (fixed_point_type(3.14159))  << '\n' << std::cosh(3.14159)       << std::endl;
  std::cout << "Tanh: "  << '\n' << tanh (fixed_point_type(3.14159))  << '\n' << std::tanh(3.14159)       << std::endl;
//  std::cout << "ASinh: " << '\n' << asinh(fixed_point_type(3.14159))  << '\n' << std::asinh(3.14159)      << std::endl;
//  std::cout << "ACosh: " << '\n' << acosh(fixed_point_type(3.14159))  << '\n' << std::acosh(3.14159)      << std::endl;
//  std::cout << "ATanh: " << '\n' << atanh(fixed_point_type(0.5))      << '\n' << std::atanh(1.0 / 2.0)    << std::endl;

  std::cout << '\n' << "Test limits" << '\n' << std::endl;
  std::cout << "Eps:     "  << std::numeric_limits<fixed_point_type>::epsilon() << std::endl;
  std::cout << "1 - Eps: "  << fixed_point_type(1) - std::numeric_limits<fixed_point_type>::epsilon() << std::endl;
  std::cout << "Min:     "  << (std::numeric_limits<fixed_point_type>::min)()   << std::endl;
  std::cout << "Max:     "  << (std::numeric_limits<fixed_point_type>::max)()   << std::endl;
}
