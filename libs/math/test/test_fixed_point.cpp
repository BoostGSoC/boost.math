#include <iomanip>
#include <iostream>
#include <limits>
#include <boost/math/fixed_point/fixed_point.hpp>

namespace
{
  struct fixed_point_digits
  {
    static const int range      = +256;
    static const int resolution = -240;
  };

  typedef math::fixed_point::negatable<fixed_point_digits::range,
                                       fixed_point_digits::resolution> fixed_point_type;
}

int main()
{
  const fixed_point_type x(fixed_point_type(2) / 3);

  std::cout << std::setprecision(std::numeric_limits<fixed_point_type>::digits10) << x << std::endl;
}
