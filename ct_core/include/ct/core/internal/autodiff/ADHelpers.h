/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

#pragma once

namespace ct {
namespace core {
namespace ADHelperFunctions {

template <typename SCALAR, typename T>
SCALAR convert_to_scalar(const T& x)
{
    return SCALAR(x);
}

#ifdef CPPAD
// Specialization for Cpp::AD that has no implicit conversion to scalar
template <typename SCALAR>
SCALAR convert_to_scalar(const CppAD::AD<SCALAR>& ad_scalar)
{
    return CppAD::Value(CppAD::Var2Par(ad_scalar));
}
#endif


}  // namespace ADHelperFunctions

class ADHelpers
{
public:
private:
};
}  // namespace core
}  // namespace ct
