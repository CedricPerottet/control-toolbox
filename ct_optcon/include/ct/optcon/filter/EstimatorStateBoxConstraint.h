/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

#pragma once

namespace ct {
namespace optcon {

/*!
 * \ingroup Filter
 *
 * \brief Simple box constraint for estimated state.
 *  
 * @tparam     STATE_DIM  The state dimension
 * @tparam     SCALAR     The Scalar type
 */
template <size_t STATE_DIM, typename SCALAR = double>
class EstimatorStateBoxConstraint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using state_vector_t = ct::core::StateVector<STATE_DIM, SCALAR>;
    using VectorXs = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>;
    using SparsityVector = Eigen::Matrix<int, STATE_DIM, 1>;

    EstimatorStateBoxConstraint();

    /**
	 * @brief      Constructor taking lower and upper state bounds directly. Assumes the box constraint is dense.
	 *
	 * @param[in]  lb  The full lower bound
	 * @param[in]  ub  The full upper bound
	 */
    EstimatorStateBoxConstraint(const state_vector_t& lb, const state_vector_t& ub);

    /**
     * @brief 	  Constructor for sparse box constraint. Takes bounds and sparsity pattern.
     * @param lb  Lower boundary values
     * @param ub  Upper boundary values
     * @param sparsity_vec Box constraint sparsity pattern as a vector
     */
    EstimatorStateBoxConstraint(const VectorXs& lb, const VectorXs& ub, const SparsityVector& sparsity_vec);


    state_vector_t stateClipping(const state_vector_t& x);

    bool not_empty() { return has_constraints_; }

private:
    bool has_constraints_;
    VectorXs lb_;
    VectorXs ub_;
    SparsityVector sparsity_vec_;
    void sanityCheck(const size_t& n_constraints, const VectorXs& lb, const VectorXs& ub) const;
};

}  // namespace optcon
}  // namespace ct