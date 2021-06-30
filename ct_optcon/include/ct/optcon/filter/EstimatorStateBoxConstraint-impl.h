/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

#pragma once

namespace ct {
namespace optcon {


template <size_t STATE_DIM, typename SCALAR>
EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::EstimatorStateBoxConstraint()
{
    has_constraints_ = false;
}

template <size_t STATE_DIM, typename SCALAR>
EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::EstimatorStateBoxConstraint(const state_vector_t& lb,
    const state_vector_t& ub)
{
    lb_ = lb;
    ub_ = ub;
    sparsity_vec_ = SparsityVector::Constant(1);
    has_constraints_ = true;
    sanityCheck(STATE_DIM, lb_, ub_);
}

template <size_t STATE_DIM, typename SCALAR>
EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::EstimatorStateBoxConstraint(const VectorXs& lb,
    const VectorXs& ub,
    const SparsityVector& sparsity_vec)
{
    lb_ = lb;
    ub_ = ub;
    sparsity_vec_ = sparsity_vec;
    has_constraints_ = true;
    // make sure the provided sparsity pattern is correct and consists only of ones and zeros
    assert(sparsity_vec.maxCoeff() <= 1);
    assert(sparsity_vec.minCoeff() >= 0);

    size_t n_constraints = (size_t)sparsity_vec.sum();

    sanityCheck(n_constraints, lb_, ub_);
}


template <size_t STATE_DIM, typename SCALAR>
typename EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::state_vector_t
EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::stateClipping(const state_vector_t& x)
{
    if (!has_constraints_)
    {
        return x;
    }
    else
    {
        state_vector_t x_clipped = x;
        size_t i_constraint = 0;
        for (size_t i_state_dim = 0; i_state_dim < STATE_DIM; i_state_dim++)
        {
            if (sparsity_vec_[i_state_dim] != 0)
            {
                x_clipped[i_state_dim] = std::min(x_clipped[i_state_dim], ub_[i_constraint]);
                x_clipped[i_state_dim] = std::max(x_clipped[i_state_dim], lb_[i_constraint]);
                i_constraint++;
            }
        }
        return x_clipped;
    }
}

template <size_t STATE_DIM, typename SCALAR>
void EstimatorStateBoxConstraint<STATE_DIM, SCALAR>::sanityCheck(const size_t& n_constraints,
    const VectorXs& lb,
    const VectorXs& ub) const
{
    // assert that the size of constraint vectors is equal to the computed/given number of constraints
    if ((lb.rows() != static_cast<int>(n_constraints)) | (ub.rows() != static_cast<int>(n_constraints)))
    {
        std::cout << "no. Constraints: " << n_constraints << std::endl;
        std::cout << "EstimatorStateBoxConstraint: lb " << lb.transpose() << std::endl;
        std::cout << "EstimatorStateBoxConstraint: ub " << ub.transpose() << std::endl;
        throw std::runtime_error("EstimatorStateBoxConstraint: wrong constraint sizes in StateConstraint");
    }

    // assert that the boundaries are meaningful
    for (size_t i = 0; i < n_constraints; i++)
    {
        if (lb(i) > ub(i))
        {
            std::cout << "EstimatorStateBoxConstraint: lb " << lb.transpose() << std::endl;
            std::cout << "EstimatorStateBoxConstraint: ub " << ub.transpose() << std::endl;
            throw std::runtime_error("EstimatorStateBoxConstraint: wrong boundaries: lb > ub");
        }
    }
}

}  // namespace optcon
}  // namespace ct