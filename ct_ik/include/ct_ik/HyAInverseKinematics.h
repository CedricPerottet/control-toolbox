#pragma once

#define IKFAST_HAS_LIBRARY
#define IKFAST_NAMESPACE hya_ik
#include <ikfast.h>
#include <Eigen/Dense>

#include <ct/rbd/robot/kinematics/InverseKinematicsBase.h>

using namespace ikfast;

namespace ct {
namespace rbd {

template <typename SCALAR = double>
class HyAInverseKinematics : InverseKinematicsBase<6, SCALAR>
{
public:
    virtual std::vector<typename tpl::JointState<6, SCALAR>::Position> computeInverseKinematics(
        const tpl::RigidBodyPose<SCALAR>& eeBasePose,
        const std::vector<SCALAR>& freeJoints = std::vector<SCALAR>()) const
    {
        // TODO: Check for valid solutions.
        IkSolutionList<double> solutions;

        if (size_t(hya_ik::GetNumFreeParameters()) != freeJoints.size()) throw "Error";

        hya_ik::ComputeIk(eeBasePose.position().toImplementation().data(),
            eeBasePose.getRotationMatrix().toImplementation().data(),
            freeJoints.size() > 0 ? freeJoints.data() : nullptr, solutions);

        size_t num_solutions = solutions.GetNumSolutions();
        std::vector<typename tpl::JointState<6, SCALAR>::Position> res(solutions.GetNumSolutions());

        for (size_t i = 0u; i < num_solutions; ++i)
        {
            const IkSolutionBase<double>& solution = solutions.GetSolution(i);
            solution.GetSolution(res[i].data(), freeJoints.size() > 0 ? freeJoints.data() : nullptr);
        }

        return res;
    }

    virtual std::vector<typename tpl::JointState<6, SCALAR>::Position> computeInverseKinematics(
        const tpl::RigidBodyPose<SCALAR>& eeWorldPose,
        const tpl::RigidBodyPose<SCALAR>& baseWorldPose,
        const std::vector<SCALAR>& freeJoints = std::vector<SCALAR>()) const
    {
        return computeInverseKinematics(eeWorldPose.inReferenceFrame(baseWorldPose), freeJoints);
    }
};
} /* namespace rbd */
} /* namespace ct */
