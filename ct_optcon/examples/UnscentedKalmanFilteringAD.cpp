/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

/*!
 * \example KalmanFiltering.cpp
 *
 * This example shows how to use the Kalman Filter to estimate the state of a simple oscillator.
 *
 */

#include <ct/optcon/optcon.h>
#include <ct/core/examples/CustomController.h>  // Using the custom controller from ct_core examples.
#include "exampleDir.h"

int main(int argc, char** argv)
{
    // file with kalman weights
    using SCALAR_AD = CppAD::AD<double>;
    typedef ct::core::tpl::SecondOrderSystem<SCALAR_AD> SecondOrderSystemAD;
    std::string settingsFile = ct::optcon::exampleDir + "/ukfWeights.info";

    // a damped oscillator has two states, position and velocity
    const size_t state_dim = SecondOrderSystemAD::STATE_DIM;      // = 2
    const size_t control_dim = SecondOrderSystemAD::CONTROL_DIM;  // = 1
    const size_t output_dim = 2;  // we assume we observe the full state (however with noise)


    // create an initial state: we initialize it at a point with unit deflection and zero velocity
    ct::core::StateVector<state_dim, SCALAR_AD> x;
    x(0) = 1.0;
    x(1) = 0.0;

    // create an oscillator
    SCALAR_AD w_n = 50;
    std::shared_ptr<SecondOrderSystemAD> oscillator(new SecondOrderSystemAD(w_n));

    // create an integrator for "simulating" the measured data
    ct::core::Integrator<state_dim, SCALAR_AD> integrator(oscillator, ct::core::IntegrationType::EULERCT);

    ct::core::StateVectorArray<state_dim, SCALAR_AD> states;
    ct::core::ControlVectorArray<control_dim, SCALAR_AD> controls;
    ct::core::tpl::TimeArray<SCALAR_AD> times;

    // simulate 100 steps
    double dt = 0.001;
    SCALAR_AD dt_AD = dt;
    size_t nSteps = 100;
    states.push_back(x);
    for (size_t i = 0; i < nSteps; i++)
    {
        // compute control (needed for filter later)
        ct::core::ControlVector<control_dim, SCALAR_AD> u_temp;
        u_temp.setRandom();
        controls.push_back(u_temp);

        integrator.integrate_n_steps(x, i * dt, 1, dt);

        states.push_back(x);
        times.push_back(i * dt);
    }

    // create system observation matrix C: we measure both position and velocity
    ct::core::OutputStateMatrix<output_dim, state_dim, SCALAR_AD> C;
    C.setIdentity();

    // load Kalman Filter weighting matrices from file
    ct::core::StateMatrix<state_dim, SCALAR_AD> Q, dFdv;
    ct::core::OutputMatrix<output_dim, SCALAR_AD> R;
    dFdv.setIdentity();
    ct::core::loadMatrix(settingsFile, "kalman_weights.Q", Q);
    ct::core::loadMatrix(settingsFile, "kalman_weights.R", R);
    std::cout << "Loaded Kalman R as " << std::endl << R << std::endl;
    std::cout << "Loaded Kalman Q as " << std::endl << Q << std::endl;

    // create a sensitivity approximator to compute A and B matrices
    std::shared_ptr<ct::core::SystemLinearizer<state_dim, control_dim, SCALAR_AD>> linearizer(
        new ct::core::SystemLinearizer<state_dim, control_dim, SCALAR_AD>(oscillator));

    std::shared_ptr<ct::core::SensitivityApproximation<state_dim, control_dim, state_dim / 2, state_dim / 2, SCALAR_AD>>
        sensApprox(
            new ct::core::SensitivityApproximation<state_dim, control_dim, state_dim / 2, state_dim / 2, SCALAR_AD>(
                dt_AD, linearizer));


    // the observer is supplied with a dynamic model identical to the one used above for data generation
    std::shared_ptr<SecondOrderSystemAD> oscillator_observer_model(new SecondOrderSystemAD(w_n));

    // set up the system model
    std::shared_ptr<ct::optcon::CTSystemModel<state_dim, control_dim, SCALAR_AD>> sysModel(
        new ct::optcon::CTSystemModel<state_dim, control_dim, SCALAR_AD>(oscillator_observer_model, sensApprox, dFdv));

    // set up the measurement model
    ct::core::OutputMatrix<output_dim, SCALAR_AD> dHdw;
    dHdw.setIdentity();
    std::shared_ptr<ct::optcon::LinearMeasurementModel<output_dim, state_dim, SCALAR_AD>> measModel(
        new ct::optcon::LTIMeasurementModel<output_dim, state_dim, SCALAR_AD>(C, dHdw));

    // set up state constraint
    ct::core::StateVector<state_dim, SCALAR_AD> lb, ub;
    lb << -24.0, -24.0;
    ub << 25.0, 25.0;
    ct::optcon::EstimatorStateBoxConstraint<state_dim, SCALAR_AD> box_constraint(lb, ub);

    // set up Filter : Unscented Kalman filter
    ct::optcon::UnscentedKalmanFilter<state_dim, control_dim, output_dim, SCALAR_AD> filter(
        sysModel, measModel, states[0], box_constraint);


    ct::core::StateMatrix<state_dim> meas_var;
    ct::core::loadMatrix(settingsFile, "measurement_noise.measurement_var", meas_var);

    ct::core::StateVectorArray<state_dim, SCALAR_AD> states_est(states.size());
    ct::core::StateVectorArray<state_dim, SCALAR_AD> states_meas(states.size());
    states_est[0] = states[0];
    states_meas[0] = states[0];

    // run the filter over the simulated data
    for (size_t i = 1; i < states.size(); ++i)
    {
        // compute an observation
        states_meas[i] = states[i];

        ct::core::OutputVector<output_dim, SCALAR_AD> y = C * states_meas[i];

        // Kalman filter prediction step
        filter.predict(controls[i], dt, dt * i);

        // Kalman filter estimation step
        ct::core::StateVector<state_dim, SCALAR_AD> x_est = filter.update(y, dt, dt * i);

        // and log for printing
        states_est[i] = x_est;
    }


    for (size_t i = 0; i < states_est.size(); ++i)
        std::cout << "State\t\tState_est\n"
                  << std::fixed << std::setprecision(6) << states[i][0] << "\t" << states_est[i][0] << std::endl
                  << states[i][1] << "\t" << states_est[i][1] << std::endl
                  << std::endl;
    return 0;
}
