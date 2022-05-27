#pragma once

/* DEBUGGING TOOLS */
// a number of safety tests is performed:
// - nullspace test
// - dLagrange multiplier test
#ifndef SAFEGUARD
#define SAFEGUARD false
#endif
// solver timings
#ifndef TIMEMEASUREMENTS
#define TIMEMEASUREMENTS true
#endif

namespace nipmhlsp
{
    struct options
    {
        // Parameters
        const int maxIter = 10; // maximum number of Newton iterations per priority level
        const double dualInitThres = 1.; // minimum initialization value for dual, only for lam_inact, w_inact is explicit
        const double linearDependency = 1e-12; // linear dependency in rank revealing qr decomposition
        const double alphaThres = 1e-12; // linear search stuck threshold
        const double KKT_thres = 1e-12; // convergence of KKT conditions
        const double mu_thres = 1e-3; // threshold for interior point convergence
        const double numThres = 1e-8; // non linear dual contribution to KKT: too large - line search stuck, too small - numerical issues FIXME: adaptive? or at least Eigen style thresholding 
        const double posThres = 0.; // line search
        const double activationThres = 1e-8; // activate from inactive constraints A_inact if corresponding Lagrange multiplier larger than this threshold
        const double LagrangeTestThres = 1e-4;

        // 0: exit Newton only if KKTnorm so small that maximum nr of iteration is reached or other numerical issues arise
        // 1: calculate KKT in each Newton iteration
        // 2: calculate dual only once ||KKT_2, ..., KKT_7||^2 < epsilon, then convergence check with full KKT // leads to slower convergence, why? Is LagrangeMultipliers_N not accurate? // NIY
        // 3: same as 2 but also check ||NTr1|| // NIY
        const int convergenceTest = 1;

        // Mehrotra's predictor corrector algorithm
        // the decomposition is only done once but the primal and dual need to be calculated twice (high ~n^2), can have significant influence on computation time
        const bool predCor = true; // FIXME: there might be a bug here since sometimes bad behaviour (kkt increase in early stages of Newton)
        const int predFactor = 3; // Nocedal: 3
        // if false, how to calculate mu
        // 0: ... / (nVar + m_inact)
        // 1: ... / (nVar + m_inact^2)
        // 2: ... / (nVar + m_inact)^2
        const int muCalc = 2;

        // Newton's method in least squares form
        const bool leastSquares = false;

        // 0: NONE : none
        // 1: CONV : basic solver info
        // 2: SOLVE : more detailed solver info
        // 3: MAT : matrices
        const int verbose = NONE;

        void print()
        {
            cout << "--- RUN NIPM_HLSP WITH FOLLOWING CONFIGURATION" << endl;
            cout << "maxIter:                " << maxIter << endl;
            cout << "dualInitThres:          " << dualInitThres << endl;
            cout << "linearDependency:       " << linearDependency << endl; 
            cout << "alphaThres:             " << alphaThres << endl;
            cout << "KKT_thres:              " << KKT_thres << endl;
            cout << "mu_thres:               " << mu_thres << endl; 
            cout << "numThres:               " << numThres << endl;
            cout << "posThres:               " << posThres << endl;
            cout << "activationThres:        " << activationThres << endl;
            cout << "LagrangeTestThres:      " << LagrangeTestThres << endl;
            cout << "convergenceTest:        " << convergenceTest << endl;
            cout << "predCor:                " << predCor << endl;
            cout << "predFactor:             " << predFactor << endl;
            cout << "muCalc:                 " << muCalc << endl;
            cout << "leastSquares:           " << leastSquares << endl;
            cout << "verbose:                " << verbose << endl;
            cout << "---" << endl;
        }
    };
}
