#ifndef _LAGRANGIAN_
#define _LAGRANGIAN_

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{
    struct hlsp;
    struct level;
    struct levelShared;
    struct data;

    class lagrangian
    {
        // The Lagrangian of the constrained optimization problem
        // does not contain vital hlsp information and can be deleted anytime

        public:
            lagrangian(int _l, shared_ptr<hlsp> _hp);

            void init(const vec& x);

            void computeHg(const vec& x, bool centered=false);
            void computeStep(const vec& x, bool centered=false);
            void makeStep(vec& x);

            double lineSearch(bool reinitDual=false);

            status convergenceTest(const vec& x);
            void KKT(const vec& x, bool full=false);

            void printStatus(const vec& x);

            shared_ptr<hlsp> hp;
            shared_ptr<data> ws;

            int l;
            int m;
            int me;
            int mi;
            int m_act;
            int m_act_l;
            int m_inact;

            void printVars(const vec& x);

        private:

            void computeH(const vec& x); // compute Lagrangian Hessian
            void computeg(const vec& x, bool centered); // compute 

            void calcMu();
            void calcSigma(bool centered=true);
            
            // calculate the increment of the Lagrange multipliers
            // works with data from A_actN
            // recursively after Dimitrov2015
            void dLagrangeMultipliers(const vec& x);
            Eigen::ColPivHouseholderQR<MatrixXd> qrd;

            double convergenceNorm = 0;
            double maxKKTnorm = 1e13;

            level& lvl;
            levelShared& lvlSh;

            double alpha = 0;

            double mu;
            double mu_l;
            double mu_aff;
            double mu_aff_l;
            double sigma = 1.;
            double sigma_l = 1.;
            double numThresCur;

            options opt;
    };
}

#endif // _LEVEL_
