#ifndef _NEWTON__
#define _NEWTON_

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{

    class lagrangian;
    struct data;

    /* Newton's method to solve the non-linear Lagrangian optimality conditions */
    class newton
    {
        public:
            newton(vec& x0, shared_ptr<lagrangian> _ip);

            int run(int nrSteps = 20);

            void printConvergenceMessage(status s=OK);

            shared_ptr<lagrangian> lag;
        private:

            bool step(bool centered=true);

            bool solve(bool centered);

            double lineSearch();

            vec& x;
            shared_ptr<data> ws;

            // saved in place (lag->H) QR data
            int rank;
		    vec hCoeffs;
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;

            int iter = 0;

            options opt;
    };
}

#endif // _NEWTON_
