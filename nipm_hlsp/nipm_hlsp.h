#ifndef _NIPMHLSP_
#define _NIPMHLSP_

/* This solver solves HLSP's with the interior point method */
// Avoids:
// - avoid constraining the same variable twice on the same level (i.e. turn -2 < x < 2 and -1 < x < 1 into one constraint -1 < x < 1
// - avoid open bounds (i.e turn -1e153(-ifty) < x < 2 into x < 2 since your robotic framework may have this convention)

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{
    struct hlsp;
    struct newton;
    struct lagrangian;
    struct level;
    struct levelShared;

    class NIpmHLSP
    {
    public:
        // default constructor
        NIpmHLSP();
        
        // initializing constructor
        NIpmHLSP(const int p_, const int nVar_);

        // set the data of the problem
    	void setData(const int l_, const MatrixXd& Ae, const VectorXd& be, const MatrixXd& Ai, const VectorXd& bi);

        // solves the IPM-HQP by running the hierarchical Newton's method on the cascade of priority levels
        bool solve();

        vec get_x() { return x; };

        double KKT = 0;
        double iter = 0;
        double time = 0;

    private:
        
        shared_ptr<hlsp> hp;
        std::vector<shared_ptr<lagrangian> > lag;
        newton* nwt;

        VectorXd x; // primal

        options opt;

    };
} // namespace solver

#endif // _NIPMHLSP
