#ifndef _LEVELSHARED_
#define _LEVELSHARED_

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{
    struct hlsp;
    struct level;
    class lagrangian;
    struct data;

    /* tracks the active and inactive constraints from all levels
     * Requires the dual data from the lagrangian and levels */
    struct levelShared
    {
        levelShared(shared_ptr<hlsp> _hp, int p, int n);

        void wact(const vec& x, shared_ptr<data> ws);
        void winact(const vec& x, shared_ptr<data> ws);
        void winactN(const vec& z, shared_ptr<data> ws);
        void _winact(const vec& x, shared_ptr<data> ws);
        void _winactN(const vec& z, shared_ptr<data> ws);
        
        // apply the nullspace of the active constraints on the left of a vector
        void applyNSOnTheLeftOnDx(int l, shared_ptr<data> ws); 
        void applyNSOnTheLeftOnx(int l, vec& x); 
        // after a newton run, activate constraints from the set of inactive ones
        void addActivatedConstraints_inact(const int l, shared_ptr<data> ws, vec& x);
        // after a newton run, activate ineq constraints from the current level
        void addActivatedConstraints_l(const int l, shared_ptr<data> ws, vec& x);
        // do the decomposition and projection with the updated active set
        void project(const int l, const bool virtualLvl = false);
        
        shared_ptr<hlsp> hp;

        int m;
        int mi;
        int me;
        veci mp;

        veci rank_p;
        int m_act = 0;
        veci m_act_p;
        mati idx_act;
        int m_inact = 0;
        veci m_inact_p;
        mati idx_inact;
        int m_activated; // inequalities that got activcated from m_inact 

        mat A_act, A_actN, A_inact, A_inactN;
        vec b_act, b_inact;

        mat Lambda_act;

        std::vector<vec> hCoeff;
        std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> > permLvl;

        options opt;
    };
}

#endif // _LEVEL_
