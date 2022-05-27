#ifndef _HLSP_
#define _HLSP_

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{
    struct level;
    struct levelShared; 
    struct data;

    struct hlsp
    {
        // contains all information about the HLSP
        public:
            hlsp(int _p, int _n);

            bool setData(int l, const mat& Ae, const vec& be, const mat& Ai, const vec& bi);

            // constructor variables
            int p; // number of priority levels; note that this number is twice the number of levels of the original problem due to virtual priority levels handling activated constraints from A_inact; l0: virtual, l1: normal, l2: virtual and so forth
            int n;
            int nr;
            int m;
            int me;
            int mi;

            std::vector<shared_ptr<level> > lvls;
            shared_ptr<levelShared> lvlSh;
            shared_ptr<data> ws;

            void print()
            {
                cout << "--- SOLVE HLSP WITH" << endl;
                cout << "p:  " << p << endl;
                cout << "n:  " << n << endl;
                cout << "m:  " << m << endl;
                cout << "me: " << me << endl;
                cout << "mi: " << mi << endl;
                cout << "---" << endl;
            };
    };
}

#endif // _HLSP_
