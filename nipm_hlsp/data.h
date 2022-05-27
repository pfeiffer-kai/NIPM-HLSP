#ifndef _DATA_
#define _DATA_

#pragma once
#include "typedefs.h"
#include "options.h"

namespace nipmhlsp
{
    struct data
    {
        // is this fast?
        
        // contains all information about the HLSP
        public:
            data(int _r, int _c);
            void reRef(int _n, int _nr, int _m_act, int _m_inact, int _me, int _mi); 

            //linear system
            Ref<mat> H();
            Ref<vec> g();
            // KKT
            Ref<vec> K(const string index);
            // --- dual
            // act is in levelShared since it is permanent data
            // inact
            Ref<vec> w_inact();
            Ref<vec> _w_inact();
            Ref<vec> w_inactd(); // the dual one
            Ref<vec> lam_inact();
            // eq
            Ref<vec> we(); // the one updated in the Newton method
            Ref<vec> _we(); // the explicit one

            Ref<vec> wed(); // the dual one
            // ineq
            Ref<vec> wi();
            Ref<vec> _wi();
            Ref<vec> wid(); // the dual one
            Ref<vec> omega();
            // dual scales
            Ref<vec> invW_inact();
            Ref<vec> sqInvWLam_inact();
            Ref<vec> sqWInvLam_inact();
            Ref<vec> invWO();
            Ref<vec> IpObWO();
            Ref<vec> sqIpObWO();
            Ref<vec> invSqIpObWO();
            // --- steps
            // primal
            Ref<vec> dx();
            // active
            Ref<vec> dlam_act();
            // inactive
            Ref<vec> dw_inact();
            Ref<vec> dlam_inact();
            // eq
            Ref<vec> dwe();
            // ineq
            Ref<vec> dwi();
            Ref<vec> domega();

            void print()
            {
                cout << "--- DATA WITH" << endl;
                cout << "r:  " << r << endl;
                cout << "c:  " << c << endl;
                cout << "n:  " << n << endl;
                cout << "nr:  " << nr << endl;
                cout << "m:  " << m << endl;
                cout << "m_act:  " << m_act << endl;
                cout << "m_inact:  " << m_inact << endl;
                cout << "me:  " << me << endl;
                cout << "mi:  " << mi << endl;
                cout << "kktrg:  " << kktrg << endl;
                cout << "inc " << inc.transpose().head(27) << endl;

                cout << "wsv\n" << wsv.transpose() << endl;
                // if (opt.verbose >= MAT) cout << "workspace\n" << wsm << endl;
                cout << "---" << endl;
            };

        private:
            int r;
            int c;

            mat wsm;
            vec wsv;

            veci inc;

            // the same data is also stored in hp, use that instead
            int n, nr, m, m_act, m_inact, me, mi, kktrg;

            options opt;

    };
}

#endif // _HLSP_
