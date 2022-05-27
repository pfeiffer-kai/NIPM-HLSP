#include "hlsp.h"
#include "level.h"
#include "levelShared.h"
#include "data.h"

namespace nipmhlsp
{
    hlsp::hlsp(int _p, int _n) : p(_p), n(_n), nr(_n) 
    { 
        lvls.resize(p);
    }

    bool hlsp::setData(int l, const mat& Ae, const vec& be, const mat& Ai, const vec& bi)
    {
        if (!lvls[l])
        {
            lvls[l] = make_shared<level>();
            lvls[l+1] = make_shared<level>(); // virtual level
        }
        lvls[l]->setData(Ae, be, Ai, bi);

        bool allLvlSet = true;
        m = 0;
        me = 0;
        mi = 0;
        for (shared_ptr<level> lvl : lvls) 
        {
            if (!lvl)
            { 
                allLvlSet = false;
                break; 
            }
            else
            {
                m += lvl->m;
                me += lvl->me;
                mi += lvl->mi;
            }
        }
        if (allLvlSet) 
        { 
            lvlSh = make_shared<levelShared>(shared_ptr<hlsp>(this), p, n);
            // vector and matrix workspace for lagrangian
            ws = make_shared<data>(m,n);
        }


        return true;
    }
}
