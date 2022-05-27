#ifndef _LEVEL_
#define _LEVEL_

#include "typedefs.h"
#include "data.h"

namespace nipmhlsp
{
    struct data;

    /* contains the raw data and its projection into the NS of the active constraints */
    struct level
    {
        int m = 0;
        int me = 0;
        int mi = 0;

        mat Ae, Ai;
        mat AeN, AiN;
        vec be, bi;

        int rank;

        bool setData(const mat& Ae, const vec& be, const mat& Ai, const vec& bi);

        void we(const vec& x, shared_ptr<data> ws);
        void wi(const vec& x, shared_ptr<data> ws);
        void _we(const vec& x, shared_ptr<data> ws);
        void _wi(const vec& x, shared_ptr<data> ws);

    };
}

#endif // _LEVEL_
