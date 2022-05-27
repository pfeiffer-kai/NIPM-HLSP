#include "level.h"
#include "data.h"

namespace nipmhlsp
{
    bool level::setData(const mat& _Ae, const vec& _be, const mat& _Ai, const vec& _bi)
    {
        Ae = _Ae;
        be = _be;
        Ai = _Ai;
        bi = _bi;
        AeN = _Ae;
        AiN = _Ai;

        me = Ae.rows();
        mi = Ai.rows();
        m = me + mi;

        return true;
    }

    void level::we(const vec& x, shared_ptr<data> ws) { ws->we() = Ae * x - be; }
    void level::wi(const vec& x, shared_ptr<data> ws) { ws->wi() = Ai * x - bi; }
    void level::_we(const vec& x, shared_ptr<data> ws) { ws->_we() = Ae * x - be; }
    void level::_wi(const vec& x, shared_ptr<data> ws) { ws->_wi() = Ai * x - bi; }
}
