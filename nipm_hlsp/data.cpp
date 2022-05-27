#include "data.h"

namespace nipmhlsp
{
    data::data(int _r, int _c) : r(30 * _r), c(_c)
    {
        wsm.resize(r, c); wsm.setZero();
        wsv.resize(r); wsv.setZero();
    }

    void data::reRef(int _n, int _nr, int _m_act, int _m_inact, int _me, int _mi)
    {
        n = _n;
        nr = _nr;
        m_act = _m_act;
        m_inact = _m_inact;
        me = _me;
        mi = _mi;

        if (opt.leastSquares) m = m_inact + mi + me;
        else m = nr;

        kktrg = n + me + 2*mi + m_act + 2*m_inact;

        inc.resize(1000); inc.setZero();
        // g
        inc[0] = 0;
        // KKT
        inc[1] = max(m, nr);
        // --- dual
        // inact 
        inc[2] = inc[1] + kktrg;
        inc[3] = inc[2] + m_inact;
        inc[4] = inc[3] + m_inact;
        inc[5] = inc[4] + m_inact;
        // eq 
        inc[6] = inc[5] + m_inact;
        inc[7] = inc[6] + me;
        inc[8] = inc[7] + me;
        // ineq
        inc[9] = inc[8] + me;
        inc[10] = inc[9] + mi;
        inc[11] = inc[10] + mi;
        inc[12] = inc[11] + mi;
        // dual scales
        inc[13] = inc[12] + mi;
        inc[14] = inc[13] + m_inact;
        inc[15] = inc[14] + m_inact;
        inc[16] = inc[15] + m_inact;
        inc[17] = inc[16] + mi;
        inc[18] = inc[17] + mi;
        inc[19] = inc[18] + mi;
        // --- steps
        // primal
        inc[20] = inc[19] + mi;
        // active
        inc[21] = inc[20] + n;
        // inactive
        inc[22] = inc[21] + m_act;
        inc[23] = inc[22] + m_inact;
        // eq
        inc[24] = inc[23] + m_inact;
        // ineq
        inc[25] = inc[24] + me;
        inc[26] = inc[25] + mi;

    }

    Ref<mat> data::H() { return wsm.topLeftCorner(m, nr); }

    // vectors 
    Ref<vec> data::g() { return wsv.head(max(m, nr)); }

    Ref<vec> data::K(const string index)
    {
        if (index == "all") return wsv.segment(inc[1], kktrg);
        else if (index == "all_inact") return wsv.segment(inc[1] + n + me + m_act, 2 * m_inact + 2 * mi);
        else if (index == "x") return wsv.segment(inc[1], n);
        else if (index == "lambda_e") return wsv.segment(inc[1] + n, me);
        else if (index == "lambda_act") return wsv.segment(inc[1] + n + me, m_act);
        else if (index == "lambda_inact") return wsv.segment(inc[1] + n + me + m_act, m_inact);
        else if (index == "w_inact") return wsv.segment(inc[1] + n + me + m_act + m_inact, m_inact);
        else if (index == "lambda_i") return wsv.segment(inc[1] + n + me + m_act + 2 * m_inact, mi);
        else if (index == "omega_i") return wsv.segment(inc[1] + n + me + m_act + 2 * m_inact + mi, mi);
        else { cout << "ERROR: K index not available" << endl; throw; }
    }

    // inact
    Ref<vec> data::w_inact() { return wsv.segment(inc[2], m_inact); }
    Ref<vec> data::_w_inact() { return wsv.segment(inc[3], m_inact); }
    Ref<vec> data::w_inactd() { return wsv.segment(inc[4], m_inact); }
    Ref<vec> data::lam_inact() { return wsv.segment(inc[5], m_inact); }
    // eq
    Ref<vec> data::we() { return wsv.segment(inc[6], me); }
    Ref<vec> data::_we() { return wsv.segment(inc[7], me); }
    Ref<vec> data::wed() { return wsv.segment(inc[8], me); }
    // ineq
    Ref<vec> data::wi() { return wsv.segment(inc[9], mi); }
    Ref<vec> data::_wi() { return wsv.segment(inc[10], mi); }
    Ref<vec> data::wid() { return wsv.segment(inc[11], mi); }
    Ref<vec> data::omega() { return wsv.segment(inc[12], mi); }
    // dual scales
    Ref<vec> data::invW_inact() { return wsv.segment(inc[13], m_inact); }
    Ref<vec> data::sqInvWLam_inact() { return wsv.segment(inc[14], m_inact); }
    Ref<vec> data::sqWInvLam_inact() { return wsv.segment(inc[15], m_inact); }
    Ref<vec> data::invWO() { return wsv.segment(inc[16], mi); }
    Ref<vec> data::IpObWO() { return wsv.segment(inc[17], mi); }
    Ref<vec> data::sqIpObWO() { return wsv.segment(inc[18], mi); }
    Ref<vec> data::invSqIpObWO() { return wsv.segment(inc[19], mi); }
    // --- steps
    // primal
    Ref<vec> data::dx() { return wsv.segment(inc[20], n); }
    // active
    Ref<vec> data::dlam_act() { return wsv.segment(inc[21], m_act); }
    // inactive
    Ref<vec> data::dw_inact() { return wsv.segment(inc[22], m_inact); }
    Ref<vec> data::dlam_inact() { return wsv.segment(inc[23], m_inact); }
    // eq
    Ref<vec> data::dwe() { return wsv.segment(inc[24], me); }
    // ineq
    Ref<vec> data::dwi() { return wsv.segment(inc[25], mi); }
    Ref<vec> data::domega() { return wsv.segment(inc[26], mi); }


}
