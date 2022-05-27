#include "levelShared.h"

#include "lagrangian.h"
#include "hlsp.h"
#include "level.h"
#include "data.h"

namespace nipmhlsp
{
    levelShared::levelShared(shared_ptr<hlsp> _hp, int p, int n) :
        hp(_hp)
    {
        m_act_p.resize(p+1);
        m_act_p.setZero();
        
        permLvl.resize(p);
        rank_p.resize(p); rank_p.setZero();
        hCoeff.resize(p);

        m_activated = 0;

        int i=0;
        m = 0;
        me = 0;
        mi = 0;
        mp.resize(hp->lvls.size() + 1);
        mp.setZero();
        for (shared_ptr<level> lv : hp->lvls)
        {
            if (!lv)
            {
                cout << "levelShared:: level not initialized" << endl;
                throw;
            }
            m += lv->m;
            me += lv->me;
            mi += lv->mi;
            mp[i+1] = mp[i] + m;
            i++;
        }

        A_act.resize(m, n);
        A_actN.resize(m, n);
        b_act.resize(m);
        A_inact.resize(mi, n);
        A_inactN.resize(mi, n);
        b_inact.resize(mi);

        Lambda_act.resize(m, p); Lambda_act.setZero();

        idx_act = mati::Zero(mi, 4); // level, row in A, virtual priority level, row in A_act
        idx_inact = mati::Zero(mi, 2); // level and row in A
    }

    void levelShared::applyNSOnTheLeftOnDx(int l, shared_ptr<data> ws)
    {
        int curRank = 0;

        if (l > 0) curRank = rank_p.head(l - 1).sum();

        for (int ll = l - 1; ll >= 0; ll--)
        {
            int rll = rank_p[ll];
            int m_act_ll = m_act_p[ll];
            // T_ll * x
            ws->dx().segment(curRank, rll) = -A_actN.block(m_act_ll, curRank + rll, rll, hp->n - curRank - rll) * ws->dx().tail(hp->n - curRank - rll);
            // R_ll^-1 * x
            A_actN.block(m_act_ll, curRank, rll, rll)
                .triangularView<Eigen::Upper>()
                .solveInPlace<Eigen::OnTheLeft>(ws->dx().segment(curRank, rll));
            // P_ll * dx
            ws->dx().tail(hp->n - curRank).applyOnTheLeft(permLvl[ll]);
            if (ll > 0)
            {
                curRank -= rank_p[ll - 1];
            }
        }
    }

    void levelShared::applyNSOnTheLeftOnx(int l, vec& x)
    {
        int curRank = 0;

        if (l > 0) curRank = rank_p.head(l - 1).sum();

        for (int ll = l - 1; ll >= 0; ll--)
        {
            int rll = rank_p[ll];
            int m_act_ll = m_act_p[ll];
            // T_ll * x
            x.segment(curRank, rll) = -A_actN.block(m_act_ll, curRank + rll, rll, hp->n - curRank - rll) * x.tail(hp->n - curRank - rll);
            // R_ll^-1 * x
            A_actN.block(m_act_ll, curRank, rll, rll)
                .triangularView<Eigen::Upper>()
                .solveInPlace<Eigen::OnTheLeft>(x.segment(curRank, rll));
            // P_ll * x
            x.tail(hp->n - curRank).applyOnTheLeft(permLvl[ll]);
            if (ll > 0)
            {
                curRank -= rank_p[ll - 1];
            }
        }
    }

    void levelShared::addActivatedConstraints_inact(const int l, shared_ptr<data> ws, vec& x)
    {
        // add activated constraints from inactive ones
        // adds a virtual priority level
        
        // handle activated inequality constraints from lower priority levels
        m_act_p[l + 1] = m_act_p[l];

        // Check whether inequality constraints from higher priority levels got activated
        int m_inact_ = m_inact;
        int i_inact = 0;

        // improves numerical stability
        winact(x, ws);

        while (i_inact < m_inact_)
        {
            // how to achieve highest accuracy with fewest number of Newton iterations?
            if (ws->w_inact()[i_inact] < opt.activationThres && ws->lam_inact()[i_inact] > opt.activationThres)
            {
                if (opt.verbose >= SOLVE) cout << "activate inactive constraint " << i_inact << " from inact (constraint from level " << idx_act(i_inact,0) << ", row in A: " << idx_inact(i_inact,1) << ") with w_inact[i_inact] " << ws->w_inact()[i_inact] << " and lam_inact[i_inact] " << ws->lam_inact()[i_inact] << std::endl;
                A_act.row(m_act_p[l + 1]) = A_inact.row(i_inact); // copy including the block operation
                A_actN.row(m_act_p[l + 1]) = A_inactN.row(i_inact); // copy including the block operation
                b_act(m_act_p[l + 1]) = b_inact(i_inact) + ws->_w_inact()(i_inact);
                A_inact.row(i_inact).swap(A_inact.row(m_inact_-1));
                A_inactN.row(i_inact).swap(A_inactN.row(m_inact_-1));
                Lambda_act(m_act_p[l + 1], l / 2);
                ws->lam_inact()[i_inact];
                Lambda_act(m_act_p[l + 1], l / 2) = ws->lam_inact()[i_inact];
                std::swap(b_inact(i_inact), b_inact(m_inact_ - 1));
                std::swap(ws->_w_inact()(i_inact), ws->_w_inact()(m_inact_ - 1));
                std::swap(ws->lam_inact()(i_inact), ws->lam_inact()(m_inact_ - 1));
                idx_act.row(m_activated).head(2) = idx_inact.row(i_inact);
                idx_act(m_activated, 2) = l;
                idx_act(m_activated, 3) = m_act_p[l + 1];
                idx_inact.row(i_inact) = idx_inact.row(m_inact_ - 1);
                idx_inact.row(m_inact_ - 1).setZero();
                m_act_p[l + 1]++;
                m_act++;
                m_activated++;
                m_inact_--;
            }
            else
                i_inact++;
        }

        m_inact = m_inact_;

        // qr decomposition of said constraints and projection of lower priority levels
        // the Lagrange multipliers are zero because for the activated constraints w_inact = 0 so need not to be recalculated
        project(l, true); 
    }

    void levelShared::addActivatedConstraints_l(const int l, shared_ptr<data> ws, vec& x)
    {
        // continue active set for THIS level
        
        level& lvl = *hp->lvls[l-1];

        m_act_p[l + 1] = m_act_p[l];
        // equality constraints
        A_act.middleRows(m_act_p[l + 1], lvl.me) = lvl.Ae;
        A_actN.middleRows(m_act_p[l + 1], lvl.me) = lvl.AeN;
        b_act.segment(m_act_p[l + 1], lvl.me) = lvl.be;
        Lambda_act.col((l - 1) / 2).segment(m_act_p[l + 1], lvl.me) = ws->we();
        Lambda_act.col(0).segment(m_act_p[l + 1], lvl.me) = ws->we();
        // calculate here so we activate the correct constraints
        m_act_p[l + 1] += lvl.me;
        m_act += lvl.me;
        if (opt.verbose >= SOLVE) cout << "activate " << lvl.me << " eq constraints from level " << l << " with we " << ws->we().transpose() << std::endl;
        // check the activity of inequality constraints on this level
        // recalculate wi
        lvl.wi(x, ws);
        // FIXME: is the amount of copying avoidable
        for (int i = 0; i < lvl.mi; i++)
        {
                if (ws->_wi()[i] < -opt.activationThres && ws->omega()[i] < opt.activationThres)
                // if (w_l[me[l - 1] + i] < -opt.activationThres)
                {
                    if (opt.verbose >= SOLVE) cout << "activate constraint " << i << " from level " << l << " with wi " << ws->wi()[i] << " and omega " << ws->omega()[i] << std::endl;
                    A_act.row(m_act_p[l + 1]) = lvl.Ai.row(i);
                    A_actN.row(m_act_p[l + 1]) = lvl.AiN.row(i);
                    b_act[m_act_p[l + 1]] = lvl.bi[i] + ws->_wi()[i];
                    Lambda_act(m_act_p[l + 1], (l - 1) / 2) = ws->_wi()[i]; // or wi
                    Lambda_act(m_act_p[l + 1], 0) = ws->_wi()[i];
                    idx_act(m_activated, 0) = (l - 1) / 2;
                    idx_act(m_activated, 1) = mp[l - 1] + lvl.me + i;
                    idx_act(m_activated, 2) = l;
                    idx_act(m_activated, 3) = m_act_p[l + 1];
                    m_activated++;
                    m_act_p[l + 1]++;
                    m_act++;
                }
                else
                {
                    if (opt.verbose >= SOLVE) cout << "add constraint " << i << " from this level to inact with wi " << ws->wi()[i] << " and omega " << ws->omega()[i] << std::endl;
                    A_inact.row(m_inact) = lvl.Ai.row(i);
                    A_inactN.row(m_inact) = lvl.AiN.row(i);
                    b_inact(m_inact) = lvl.bi[i];
                    idx_inact(m_inact, 0) = (l - 1) / 2;
                    idx_inact(m_inact, 1) = mp[l - 1] + lvl.me + i;
                    m_inact++;
                }
        }
        // qr decomposition of said constraints, projection of lower priority levels and recalculation of Lagrange multipliers
        if (hp->nr > 0)
            project(l);
    }

    void levelShared::project(const int l, const bool virtualLvl)
    {
        // Nullspace calculation of active constraints and projection of lower priority matrices into it
        // for safety here
        permLvl[l].resize(hp->nr);
        permLvl[l].setIdentity();
        if (m_act_p[l + 1] - m_act_p[l] > 0)
        {
            // QR decomposition of A_actN; theoretically, it is already available from the Newton's method but with an error depending on the magnitude of the KKT conditions when the Newton's method is stopped
            mat mat2save = A_actN.block(m_act_p[l], hp->n - hp->nr, m_act_p[l + 1] - m_act_p[l], hp->nr);
            Eigen::Ref<mat> mat2dec = A_actN.block(m_act_p[l], hp->n - hp->nr, m_act_p[l + 1] - m_act_p[l], hp->nr);
            Eigen::ColPivHouseholderQR<Eigen::Ref<mat> > qr(mat2dec);

            // save the data for dual computation
            permLvl[l] = qr.colsPermutation();
            rank_p[l] = qr.rank();
            hCoeff[l] = qr.hCoeffs();

#if SAFEGUARD
            // nullspace test
            mat2save.rightCols(hp->nr).applyOnTheRight(permLvl[l]);
            mat2dec.topLeftCorner(rank_p[l], rank_p[l])
                .triangularView<Eigen::Upper>()
                .solveInPlace<Eigen::OnTheRight>(mat2save.leftCols(rank_p[l]));
            mat2save.rightCols(hp->nr - rank_p[l]) -= mat2save.leftCols(rank_p[l]) *
                mat2dec.block(0, rank_p[l], rank_p[l], hp->nr - rank_p[l]);
            if (mat2save.rightCols(hp->nr - rank_p[l]).norm() > 1e-12)
            {
                cout << "lvlShared::project: NS corrupted with error norm " << mat2save.rightCols(hp->nr - rank_p[l]).norm() << endl;
                throw;
            }
#endif // SAFEGUARD

            // inactive constraints
            A_inactN.rightCols(hp->nr).topRows(m_inact).applyOnTheRight(permLvl[l]);
            mat2dec.topLeftCorner(rank_p[l], rank_p[l])
                .triangularView<Eigen::Upper>()
                .solveInPlace<Eigen::OnTheRight>(A_inactN.block(0, hp->n - hp->nr, m_inact, rank_p[l]));
            A_inactN.block(0, hp->n-hp->nr + rank_p[l], m_inact, hp->nr - rank_p[l]) -= A_inactN.block(0, hp->n-hp->nr, m_inact, rank_p[l]) *
                mat2dec.block(0, rank_p[l], rank_p[l], hp->nr - rank_p[l]);
            
            for (int i = l + 1; i < hp->p; i++)
            {
                if (i % 2 == 0)
                {
                    level& lvl = *hp->lvls[i];
                    // lower level equality constraints
                    lvl.AeN.rightCols(hp->nr).applyOnTheRight(permLvl[l]);
                    mat2dec.topLeftCorner(rank_p[l], rank_p[l])
                        .triangularView<Eigen::Upper>()
                        .solveInPlace<Eigen::OnTheRight>(lvl.AeN.middleCols(hp->n - hp->nr, rank_p[l]));
                    lvl.AeN.middleCols(hp->n-hp->nr + rank_p[l], hp->nr - rank_p[l]) -= lvl.AeN.middleCols(hp->n-hp->nr, rank_p[l]) *
                        mat2dec.block(0, rank_p[l], rank_p[l], hp->nr - rank_p[l]);

                    // lower level iequality constraints
                    lvl.AiN.rightCols(hp->nr).applyOnTheRight(permLvl[l]);
                    mat2dec.topLeftCorner(rank_p[l], rank_p[l])
                        .triangularView<Eigen::Upper>()
                        .solveInPlace<Eigen::OnTheRight>(lvl.AiN.middleCols(hp->n - hp->nr, rank_p[l]));
                    lvl.AiN.middleCols(hp->n-hp->nr + rank_p[l], hp->nr - rank_p[l]) -= lvl.AiN.middleCols(hp->n-hp->nr, rank_p[l]) *
                        mat2dec.block(0, rank_p[l], rank_p[l], hp->nr - rank_p[l]);
                }
            }
        }
        hp->nr -= rank_p[l];
    }

    void levelShared::wact(const vec& x, shared_ptr<data> ws) { if (m_act > 0) Lambda_act.col(0).head(m_act) = A_act.topRows(m_act) * x - b_act.head(m_act); }
    void levelShared::winact(const vec& x, shared_ptr<data> ws) { if (m_inact > 0) ws->w_inact() = A_inact.topRows(m_inact) * x - b_inact.head(m_inact); }
    void levelShared::winactN(const vec& z, shared_ptr<data> ws) { if (m_inact > 0) ws->w_inact() = A_inactN.topRows(m_inact) * z - b_inact.head(m_inact); }
    void levelShared::_winact(const vec& x, shared_ptr<data> ws) { if (m_inact > 0) ws->_w_inact() = A_inact.topRows(m_inact) * x - b_inact.head(m_inact); }
    void levelShared::_winactN(const vec& z, shared_ptr<data> ws) { if (m_inact > 0) ws->_w_inact() = A_inactN.topRows(m_inact) * z - b_inact.head(m_inact); }
}
