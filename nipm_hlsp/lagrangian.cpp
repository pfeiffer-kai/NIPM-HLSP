#include "lagrangian.h"

#include "hlsp.h"
#include "level.h"
#include "levelShared.h"
#include "data.h"

namespace nipmhlsp
{
    lagrangian::lagrangian(int _l, shared_ptr<hlsp> _hp) :
        hp(_hp),
        ws(_hp->ws),
        l(_l),
        m(hp->lvls[l]->m),
        me(hp->lvls[l]->me),
        mi(hp->lvls[l]->mi),
        m_act(hp->lvlSh->m_act),
        m_act_l(hp->lvlSh->m_act_p[l]),
        m_inact(hp->lvlSh->m_inact),
        lvl(*_hp->lvls[l]),
        lvlSh(*_hp->lvlSh)
    {

        // re-reference the workspace
        ws->reRef(hp->n, hp->nr, m_act_l, m_inact, me, mi); 

        numThresCur = opt.numThres;
    }

    void lagrangian::init(const vec& x)
    {
        // initialize dual
        level& lvl = *hp->lvls[l];
        levelShared& lvlSh = *hp->lvlSh;

        lvl.we(x, ws);
        lvl.wi(x, ws);
        ws->wi() = ws->wi().cwiseMin(-opt.dualInitThres);
        ws->omega().setConstant(opt.dualInitThres);
        ws->invWO() = (ws->wi() - ws->omega()).cwiseMin(-numThresCur).cwiseInverse();
        ws->IpObWO() = vec::Ones(mi) + ws->invWO().cwiseProduct(ws->omega());
        for (int i = 0; i < mi; i++)
        {
            if (1 + ws->invWO()(i) * ws->omega()(i) > opt.linearDependency)
            {
                ws->sqIpObWO()(i) = sqrt(1 + ws->invWO()(i) * ws->omega()(i));
                ws->invSqIpObWO()(i) = 1 / ws->sqIpObWO()(i);
            }
            else
            {
                ws->sqIpObWO()(i) = 0;
                ws->invSqIpObWO()(i) = 0;
            }
        }
        // inact
        lvlSh._winact(x, ws); 
        ws->w_inact() = ws->_w_inact().cwiseMax(opt.dualInitThres);
        // ws->lam_inact().head(m_inact).setConstant(opt.dualInitThres);
        ws->lam_inact() = ws->_w_inact().cwiseInverse();

        //
        m_act_l = lvlSh.m_act_p[l];
        lvlSh.Lambda_act.col(l / 2).head(m_act_l).setZero();
    }

    void lagrangian::calcMu()
    {
        if (opt.predCor)
        {
            if (mi > 0) mu_l = (double)(-ws->omega().transpose() * ws->wi()) / (double)(mi);
            else mu_l = 1;
            if (m_inact > 0) mu = (double)(ws->lam_inact().transpose() * ws->w_inact()) / (double)(m_inact);
            else mu = 1;
        }
        else
        {
            if (opt.muCalc == 0)
            {
                mu_l = (double)(-ws->omega().transpose() * ws->wi()) / (double)(hp->n + mi);
                mu = (double)(ws->lam_inact().transpose() * ws->w_inact()) / (double)(hp->n + m_inact);
            }
            else if (opt.muCalc == 1)
            {
                mu_l = (double)(-ws->omega().transpose() * ws->wi()) / (double)(hp->n + pow(mi, 1));
                mu = (double)(ws->lam_inact().transpose() * ws->w_inact()) / (double)(hp->n + pow(m_inact, 1));
            }
            else if (opt.muCalc == 2)
            {
                mu_l = (double)(-ws->omega().transpose() * ws->wi()) / (double)(pow(hp->n + mi, 2));
                mu = (double)(ws->lam_inact().transpose() * ws->w_inact()) / (double)(pow(hp->n + m_inact, 2));
            }
        }
    }

    void lagrangian::calcSigma(bool centered)
    {
        if (opt.predCor)
        {
            if (centered)
            {
                sigma = 0;
                sigma_l = 0;
            }
            else // corrector step
            {
                // we have previously calculated the affine step ws->wi()th sigma = 0
                // now the corrector step
                alpha = lineSearch(false);
                if (mi > 0)
                {
                    mu_aff_l = -1. * (double)((ws->omega() + alpha * ws->domega()).transpose() * (ws->wi() + alpha * ws->dwi())) / (double)(mi);
                    // centering parameter
                    sigma_l = min(1., pow(mu_aff_l / max(1e-24, mu_l), opt.predFactor));
                }
                else
                {
                    sigma_l = 0;
                }
                ws->K("omega_i") = ws->omega().cwiseProduct(ws->wi()) + ws->domega().cwiseProduct(ws->dwi()) + sigma_l * mu_l * vec::Ones(mi);
                if (m_inact > 0)
                {
                    mu_aff = (double)((ws->lam_inact() + alpha * ws->dlam_inact()).transpose() * (ws->w_inact() + alpha * ws->dw_inact())) / (double)(m_inact);
                    // centering parameter
                    sigma = min(1., pow(mu_aff / max(1e-24, mu), opt.predFactor));
                }
                else
                {
                    sigma = 0;
                }
                ws->K("w_inact") = ws->lam_inact().cwiseProduct(ws->w_inact()) + ws->dlam_inact().cwiseProduct(ws->dw_inact()) - sigma * mu * vec::Ones(m_inact);
            }
        }
        else
        {
            sigma = 1;
            sigma_l = 1;
        }
    }


    void lagrangian::computeHg(const vec& x, bool centered)
    {
        // compute common variables
        // dual variables
        if (centered)
        {
            if (!opt.predCor) calcMu(); else { mu=0; mu_l=0; };
            calcSigma(centered);
            KKT(x);
        }
        else
        {
            calcMu();
            calcSigma(centered);
        }

        // explicit slacks
        lvl._we(x, ws);
        lvl._wi(x, ws);
        lvlSh._winact(x, ws);

        // auxiliary variables
        ws->invWO() = (ws->wi() - ws->omega()).cwiseMin(-numThresCur).cwiseInverse();
        ws->IpObWO() = vec::Ones(mi) + ws->invWO().cwiseProduct(ws->omega());
        // ws->sqIpObWO() = ws->IpObWO().cwiseSqrt().cwiseMax(numThresCur);
        // ws->invSqIpObWO() = ws->sqIpObWO().cwiseInverse();
        // thresholding
        for (int i = 0; i < mi; i++)
        {
            if (1 + ws->invWO()(i) * ws->omega()(i) > opt.linearDependency)
            {
                ws->sqIpObWO()(i) = sqrt(1 + ws->invWO()(i) * ws->omega()(i));
                ws->invSqIpObWO()(i) = 1 / ws->sqIpObWO()(i);
            }
            else
            {
                ws->sqIpObWO()(i) = 0;
                ws->invSqIpObWO()(i) = 1e12;
            }
        }
        ws->invW_inact() = ws->w_inact().cwiseMax(numThresCur).cwiseInverse();

        if (centered)
        {
            computeH(x);
            computeg(x, centered);
        }
        else
            computeg(x, centered);
    }

    void lagrangian::computeH(const vec& x)
    {
        ws->H().setZero();
        if (opt.leastSquares)
        {
            // inactive constraints
            ws->sqInvWLam_inact() = (ws->lam_inact().cwiseMax(numThresCur).cwiseProduct(ws->w_inact().cwiseMax(numThresCur).cwiseInverse())).cwiseSqrt();
            ws->H().topRows(m_inact) = lvlSh.A_inactN.topRightCorner(m_inact, hp->nr);
            ws->H().topRows(m_inact).applyOnTheLeft(ws->sqInvWLam_inact().asDiagonal());

            // inequality constraints of l
            ws->H().middleRows(m_inact, mi) = lvl.AiN.rightCols(hp->nr);
            ws->H().middleRows(m_inact, mi).applyOnTheLeft(ws->sqIpObWO().asDiagonal());

            // equality constraints of l
            ws->H().bottomRows(me) = lvl.AeN.rightCols(hp->nr);
        }
        else
        {
            if (m_inact > 0) ws->H().triangularView<Eigen::Upper>() +=
                (lvlSh.A_inactN.topRightCorner(m_inact, hp->nr).transpose() * (ws->invW_inact().cwiseProduct(ws->lam_inact())).asDiagonal() * lvlSh.A_inactN.topRightCorner(m_inact, hp->nr));
            if (mi > 0) ws->H().triangularView<Eigen::Upper>() += 
                lvl.AiN.rightCols(hp->nr).transpose() * (vec::Ones(mi) + ws->invWO().cwiseProduct(ws->omega())).asDiagonal() * lvl.AiN.rightCols(hp->nr);
            if (me > 0) ws->H().triangularView<Eigen::Upper>() +=
                lvl.AeN.rightCols(hp->nr).transpose() * lvl.AeN.rightCols(hp->nr);
            // copy for symmetric matrix
            ws->H().triangularView<Eigen::Lower>() = ws->H().triangularView<Eigen::Upper>().transpose();
        }
    }

    void lagrangian::computeg(const vec& x, bool centered)
    {
        // compute the rhs in the normal equations and least-squares case

        ws->g().setZero();

        if (centered)
        {
            if (opt.leastSquares)
            {
                // inactive constraints
                ws->sqInvWLam_inact() = (ws->w_inact().cwiseMax(numThresCur).cwiseProduct(ws->lam_inact().cwiseMax(numThresCur).cwiseInverse())).cwiseSqrt();
                ws->g().head(m_inact) = ws->sqInvWLam_inact().cwiseProduct(ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(ws->_w_inact()) + sigma * mu * vec::Ones(m_inact)));
                // inequality constraints this level
                ws->g().segment(m_inact, mi) = ws->invSqIpObWO().cwiseProduct(-ws->_wi() + ws->omega() - ws->invWO().cwiseProduct(sigma_l * mu_l * vec::Ones(mi) + ws->omega().cwiseProduct(ws->_wi() - ws->omega())));
                // equality constraints this level
                ws->g().segment(m_inact + mi, me) = -ws->_we();
            }
            else
            {
                // normal form

                if (m_inact > 0)
                    ws->g() += lvlSh.A_inactN.block(0, hp->n - hp->nr, m_inact, hp->nr).transpose() * (ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(ws->_w_inact()) + sigma * mu * vec::Ones(m_inact)));
                if (mi > 0)
                    ws->g() += lvl.AiN.rightCols(hp->nr).transpose() * (-ws->_wi() + ws->omega() - ws->invWO().cwiseProduct(sigma_l * mu_l * vec::Ones(mi) + ws->omega().cwiseProduct(ws->_wi() - ws->omega())));
                if (me > 0)
                    ws->g() += lvl.AeN.rightCols(hp->nr).transpose() * (-ws->_we());
            }
        }
        else
        {
            if (opt.leastSquares)
            {
                // inactive constraints
                ws->sqInvWLam_inact() = (ws->w_inact().cwiseMax(numThresCur).cwiseProduct(ws->lam_inact().cwiseMax(numThresCur).cwiseInverse())).cwiseSqrt();
                ws->g().head(m_inact) = ws->sqInvWLam_inact().cwiseProduct(ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(ws->_w_inact()) - ws->dlam_inact().cwiseProduct(ws->dw_inact()) + sigma * mu * vec::Ones(m_inact)));
                // inequality constraints this level
                ws->g().segment(m_inact, mi) = ws->invSqIpObWO().cwiseProduct(-ws->_wi() + ws->omega() - ws->invWO().cwiseProduct(sigma_l * mu_l * vec::Ones(mi) + ws->domega().cwiseProduct(ws->dwi()) + ws->omega().cwiseProduct(ws->_wi() - ws->omega())));
                // equality constraints this level
                ws->g().segment(m_inact + mi, me) = -ws->_we();
            }
            else
            {
                // normal form
                if (m_inact > 0)
                    ws->g() += lvlSh.A_inactN.block(0, hp->n - hp->nr, m_inact, hp->nr).transpose() * (ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(ws->_w_inact()) - ws->dlam_inact().cwiseProduct(ws->dw_inact()) + sigma * mu * vec::Ones(m_inact)));
                if (mi > 0)
                    ws->g() += lvl.AiN.rightCols(hp->nr).transpose() * (-ws->_wi() + ws->omega() - ws->invWO().cwiseProduct(sigma_l * mu_l * VectorXd::Ones(mi) + ws->domega().cwiseProduct(ws->dwi()) + ws->omega().cwiseProduct(ws->_wi() - ws->omega())));
                if (me > 0)
                    ws->g() += lvl.AeN.rightCols(hp->nr).transpose() * (-ws->_we());
            }
        }

    }

    void lagrangian::computeStep(const vec& x, bool centered)
    {
        // turn into format dz = [0, P * [dz 0]] for dx = N dz in lagrangian
        ws->dx().head(hp->n - hp->nr).setZero();
        ws->dx().tail(hp->nr) = ws->g().head(hp->nr);

        lvlSh.applyNSOnTheLeftOnDx(l, ws);
        // if the previous level converged badly then this component ws->will never be below threshold
        // K.segment(nVar + me[l] + mi[l], m_act[l]) = b_act.topRows(m_act[l]) - A_act.topRows(m_act[l]) * x + Lambda_act.col(0).head(m_act[l]);

        // equalities
        if (!opt.predCor || !centered)
        {
            ws->dwe() = lvl.AeN.rightCols(hp->nr) * ws->g().head(hp->nr) - ws->K("lambda_e");
            if (opt.convergenceTest == 1 && m_act_l > 0)
            {
                // calculate lam_act + ws->dlam_act()
                // LagrangeMultipliers_N(l); // leads to worse behaviour somehow
                ws->dlam_act().head(m_act_l).setZero();
                dLagrangeMultipliers(x);
                ws->dlam_act() -= lvlSh.Lambda_act.col(l / 2).head(m_act_l);
            }
        }

        // inequalities, necessary for line search
        ws->domega() = ws->invWO().asDiagonal() * (ws->omega().asDiagonal() * (-lvl.Ai * ws->dx() + ws->K("lambda_i")) - ws->K("omega_i"));
        ws->dwi() = lvl.AiN.rightCols(hp->nr) * ws->g().head(hp->nr) - ws->domega() - ws->K("lambda_i"); 
        ws->dlam_inact() = ws->invW_inact().cwiseProduct(ws->lam_inact().cwiseProduct(ws->K("lambda_inact")) - ws->K("w_inact")
                - ws->lam_inact().cwiseProduct(lvlSh.A_inactN.topRows(m_inact).rightCols(hp->nr) * ws->g().head(hp->nr)));
        ws->dw_inact() = lvlSh.A_inactN.topRows(m_inact).rightCols(hp->nr) * ws->g().head(hp->nr) - ws->K("lambda_inact");

#if SAFEGUARD
        // NaN check
        for (int v = 0; v < hp->nr; v++) if (isnan(ws->g()[v])) { cout << "ERROR: dz NaN "<<endl; throw; }
#endif
    }

    void lagrangian::makeStep(vec& x)
    {
        alpha = lineSearch();

        if (abs(alpha) < opt.alphaThres)
        {
            // should not happen due to the modified line search
#if SAFEGUARD
            cout << "LFLFLFLF ERROR: line search stuck at level " << l  << " with KKTres " << ws->K("all").norm() << ", reset dual and swwitch from opt.muCalc = " << opt.muCalc << " to opt.muCalc = " << max(0, opt.muCalc - 1) << "!!!!!!!!!!!!!!!!!!" << endl;
            throw;
#endif // SAFEGUARD
        }

        // make step
        x += alpha * ws->dx();
        lvlSh.Lambda_act.col(l / 2).head(m_act_l) += alpha * ws->dlam_act().head(m_act_l);
        ws->w_inact() += alpha * ws->dw_inact();
        ws->lam_inact() += alpha * ws->dlam_inact();
        ws->omega() += alpha * ws->domega();
        ws->we() += alpha * ws->dwe();
        ws->wi() += alpha * ws->dwi();
    }

    void lagrangian::dLagrangeMultipliers(const vec& x)
    {
        // actually computes lam_act + ws->dlam_act()
        // FIXME: is there potential for further improvement?

        lvl._we(x + ws->dx(), ws);
        // mi
        if (lvl.mi > 0)
        {
            if (opt.predCor)
            {
                if (sigma == 0)
                {
                    ws->wid() = ws->IpObWO().asDiagonal() * lvl.Ai * ws->dx() + ws->_wi() - ws->omega() + ws->invWO().cwiseProduct(ws->omega().cwiseProduct(ws->_wi() - ws->omega()));
                }
                else
                {
                    ws->wid() = ws->IpObWO().asDiagonal() * lvl.Ai * ws->dx() + ws->_wi() - ws->omega() + ws->invWO().cwiseProduct(sigma_l * mu_l * vec::Ones(mi) + ws->dwi().cwiseProduct(ws->domega()) + ws->omega().cwiseProduct(ws->_wi() - ws->omega()));
                }
            }
            else
            {
                ws->wid() = ws->IpObWO().asDiagonal() * lvl.Ai * ws->dx() + ws->_wi() - ws->omega() + ws->invWO().cwiseProduct(sigma_l * mu_l * vec::Ones(mi) + ws->omega().cwiseProduct(ws->_wi() - ws->omega()));
            }
        }
        // m_inact
        if (m_inact > 0)
        {
            if (opt.predCor)
            {
                if (sigma == 0)
                {
                    ws->w_inactd() = -(ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(lvlSh.A_inact.topRows(m_inact) * (x + ws->dx()) - lvlSh.b_inact.head(m_inact))));
                }
                else
                {
                    ws->w_inactd() = -(ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(lvlSh.A_inact.topRows(m_inact) * (x + ws->dx()) - lvlSh.b_inact.head(m_inact)) - ws->dlam_inact().cwiseProduct(ws->dw_inact()) + sigma * mu * vec::Ones(m_inact)));
                }
            }
            else
            {
                ws->w_inactd() = -(ws->lam_inact() + ws->invW_inact().cwiseProduct(-ws->lam_inact().cwiseProduct(lvlSh.A_inact.topRows(m_inact) * (x + ws->dx()) - lvlSh.b_inact.head(m_inact)) + sigma * mu * vec::Ones(m_inact)));
            }
        }

        // dual, Dimitrov, only Matrix Vector multiplications
        int curRank = 0;
        if (l > 0)
            curRank = lvlSh.rank_p.head(l - 1).sum();
        for (int ll = l - 1; ll >= 1; --ll)
        {
            if (lvlSh.m_act_p[ll + 1] - lvlSh.m_act_p[ll] > 0)
            {
                ws->dlam_act().segment(lvlSh.m_act_p[ll], lvlSh.rank_p[ll]) = lvl.AeN.middleCols(curRank, lvlSh.rank_p[ll]).transpose() * ws->_we()
                    + lvl.AiN.middleCols(curRank, lvlSh.rank_p[ll]).transpose() * ws->wid()
                    + lvlSh.A_inactN.block(0, curRank, m_inact, lvlSh.rank_p[ll]).transpose() * ws->w_inactd();
                ws->dlam_act().segment(lvlSh.m_act_p[ll], lvlSh.rank_p[ll]) -= lvlSh.A_actN.block(lvlSh.m_act_p[ll + 1], curRank, lvlSh.m_act_p[l] - lvlSh.m_act_p[ll + 1], lvlSh.rank_p[ll]).transpose() *
                    ws->dlam_act().segment(lvlSh.m_act_p[ll + 1], lvlSh.m_act_p[l] - lvlSh.m_act_p[ll + 1]);
                ws->dlam_act().segment(lvlSh.m_act_p[ll] + lvlSh.rank_p[ll], lvlSh.m_act_p[ll + 1] - lvlSh.m_act_p[ll] - lvlSh.rank_p[ll]).setZero();
                // apply QT
                ws->dlam_act().segment(lvlSh.m_act_p[ll], lvlSh.m_act_p[ll + 1] - lvlSh.m_act_p[ll]).applyOnTheLeft(householderSequence(lvlSh.A_actN.block(lvlSh.m_act_p[ll], curRank, lvlSh.m_act_p[ll + 1] - lvlSh.m_act_p[ll], lvlSh.rank_p[ll]), lvlSh.hCoeff[ll]));
            }
            if (ll - 1 >= 0)
                curRank -= lvlSh.rank_p[ll - 1];
        }

#if SAFEGUARD
        // test dLagrange multipliers
        VectorXd dLagrangeTest = lvl.Ae.transpose() * ws->_we()
            + lvl.Ai.transpose() * ws->wid()
            + lvlSh.A_inact.topRows(m_inact).transpose() * ws->w_inactd();
        dLagrangeTest -= lvlSh.A_act.topRows(lvlSh.m_act_p[l]).transpose() * ws->dlam_act().head(lvlSh.m_act_p[l]);

        if (dLagrangeTest.norm() > opt.LagrangeTestThres)
        {
            if (opt.verbose > NONE) cout << "MFMFMFMF dLagrange multiplier test FAILED; error norm: " << dLagrangeTest.norm() << endl; // << dLagrangeTest << endl;
            // throw;
        }
#endif
    }


    void lagrangian::KKT(const vec& x, bool full)
    {
        if (full)
            ws->K("x") = lvl.Ae.transpose() * ws->we() + lvl.Ai.transpose() * ws->wi() - lvlSh.A_act.topRows(m_act_l).transpose() * lvlSh.Lambda_act.col(l / 2).head(m_act_l) - lvlSh.A_inact.topRows(m_inact).transpose() * ws->lam_inact().head(m_inact);
        ws->K("lambda_e") = lvl.be - lvl.Ae * x + ws->we();
        ws->K("lambda_act").setZero();
        // if the previous level converged badly then this component will never be below threshold because of the disturbed x (not because of not computed explicit ws->_wi())
        // ws->K("lambda_act") = lvlSh.b_act.head(m_act_l) - lvlSh.A_act.topRows(m_act_l) * x + Lambda_act.col(0).head(m_act_l);
        ws->K("lambda_inact") = lvlSh.b_inact.head(m_inact) - lvlSh.A_inact.topRows(m_inact) * x + ws->w_inact();
        ws->K("w_inact") = ws->lam_inact().cwiseProduct(ws->w_inact()) - sigma * mu * vec::Ones(m_inact);
        ws->K("lambda_i") = lvl.bi - lvl.Ai * x + ws->wi() + ws->omega().head(mi);
        ws->K("omega_i") = ws->omega().cwiseProduct(ws->wi()) + sigma_l * mu_l * vec::Ones(mi);
    }


    double lagrangian::lineSearch(bool reinitDual)
    {
        double _alpha = 1.;
        double _alpha_ = 1.;
        // on l
        for (int i = 0; i < mi; i++)
        {
            if (ws->omega()[i] + ws->domega()[i] < opt.posThres) // if the delta leads to a negative value
            {
                _alpha_ = min(_alpha, (opt.posThres - ws->omega()[i]) / ws->domega()[i]);
                if (_alpha_ < opt.alphaThres)
                {
                    // reset the dual step if the would lead to a stuck line search
                    if (reinitDual)
                    {
                        ws->domega()[i] = opt.dualInitThres;
                    }
                    else
                    {
                        // cout << "reset ws->domega() of ctr " << i << endl;
                        ws->domega()[i] = 0;
                        // omega_l[i] = 0;
                    }
                }
                else
                {
                    _alpha = _alpha_;
                }
            }
            if (ws->wi()[i] + ws->dwi()[i] > opt.posThres) // if the delta leads to positive value
            {
                _alpha_ = min(_alpha, (opt.posThres - ws->wi()[i]) / ws->dwi()[i]);
                if (_alpha_ < opt.alphaThres)
                {
                    if (reinitDual)
                    {
                        ws->dwi()[i] = -opt.dualInitThres;
                    }
                    else
                    {
                        // cout << "reset dw of ctr " << i << " dw " << dw_l[me[l] + i] << " w_l " << w_l[me[l] + i] << endl;
                        ws->dwi()[i] = 0;
                        // w_l[me[l] + i] = 0;
                    }
                }
                else
                {
                    _alpha = _alpha_;
                }
            }
        }

        // inact
        for (int i = 0; i < m_inact; i++)
        {
            if (ws->lam_inact()[i] + ws->dlam_inact()[i] < opt.posThres) // if the delta leads to a negative value
            {
                _alpha_ = min(_alpha, (opt.posThres - ws->lam_inact()[i]) / ws->dlam_inact()[i]);
                if (_alpha_ < opt.alphaThres)
                {
                    // reset the dual step if the would lead to a stuck line search
                    if (reinitDual)
                    {
                        ws->dlam_inact()[i] = opt.dualInitThres;
                    }
                    else
                    {
                        ws->dlam_inact()[i] = 0;
                    }
                }
                else
                {
                    _alpha = _alpha_;
                }
            }
            if (ws->w_inact()[i] + ws->dw_inact()[i] < opt.posThres)
            {
                _alpha_ = min(_alpha, (opt.posThres - ws->w_inact()[i]) / ws->dw_inact()[i]);
                if (_alpha_ < opt.alphaThres)
                {
                    if (reinitDual)
                    {
                        ws->dw_inact()[i] = opt.dualInitThres;
                    }
                    else
                    {
                        ws->dw_inact()[i] = 0;
                    }
                }
                else
                {
                    _alpha = _alpha_;
                }
            }
        }
        // cout << "alpha " << _alpha << endl;
        return _alpha;
    }


    status lagrangian::convergenceTest(const vec& x)
    {
        // this Newton iteration is finished, check for convergence
        if (opt.convergenceTest == 1)
        {
            calcMu();
            KKT(x, true);
            if (opt.verbose > NONE) printVars(x);
            if (ws->K("all").norm() < opt.KKT_thres)
            {
                convergenceNorm = ws->K("all").norm();
                if (convergenceNorm > maxKKTnorm)
                    maxKKTnorm = convergenceNorm;
                return SUCCESS;
            }
        }
        else if (opt.convergenceTest == 2 || opt.convergenceTest == 3)
        {
            cout << "convergenceTest 2, 3 NIY" << endl; throw;
            if (ws->K("all_inact").norm() < 1e-6)
            {
                if (m_act_l > 0)
                {
                    // LagrangeMultipliers_N(l); // NIY
                    calcMu();
                    KKT(x, true); 
                }
                if (ws->K("all").norm() < opt.KKT_thres)
                {
                    convergenceNorm = ws->K("all").norm();
                    if (convergenceNorm > maxKKTnorm)
                        maxKKTnorm = convergenceNorm;
                    return SUCCESS;
                }
            }
        }
        return OK;
    }

    void lagrangian::printStatus(const vec& x)
    {
        cout << "--- Lagrangian status: level: " << l << ", remaining variables: " << hp->nr << " of " << hp->n << ", inactive constraints: " << m_inact << ", active constraints: " << m_act_l << ", inequality constraints: " << mi << ", equality constraints: " << mi << ", KKT residual: " << ws->K("all").norm();
        cout << " ---" << endl;
    }

    void lagrangian::printVars(const vec& x)
    {
        // // calculate full KKT
        // KKT(x, true);
        // calculate the corresponding KKT residual
        if (opt.verbose >= SOLVE)
        {
            cout << " - - - lagrangian::printVars - - -\n";
            cout << "x " << x.transpose() << endl;
            cout << "-> dx " << ws->dx().transpose() << endl;
            if (m_act_l > 0)
            {
                cout << "ACT Lambda_act " << lvlSh.Lambda_act.col(l / 2).head(lvlSh.m_act_p[l]).transpose() << endl;
                cout << "ACT -> dlam_act " << ws->dlam_act().transpose() << endl;
                // cout << "w_act " << lvlSh.Lambda_act.col(0).head(m_act[l]).transpose() << endl;
                // cout << "w_act as A_actx - b_act " << (A_act.topRows(m_act[l]) * x - b_act.topRows(m_act[l])).transpose() << endl;
            }
            if (m_inact > 0)
            {
                cout << "INACT lam_inact " << ws->lam_inact().transpose() << endl;
                cout << "INACT -> dlam_inact " << ws->dlam_inact().transpose() << endl;
                cout << "INACT w_inact " << ws->w_inact().transpose() << endl;
                cout << "INACT ws->dw_inact " << ws->dw_inact().transpose() << endl;
                cout << "INACT mu " << mu << "sigma " << sigma << endl;
            }
            if (me > 0)
            {
                cout << "EQ we " << ws->we().transpose() << endl;
                cout << "EQ -> dwe " << ws->dwe().transpose() << endl;
            }
            if (mi > 0)
            {
                cout << "IQ wi " << ws->wi().transpose() << endl;
                cout << "IQ dwi " << ws->dwi().transpose() << endl;
                cout << "IQ omega " << ws->omega().transpose() << endl;
                cout << "IQ domega " << ws->domega().transpose() << endl;
                cout << "IQ mu_l " << mu_l << " sigma_l " << sigma_l << endl;
            }
            cout << "LINESEARCH: v -> alpha dv with alpha = " << alpha << endl;
        }

        if (opt.verbose > NONE)
        {
            cout << "KKT ||KKT||_2: " << ws->K("all").norm() << endl;
            if (opt.verbose >= SOLVE)
            {
                cout << "KKT ||x||_2: " << ws->K("x").norm() << endl;
                cout << "KKT ||lambda_e||_2: " << ws->K("lambda_e").norm() << endl;
                cout << "KKT ||lambda_i||_2: " << ws->K("lambda_i").norm() << endl;
                cout << "KKT ||omega_i||_2: " << ws->K("omega_i").norm() << endl;
                cout << "KKT ||lambda_act||_2: " << ws->K("lambda_act").norm() << endl;
                cout << "KKT ||lambda_inact||_2: " << ws->K("lambda_inact").norm() << endl;
                cout << "KKT ||w_inact()||_2: " << ws->K("w_inact").norm() << endl;
            }
        }
        if (opt.verbose >= MAT) ws->print();
        if (opt.verbose >= SOLVE) cout << " - - -   - - -   - - -    - - -\n";
    }
} // namespace nipmhlsp
