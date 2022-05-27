#include "newton.h"

#include "lagrangian.h"
#include "hlsp.h"
#include "data.h"

namespace nipmhlsp
{
    newton::newton(vec& x0, shared_ptr<lagrangian> _lag) : lag(_lag), x(x0), ws(_lag->hp->ws) 
    { 
        lag->init(x);
    }

    int newton::run(int nrSteps)
    {
        for (iter=0; iter<nrSteps; iter++)
        {
            if (opt.verbose >= SOLVE) cout << "==== Newton iteration " << iter << " with current x\n" << x.transpose() << endl;

            if (opt.predCor)
            {
                // centered step
                step(true);
                // corrector step
                step(false);
            }
            else
                step();

            if (opt.verbose > NONE) { cout << "l " << lag->l << " newton iter " << iter << ": "; };
            status s = lag->convergenceTest(x);

            if (s == SUCCESS)
            {
                if (opt.verbose >= CONV) printConvergenceMessage(s);
                iter++;
                return iter;
            }
        }

        if (opt.verbose >= CONV) printConvergenceMessage(EXCEED);

        iter++;
        return iter;
    }

    bool newton::step(bool centered)
    {
        // compute the gradient and Hessian of the Lagrangian
        lag->computeHg(x, centered);

        // solve the linear system
        solve(centered);

        // compute and make step with the calculated NS increment
        lag->computeStep(x, centered);
        if (!centered || !opt.predCor)
            lag->makeStep(x);

        return true;
    }

    bool newton::solve(bool centered)
    {
        // solve rank deficient normal equations with rank revealing QR decomposition in place
        if (opt.verbose >= MAT) cout << "g\n" << ws->g().transpose() << endl;

        int rg = lag->hp->nr;
        if (opt.leastSquares)
            rg = lag->m_inact+lag->m;

#if SAFEGUARD
        mat _H = ws->H();
        vec _g = ws->g().head(rg);
#endif

        Eigen::Ref<mat> lhs = ws->H();
        if (centered)
        {
            if(opt.verbose >= MAT) cout << "H\n" << lhs << endl;
            Eigen::ColPivHouseholderQR<Eigen::Ref<mat> > qr(lhs);
            // if (opt.predCor && lag->mi + lag->m_inact > 0)
            // {
            // save qr data for corrector solve
            rank = qr.rank();
            hCoeffs = qr.hCoeffs();
            P = qr.colsPermutation();
            // }
            // else
            // {
            //     // apply Q
            //     lag->g.head(lag->hp->nr).applyOnTheLeft(householderSequence(lhs.topLeftCorner(lag->hp->nr, lag->hp->nr), qr.hCoeffs()).transpose());
            //     // solve 
            //     lhs.topLeftCorner(qr.rank(), qr.rank())
            //         .triangularView<Eigen::Upper>()
            //         .solveInPlace<Eigen::OnTheLeft>(lag->g.head(qr.rank()));
            //     // apply P
            //     lag->g.segment(qr.rank(), lag->hp->nr - qr.rank()).setZero();
            //     lag->g.head(lag->hp->nr).applyOnTheLeft(qr.colsPermutation());
            // }
        }

        // corrector step
        // if (opt.predCor && lag->mi + lag->m_inact > 0) 
        // {
        // apply Q
        ws->g().head(rg).applyOnTheLeft(householderSequence(ws->H().topRows(rg), hCoeffs).transpose());
        // solve 
        ws->H().topLeftCorner(rank, rank)
            .triangularView<Eigen::Upper>()
            .solveInPlace<Eigen::OnTheLeft>(ws->g().head(rank));
        // apply P
        ws->g().segment(rank, lag->hp->nr - rank).setZero();
        ws->g().head(lag->hp->nr).applyOnTheLeft(P);
        // }

        if (opt.verbose >= SOLVE) cout << "NEWTON dz\n" << ws->g().head(lag->hp->nr).transpose() << endl;
#if SAFEGUARD
        if (opt.verbose >= SOLVE) cout << "NEWTON: solve test " << (_H * ws->g().head(lag->hp->nr) - _g).norm() << endl;
#endif

        return true;
    }

    void newton::printConvergenceMessage(status s)
    {
        if (s == FAILURE) cout << "++++++++++++++++ NEWTON FAILURE at iteration " << iter << ": unknown failure with ";
        else if (s == EXCEED) cout << "++++++++++++++++ NEWTON EXCEED at iteration " << iter << ": maximum iteration exceeded with ";
        else if (s == OK) cout << "++++++++++++++++ NEWTON OK at iteration " << iter << ": everything ok with ";
        else if (s == SUCCESS) cout << "++++++++++++++++ NEWTON SUCCESS at iteration " << iter << ": convergence with ";
        lag->printStatus(x);
    }
} // namespace nipmhlsp
