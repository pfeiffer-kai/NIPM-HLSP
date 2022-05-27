#include "nipm_hlsp.h"

#include "hlsp.h"
#include "newton.h"
#include "lagrangian.h"
#include "level.h"
#include "levelShared.h"
#include "data.h"

namespace nipmhlsp
{
    NIpmHLSP::NIpmHLSP()
    {}

    NIpmHLSP::NIpmHLSP(const int p, const int n)
    {
        hp = make_shared<hlsp>(2 * p, n);
        lag.resize(2 * p);
        x = VectorXd::Zero(n);
    }

    void NIpmHLSP::setData(int l_, const MatrixXd& Ae, const VectorXd& be, const MatrixXd& Ai, const VectorXd& bi)
    {
        if (opt.verbose >= MAT)
        {
            std::cout << "setdata NIpmHLSP" << std::endl;
            std::cout << "Ae\n" << Ae << std::endl;
            std::cout << "be\n" << be.transpose() << std::endl;
            std::cout << "Ai\n" << Ai << std::endl;
            std::cout << "bi\n" << bi.transpose() << std::endl;
        }
        int l = 2 * l_;
        if (!hp)
        {
            cout << "Run nipm_hlsp resize first" << endl;
            throw;
        }
        hp->setData(l, Ae, be, Ai, bi);
    }

    /* ================================
     * SOLVE THE HLSP
     * ================================ */
    bool NIpmHLSP::solve()
    {
        if (opt.verbose > NONE) { opt.print(); hp->print(); }

#if TIMEMEASUREMENTS
        auto startHLSP = std::chrono::steady_clock::now();
#endif // TIMEMEASUREMENTS
        iter = 0;
        for (int l = 0; l < hp->p; l++)
        {
#if TIMEMEASUREMENTS
            auto startnewton = std::chrono::steady_clock::now();
#endif // TIMEMEASUREMENTS
            lag[l] = make_shared<lagrangian>(l, hp); // FIXME: there should only be one instance over runtime
            nwt = new newton(x, lag[l]);
            int nrIter = nwt->run(opt.maxIter);
            // virtual priority level
            hp->lvlSh->addActivatedConstraints_inact(l, hp->ws, x);
            // this level
            hp->lvlSh->addActivatedConstraints_l(l+1, hp->ws, x);

#if TIMEMEASUREMENTS
            auto endnewton = std::chrono::steady_clock::now();
            std::chrono::duration<double> differencenewton = endnewton - startnewton;
            if (opt.verbose > NONE) std::cout << "newton Level " << l << " with " << hp->nr << " variables and " << lag[l]->m_inact << " inactive constraints with running time " << differencenewton.count() << "[s] and " << nrIter << " Newton iterations +++" << std::endl;
#endif // TIMEMEASUREMENTS

            iter += nrIter;
            KKT += hp->ws->K("all").norm();

            l++;
            if (hp->nr == 0)
                break;
        }
#if TIMEMEASUREMENTS
        auto endHLSP = std::chrono::steady_clock::now();
        std::chrono::duration<double> differenceHLSP = endHLSP - startHLSP;
        time = differenceHLSP.count();
        if (opt.verbose > NONE) std::cout << "++++++++++++++++++ HLSP solved in " << iter << " Newton iterations and " << differenceHLSP.count() << " [s] ++++++++++++++++++" << std::endl;
#endif // TIMEMEASUREMENTS

        return true;
    }
}

