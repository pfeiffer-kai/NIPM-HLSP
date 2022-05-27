#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

#include <nipm_hlsp/nipm_hlsp.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1> >::Index Index;

int main()
{
    int nVar = 0;
    int p = 0;
    bool randomMatrix = false;
    VectorXi me;
    VectorXi mi;

    for (int testnr = 0; testnr < 17; testnr++)
    {
        cout << "\n============= RUNNING TEST " << testnr << " ... " << endl;
        if (testnr == 0)
        {
            nVar = 2;
            p = 1;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            randomMatrix = true;
        }
        else if (testnr == 1)
        {
            nVar = 4;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            me[1] = 4;
            randomMatrix = true;
        }
        else if (testnr == 2)
        {
            nVar = 1;
            p = 1;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            mi[0] = 1;
            randomMatrix = false;
        }
        else if (testnr == 3)
        {
            nVar = 1;
            p = 1;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            mi[0] = 2;
            randomMatrix = false;
        }
        else if (testnr == 4)
        {
            nVar = 1;
            p = 1;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 1;
            mi[0] = 1;
            randomMatrix = false;
        }
        else if (testnr == 5)
        {
            nVar = 2;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[1] = 2;
            mi[0] = 1;
            randomMatrix = false;
        }
        else if (testnr == 6)
        {
            nVar = 4;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[1] = 4;
            mi[0] = 2;
            randomMatrix = false;
        }
        else if (testnr == 7)
        {
            nVar = 4;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 1;
            me[1] = 4;
            mi[0] = 2;
            mi[1] = 1;
            randomMatrix = false;
        }
        else if (testnr == 8)
        {
            nVar = 3;
            p = 3;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 1;
            me[1] = 2;
            me[2] = 3;
            randomMatrix = false;
        }
        else if (testnr == 9)
        {
            nVar = 6;
            p = 3;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            me[1] = 4;
            me[2] = 6;
            randomMatrix = false;
        }
        else if (testnr == 10)
        {
            nVar = 6;
            p = 3;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 1;
            me[1] = 3;
            me[2] = 6;
            mi[0] = 2;
            mi[1] = 4;
            mi[2] = 6;
            randomMatrix = false;
        }
        else if (testnr == 11)
        {
            nVar = 2;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[1] = 2;
            randomMatrix = false;
        }
        else if (testnr == 12)
        {
            nVar = 4;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            me[1] = 4;
            randomMatrix = true;
        }
        else if (testnr == 13)
        {
            nVar = 6;
            p = 3;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            me[1] = 4;
            me[2] = 6;
            randomMatrix = true;
        }
        else if (testnr == 14)
        {
            nVar = 2;
            p = 1;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            me[0] = 2;
            mi[0] = 1;
            randomMatrix = true;
        }
        else if (testnr == 15)
        {
            nVar = 4;
            p = 2;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            mi[0] = 2;
            me[0] = 1;
            mi[1] = 2;
            me[1] = 4;
            randomMatrix = true;
        }
        else if (testnr == 16)
        {
            nVar = 9;
            p = 3;
            me = VectorXi::Zero(p);
            mi = VectorXi::Zero(p);
            mi[0] = 2;
            me[0] = 1;
            mi[1] = 2;
            me[1] = 3;
            mi[2] = 2;
            me[2] = 9;
            randomMatrix = true;
        }
        

        int m_all = 0;
        VectorXi ml = VectorXi::Zero(p+1);
        for (int l = 0; l < p; l++)
        {
            ml[l + 1] = ml[l] + me[l] + mi[l];
            m_all += me[l] + mi[l];
        }

        MatrixXd A = MatrixXd::Zero(m_all, nVar);
        VectorXd b = VectorXd::Zero(m_all);

        if (!randomMatrix)
        {
            if (p > 0)
            {
                if (testnr == 3)
                {
                    A(0,0) = -1; b[0] = 0;
                    A(1,0) = 1; b[1] = -0;
                }
                else
                {
                    A.block(0, 0, me[0], me[0]) = MatrixXd::Identity(me[0], me[0]);
                    b.head(me[0]) = 1 * VectorXd::Ones(me[0]);
                    A.block(me[0], 0, mi[0], mi[0]) = -MatrixXd::Identity(mi[0], mi[0]);
                    b.segment(me[0],mi[0]) = -1 * VectorXd::Ones(mi[0]);
                }
            }

            if (p > 1)
            {
                A.block(ml[1], 0, me[1], me[1]) = MatrixXd::Identity(me[1], me[1]);
                b.segment(ml[1],me[1]) = 10 * VectorXd::Ones(me[1]);
                A.block(ml[1] + me[1], 0, mi[1], mi[1]) = -MatrixXd::Identity(mi[1], mi[1]);
                b.segment(ml[1] + me[1],mi[1]) = -2 * VectorXd::Ones(mi[1]);
            }

            if (p > 2)
            {
                A.block(ml[2], 0, me[2], me[2]) = MatrixXd::Identity(me[2], me[2]);
                b.segment(ml[2],me[2]) = 20 * VectorXd::Ones(me[2]);
                A.block(ml[2] + me[2], 0, mi[2], mi[2]) = -MatrixXd::Identity(mi[2], mi[2]);
                b.segment(ml[2] + me[2],mi[2]) = -3 * VectorXd::Ones(mi[2]);
            }
        }
        else
        {
            if (p > 0)
            {
                A.block(0, 0, me[0], nVar) = MatrixXd::Random(me[0], nVar);
                b.head(me[0]) = 10 * VectorXd::Ones(me[0]);
                A.block(me[0], 0, mi[0], nVar) = -MatrixXd::Random(mi[0], nVar);
                b.segment(me[0],mi[0]) = -2 * VectorXd::Ones(mi[0]);
            }

            if (p > 1)
            {
                A.block(ml[1], 0, me[1], nVar) = MatrixXd::Random(me[1], nVar);
                b.segment(ml[1],me[1]) = 20 * VectorXd::Ones(me[1]);
                A.block(ml[1] + me[1], 0, mi[1], nVar) = MatrixXd::Random(mi[1], nVar);
                b.segment(ml[1] + me[1],mi[1]) = -3 * VectorXd::Ones(mi[1]);
            }

            if (p > 2)
            {
                A.block(ml[2], 0, me[2], nVar) = MatrixXd::Random(me[2], nVar);
                b.segment(ml[2],me[2]) = 30 * VectorXd::Ones(me[2]);
                A.block(ml[2] + me[2], 0, mi[2], nVar) = MatrixXd::Random(mi[2], nVar);
                b.segment(ml[2] + me[2],mi[2]) = -4 * VectorXd::Ones(mi[2]);
            }
        }

        // VectorXd x_prev = VectorXd::Zero(nVar);
        nipmhlsp::NIpmHLSP solver(p, nVar);
        for (Index l=0; l < p; l++)
                solver.setData(l, A.middleRows(ml[l],me[l]), b.segment(ml[l],me[l]), A.middleRows(ml[l]+me[l],mi[l]), b.segment(ml[l]+me[l],mi[l]));
        solver.solve();

        std::cout << "=============== Test " << testnr << " with " << nVar << " variables and " << p << " levels finished with KKT " << solver.KKT << " in " << solver.iter << " iterations and " << solver.time << " [s] with primal x: " << solver.get_x().transpose() << std::endl;
        if (solver.KKT > 1e-3) { cout << "ERROR: KKT norm too high, something's wrong" << endl; throw; }
    }

    return 0;
}


