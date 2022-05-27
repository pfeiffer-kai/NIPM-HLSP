#pragma once

#include <Eigen/Dense>
#include <Eigen/Jacobi>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <memory>

using namespace Eigen;

typedef Eigen::Block<MatrixXd, Eigen::Dynamic, Eigen::Dynamic> dBlockType;
typedef Eigen::VectorBlock<VectorXd, Eigen::Dynamic>            dVectorBlockType;
typedef Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1> >::Index Index;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mati;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;

typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi veci;

using namespace std;

enum status {FAILURE=-2, EXCEED, OK, SUCCESS};
enum verboseLevel {NONE = 0, CONV, SOLVE, MAT};
