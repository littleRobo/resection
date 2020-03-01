#ifndef PHOTOGRAMMETRY_H
#define PHOTOGRAMMETRY_H

#include <iostream>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Geometry"

using namespace std;
using namespace Eigen;

///===  Photogrammetry.2nd by Zhang. Pan. Wang ====
///

struct R_abc{
    double a1;
    double a2;
    double a3;
    double b1;
    double b2;
    double b3;
    double c1;
    double c2;
    double c3;
};

//--- R = [a1, a2, a3]---[x] = [a1, b1, c1][X-Xs]
//--------[b1, b2, b3]---[y]   [a2, b2, c2][Y-Ys]
//--------[c1, c2, c3]---[-f]  [a3, b3, c3][Z-Zs]

struct EO{
    double Xs;
    double Ys;
    double Zs;
    double phi;
    double omg;
    double kaf;
};

struct IO{
    double f;
    double x0;
    double y0;
};

struct coorW{
    double X;
    double Y;
    double Z;
};//coordinate in world frame

struct coorC{
    double x;
    double y;
};

struct dcoorC_dEO{
    double a11;
    double a12;
    double a13;
    double a21;
    double a22;
    double a23;
    double a14;
    double a15;
    double a16;
    double a24;
    double a25;
    double a26;
};

struct regisEle{
    double a1;
    double a2;
    double a3;
    double b1;
    double b2;
    double b3;
};
//---- x2 = a1*x1 + a2*y1 + a3
//---- y2 = b1*x1 + b2*y1 + b3

class photogrammetry
{
public:
    photogrammetry();
    void angle2Rabc(double phi, double omg, double kaf, R_abc &rabc);//phi,omiga,kafa(2-3-3)
    void Rabc2R(R_abc rabc, Matrix3d &R);
    void Rabc2Rinv(R_abc rabc, Matrix3d &Rinv);
    void getdcoorC_dEO(coorC coc, coorW cow, IO io, EO eo, dcoorC_dEO &dcdeo);
    void coorW2coorC(EO eo, IO io, coorW cow, coorC &coc);//(2-3-5)

    int resection(vector<coorC> cocs, vector<coorW> cows, IO io, EO &eo);

    int registration(vector<coorC> cocs1, vector<coorC> cocs2, regisEle &reg);

};

#endif // PHOTOGRAMMETRY_H
