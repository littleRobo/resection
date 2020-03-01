#include "photogrammetry.h"

photogrammetry::photogrammetry()
{

}

void photogrammetry::angle2Rabc(double phi, double omg, double kaf, R_abc &rabc)
{
    rabc.a1 = cos(phi)*cos(kaf) - sin(phi)*sin(omg)*sin(kaf);
    rabc.a2 = -1*cos(phi)*sin(kaf) - sin(phi)*sin(omg)*cos(kaf);
    rabc.a3 = -1*sin(phi)*cos(omg);
    rabc.b1 = cos(omg)*sin(kaf);
    rabc.b2 = cos(omg)*cos(kaf);
    rabc.b3 = -1*sin(omg);
    rabc.c1 = sin(phi)*cos(kaf) + cos(phi)*sin(omg)*sin(kaf);
    rabc.c2 = -1*sin(phi)*sin(kaf) + cos(phi)*sin(omg)*cos(kaf);
    rabc.c3 = cos(phi)*cos(omg);
}

void photogrammetry::Rabc2R(R_abc rabc, Matrix3d &R)
{
    R << rabc.a1, rabc.a2, rabc.a3,
            rabc.b1, rabc.b2, rabc.b3,
            rabc.c1, rabc.c2, rabc.c3;
}

void photogrammetry::Rabc2Rinv(R_abc rabc, Matrix3d &Rinv)
{
    Rinv << rabc.a1, rabc.b1, rabc.c1,
            rabc.a2, rabc.b2, rabc.c2,
            rabc.a3, rabc.b3, rabc.c3;
}

void photogrammetry::getdcoorC_dEO(coorC coc, coorW cow, IO io, EO eo, dcoorC_dEO &dcdeo)
{
    R_abc rabc;
    angle2Rabc(eo.phi, eo.omg, eo.kaf, rabc);
    Matrix3d Rinv;
    Rabc2Rinv(rabc, Rinv);
    Vector3d XYZ;
    XYZ << cow.X-eo.Xs, cow.Y-eo.Ys, cow.Z-eo.Zs;
    Vector3d XYZ_;
    XYZ_ = Rinv * XYZ;//(2.5.7)
    double X_, Y_, Z_;
    X_ = XYZ_(0);
    Y_ = XYZ_(1);
    Z_ = XYZ_(2);
    //(2.5.9)
    dcdeo.a11 = 1.0/Z_*(rabc.a1*io.f+rabc.a3*(coc.x-io.x0));
    dcdeo.a12 = 1.0/Z_*(rabc.b1*io.f+rabc.b3*(coc.x-io.x0));
    dcdeo.a13 = 1.0/Z_*(rabc.c1*io.f+rabc.c3*(coc.x-io.x0));

    dcdeo.a21 = 1.0/Z_*(rabc.a2*io.f+rabc.a3*(coc.y-io.y0));
    dcdeo.a22 = 1.0/Z_*(rabc.b2*io.f+rabc.b3*(coc.y-io.y0));
    dcdeo.a23 = 1.0/Z_*(rabc.c2*io.f+rabc.c3*(coc.y-io.y0));

    dcdeo.a14 = (coc.y-io.y0)*sin(eo.omg)-
            ((coc.x-io.x0)/io.f*((coc.x-io.x0)*cos(eo.kaf)-(coc.y-io.y0)*sin(eo.kaf))+io.f*cos(eo.kaf))*cos(eo.omg);
    dcdeo.a15 = -1*io.f*sin(eo.kaf)-(coc.x-io.x0)/io.f*((coc.x-io.x0)*sin(eo.kaf)+(coc.y-io.y0)*cos(eo.kaf));
    dcdeo.a16 = coc.y-io.y0;

    dcdeo.a24 = -1*(coc.x-io.x0)*sin(eo.omg)-
            ((coc.y-io.y0)/io.f*((coc.x-io.x0)*cos(eo.kaf)-(coc.y-io.y0)*sin(eo.kaf))-io.f*sin(eo.kaf))*cos(eo.omg);
    dcdeo.a25 = -1*io.f*cos(eo.kaf)-(coc.y-io.y0)/io.f*((coc.x-io.x0)*sin(eo.kaf)+(coc.y-io.y0)*cos(eo.kaf));
    dcdeo.a26 = io.x0-coc.x;
}

void photogrammetry::coorW2coorC(EO eo, IO io, coorW cow, coorC &coc)
{
    R_abc rabc;
    angle2Rabc(eo.phi, eo.omg, eo.kaf, rabc);
    coc.x = io.x0 - io.f*(rabc.a1*(cow.X-eo.Xs)+rabc.b1*(cow.Y-eo.Ys)+rabc.c1*(cow.Z-eo.Zs))
            /(rabc.a3*(cow.X-eo.Xs)+rabc.b3*(cow.Y-eo.Ys)+rabc.c3*(cow.Z-eo.Zs));
    coc.y = io.y0 - io.f*(rabc.a2*(cow.X-eo.Xs)+rabc.b2*(cow.Y-eo.Ys)+rabc.c2*(cow.Z-eo.Zs))
            /(rabc.a3*(cow.X-eo.Xs)+rabc.b3*(cow.Y-eo.Ys)+rabc.c3*(cow.Z-eo.Zs));
}


bool isGoodX(Matrix<double,6,1> X, double thresholdAngle, double thresholdDist){
    bool flag = true;
    for(int i=0;i<3;i++){
        if(abs(X(i,0))>thresholdDist) flag=false;
        if(abs(X(i+3,0))>thresholdAngle) flag=false;
    }
    return flag;
}
int photogrammetry::resection(vector<coorC> cocs, vector<coorW> cows, IO io, EO &eo)
{
    int lenc = cocs.size();
    int lenw = cows.size();
    if(lenc!=lenw) return 0;
    if(lenc < 3) return 0;

    //initial EO
    eo.Xs = 0;
    eo.Ys = 0;
    eo.Zs = 0;
    eo.phi = 0;
    eo.omg = 0;
    eo.kaf = 0;
    for(int i=0; i<lenw; i++){
        eo.Xs += cows[i].X;
        eo.Ys += cows[i].Y;
        eo.Zs += cows[i].Z;
    }
    eo.Xs /= lenw;
    eo.Ys /= lenw;
    eo.Zs /= lenw;
    eo.Zs *=2;//for right result

    //LS solver(2-5-12)
    MatrixXd A, L;
    A.resize(2*lenc,6);
    L.resize(2*lenc,1);

    Matrix<double,6,1> X;//Xs,Ys,Zs,phi,omg,kaf
    X <<1,1,1,1,1,1;

    double thresholdDist = 1e-4;//0.0001m
    double thresholdAngle = 4.8e-7;//0.1/206265rad

    bool isGoodResult =true;
    int stepnum = 20;
    int step =0;
    while(!isGoodX(X,thresholdAngle,thresholdDist) && step<stepnum){
        for(int i=0; i<lenc; i++){
            coorC cocpro;
            coorW2coorC(eo, io, cows[i], cocpro);
            L(2*i,0) = cocs[i].x - cocpro.x;
            L(2*i+1,0) = cocs[i].y - cocpro.y;

            dcoorC_dEO dcdeoi;
            getdcoorC_dEO(cocs[i], cows[i], io, eo, dcdeoi);
            A(2*i,0) = dcdeoi.a11;
            A(2*i,1) = dcdeoi.a12;
            A(2*i,2) = dcdeoi.a13;
            A(2*i,3) = dcdeoi.a14;
            A(2*i,4) = dcdeoi.a15;
            A(2*i,5) = dcdeoi.a16;

            A(2*i+1,0) = dcdeoi.a21;
            A(2*i+1,1) = dcdeoi.a22;
            A(2*i+1,2) = dcdeoi.a23;
            A(2*i+1,3) = dcdeoi.a24;
            A(2*i+1,4) = dcdeoi.a25;
            A(2*i+1,5) = dcdeoi.a26;
        }

        Matrix<double,6,6> N;
        Matrix<double,6,1> U;

        N = A.transpose()*A;
        U = A.transpose()*L;

        X = N.inverse()*U;

        //update
        eo.Xs+=X(0,0);
        eo.Ys+=X(1,0);
        eo.Zs+=X(2,0);
        eo.phi+=X(3,0);
        eo.omg+=X(4,0);
        eo.kaf+=X(5,0);

        step++;
        if(step==stepnum) isGoodResult =false;
    }

    if(!isGoodResult) return 0;
    return 1;
}

int photogrammetry::registration(vector<coorC> cocs1, vector<coorC> cocs2, regisEle &reg)
{
    int lenc1 = cocs1.size();
    int lenc2 = cocs2.size();
    if(lenc1!=lenc2) return 0;
    if(lenc1 < 3) return 0;

    MatrixXd A, L;
    A.resize(2*lenc1,6);
    L.resize(2*lenc1,1);
    Matrix<double,6,1> X;//a1,a2,a3,b1,b2,b3

    for(int i=0; i<lenc1; i++){
        L(2*i,0) = cocs2[i].x;
        L(2*i+1,0) = cocs2[i].y;

        for(int j=0; j<6; j++){
            A(2*i,j) = 0;
            A(2*i+1,j) = 0;
        }

        A(2*i,0) = cocs1[i].x;
        A(2*i,1) = cocs1[i].y;
        A(2*i,2) = 1;
        A(2*i+1,3) = cocs1[i].x;
        A(2*i+1,4) = cocs1[i].y;
        A(2*i+1,5) = 1;
    }
    X = (A.transpose()*A).inverse()*A.transpose()*L;

    reg.a1 = X(0,0);
    reg.a2 = X(1,0);
    reg.a3 = X(2,0);
    reg.b1 = X(3,0);
    reg.b2 = X(4,0);
    reg.b3 = X(5,0);

}
