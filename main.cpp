#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>

#include "variables.h"
#include "objects.h"
#include "initialize.h"
#include "readmesh.h"
//

float ul_, ur_, cl_, cr_, uf_, mf_, rf_, ruf_, rvf_, ref_, pf_; // upwind

float vcell, vnorm;

float res12, res120;

void flux(int i){

    for (int j = 0; j<3; j++){

        if (cells[i].neighID[j] - 1 > i || cells[i].neighID[j] < 0){

            ql_.rho          =   q[i].rho;
            ql_.rhou         =   q[i].rhou * cells[i].surfaceC[j] + q[i].rhov * cells[i].surfaceS[j];
            ql_.rhov         = - q[i].rhou * cells[i].surfaceS[j] + q[i].rhov * cells[i].surfaceC[j];
            ql_.e            =   q[i].e;
            ql_.p            =   q[i].p;
            ql_.T            =   q[i].T;

            if (cells[i].neighID[j] > 0){

                qr_.rho          =   q[cells[i].neighID[j] - 1].rho;
                qr_.rhou         =   q[cells[i].neighID[j] - 1].rhou * cells[i].surfaceC[j] + q[cells[i].neighID[j] - 1].rhov * cells[i].surfaceS[j];
                qr_.rhov         = - q[cells[i].neighID[j] - 1].rhou * cells[i].surfaceS[j] + q[cells[i].neighID[j] - 1].rhov * cells[i].surfaceC[j];
                qr_.e            =   q[cells[i].neighID[j] - 1].e;
                qr_.p            =   q[cells[i].neighID[j] - 1].p;
                qr_.T            =   q[cells[i].neighID[j] - 1].T;

            }else if(cells[i].neighID[j] == -1){// BC FF

                qr_.rho          =   rho_;
                qr_.rhou         =   rho_ * u_ * cells[i].surfaceC[j] + rho_ * v_ * cells[i].surfaceS[j];
                qr_.rhov         = - rho_ * u_ * cells[i].surfaceS[j] + rho_ * v_ * cells[i].surfaceC[j];
                qr_.e            =   e_;
                qr_.p            =   p_;
                qr_.T            =   T_;

            }else{// BC EulerWall

                qr_.rho          =   ql_.rho;
                qr_.rhou         = - ql_.rhou;
                qr_.rhov         =   ql_.rhov;
                qr_.e            =   ql_.e;
                qr_.p            =   ql_.p;
                qr_.T            =   ql_.T;

            }

            // UPWIND

            ul_         = ql_.rhou / ql_.rho;
            ur_         = qr_.rhou / qr_.rho;

            cl_         = sqrt(gamma_ * ql_.p / ql_.rho);
            cr_         = sqrt(gamma_ * qr_.p / qr_.rho);

            uf_         = 0.5 * (ul_ + ur_);
            mf_         = 0.5 * (ul_ / cl_ + ur_ / cr_);

            if (mf_ > 0){

                rf_     = ql_.rho;
                ruf_    = ql_.rhou;
                rvf_    = ql_.rhov;
                ref_    = ql_.e;
                pf_     = ql_.p;

            }else{

                rf_     = qr_.rho;
                ruf_    = qr_.rhou;
                rvf_    = qr_.rhov;
                ref_    = qr_.e;
                pf_     = qr_.p;

            }

            if (sqrt(pow(mf_, 2.)) < 1.){ pf_ = 0.5 * (ql_.p + qr_.p);}

            flux_.f1    = uf_ * rf_;
            flux_.f2    = uf_ * ruf_ + pf_;
            flux_.f3    = uf_ * rvf_;
            flux_.f4    = uf_ * (ref_ + pf_);

            flux_.f1        = flux_.f1 * cells[i].surfaceN[j];
            flux_.f2        = flux_.f2 * cells[i].surfaceN[j];
            flux_.f3        = flux_.f3 * cells[i].surfaceN[j];
            flux_.f4        = flux_.f4 * cells[i].surfaceN[j];

            flux__.f2       = flux_.f2;

            flux_.f2        = flux__.f2 * cells[i].surfaceC[j] - flux_.f3 * cells[i].surfaceS[j];
            flux_.f3        = flux__.f2 * cells[i].surfaceS[j] + flux_.f3 * cells[i].surfaceC[j];

            tflux[i].f1     = tflux[i].f1 + flux_.f1;
            tflux[i].f2     = tflux[i].f2 + flux_.f2;
            tflux[i].f3     = tflux[i].f3 + flux_.f3;
            tflux[i].f4     = tflux[i].f4 + flux_.f4;

            if (cells[i].neighID[j] > 0){

                tflux[cells[i].neighID[j] - 1].f1     = tflux[cells[i].neighID[j] - 1].f1 - flux_.f1;
                tflux[cells[i].neighID[j] - 1].f2     = tflux[cells[i].neighID[j] - 1].f2 - flux_.f2;
                tflux[cells[i].neighID[j] - 1].f3     = tflux[cells[i].neighID[j] - 1].f3 - flux_.f3;
                tflux[cells[i].neighID[j] - 1].f4     = tflux[cells[i].neighID[j] - 1].f4 - flux_.f4;

            }

        }
    }
}

void localdt(){

    for (int i=0; i<cellNumber; i++){

        vcell       = sqrt(pow((q[i].rhou / q[i].rho), 2.) + pow((q[i].rhov / q[i].rho), 2.));
        vnorm       = sqrt(gamma_ * q[i].p / q[i].rho) + vcell;
        dtmin[i]    = cfl * cells[i].dxmin / vnorm;

    }
}

void solve(){

    for (int i=0; i<cellNumber;i++){

        q_[i].rho   = q[i].rho;
        q_[i].rhou  = q[i].rhou;
        q_[i].rhov  = q[i].rhov;
        q_[i].e     = q[i].e;

    }
    for (int rk3=0; rk3<3; rk3++){

        for (int i=0; i<cellNumber;i++){

            tflux[i].f1     = 0.;
            tflux[i].f2     = 0.;
            tflux[i].f3     = 0.;
            tflux[i].f4     = 0.;

        }

        for (int i=0; i<cellNumber; i++){

            flux(i);

            dq[i].rho   = - dtmin[i] * tflux[i].f1 / cells[i].area;
            dq[i].rhou  = - dtmin[i] * tflux[i].f2 / cells[i].area;
            dq[i].rhov  = - dtmin[i] * tflux[i].f3 / cells[i].area;
            dq[i].e     = - dtmin[i] * tflux[i].f4 / cells[i].area;

            q[i].rho    = q_[i].rho  + rk3cf[rk3] * dq[i].rho;
            q[i].rhou   = q_[i].rhou + rk3cf[rk3] * dq[i].rhou;
            q[i].rhov   = q_[i].rhov + rk3cf[rk3] * dq[i].rhov;
            q[i].e      = q_[i].e    + rk3cf[rk3] * dq[i].e;

            q[i].T      = gamma_ * gamm1_ * (q[i].e / q[i].rho - 0.5 * (pow(q[i].rhou, 2.) + pow(q[i].rhov, 2.)) / q[i].rho / q[i].rho);
            q[i].p      = q[i].rho * q[i].T / gamma_;

        }
    }

    res12       = 0.;

    for (int i=0; i<cellNumber; i++){

        res12   = res12 + dq[i].rho * dq[i].rho;

    }

    res12       = sqrt(res12 / cellNumber);

    if(timestep == 1){res120    = res12;};

    res12       = log10(res12 / res120);

    std::cout << "LOG(res_L2): " << res12 << "Timestep: " << timestep << std::endl;

}

void postProcess(){

    float sumr1[cellNumber], r1, rrho, ru, rv, rp, rmach, rt;   // post-process

    qnode   = qnode_;

    for (int i=0; i<cellNumber; i++){

        for (int j=0; j<3; j++){

            r1              = 1. / sqrt(pow((nodes[cells[i].nodeIDs[j] - 1].x - cells[i].centerX), 2.) + pow((nodes[cells[i].nodeIDs[j] - 1].y - cells[i].centerY), 2.));

            qnode[cells[i].nodeIDs[j] - 1].rho      = qnode[cells[i].nodeIDs[j] - 1].rho  + q[i].rho  * r1;
            qnode[cells[i].nodeIDs[j] - 1].rhou     = qnode[cells[i].nodeIDs[j] - 1].rhou + q[i].rhou * r1;
            qnode[cells[i].nodeIDs[j] - 1].rhov     = qnode[cells[i].nodeIDs[j] - 1].rhov + q[i].rhov * r1;
            qnode[cells[i].nodeIDs[j] - 1].e        = qnode[cells[i].nodeIDs[j] - 1].e    + q[i].e    * r1;
            qnode[cells[i].nodeIDs[j] - 1].p        = qnode[cells[i].nodeIDs[j] - 1].rho  + q[i].p    * r1;
            qnode[cells[i].nodeIDs[j] - 1].T        = qnode[cells[i].nodeIDs[j] - 1].rho  + q[i].T    * r1;

            sumr1[cells[i].nodeIDs[j] - 1]          = sumr1[cells[i].nodeIDs[j] - 1] + r1;
        }

    }

    for (int i=0; i<nodeNumber; i++){

        qnode[i].rho        = qnode[i].rho  / sumr1[i];
        qnode[i].rhou       = qnode[i].rhou / sumr1[i];
        qnode[i].rhov       = qnode[i].rhov / sumr1[i];
        qnode[i].e          = qnode[i].e    / sumr1[i];
        qnode[i].p          = qnode[i].p    / sumr1[i];
        qnode[i].T          = qnode[i].T    / sumr1[i];

    }

    std::ofstream   solutionFile;

    solutionFile.open("solution.dat");

    solutionFile << "VARIABLES= \"x\",\"y\",\"Mach\",\"Temp\",\"Pressure\",\"Rho\",\"u\", \"v\" ZONE N="<< nodeNumber <<" E= "<< cellNumber <<" F=FEPOINT  ET=triangle\n";

    for (int i=0; i<nodeNumber; i++){

        rrho    = qnode[i].rho;
        ru      = qnode[i].rhou / rrho;
        rv      = qnode[i].rhov / rrho;
        rp      = gamm1_ * (qnode[i].e - 0.5 * rrho * (ru*ru + rv*rv));
        rmach   = sqrt(ru*ru + rv*rv) / (gamma_ * rp * rrho);
        rt      = qnode[i].T;

        solutionFile << nodes[i].x << "\t" << nodes[i].y << "\t" << rmach << "\t" << rt << "\t" << rp << "\t" << rrho << "\t" << ru << "\t" << rv << "\n";

    }

    for (int i=0; i<cellNumber; i++){

        solutionFile << cells[i].nodeIDs[0] << "\t" << cells[i].nodeIDs[1] << "\t" << cells[i].nodeIDs[2] << "\n";

    }

    solutionFile.close();

}

int main(){

    readMesh();

    initialise();

    while(timestep <= maxstep){

        if(timestep%100 == 1){localdt();}

        solve();

        timestep++;

    }

    postProcess();

    delete cells, nodes, q, q_, dq, qnode, qnode_, tflux, tflux_; 

    return 0;
}

