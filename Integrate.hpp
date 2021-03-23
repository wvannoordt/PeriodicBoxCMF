#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP
#include "cmf.h"

double Enstrophy(cmf::CartesianMeshArray& prims)
{
    double integratedEnstrophyLocal = 0.0;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockInfo info = prims.Mesh()->GetBlockInfo(lb);
        double dx = info.dx[0];
        double dy = info.dx[1];
        double dz = info.dx[2];
        for (cmf::cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
        {
            for (cmf::cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
            {
                for (cmf::cell_t i = primsLb.imin; i < primsLb.imax; i++)
                {
                    double du_dx = (primsLb(2, i+1, j, k)-primsLb(2, i-1, j, k))/dx;
                    double dv_dx = (primsLb(3, i+1, j, k)-primsLb(3, i-1, j, k))/dx;
                    double dw_dx = (primsLb(4, i+1, j, k)-primsLb(4, i-1, j, k))/dx;
                    
                    double du_dy = (primsLb(2, i, j+1, k)-primsLb(2, i, j-1, k))/dy;
                    double dv_dy = (primsLb(3, i, j+1, k)-primsLb(3, i, j-1, k))/dy;
                    double dw_dy = (primsLb(4, i, j+1, k)-primsLb(4, i, j-1, k))/dy;
                    
                    double du_dz = (primsLb(2, i, j, k+1)-primsLb(2, i, j, k-1))/dz;
                    double dv_dz = (primsLb(3, i, j, k+1)-primsLb(3, i, j, k-1))/dz;
                    double dw_dz = (primsLb(4, i, j, k+1)-primsLb(4, i, j, k-1))/dz;
                    
                    double omega1 = dw_dy - dv_dz;
                    double omega2 = du_dz - dw_dx;
                    double omega3 = dv_dx - du_dy;
                    
                    integratedEnstrophyLocal += (omega1*omega1 + omega2*omega2 + omega3*omega3)*dx*dy*dz;
                }
            }
        }
    }
    double integratedEnstrophyGlobal = cmf::globalGroup.Sum(integratedEnstrophyLocal);
    return integratedEnstrophyGlobal;
}

#endif