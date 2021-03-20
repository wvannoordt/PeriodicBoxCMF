#ifndef INITIAL_CONDITION_HPP
#define INITIAL_CONDITION_HPP
#include "cmf.h"
#include "InputParams.h"
void InitialCondition(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, InputParams& params)
{
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> rhsLb = rhs[lb];
        cmf::BlockInfo binfo = prims.Mesh()->GetBlockInfo(lb);
        double xmin = binfo.blockBounds[0];
        double ymin = binfo.blockBounds[2];
        double zmin = binfo.blockBounds[4];
        double dx = binfo.dx[0];
        double dy = binfo.dx[1];
        double dz = binfo.dx[2];
        for (cmf::cell_t k = primsLb.kmin - primsLb.exchangeK; k < primsLb.kmax + primsLb.exchangeK; k++)
        {
            for (cmf::cell_t j = primsLb.jmin - primsLb.exchangeJ; j < primsLb.jmax + primsLb.exchangeJ; j++)
            {
                for (cmf::cell_t i = primsLb.imin - primsLb.exchangeI; i < primsLb.imax + primsLb.exchangeI; i++)
                {
                    double xi = xmin + dx*((double)i+0.5);
                    double yi = ymin + dy*((double)j+0.5);
                    double zi = zmin + dz*((double)k+0.5);
                    for (int v = 0; v < 5; v++) rhsLb(v, i, j, k) = 0.0;
                    double gm1  = params.gamma - 1.0;
                    double cp   = params.cp;
                    double Rgas = cp*gm1/params.gamma;
                    double L    = 0.159154943091895;
                    double p0   = 10.0;
                    double T0   = 1 / params.gamma;
                    double rho0 = p0 / Rgas*T0;
                    double c0   = sqrt(params.gamma*Rgas*T0);
                    double M0   = 0.1;
                    double u0   = M0*c0;
                    primsLb(0, i, j, k) =  p0 + (rho0*u0*u0/16.0)*(cos(2.0*xi/L) + cos(2.0*yi/L))*(cos(2.0*zi/L)+2.0);
                    primsLb(1, i, j, k) =  T0;
                    primsLb(2, i, j, k) =  u0*sin(xi/L)*cos(yi/L)*cos(zi/L);
                    primsLb(3, i, j, k) = -u0*cos(xi/L)*sin(yi/L)*cos(zi/L);
                    primsLb(4, i, j, k) =  0.0;
                }
            }
        }
    }
}

#endif