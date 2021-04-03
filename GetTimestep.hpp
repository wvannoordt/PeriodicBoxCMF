#ifndef GET_TIMESTEP_HPP
#define GET_TIMESTEP_HPP
#include "cmf.h"
#include "Rhs.hpp"
#define m_min(a, b) (((a)<(b))?(a):(b))
#define m_max(a, b) (((a)>(b))?(a):(b))
using cmf::print;
double GetTimestep(cmf::CartesianMeshArray& prims, InputParams& params)
{
    double dxmin = 1e20;
    double umax = UMax(prims);
    double sosMax = -1;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockInfo binfo = prims.Mesh()->GetBlockInfo(lb);
        double xmin = binfo.blockBounds[0];
        double ymin = binfo.blockBounds[2];
        double zmin = binfo.blockBounds[4];
        double dx = binfo.dx[0];
        double dy = binfo.dx[1];
        double dz = binfo.dx[2];
        dxmin = m_min(dxmin, dx);
        dxmin = m_min(dxmin, dy);
        dxmin = m_min(dxmin, dz);
        for (cmf::cell_t k = primsLb.kmin - primsLb.exchangeK; k < primsLb.kmax + primsLb.exchangeK; k++)
        {
            for (cmf::cell_t j = primsLb.jmin - primsLb.exchangeJ; j < primsLb.jmax + primsLb.exchangeJ; j++)
            {
                for (cmf::cell_t i = primsLb.imin - primsLb.exchangeI; i < primsLb.imax + primsLb.exchangeI; i++)
                {
                    double gm1  = params.gamma - 1.0;
                    double cp   = params.cp;
                    double gamma = params.gamma;
                    double Rgas = cp*gm1/gamma;
                    double temperature = primsLb(1, i, j, k);
                    double sos = sqrt(gamma*Rgas*temperature);
                    sosMax = m_max(sosMax, sos);
                }
            }
        }
    }
    double sosMaxGlob = cmf::globalGroup.Max(sosMax);
    return params.CFL*dxmin/(umax+sosMaxGlob);
}
#endif