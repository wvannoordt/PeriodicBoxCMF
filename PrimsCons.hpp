#ifndef PRIMS_CONS_HPP
#define PRIMS_CONS_HPP
#include "cmf.h"
#include "InputParams.h"
void ConsToPrims(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, InputParams& params)
{
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> consLb  = cons[lb];
        for (cmf::cell_t k = primsLb.kmin - primsLb.exchangeK; k < primsLb.kmax + primsLb.exchangeK; k++)
        {
            for (cmf::cell_t j = primsLb.jmin - primsLb.exchangeJ; j < primsLb.jmax + primsLb.exchangeJ; j++)
            {
                for (cmf::cell_t i = primsLb.imin - primsLb.exchangeI; i < primsLb.imax + primsLb.exchangeI; i++)
                {
                    double rho = consLb(0, i, j ,k);
                    double invrho = 1.0/rho;
                    double u = invrho*consLb(2, i, j, k);
                    double v = invrho*consLb(3, i, j, k);
                    double w = invrho*consLb(4, i, j, k);
                    double rhoU2 = rho*(u*u+v*v+w*w);
                    double p = (params.gamma - 1.0)*(consLb(1, i, j ,k) - 0.5*rhoU2);
                    double Rgas = params.cp*(params.gamma - 1.0)/params.gamma;
                    double T = p/(Rgas*rho);
                    primsLb(0, i, j, k) = p;
                    primsLb(1, i, j, k) = T;
                    primsLb(2, i, j, k) = u;
                    primsLb(3, i, j, k) = v;
                    primsLb(4, i, j, k) = w;
                }
            }
        }
    }
}

void PrimsToCons(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, InputParams& params)
{
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> consLb  = cons[lb];
        for (cmf::cell_t k = primsLb.kmin - primsLb.exchangeK; k < primsLb.kmax + primsLb.exchangeK; k++)
        {
            for (cmf::cell_t j = primsLb.jmin - primsLb.exchangeJ; j < primsLb.jmax + primsLb.exchangeJ; j++)
            {
                for (cmf::cell_t i = primsLb.imin - primsLb.exchangeI; i < primsLb.imax + primsLb.exchangeI; i++)
                {
                    double p = primsLb(0, i, j, k);
                    double T = primsLb(1, i, j, k);
                    double u = primsLb(2, i, j, k);
                    double v = primsLb(3, i, j, k);
                    double w = primsLb(4, i, j, k);
                    double Rgas = params.cp*(params.gamma - 1.0)/params.gamma;
                    double rho = p / (Rgas*T);
                    double rhoU2 = rho*(u*u+v*v+w*w);
                    double rhoE = rhoU2 + (p/((params.gamma - 1.0)));
                    double rhoU = rho*u;
                    double rhoV = rho*v;
                    double rhoW = rho*w;
                    consLb(0, i, j, k) = rho;
                    consLb(1, i, j, k) = rhoE;
                    consLb(2, i, j, k) = rhoU;
                    consLb(3, i, j, k) = rhoV;
                    consLb(4, i, j, k) = rhoW;
                }
            }
        }
    }
}

#endif