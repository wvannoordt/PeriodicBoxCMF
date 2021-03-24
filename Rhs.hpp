#ifndef RHS_HPP
#define RHS_HPP
#include "cmf.h"
#include "InputParams.h"
#include "util.h"
#define stencilIdx(v,j) ((v)+(5+3)*(j))
#define f_DivSplit(q,j,l,v1)         (0.500*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))]))
#define fg_QuadSplit(q,j,l,v1,v2)    (0.250*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))]))
#define fg_CubeSplit(q,j,l,v1,v2,v3) (0.125*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))])*(q[stencilIdx((v3),(j))] + q[stencilIdx((v3),(j)+(l))]))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*((q[stencilIdx((v1),(j)+(l))]*q[stencilIdx((v2),(j))]) + (q[stencilIdx((v1),(j))]*q[stencilIdx((v2),(j)+(l))])))
void ComputeRhs(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, cmf::CartesianMeshArray& rhs, InputParams& params)
{
    double centerCoef[4] = {0.0};
    switch (params.centOrder)
    {
        case 2: {centerCoef[0] = 1.0/2.0; break;}
        case 4: {centerCoef[0] = 2.0/3.0; centerCoef[1] = -1.0/12.0; break;}
        case 6: {centerCoef[0] = 3.0/4.0; centerCoef[1] = -3.0/20.0; centerCoef[2] = 1.0/60.0; break;}
        case 8: {centerCoef[0] = 4.0/5.0; centerCoef[1] = -1.0/5.0 ; centerCoef[2] = 4.0/105 ; centerCoef[3] = -1.0/280.0; break;}
        default: {std::cout << "Bad central scheme order." << std::endl; abort();}
    }
    double Rgas = params.cp*(params.gamma - 1.0)/params.gamma;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> consLb  = cons[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        cmf::BlockInfo info = rhs.Mesh()->GetBlockInfo(lb);
        int stencilWid = params.centOrder/2;
        for (int idir = 0; idir < 3; idir++)
        {
            int dijk[3] = {0};
            dijk[idir] = 1;
            for (cmf::cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
            {
                for (cmf::cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
                {
                    for (cmf::cell_t i = primsLb.imin; i < primsLb.imax; i++)
                    {
                        
                        double stencilData[9*(8)]; //ie,ke,T,P,rho,u,v,w
                        for (int n = 0; n < params.centOrder + 1; n++)
                        {
                            for (int v = 3; v < (5+3); v++)
                            {
                                int ii = i+dijk[0]*(n-stencilWid);
                                int jj = j+dijk[1]*(n-stencilWid);
                                int kk = k+dijk[2]*(n-stencilWid);
                                stencilData[stencilIdx(v,n)] = primsLb(v-3, ii, jj, kk);
                            }
                            // stencilData = ? ? ? P T u v w
                            // T
                            stencilData[stencilIdx(2,n)] = stencilData[stencilIdx(4,n)];
                            
                            // stencilData = ? ? T P T u v w
                            //rho
                            stencilData[stencilIdx(4,n)] = stencilData[stencilIdx(3,n)]/(Rgas*stencilData[stencilIdx(2,n)]);//p = rho r t -> rho = p/RT
                            
                            // IE = P/(rho*(gamma - 1))
                            stencilData[stencilIdx(0,n)] = stencilData[stencilIdx(3,n)]/(stencilData[stencilIdx(4,n)]*(params.gamma - 1.0));
                            
                            // ke (don't care)
                            stencilData[stencilIdx(1,n)] = 0.0;

                            // Not needed per se starts
                            for (int vel_comp = 0; vel_comp < 3; vel_comp ++)
                            {
                                stencilData[stencilIdx(1,n)] += 0.5*stencilData[stencilIdx(5+vel_comp,n)]*stencilData[stencilIdx(5+vel_comp,n)];
                            }
                            // Not needed per se ends
                        }
                        
                        double C[2]     = {0.0};
                        double M[6]     = {0.0};
                        double PGRAD[2] = {0.0};
                        double KE[2]    = {0.0};
                        double IE[2]    = {0.0};
                        double PDIFF[2] = {0.0};
                        
                        for (int l = 1; l <= stencilWid; l++)
                        {
                            double al = centerCoef[l-1];
                            int jf = stencilWid;
                            for (int m = 0; m <= (l-1); m++)
                            {
                                C[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,4,5+idir);
                                C[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir);
                                for (int idir_mom = 0; idir_mom < 3; idir_mom++)
                                {
                                    M[idir_mom    ] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,4,5+idir,5+idir_mom);
                                    M[idir_mom + 3] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,4,5+idir,5+idir_mom);
                                }

                                PGRAD[1] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf-m, l,3);
                                PGRAD[0] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf+m,-l,3);

                                for (int vel_comp = 0;  vel_comp < 3; vel_comp ++)
                                {
                                    KE[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,4,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf-m)]*stencilData[stencilIdx(5+vel_comp,jf-m+l)]);
                                    KE[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf+m)]*stencilData[stencilIdx(5+vel_comp,jf+m-l)]);
                                }

                                IE[1] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,4,0,5+idir);
                                IE[0] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,4,0,5+idir);

                                PDIFF[1] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                                PDIFF[0] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf+m,-l,5+idir,3);
                            }
                        }
                        
                        rhsLb(0, i, j, k)      -= info.dxInv[idir]*(C[1] - C[0]);
                        rhsLb(1, i, j, k)      -= info.dxInv[idir]*(IE[1] + KE[1] + PDIFF[1] - IE[0] - KE[0] - PDIFF[0]);
                        rhsLb(2, i, j, k)      -= info.dxInv[idir]*(M[0] - M[3]);
                        rhsLb(3, i, j, k)      -= info.dxInv[idir]*(M[1] - M[4]);
                        rhsLb(4, i, j, k)      -= info.dxInv[idir]*(M[2] - M[5]);
                        rhsLb(2+idir, i, j, k) -= info.dxInv[idir]*(PGRAD[1] - PGRAD[0]);
                    }
                }
            }
            dijk[idir] = 0;
        }
    }
}

void Advance(cmf::CartesianMeshArray& cons, cmf::CartesianMeshArray& rhs, cmf::CartesianMeshArray& consRKTemp, cmf::CartesianMeshArray& rhsRKTemp, InputParams& params)
{
    for (auto lb: cons)
    {
        cmf::BlockArray<double, 1> consLb  = cons[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        for (cmf::cell_t k = consLb.kmin; k < consLb.kmax; k++)
        {
            for (cmf::cell_t j = consLb.jmin; j < consLb.jmax; j++)
            {
                for (cmf::cell_t i = consLb.imin; i < consLb.imax; i++)
                {
                    consLb(0, i, j, k) += params.deltaT * rhsLb(0, i, j, k);
                    consLb(1, i, j, k) += params.deltaT * rhsLb(1, i, j, k);
                    consLb(2, i, j, k) += params.deltaT * rhsLb(2, i, j, k);
                    consLb(3, i, j, k) += params.deltaT * rhsLb(3, i, j, k);
                    consLb(4, i, j, k) += params.deltaT * rhsLb(4, i, j, k);
                }
            }
        }
    }
}

void AEqualsBXPlusC(cmf::CartesianMeshArray& a, double b, cmf::CartesianMeshArray& x, cmf::CartesianMeshArray& c)
{
    for (auto lb: a)
    {
        cmf::BlockArray<double, 1> aLb = a[lb];
        cmf::BlockArray<double, 1> xLb = x[lb];
        cmf::BlockArray<double, 1> cLb = c[lb];
        for (cmf::cell_t k = cLb.kmin; k < cLb.kmax; k++)
        {
            for (cmf::cell_t j = cLb.jmin; j < cLb.jmax; j++)
            {
                for (cmf::cell_t i = cLb.imin; i < cLb.imax; i++)
                {
                    aLb(0, i, j, k) = b*xLb(0, i, j, k) + cLb(0, i, j, k);
                    aLb(1, i, j, k) = b*xLb(1, i, j, k) + cLb(1, i, j, k);
                    aLb(2, i, j, k) = b*xLb(2, i, j, k) + cLb(2, i, j, k);
                    aLb(3, i, j, k) = b*xLb(3, i, j, k) + cLb(3, i, j, k);
                    aLb(4, i, j, k) = b*xLb(4, i, j, k) + cLb(4, i, j, k);
                }
            }
        }
    }
}

void PlusEqualsKX(cmf::CartesianMeshArray& a, double kd, cmf::CartesianMeshArray& x)
{
    for (auto lb: a)
    {
        cmf::BlockArray<double, 1> aLb = a[lb];
        cmf::BlockArray<double, 1> xLb = x[lb];
        for (cmf::cell_t k = aLb.kmin; k < aLb.kmax; k++)
        {
            for (cmf::cell_t j = aLb.jmin; j < aLb.jmax; j++)
            {
                for (cmf::cell_t i = aLb.imin; i < aLb.imax; i++)
                {
                    aLb(0, i, j, k) += kd*xLb(0, i, j, k);
                    aLb(1, i, j, k) += kd*xLb(1, i, j, k);
                    aLb(2, i, j, k) += kd*xLb(2, i, j, k);
                    aLb(3, i, j, k) += kd*xLb(3, i, j, k);
                    aLb(4, i, j, k) += kd*xLb(4, i, j, k);
                }
            }
        }
    }    
}

void ZeroRhs(cmf::CartesianMeshArray& rhs)
{
    for (auto lb: rhs)
    {
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        for (cmf::cell_t k = rhsLb.kmin; k < rhsLb.kmax; k++)
        {
            for (cmf::cell_t j = rhsLb.jmin; j < rhsLb.jmax; j++)
            {
                for (cmf::cell_t i = rhsLb.imin; i < rhsLb.imax; i++)
                {
                    rhsLb(0, i, j, k) = 0.0;
                    rhsLb(1, i, j, k) = 0.0;
                    rhsLb(2, i, j, k) = 0.0;
                    rhsLb(3, i, j, k) = 0.0;
                    rhsLb(4, i, j, k) = 0.0;
                }
            }
        }
    }
}
double UMax(cmf::CartesianMeshArray& prims)
{
    double umax = -1;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb  = prims[lb];
        for (cmf::cell_t k = primsLb.kmin - primsLb.exchangeK; k < primsLb.kmax + primsLb.exchangeK; k++)
        {
            for (cmf::cell_t j = primsLb.jmin - primsLb.exchangeJ; j < primsLb.jmax + primsLb.exchangeJ; j++)
            {
                for (cmf::cell_t i = primsLb.imin - primsLb.exchangeI; i < primsLb.imax + primsLb.exchangeI; i++)
                {
                    double u = primsLb(2, i, j, k);
                    double v = primsLb(3, i, j, k);
                    double w = primsLb(4, i, j, k);
                    double u2 = u*u+v*v+w*w;
                    umax = (umax<u2)?u2:umax;
                }
            }
        }
    }
    return sqrt(umax);
}

#endif