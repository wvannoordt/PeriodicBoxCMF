#ifndef RHS_HPP
#define RHS_HPP
#include "cmf.h"
#include "InputParams.h"
#include "util.h"

using cmf::face_t;
using cmf::cell_t;

#define stencilIdx(v,j) ((v)+(5+3)*(j))
#define f_DivSplit(q,j,l,v1)         (0.500*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))]))
#define fg_QuadSplit(q,j,l,v1,v2)    (0.250*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))]))
#define fg_CubeSplit(q,j,l,v1,v2,v3) (0.125*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))])*(q[stencilIdx((v3),(j))] + q[stencilIdx((v3),(j)+(l))]))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*((q[stencilIdx((v1),(j)+(l))]*q[stencilIdx((v2),(j))]) + (q[stencilIdx((v1),(j))]*q[stencilIdx((v2),(j)+(l))])))
void ComputeConv(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, cmf::CartesianMeshArray& rhs, InputParams& params)
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

void ComputeVisc(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, cmf::CartesianMeshArray& rhs, InputParams& params)
{
    double beta = 0.0;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> consLb  = cons[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        
        cmf::BlockInfo info = rhs.Mesh()->GetBlockInfo(lb);
        
        double fluxLRAr[10];
        cmf::MdArray<double, 1> fluxLeft (&fluxLRAr[0], 5);
        cmf::MdArray<double, 1> fluxRight(&fluxLRAr[5], 5);
        
        double velGradArr[18];
        cmf::MdArray<double, 2> velGradLeft (&velGradArr[0], 3, 3);
        cmf::MdArray<double, 2> velGradRight(&velGradArr[9], 3, 3);
        
        for (int idir = 0; idir < 3; idir++)
        {
            int dir1 = (idir+1) % 3;
            int dir2 = (idir+2) % 3;
            cell_t ijkCell[3] = {0};
            cell_t ijkCellR[3] = {0};
            cell_t ijkCellL[3] = {0};
            int dijk1[3] = {0};
            int dijk2[3] = {0};
            dijk1[dir1] = 1;
            dijk2[dir2] = 1;
            for (cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
            {
                for (cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
                {
                    for (cell_t i = primsLb.imin; i < primsLb.imax; i++)
                    {
                        ijkCell[0]  = i;
                        ijkCellR[0] = i;
                        ijkCellL[0] = i;
                        ijkCell[1]  = j;
                        ijkCellR[1] = j;
                        ijkCellL[1] = j;
                        ijkCell[2]  = k;
                        ijkCellR[2] = k;
                        ijkCellL[2] = k;
                        
                        ijkCellR[idir]++;
                        ijkCellL[idir]--;
                        
                        double u_UR, u_UL, u_LR, u_LL;
                        for (int uvw = 0; uvw < 3; uvw++)
                        {
                            velGradLeft (uvw, idir) = info.dxInv[idir]*(primsLb(2+uvw, ijkCell[0],  ijkCell[1],  ijkCell[2]) -primsLb(2+uvw, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                            velGradRight(uvw, idir) = info.dxInv[idir]*(primsLb(2+uvw, ijkCellR[0], ijkCellR[1], ijkCellR[2])-primsLb(2+uvw, ijkCell[0],  ijkCell[1],  ijkCell[2]));
                            
                            //left, dir1
                            u_UR = primsLb(2+uvw, ijkCell [0] + dijk1[0],  ijkCell [1] + dijk1[1],  ijkCell [2] + dijk1[2]);
                            u_UL = primsLb(2+uvw, ijkCellL[0] + dijk1[0],  ijkCellL[1] + dijk1[1],  ijkCellL[2] + dijk1[2]);
                            u_LR = primsLb(2+uvw, ijkCell [0] - dijk1[0],  ijkCell [1] - dijk1[1],  ijkCell [2] - dijk1[2]);
                            u_LL = primsLb(2+uvw, ijkCellL[0] - dijk1[0],  ijkCellL[1] - dijk1[1],  ijkCellL[2] - dijk1[2]);
                            velGradLeft (uvw, dir1) = 0.25*info.dxInv[dir1]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //left, dir2
                            u_UR = primsLb(2+uvw, ijkCell [0] + dijk2[0],  ijkCell [1] + dijk2[1],  ijkCell [2] + dijk2[2]);
                            u_UL = primsLb(2+uvw, ijkCellL[0] + dijk2[0],  ijkCellL[1] + dijk2[1],  ijkCellL[2] + dijk2[2]);
                            u_LR = primsLb(2+uvw, ijkCell [0] - dijk2[0],  ijkCell [1] - dijk2[1],  ijkCell [2] - dijk2[2]);
                            u_LL = primsLb(2+uvw, ijkCellL[0] - dijk2[0],  ijkCellL[1] - dijk2[1],  ijkCellL[2] - dijk2[2]);
                            velGradLeft (uvw, dir2) = 0.25*info.dxInv[dir2]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //right, dir1
                            u_UR = primsLb(2+uvw, ijkCellR[0] + dijk1[0],  ijkCellR[1] + dijk1[1],  ijkCellR[2] + dijk1[2]);
                            u_UL = primsLb(2+uvw, ijkCell [0] + dijk1[0],  ijkCell [1] + dijk1[1],  ijkCell [2] + dijk1[2]);
                            u_LR = primsLb(2+uvw, ijkCellR[0] - dijk1[0],  ijkCellR[1] - dijk1[1],  ijkCellR[2] - dijk1[2]);
                            u_LL = primsLb(2+uvw, ijkCell [0] - dijk1[0],  ijkCell [1] - dijk1[1],  ijkCell [2] - dijk1[2]);
                            velGradRight(uvw, dir1) = 0.25*info.dxInv[dir1]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //right, dir2
                            u_UR = primsLb(2+uvw, ijkCellR[0] + dijk2[0],  ijkCellR[1] + dijk2[1],  ijkCellR[2] + dijk2[2]);
                            u_UL = primsLb(2+uvw, ijkCell [0] + dijk2[0],  ijkCell [1] + dijk2[1],  ijkCell [2] + dijk2[2]);
                            u_LR = primsLb(2+uvw, ijkCellR[0] - dijk2[0],  ijkCellR[1] - dijk2[1],  ijkCellR[2] - dijk2[2]);
                            u_LL = primsLb(2+uvw, ijkCell [0] - dijk2[0],  ijkCell [1] - dijk2[1],  ijkCell [2] - dijk2[2]);
                            velGradRight(uvw, dir2) = 0.25*info.dxInv[dir2]*(u_UR + u_LR - u_UL - u_LL);
                        }
                        
                        double tempGradLeft  = info.dxInv[idir]*(primsLb(1, ijkCell [0], ijkCell [1], ijkCell [2]) - primsLb(1, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                        double tempGradRight = info.dxInv[idir]*(primsLb(1, ijkCellR[0], ijkCellR[1], ijkCellR[2]) - primsLb(1, ijkCell [0], ijkCell [1], ijkCell[2]));
                        double uLeft[3] = {0};
                        double uRight[3] = {0};
                        double divLeft = 0;
                        double divRight = 0;
                        for (int d = 0; d < 3; d++)
                        {
                            uLeft[d]  = 0.5*(primsLb(2+d, ijkCell [0], ijkCell [1], ijkCell [2]) + primsLb(2+d, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                            uRight[d] = 0.5*(primsLb(2+d, ijkCellR[0], ijkCellR[1], ijkCellR[2]) + primsLb(2+d, ijkCell [0], ijkCell [1], ijkCell[2]));
                            divLeft  += velGradLeft (d, d);
                            divRight += velGradRight(d, d);
                        }
                        
                        fluxLeft (0) = 0.0;
                        fluxLeft (1) = (params.viscosity / params.prandtl) * tempGradLeft;
                        fluxLeft (2) = params.viscosity*(velGradLeft(1, idir) + velGradLeft(idir, 1));
                        fluxLeft (3) = params.viscosity*(velGradLeft(2, idir) + velGradLeft(idir, 2));
                        fluxLeft (4) = params.viscosity*(velGradLeft(3, idir) + velGradLeft(idir, 3));
                        fluxLeft (2+idir) -= 0.666666666667*params.viscosity*divLeft;
                        
                        
                        fluxRight(0) = 0.0;
                        fluxRight(1) = (params.viscosity / params.prandtl) * tempGradRight;
                        fluxRight(2) = params.viscosity*(velGradRight(1, idir) + velGradRight(idir, 1));
                        fluxRight(3) = params.viscosity*(velGradRight(2, idir) + velGradRight(idir, 2));
                        fluxRight(4) = params.viscosity*(velGradRight(3, idir) + velGradRight(idir, 3));
                        fluxRight(2+idir) -= 0.666666666667*params.viscosity*divRight;
                        
                        for (int d = 0; d < 3; d++)
                        {
                            fluxLeft (1) +=uLeft[d] *fluxLeft (2+d);
                            fluxRight(1) +=uRight[d]*fluxRight(2+d);
                        }
                        for (int v = 0; v < 5; v++)
                        {
                            rhsLb(v, i, j, k) += info.dxInv[idir]*(fluxRight(v) - fluxLeft(v));
                        }
                    }                    
                }
            }
        }
    }
}

void ComputeRhs(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& cons, cmf::CartesianMeshArray& rhs, InputParams& params)
{
    ComputeConv(prims, cons, rhs, params);
    if (params.viscous)
    {
        ComputeVisc(prims, cons, rhs, params);
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
    double um = sqrt(umax);
    return cmf::globalGroup.Max(um);
}

#endif