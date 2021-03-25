#ifndef OUTPUT_DATA_HPP
#define OUTPUT_DATA_HPP

#include "cmf.h"
#include "InputParams.h"
using cmf::ZFill;
using cmf::strformat;
using cmf::print;
void OutputData(int nt, cmf::CartesianMeshArray& array, InputParams& params)
{
    int iGuards = params.guardOutput?1:0;
    double multGuard = 1.0 + iGuards;
    for (int process = 0; process < cmf::globalGroup.Size(); process++)
    {
        if (cmf::globalGroup.Rank()==process)
        {
            print("Output solution, rank", process);
            int p = 0;
            for (auto lb: array)
            {
                print(" >> block ", p);
                auto info = array.GetBlockInfo(lb);
                cmf::BlockArray<double, 1> block = array[lb];
                std::string filename = strformat("output/block_{}_proc_{}_nt_{}.vtk", ZFill(p, 4), ZFill(cmf::globalGroup.Rank(), 4), ZFill(nt, 7));
                std::ofstream myfile;
                myfile.open(filename.c_str());
                myfile << "# vtk DataFile Version 3.0" << std::endl;
                myfile << "vtk output" << std::endl;
                myfile << "ASCII" << std::endl;
                myfile << "DATASET STRUCTURED_POINTS" << std::endl;
                myfile << "DIMENSIONS ";
                int ni = block.imax - block.imin + 2*iGuards*block.exchangeI;
                int nj = block.jmax - block.jmin + 2*iGuards*block.exchangeJ;
                int nk = block.kmax - block.kmin + 2*iGuards*block.exchangeK;
                myfile << (ni+1) << " ";
                myfile << (nj+1) << " ";
                myfile << (nk+1) << std::endl;
                
                double origin[3] = {0.0};
                double spacing[3] = {1.0};
                for (int i = 0; i < 3; i++) origin[i] = multGuard*info.blockBounds[2*i];
                for (int i = 0; i < 3; i++) spacing[i] = info.dx[i];
                
                myfile << "ORIGIN "  << origin[0]  << " " << origin[1]  << " " << origin[2]  << std::endl;
                myfile << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
                myfile << "CELL_DATA " << ni*nj*nk << std::endl;
                for (int var = 0; var < 5; var++)
                {
                    myfile << "SCALARS " << array.ComponentName({var}) << " double"  << std::endl;
                    myfile << "LOOKUP_TABLE default" << std::endl;
                    for (cmf::cell_t k = block.kmin-iGuards*block.exchangeK; k < block.kmax+iGuards*block.exchangeK; k++)
                    {
                        for (cmf::cell_t j = block.jmin-iGuards*block.exchangeJ; j < block.jmax+iGuards*block.exchangeJ; j++)
                        {
                            for (cmf::cell_t i = block.imin-iGuards*block.exchangeI; i < block.imax+iGuards*block.exchangeI; i++)
                            {
                                myfile << block(var, i, j, k) << std::endl;
                            }
                        }
                    }
                }
                myfile.close();
                p++;
            }
        }
        cmf::globalGroup.Synchronize();
    }
}

#endif