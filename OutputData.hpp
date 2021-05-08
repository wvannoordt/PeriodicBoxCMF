#ifndef OUTPUT_DATA_HPP
#define OUTPUT_DATA_HPP

#include "cmf.h"
#include "InputParams.h"
using cmf::ZFill;
using cmf::strformat;
using cmf::print;
std::string spaces(int i)
{
    std::string output = "";
    if (i==0) return output;
    while ((output+= " ").length()<i) {}
    return output;
}
void OutputData(int nt, cmf::CartesianMeshArray& array, InputParams& params)
{
#if (!CMF_IS3D)
    print("not working for 2D");
    abort();
#endif
    int iGuards = params.guardOutput?1:0;
    double multGuard = 1.0 + iGuards;
    int numBlocksWritten = 0;
    std::string filenamePrototype = "output/block_{}_proc_{}_nt_{}.vtk";
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
                std::string filename = strformat(filenamePrototype, ZFill(p, 4), ZFill(cmf::globalGroup.Rank(), 4), ZFill(nt, 7));
                
                int ni = block.imax - block.imin + 2*iGuards*block.exchangeI;
                int nj = block.jmax - block.jmin + 2*iGuards*block.exchangeJ;
                int nk = block.kmax - block.kmin + 2*iGuards*block.exchangeK;
                int numCells = ni*nj*nk;
                int numPoints = (ni+1)*(nj+1)*(nk+1);
                
                std::ofstream myfile;
                myfile.open(filename.c_str());
                numBlocksWritten++;
                myfile << "<VTKFile type=\"UnstructuredGrid\">" << std::endl;
                myfile << spaces(4) <<"<UnstructuredGrid>" << std::endl;
                myfile << spaces(8) <<"<Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">" << std::endl;
                myfile << spaces(12) << "<CellData>" << "..." << "</CellData>" << std::endl;
                myfile << spaces(12) << "<Points>" << std::endl;
                myfile << spaces(16) << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
                int lw = 0;
                std::string sp = spaces(20);
                const char* spp = sp.c_str();
                for (cmf::cell_t k = block.kmin-iGuards*block.exchangeK; k <= block.kmax+iGuards*block.exchangeK; k++)
                {
                    for (cmf::cell_t j = block.jmin-iGuards*block.exchangeJ; j <= block.jmax+iGuards*block.exchangeJ; j++)
                    {
                        for (cmf::cell_t i = block.imin-iGuards*block.exchangeI; i <= block.imax+iGuards*block.exchangeI; i++)
                        {
                            double xx = info.blockBounds[0]+i*info.dx[0];
                            double yy = info.blockBounds[2]+j*info.dx[1];
                            double zz = 0.0;
#if (CMF_IS3D)
                            zz = info.blockBounds[4]+k*info.dx[2];
#endif
                            lw++;
                            myfile << spp << xx << " " << yy << " " << zz << std::endl;
                        }
                    }
                }
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(12) << "</Points>" << std::endl;
                myfile << spaces(12) << "<Cells>" << std::endl;
                myfile << spaces(16) << "<DataArray type=\"Int32\" Name==\"connectivity\" format=\"ascii\">" << std::endl;
                int cellNum = 0;
                for (cmf::cell_t k = block.kmin-iGuards*block.exchangeK; k < block.kmax+iGuards*block.exchangeK; k++)
                {
                    for (cmf::cell_t j = block.jmin-iGuards*block.exchangeJ; j < block.jmax+iGuards*block.exchangeJ; j++)
                    {
                        for (cmf::cell_t i = block.imin-iGuards*block.exchangeI; i < block.imax+iGuards*block.exchangeI; i++)
                        {
                            
                            
                        }
                    }
                }                
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(16) << "<DataArray type=\"Int32\" Name==\"offsets\" format=\"ascii\">" << std::endl;
                
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(16) << "<DataArray type=\"Uint8\" Name==\"types\" format=\"ascii\">" << std::endl;
                
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(12) << "</Cells>" << std::endl;
                myfile << spaces(8) << "</Piece>" << std::endl;
                myfile << spaces(4) << "</UnstructuredGrid>" << std::endl;
                myfile << "</VTKFile>" << std::endl;
                myfile << (ni+1) << " ";
                myfile << (nj+1) << " ";
                myfile << (nk+1) << std::endl;
                
                double origin[3] = {0.0};
                double spacing[3] = {1.0};
                for (int i = 0; i < 3; i++) origin[i] = multGuard*info.blockBounds[2*i]-iGuards*info.exchangeDim[i]*info.dx[i];
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
        
        int* allNumBlocksWritten = cmf::globalGroup.SharedValues(numBlocksWritten);
        int totalNumBlocksWritten = cmf::globalGroup.Sum(numBlocksWritten);
        
        if (cmf::globalGroup.IsRoot())
        {
            std::ofstream myfile;
            std::string metaFile = "output/db_" + ZFill(nt, 7) + ".visit";
            myfile.open(metaFile.c_str());
            myfile << "!NBLOCKS " << totalNumBlocksWritten << std::endl;
            for (int i = 0; i < cmf::globalGroup.Size(); i++)
            {
                int rank = i;
                int numBlocksLoc = allNumBlocksWritten[i];
                for (int j = 0; j < numBlocksLoc; j++)
                {
                    std::string filename = strformat(filenamePrototype, ZFill(j, 4), ZFill(i, 4), ZFill(nt, 7));
                    myfile << filename << std::endl;
                }
            }
            myfile.close();
        }
    }
}

#endif