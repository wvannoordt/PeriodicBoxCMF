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
    auto vnames = array.GetComponentNames();
    std::string varsString = "";
    for (int l = 0; l < vnames.size(); l++)
    {
        varsString += ((l==0)?"":",");
        varsString += vnames[l];
    }
    int iGuards = params.guardOutput?1:0;
    double multGuard = 1.0 + iGuards;
    int numBlocksWritten = 0;
    std::string filenamePrototype = "output/block_{}_nt{}.vtr";
    int numBlocksLoc = 0;
    for (auto lb: array)
    {
        numBlocksLoc++;
    }
    int* allBlocks = cmf::globalGroup.SharedValues(numBlocksLoc);
    int blocksBefore = 0;
    int totalNumBlocksWritten = cmf::globalGroup.Sum(numBlocksLoc);
    for (int p = 0; p < cmf::globalGroup.Rank(); p++)
    {
        blocksBefore += allBlocks[p];
    }
    int blockWrittenByMe = 0;
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
                std::string filename = strformat(filenamePrototype, ZFill(blockWrittenByMe+blocksBefore, 7), ZFill(nt, 7));
                
                int nCellsi = block.imax - block.imin;
                int nCellsj = block.jmax - block.jmin;
                int nCellsk = block.kmax - block.kmin;
                int nGuardi = block.exchangeI;
                int nGuardj = block.exchangeJ;
                int nGuardk = block.exchangeK;
                int nTotali = nCellsi + 2*nGuardi;
                int nTotalj = nCellsj + 2*nGuardj;
                int nTotalk = nCellsk + 2*nGuardk;
                
                double bds[6] = {0.0};
                for (int i = 0; i < 2*CMF_DIM; i++) bds[i] = info.blockBounds[i];
                double ghostBnds[6] = {0.0};
                for (int i = 0; i < CMF_DIM; i++)
                {
                    ghostBnds[2*i] = info.blockBounds[2*i] - info.exchangeDim[i]*info.dx[i];
                    ghostBnds[2*i+1] = info.blockBounds[2*i+1] + info.exchangeDim[i]*info.dx[i];
                }
                
                std::ofstream myfile;
                myfile.open(filename.c_str());
                std::string sp20 = spaces(20);
                const char* csp20 = sp20.c_str();
                
                numBlocksWritten++;
                myfile << "<?xml version=\"1.0\"?>\n";
                myfile << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl;
                myfile << spaces(4) << strformat("<RectilinearGrid WholeExtent=\"0 {} 0 {} 0 {}\">", nTotali, nTotalj, nTotalk) << std::endl;
                myfile << spaces(8) << "<FieldData>" << std::endl;
                myfile << spaces(12) << "<DataArray type=\"Int32\" Name=\"avtRealDims\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
                myfile << spaces(16) << strformat("{} {} {} {} {} {}", nGuardi, nGuardj+nCellsi, nGuardj, nGuardj+nCellsj, nGuardk, nGuardk+nCellsk) << std::endl;
                myfile << spaces(12) << "</DataArray>" << std::endl;
                myfile << spaces(12) << "<DataArray type=\"Float64\" Name=\"avtOriginalBounds\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
                myfile << spaces(16) << strformat("{} {} {} {} {} {}", bds[0], bds[1], bds[2], bds[3], bds[4], bds[5]) << std::endl;
                myfile << spaces(12) << "</DataArray>" << std::endl;
                myfile << spaces(8) << "</FieldData>" << std::endl;
                myfile << spaces(8) << strformat("<Piece Extent=\"0 {} 0 {} 0 {}\">", nTotali, nTotalj, nTotalk) << std::endl;
                
                myfile << spaces(12) << "<CellData Scalars=\"" << varsString << "\">" << std::endl;
                for (int var = 0; var < vnames.size(); var++)
                {
                    myfile << spaces(16) << "<DataArray type=\"Float64\" Name=\"" << vnames[var] << "\" format=\"ascii\">" << std::endl;
                    for (cmf::cell_t k = block.kmin-block.exchangeK; k < block.kmax+block.exchangeK; k++)
                    {
                        for (cmf::cell_t j = block.jmin-block.exchangeJ; j < block.jmax+block.exchangeJ; j++)
                        {
                            for (cmf::cell_t i = block.imin-block.exchangeI; i < block.imax+block.exchangeI; i++)
                            {
                                myfile << csp20 << block(var, i, j, k) << "\n";
                            }
                        }
                    }
                    myfile << spaces(16) << "</DataArray>" << std::endl;
                }
                myfile << spaces(16) << "<DataArray type=\"UInt8\" Name=\"avtGhostZones\" format=\"ascii\">" << std::endl;
                auto isGhost = [&](int i, int j, int k) -> bool {return (i<0)||(i>=nCellsi)||(j<0)||(j>=nCellsj)||(k<0)||(k>=nCellsk);};
                for (cmf::cell_t k = block.kmin-block.exchangeK; k < block.kmax+block.exchangeK; k++)
                {
                    for (cmf::cell_t j = block.jmin-block.exchangeJ; j < block.jmax+block.exchangeJ; j++)
                    {
                        for (cmf::cell_t i = block.imin-block.exchangeI; i < block.imax+block.exchangeI; i++)
                        {
                            myfile << csp20 << (isGhost(i, j, k)?16:0) << "\n";
                        }
                    }
                }
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(12) << "</CellData>" << std::endl;
                myfile << spaces(12) << "<Coordinates>" << std::endl;
                myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[0], ghostBnds[1]) << std::endl;
                for (int i = -nGuardi; i <=nCellsi+nGuardi; i++)
                {
                    myfile << csp20 << info.blockBounds[0] + i*info.dx[0] << "\n";
                }
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[2], ghostBnds[3]) << std::endl;
                for (int j = -nGuardj; j <=nCellsj+nGuardj; j++)
                {
                    myfile << csp20 << info.blockBounds[2] + j*info.dx[1] << "\n";
                }
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(16) << strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", ghostBnds[4], ghostBnds[4]) << std::endl;
                for (int k = -nGuardk; k <=nCellsk+nGuardk; k++)
                {
                    myfile << csp20 << info.blockBounds[4] + k*info.dx[2] << "\n";
                }
                myfile << spaces(16) << "</DataArray>" << std::endl;
                myfile << spaces(12) << "</Coordinates>" << std::endl;
                myfile << spaces(8) << "</Piece>" << std::endl;
                myfile << spaces(4) << "</RectilinearGrid>" << std::endl;
                myfile << "</VTKFile>" << std::endl;
                
                
                myfile.close();
                p++;
                blockWrittenByMe++;
            }
        }
        if (cmf::globalGroup.IsRoot())
        {
            std::string filename = strformat("output/data_nt{}.vtm", ZFill(nt, 7));
            std::ofstream myfile;
            myfile.open(filename.c_str());
            myfile << "<?xml version=\"1.0\"?>\n";
            myfile << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">" << std::endl;
            myfile << spaces(4) << "<vtkMultiBlockDataSet>" << std::endl;
            myfile << spaces(8) << "<Block index =\"0\">" << std::endl;
            for (int b = 0; b < totalNumBlocksWritten; b++)
            {
                std::string blockFileName = strformat("block_{}_nt{}.vtr", ZFill(b, 7), ZFill(nt, 7));
                myfile << spaces(12) << strformat("<DataSet index=\"{}\" file=\"{}\"/>", b, blockFileName) << std::endl;
            }
            myfile << spaces(8) << "</Block>" << std::endl;
            myfile << spaces(4) << "</vtkMultiBlockDataSet>" << std::endl;
            myfile << "</VTKFile>" << std::endl;
            myfile.close();
        }
        cmf::globalGroup.Synchronize();
    }
}

#endif