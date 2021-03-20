#include "cmf.h"
#include "InputParams.h"
#include "InitialCondition.hpp"
#include "PrimsCons.hpp"
#include "Rhs.hpp"
#include <chrono>
using cmf::print;
std::string GetInputFile(int argc, char** argv)
{
	if (argc<=1) return "input.ptl";
	std::string output(argv[1]);
	return output;
}

int main(int argc, char** argv)
{
	std::string inputFile = GetInputFile(argc, argv);
	PTL::Interactive ptlInter(argc, argv, &cmf::mainInput);
	cmf::ReadInput(inputFile);
    cmf::globalSettings = cmf::GlobalSettings(cmf::mainInput["GlobalSettings"]);
    cmf::CreateParallelContext(&argc, &argv);
	
    PTL::PropertySection& inputSection = cmf::mainInput["Solver"];
	InputParams params(inputSection);
	cmf::CartesianMeshInputInfo inputInfo(cmf::mainInput["Domain"]);
    cmf::CartesianMesh domain(inputInfo);
	
	auto& prims = domain.DefineVariable("prims", sizeof(double), {5});
	auto& cons  = domain.DefineVariable("cons",  sizeof(double), {5});
	auto& rhs   = domain.DefineVariable("rhs",   sizeof(double), {5});
	
	prims.ComponentName({0}) = "P";
	prims.ComponentName({1}) = "T";
	prims.ComponentName({2}) = "U";
	prims.ComponentName({3}) = "V";
	prims.ComponentName({4}) = "W";
	
	cons.ComponentName({0}) = "Rho";
	cons.ComponentName({1}) = "RhoE";
	cons.ComponentName({2}) = "RhoU";
	cons.ComponentName({3}) = "RhoV";
	cons.ComponentName({4}) = "RhoW";
	
	InitialCondition(prims, rhs, params);
	PrimsToCons(prims, cons, params);
	double elapsedTime = 0.0;
	for (int nt = 0; nt < params.maxStep; nt++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		ComputeRhs(prims, cons, rhs, params);
		Advance(cons, rhs, params);
		ConsToPrims(prims, cons, params);
		cons.Exchange();
		prims.Exchange();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		double timeMS = 1000*elapsed.count();
		if (cmf::globalGroup.IsRoot())
		{
			print("Timestep", nt, "\nElapsed:", timeMS, "ms");
		}
        elapsedTime += timeMS;
	}
	return 0;
}