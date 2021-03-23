#include "cmf.h"
#include "InputParams.h"
#include "InitialCondition.hpp"
#include "PrimsCons.hpp"
#include "Rhs.hpp"
#include "OutputData.hpp"
#include "Integrate.hpp"
#include <chrono>

struct TimeSeries
{
	std::vector<double> times;
	std::vector<double> values;
	std::string filename;
	int outputInterval;
	int index;
	TimeSeries(std::string filename_in, int outputInterval_in=1000)
	{
		std::ofstream outfile;
		if (cmf::globalGroup.IsRoot())
		{
			outfile.open(filename.c_str());
			outfile<<"";
			outfile.close();
		}
		filename = filename_in;
		outputInterval = outputInterval_in;
		times.resize(outputInterval, 0.0);
		values.resize(outputInterval, 0.0);
		index = 0;
	}
	void Write()
	{
		if (cmf::globalGroup.IsRoot())
		{
			print("Output", filename);
			std::ofstream outfile;
			outfile.open(filename.c_str(), std::ios_base::app);
			for (int i = 0; i < index; i++)
			{
				outfile << times[i] << ", " << values[i] << std::endl;
			}
			outfile.close();
		}
		index = 0;
		cmf::globalGroup.Synchronize();
	}
	void AddEntry(double time, double val)
	{
		times[index] = time;
		values[index] = val;
		index++;
		if (index==outputInterval)
		{
			Write();
		}
	}
	~TimeSeries(void)
	{
		Write();
	}
};

using cmf::print;
using cmf::strformat;
std::string GetInputFile(int argc, char** argv)
{
	if (argc<=1) return "input.ptl";
	std::string output(argv[1]);
	return output;
}

std::string GetRemainingTime(int ntCurrent, int ntTotal, double elapsedTimeMs)
{
	int numStepsPassed = ntCurrent + 1;
	int numStepsLeft = ntTotal - ntCurrent;
	double avTimeMs = elapsedTimeMs/numStepsPassed;
	double timeLeftMs = avTimeMs*numStepsLeft;
	int hoursLeft = floor(timeLeftMs/(1000*3600));
	timeLeftMs -= ((double)hoursLeft)*(1000*3600);
	int minsLeft = floor(timeLeftMs/(1000*60));
	timeLeftMs -= ((double)minsLeft)*(1000*60);
	double secondsLeft = timeLeftMs/1000;
	return strformat("{} h, {} m, {} s", hoursLeft, minsLeft, secondsLeft);
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
	bool isRoot = cmf::globalGroup.IsRoot();
	double time = 0;
	TimeSeries enstrophySeries("series/enstrophy.csv", 50);
	for (int nt = 0; nt <= params.maxStep; nt++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		ZeroRhs(rhs);
		ComputeRhs(prims, cons, rhs, params);
		Advance(cons, rhs, params);
		double umax = UMax(prims);
		ConsToPrims(prims, cons, params);
		cmf::globalGroup.Synchronize();
		cons.Exchange();
		cmf::globalGroup.Synchronize();
		prims.Exchange();
		cmf::globalGroup.Synchronize();
		double integratedEnstrophy = Enstrophy(prims);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		double timeMS = 1000*elapsed.count();
		if (isRoot) print("Timestep", nt, "\nElapsed:", timeMS, "ms", "\nUmax:", umax, "\nRemaining time:", GetRemainingTime(nt, params.maxStep, elapsedTime), "\nEnstrophy:", integratedEnstrophy, "\nTime:", time, "\n");
		if (nt%params.outputInterval==0)
		{
			OutputData(nt, prims);
		}
        elapsedTime += timeMS;
		if (umax>params.uLimit)
		{
			KILL;
		}
		time += params.deltaT;
		enstrophySeries.AddEntry(time, integratedEnstrophy);
	}
	return 0;
}