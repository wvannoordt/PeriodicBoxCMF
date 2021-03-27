#include "cmf.h"
#include "InputParams.h"
#include "InitialCondition.hpp"
#include "PrimsCons.hpp"
#include "Rhs.hpp"
#include "OutputData.hpp"
#include "Integrate.hpp"
#include "GetTimestep.hpp"
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
	
	auto& k1 = domain.DefineVariable("k1", sizeof(double), {5});
	auto& k2 = domain.DefineVariable("k2", sizeof(double), {5});
	auto& k3 = domain.DefineVariable("k3", sizeof(double), {5});
	auto& k4 = domain.DefineVariable("k4", sizeof(double), {5});
	
	auto& c1 = domain.DefineVariable("c1", sizeof(double), {5});
	auto& c2 = domain.DefineVariable("c2", sizeof(double), {5});
	auto& c3 = domain.DefineVariable("c3", sizeof(double), {5});
	
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
	
	switch(params.cfdCase)
	{
		case CFDCase::TGV:
		{
			InitialConditionTgv(prims, rhs, params);
			break;
		}
		case CFDCase::IsentropicVortex:
		{
			InitialConditionVort(prims, rhs, params);
			break;
		}
		default: 
		{
			print("Bad case");
			KILL;
		}
	}
	PrimsToCons(prims, cons, params);
	OutputData(0, prims, params);
	double elapsedTime = 0.0;
	bool isRoot = cmf::globalGroup.IsRoot();
	
	double deltaT = GetTimestep(prims, params);
	
	
	double time = 0;
	TimeSeries enstrophySeries("series/enstrophy.csv", 50);
	TimeSeries energySeries("series/kinetic.csv", 50);
	for (int nt = 0; nt <= params.maxStep; nt++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		double umax = UMax(prims);
		ZeroRhs(rhs);
		if (params.useRK4)
		{
			ZeroRhs(c1);
			ZeroRhs(c2);
			ZeroRhs(c3);
			
			ZeroRhs(k1);
			ZeroRhs(k2);
			ZeroRhs(k3);
			ZeroRhs(k4);
			
			ComputeRhs(prims, cons, k1, params);
			PlusEqualsKX(rhs, deltaT/6.0, k1);
			PlusEqualsKX(c1, 1.0, cons);
			PlusEqualsKX(c1, 0.5*deltaT, k1);
			
			c1.Exchange();
			ConsToPrims(prims, c1, params);
			
			ComputeRhs(prims, c1, k2, params);
			PlusEqualsKX(rhs, deltaT/3.0, k2);
			PlusEqualsKX(c2, 1.0, cons);
			PlusEqualsKX(c2, 0.5*deltaT, k2);
			
			c2.Exchange();
			ConsToPrims(prims, c2, params);
			
			ComputeRhs(prims, c2, k3, params);
			PlusEqualsKX(rhs, deltaT/3.0, k3);
			PlusEqualsKX(c3, 1.0, cons);
			PlusEqualsKX(c3, 0.5*deltaT, k3);
			
			c3.Exchange();
			ConsToPrims(prims, c3, params);
			
			ComputeRhs(prims, c3, k4, params);
			PlusEqualsKX(rhs, deltaT/6.0, k4);
			
			PlusEqualsKX(cons, 1.0, rhs);
			cons.Exchange();
			ConsToPrims(prims, cons, params);
		}
		else
		{
			ComputeRhs(prims, cons, rhs, params);
			PlusEqualsKX(cons, deltaT, rhs);
			cons.Exchange();
			ConsToPrims(prims, cons, params);
			prims.Exchange();
		}
		
		double integratedEnstrophy = Enstrophy(prims);
		double integratedKE = KineticEnergy(prims, params);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		double timeMS = 1000*elapsed.count();
		double effTimeStep = GetTimestep(prims, params);
		if (isRoot)
		{
			double maxCFL = params.CFL*deltaT/effTimeStep;
			print("Timestep", nt);
			print("Elapsed:", timeMS, "ms");
			print("Umax:", umax);
			print("Remaining time:", GetRemainingTime(nt, params.maxStep, elapsedTime));
			print("Enstrophy:", integratedEnstrophy);
			print("Energy:", integratedKE);
			print("Time:", time);
			print("deltaT:", deltaT);
			print("CFL:", maxCFL);
			print("");
			if (!(params.maxCFL > maxCFL))
			{
				KILL;
			}
		}
		if (nt%params.outputInterval==0 && nt > 0)
		{
			OutputData(nt, prims, params);
		}
        elapsedTime += timeMS;
		time += deltaT;
		enstrophySeries.AddEntry(time, integratedEnstrophy);
		energySeries.AddEntry(time, integratedKE);
	}
	return 0;
}