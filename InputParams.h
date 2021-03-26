#ifndef INPUT_PARAMS_H
#define INPUT_PARAMS_H
#include "PTL.h"
namespace CFDCase
{
	enum CFDCase
	{
		TGV,
		IsentropicVortex
	};
}

std::string CaseStr(int i)
{
	switch (i)
	{
		case CFDCase::TGV:
        {
            return "TGV";
        }
		case CFDCase::IsentropicVortex:
        {
            return "IsentropicVortex";
        }
		default:
        {
            return PTL_AUTO_ENUM_TERMINATOR;
        }
	}
}
struct InputParams
{
    double gamma;
    double cp;
    double deltaT;
    int maxStep;
    int centOrder;
    int outputInterval;
    double uLimit;
    bool useRK4;
    bool guardOutput;
    double viscosity;
    int cfdCase;
    double prandtl;
    InputParams(PTL::PropertySection& section)
    {
        section["gamma"].MapTo(&gamma)                   = new PTL::PTLDouble(1.4, "Specific heat ratio");
        section["cp"].MapTo(&cp)                         = new PTL::PTLDouble(1005.0, "Specific heat");
        section["deltaT"].MapTo(&deltaT)                 = new PTL::PTLDouble(1e-4, "Specific heat");
        section["maxStep"].MapTo(&maxStep)               = new PTL::PTLInteger(1000, "Number of timesteps");
        section["centOrder"].MapTo(&centOrder)           = new PTL::PTLInteger(4, "Order of central scheme");
        section["outputInterval"].MapTo(&outputInterval) = new PTL::PTLInteger(1000, "Timesteps between data outputs");
        section["uLimit"].MapTo(&uLimit)                 = new PTL::PTLDouble(0.5, "Max velocity");
        section["useRK4"].MapTo(&useRK4)                 = new PTL::PTLBoolean(true, "Use RK4 time integration?");
        section["guardOutput"].MapTo(&guardOutput)       = new PTL::PTLBoolean(false, "output the guard cells?");
        section["viscosity"].MapTo(&viscosity)           = new PTL::PTLDouble(0.0, "Laminar viscosity");
        section["cfdCase"].MapTo(&cfdCase)               = new PTL::PTLAutoEnum(CFDCase::TGV, CaseStr, "Which case to run");
        section["prandtl"].MapTo(&prandtl)               = new PTL::PTLDouble(0.72, "Laminar viscosity");
        
        
        
        section.StrictParse();
    }
};

#endif