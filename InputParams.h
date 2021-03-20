#ifndef INPUT_PARAMS_H
#define INPUT_PARAMS_H
#include "PTL.h"
struct InputParams
{
    double gamma;
    double cp;
    double deltaT;
    int maxStep;
    int centOrder;
    InputParams(PTL::PropertySection& section)
    {
        section["gamma"].MapTo(&gamma)     = new PTL::PTLDouble(1.4, "Specific heat ratio");
        section["cp"].MapTo(&cp)           = new PTL::PTLDouble(1005.0, "Specific heat");
        section["deltaT"].MapTo(&deltaT)   = new PTL::PTLDouble(1e-4, "Specific heat");
        section["maxStep"].MapTo(&maxStep) = new PTL::PTLInteger(1000, "Number of timesteps");
        section["centOrder"].MapTo(&centOrder) = new PTL::PTLInteger(4, "Order of central scheme");
        
        
        section.StrictParse();
    }
};

#endif