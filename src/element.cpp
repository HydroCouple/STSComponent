#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "stsmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  STSModel *model)
  : id(id),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    beta(0.0),
    length(0.0),
    depth(0.0),
    width(0.0),
    relativeHumidity(0.0),
    evaporationRate(0.0),
    evaporationHeatFlux(0.0),
    saturationVaporPressureAir(0.0),
    saturationVaporPressureWater(0.0),
    vaporPressureAir(0.0),
    vaporPressureWater(0.0),
    windSpeed(0.0),
    airTemperature(0.0),
    convectionHeatFlux(0.0),
    externalHeatFluxes(0.0),
    radiationFluxes(0.0),
    externalSoluteFluxes(nullptr),
    totalHeatBalance(0.0),
    totalRadiationFluxesHeatBalance(0.0),
    totalEvaporativeHeatFluxesBalance(0.0),
    totalConvectiveHeatFluxesBalance(0.0),
    totalExternalHeatFluxesBalance(0.0),
    totalSoluteMassBalance(nullptr),
    totalExternalSoluteFluxesMassBalance(nullptr),
    model(model)
{

  initializeSolutes();

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

}

Element::~Element()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
    delete[] mainChannelSoluteConcs; mainChannelSoluteConcs = nullptr;
    delete[] soluteExchangeCoefficients;  soluteExchangeCoefficients = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
    delete[] mcSoluteExchangeFlux; mcSoluteExchangeFlux = nullptr;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  //set upstream and downstream elements

  totalHeatBalance = totalEvaporativeHeatFluxesBalance =
      totalConvectiveHeatFluxesBalance = totalRadiationFluxesHeatBalance =
      relativeHumidity = evaporationRate = saturationVaporPressureAir = saturationVaporPressureWater =
      vaporPressureAir = vaporPressureWater = windSpeed = airTemperature =
      evaporationHeatFlux = convectionHeatFlux = 0.0;

  for(int i = 0; i < numSolutes; i++)
  {
    totalSoluteMassBalance[i] = 0.0;
    totalExternalSoluteFluxesMassBalance[i] = 0.0;
  }
}

void Element::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] mainChannelSoluteConcs; mainChannelSoluteConcs = nullptr;
    delete[] soluteExchangeCoefficients;  soluteExchangeCoefficients = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
    delete[] mcSoluteExchangeFlux; mcSoluteExchangeFlux = nullptr;
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    mainChannelSoluteConcs = new double[numSolutes];
    soluteExchangeCoefficients = new double[numSolutes];
    totalSoluteMassBalance = new double[numSolutes]();
    totalExternalSoluteFluxesMassBalance = new double[numSolutes]();
    mcSoluteExchangeFlux = new double[numSolutes]();
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;

  if(xSectionArea.value >= 0)
  {
    double volume = xSectionArea.value * length;
    double rho_cp = model->m_waterDensity * model->m_cp;
    double rho_cp_vol = rho_cp * volume;
    double top_area = beta * width * length;

    double eTemp = T[index];
    mcHeatExchangeFlux = rho_cp * temperatureExchangeCoefficient * depth * length * (mainChannelTemperature - eTemp ) / (beta * width);
    DTDt += mcHeatExchangeFlux / rho_cp_vol;

    //External sources
    {
      DTDt += radiationFluxes  * length * width * beta / rho_cp_vol;

      DTDt += externalHeatFluxes / rho_cp_vol;
    }

    if(model->m_useEvaporation)
    {
      DTDt += evaporationHeatFlux * top_area / rho_cp_vol;
      if(model->m_useConvection)
      {
        DTDt += convectionHeatFlux * top_area / rho_cp_vol;
      }
    }
  }

  return DTDt;
}

void Element::computeDTDtEvaporation()
{

  double Le = 1000.0 * (2499.0 - 2.36 * temperature.value);

  saturationVaporPressureAir = 0.61275 * exp(17.27 * airTemperature / (237.3 + airTemperature));

  saturationVaporPressureWater = 0.61275 * exp(17.27 * temperature.value / (237.3 + temperature.value));

  vaporPressureAir = relativeHumidity * saturationVaporPressureAir / 100.0;

  vaporPressureWater = relativeHumidity * saturationVaporPressureWater / 100.0;

  windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);

  evaporationRate = windFunction * (saturationVaporPressureWater - vaporPressureAir);

  evaporationHeatFlux = -Le * evaporationRate * model->m_waterDensity;
}

void Element::computeDTDtConvection()
{

  double Le = 1000.0 * (2499.0 - 2.36 * temperature.value);
  windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);
  convectionHeatFlux = - model->m_waterDensity * Le * windFunction * model->m_bowensCoeff * (temperature.value - airTemperature);
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;
  double volume = xSectionArea.value * length;
  //  double rho_vol = model->m_waterDensity * volume;

  mcSoluteExchangeFlux[soluteIndex] = soluteExchangeCoefficients[soluteIndex] * depth * length * (mainChannelSoluteConcs[soluteIndex] - S[index]) / (beta * width);
  DSoluteDt += mcSoluteExchangeFlux[soluteIndex] / (volume);

  //Add external sources
  {
    DSoluteDt += externalSoluteFluxes[soluteIndex] / volume;
  }

  return DSoluteDt;
}

void Element::calculateHydraulicVariables()
{
  if(!xSectionArea.isBC)
  {
    xSectionArea.value = depth * beta * width;
  }
}

double Element::computeDispersionFactor() const
{
  double maxFactor = 2.0 * temperatureExchangeCoefficient / (beta * beta * length * length);

  for(int i = 0 ; i < numSolutes ; i++)
  {
    maxFactor = max(maxFactor , 2.0 * soluteExchangeCoefficients[i] / (beta * beta * length * length));
  }

  return maxFactor ;
}

void Element::computeHeatBalance(double timeStep)
{
  double radiationEnergy = radiationFluxes * length * width * beta * timeStep / 1000.0;
  totalRadiationFluxesHeatBalance += radiationEnergy;

  double externalEnergy = externalHeatFluxes * xSectionArea.value * length * timeStep / 1000.0;
  totalExternalHeatFluxesBalance += externalEnergy;

  double totalHeatEnergy = model->m_waterDensity * model->m_cp * xSectionArea.value * length * (temperature.value - prevTemperature.value) / 1000.0;
  totalHeatBalance +=  totalHeatEnergy;

  double totalEvaporationHeat = evaporationHeatFlux * length * width * beta * timeStep / 1000.0;
  totalEvaporativeHeatFluxesBalance += totalEvaporationHeat;

  double totalConvectiveHeat = convectionHeatFlux * length * width * beta * timeStep / 1000.0;
  totalConvectiveHeatFluxesBalance += totalConvectiveHeat;

}

void Element::computeSoluteBalance(double timeStep, int soluteIndex)
{
  double externalMass = externalSoluteFluxes[soluteIndex] * xSectionArea.value * length * timeStep;
  totalExternalSoluteFluxesMassBalance[soluteIndex] += externalMass;

  double totalMass = model->m_waterDensity * xSectionArea.value * length * (soluteConcs[soluteIndex].value - prevSoluteConcs[soluteIndex].value) ;
  totalSoluteMassBalance[soluteIndex] +=  totalMass;
}
