#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "stsmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  STSModel *model)
  : id(id),
    mainChannelTemperature(0),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    mainChannelSoluteConcs(nullptr),
    temperatureExchangeCoefficient(0.0),
    soluteExchangeCoefficients(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    length(0.0),
    depth(0.0),
    width(0.0),
    beta(0.0),
    bank(Bank::Lumped),
    externalHeatFluxes(0.0),
    radiationFluxes(0.0),
    externalSoluteFluxes(nullptr),
    model(model)
{
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
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  //set upstream and downstream elements
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
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    mainChannelSoluteConcs = new double[numSolutes];
    soluteExchangeCoefficients = new double[numSolutes];
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;
  double volume = xSectionArea * length;
  double rho_cp_vol = model->m_waterDensity * model->m_cp * volume;

  DTDt += temperatureExchangeCoefficient * depth * length * (mainChannelTemperature - T[index]) / volume;

  DTDt += radiationFluxes / depth / rho_cp_vol;

  DTDt += externalHeatFluxes / rho_cp_vol;

  return DTDt;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;
  double volume = xSectionArea * length;
  double rho_vol = model->m_waterDensity * volume;

  DSoluteDt += soluteExchangeCoefficients[soluteIndex] * depth * length * (mainChannelSoluteConcs[soluteIndex] - S[index]) / volume;

  //Add external sources
  {
    DSoluteDt += externalSoluteFluxes[soluteIndex] / rho_vol;
  }

  return DSoluteDt;
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
