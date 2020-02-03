/*!
*  \file    stmproject.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
*  either version 3 of the License, or (at your option) any later version.
*  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
*  \date 2018
*  \pre
*  \bug
*  \todo
*  \warning
*/

#include "stsmodel.h"
#include "stscomponent.h"
#include "spatial/point.h"
#include "spatial/network.h"
#include "element.h"
#include "elementjunction.h"
#include "spatial/edge.h"
#include "iboundarycondition.h"

using namespace std;

STSModel::STSModel(STSComponent *component)
  : QObject(component),
    m_timeStep(0.0001), //seconds
    m_maxTimeStep(0.5), //seconds
    m_minTimeStep(0.001), //seconds,
    m_timeStepRelaxationFactor(0.8),
    m_numInitFixedTimeSteps(2),
    m_numCurrentInitFixedTimeSteps(0),
    m_printFrequency(10),
    m_currentPrintCount(0),
    m_flushToDiskFrequency(10),
    m_currentflushToDiskCount(0),
    m_computeDispersion(false),
    m_useAdaptiveTimeStep(true),
    m_verbose(false),
    m_numHeatElementJunctions(0),
    m_odeSolver(nullptr),
    m_waterDensity(1000.0), //kg/m^3
    m_cp(4184.0), //4187.0 J/kg/C
    m_evapWindFuncCoeffA(1.505e-8),
    m_evapWindFuncCoeffB(1.600e-8),
    m_bowensCoeff(0.06266151315), //KPa / C
    #ifdef USE_NETCDF
    m_outputNetCDF(nullptr),
    #endif
    m_retrieveCouplingDataFunction(nullptr),
    m_component(component)
{
  m_odeSolver = new ODESolver(1, ODESolver::CVODE_ADAMS);
}

STSModel::~STSModel()
{

  for(Element *element : m_elements)
    delete element;

  m_elements.clear();
  m_elementsById.clear();


  for(ElementJunction *elementJunction : m_elementJunctions)
    delete elementJunction;

  m_elementJunctions.clear();
  m_elementJunctionsById.clear();

  if(m_odeSolver)
    delete m_odeSolver;


  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();
}

double STSModel::minTimeStep() const
{
  return m_minTimeStep;
}

void STSModel::setMinTimeStep(double timeStep)
{
  m_minTimeStep = timeStep;
}

double STSModel::maxTimeStep() const
{
  return m_maxTimeStep;
}

void STSModel::setMaxTimeStep(double timeStep)
{
  m_maxTimeStep = timeStep;
}

bool STSModel::useAdaptiveTimeStep() const
{
  return m_useAdaptiveTimeStep;
}

void STSModel::setUseAdaptiveTimeStep(bool use)
{
  m_useAdaptiveTimeStep = use;
}

double STSModel::timeStepRelaxationFactor() const
{
  return m_timeStepRelaxationFactor;
}

void STSModel::setTimeStepRelaxationFactor(double tStepRelaxFactor)
{
  if(tStepRelaxFactor > 0)
    m_timeStepRelaxationFactor = tStepRelaxFactor;
}

double STSModel::currentTimeStep() const
{
  return m_timeStep;
}

double STSModel::startDateTime() const
{
  return m_startDateTime;
}

void STSModel::setStartDateTime(double dateTime)
{
  m_startDateTime = dateTime;
}

double STSModel::endDateTime() const
{
  return m_endDateTime;
}

void STSModel::setEndDateTime(double dateTime)
{
  m_endDateTime = dateTime;
}

double STSModel::outputInterval() const
{
  return m_outputInterval;
}

void STSModel::setOutputInterval(double interval)
{
  m_outputInterval = interval;
}

double STSModel::currentDateTime() const
{
  return m_currentDateTime;
}

ODESolver *STSModel::solver() const
{
  return m_odeSolver;
}

double STSModel::waterDensity() const
{
  return m_waterDensity;
}

void STSModel::setWaterDensity(double value)
{
  m_waterDensity = value;
}

double STSModel::specificHeatCapacityWater() const
{
  return m_cp;
}

void STSModel::setSpecificHeatCapacityWater(double value)
{
  m_cp = value;
}

int STSModel::numSolutes() const
{
  return (int)m_solutes.size();
}

void STSModel::setNumSolutes(int numSolutes)
{
  if(numSolutes >= 0)
  {

    m_solutes.resize(numSolutes);
    m_maxSolute.resize(numSolutes);
    m_minSolute.resize(numSolutes);
    m_totalSoluteMassBalance.resize(numSolutes);
    m_totalExternalSoluteFluxMassBalance.resize(numSolutes);

    for(size_t i = 0 ; i < m_solutes.size(); i++)
    {
      m_solutes[i] = "Solute_" + std::to_string(i + 1);
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size()  ; i++)
    {
      Element *element = m_elements[i];
      element->initializeSolutes();
    }
  }
}

void STSModel::setSoluteName(int soluteIndex, const string &soluteName)
{
  m_solutes[soluteIndex] = soluteName;
}

string STSModel::solute(int soluteIndex) const
{
  return m_solutes[soluteIndex];
}

int STSModel::numElementJunctions() const
{
  return m_elementJunctions.size();
}

ElementJunction *STSModel::addElementJunction(const string &id, double x, double y, double z)
{
  if(m_elementJunctionsById.find(id) == m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = new ElementJunction(id, x, y, z, this);
    eJunction->index = m_elementJunctions.size();
    m_elementJunctions.push_back(eJunction);
    m_elementJunctionsById[id] = eJunction;
    return eJunction;
  }

  return nullptr;
}

void STSModel::deleteElementJunction(const string &id)
{
  std::unordered_map<string,ElementJunction*>::iterator eJIter =  m_elementJunctionsById.find(id) ;

  if(eJIter != m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = eJIter->second;
    m_elementJunctionsById.erase(eJIter);

    std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
    if(it != m_elementJunctions.end())
    {
      m_elementJunctions.erase(it);
    }

    delete eJunction;
  }
}

void STSModel::deleteElementJunction(int index)
{
  ElementJunction *eJunction = m_elementJunctions[index];

  m_elementJunctionsById.erase(eJunction->id);

  std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
  if(it != m_elementJunctions.end())
    m_elementJunctions.erase(it);

  delete eJunction;
}

ElementJunction *STSModel::getElementJunction(const string &id)
{
  return m_elementJunctionsById[id];
}

ElementJunction *STSModel::getElementJunction(int index)
{
  return m_elementJunctions[index];
}

int STSModel::numElements() const
{
  return m_elements.size();
}

Element *STSModel::addElement(const string &id, ElementJunction *upStream, ElementJunction *downStream)
{
  if(upStream && downStream)
  {
    Element *element = new Element(id, upStream, downStream, this);
    element->index = m_elements.size();
    m_elements.push_back(element);
    m_elementsById[id] = element;
    return element;
  }

  return nullptr;
}

void STSModel::deleteElement(const string &id)
{
  unordered_map<string,Element*>::iterator eIter = m_elementsById.find(id);

  if(eIter != m_elementsById.end())
  {
    Element *element = eIter->second;
    m_elementsById.erase(eIter);

    vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);
    if(it != m_elements.end())
      m_elements.erase(it);

    delete element;
  }
}

void STSModel::deleteElement(int index)
{
  Element *element = m_elements[index];
  m_elementJunctionsById.erase(element->id);

  vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);

  if(it != m_elements.end())
    m_elements.erase(it);

  delete element;
}

Element *STSModel::getElement(const string &id)
{
  return m_elementsById[id];
}

Element *STSModel::getElement(int index)
{
  return m_elements[index];
}

RetrieveCouplingData STSModel::retrieveCouplingDataFunction() const
{
  return m_retrieveCouplingDataFunction;
}

void STSModel::setRetrieveCouplingDataFunction(RetrieveCouplingData retrieveCouplingDataFunction)
{
  m_retrieveCouplingDataFunction = retrieveCouplingDataFunction;
}

bool STSModel::initialize(list<string> &errors)
{
  bool initialized = initializeInputFiles(errors) &&
                     initializeTimeVariables(errors) &&
                     initializeElements(errors) &&
                     initializeSolver(errors) &&
                     initializeOutputFiles(errors) &&
                     initializeBoundaryConditions(errors);


  if(initialized)
  {
    applyInitialConditions();
  }

  return initialized;
}

bool STSModel::finalize(std::list<string> &errors)
{
  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();

  return true;
}

bool STSModel::initializeTimeVariables(std::list<string> &errors)
{
  if(m_startDateTime >= m_endDateTime)
  {
    errors.push_back("End datetime must be greater than startdatetime");
    return false;
  }

  if( (m_endDateTime - m_startDateTime) *  86400.0 < m_minTimeStep )
  {
    errors.push_back("Make sure timestep is less than the simulation interval");
    return false;
  }

  if(m_minTimeStep <=  0 || m_maxTimeStep <= 0)
  {
    errors.push_back("Make sure time steps are greater 0");
    return false;
  }

  if(m_minTimeStep >= m_maxTimeStep)
  {
    errors.push_back("");
    return false;
  }

  m_numCurrentInitFixedTimeSteps = 0;

  m_currentDateTime = m_startDateTime;
  m_nextOutputTime = m_currentDateTime;

  m_currentPrintCount = 0;
  m_currentflushToDiskCount = 0;

  return true;
}

bool STSModel::initializeElements(std::list<string> &errors)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
    size_t numElements = elementJunction->incomingElements.size() + elementJunction->outgoingElements.size();

    switch (numElements)
    {
      case 0:
        elementJunction->junctionType = ElementJunction::NoElement;
        break;
      case 1:
        elementJunction->junctionType = ElementJunction::SingleElement;
        break;
      case 2:
        elementJunction->junctionType = ElementJunction::DoubleElement;
        break;
      default:
        elementJunction->junctionType = ElementJunction::MultiElement;
        break;
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->index = i;
    element->initialize();
  }

  //Set tempearture continuity junctions
  m_numHeatElementJunctions = 0;

  //Number of junctions where continuity needs to be enforced.
  m_numSoluteElementJunctions.resize(m_solutes.size(), 0);

  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
  }

  return true;
}

bool STSModel::initializeSolver(std::list<string> &errors)
{
  m_solverSize = m_elements.size();

  if(m_solutes.size())
  {
    m_solverSize += m_elements.size() * m_solverSize;
  }

  m_solverCurrentValues.resize(m_solverSize, 0.0);
  m_solverOutputValues.resize(m_solverSize, 0.0);

  m_odeSolver->setSize(m_solverSize);
  m_odeSolver->initialize();

  return true;
}

bool STSModel::initializeBoundaryConditions(std::list<string> &errors)
{
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->clear();
    boundaryCondition->findAssociatedGeometries();
    boundaryCondition->prepare();
  }

  return true;
}

bool STSModel::findProfile(Element *from, Element *to, std::vector<Element *> &profile)
{
  if(from == to)
  {
    profile.push_back(from);
    return true;
  }
  else
  {
    for(Element *outgoing : from->downstreamJunction->outgoingElements)
    {
      if(outgoing == to)
      {
        profile.push_back(from);
        profile.push_back(outgoing);
        return true;
      }
      else if(findProfile(outgoing, to, profile))
      {
        profile.insert(profile.begin(), from);
        return true;
      }
    }
  }

  return false;
}
