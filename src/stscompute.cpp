#include "stsmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "iboundarycondition.h"
#include  "stscomponent.h"


using namespace std;

void STSModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    applyBoundaryConditions(m_currentDateTime);

    if(m_component)
      m_component->applyInputValues();

    m_prevTimeStep = m_timeStep;

    m_timeStep = computeTimeStep();

    //Retrieve external data from other coupled models
    if(m_retrieveCouplingDataFunction)
    {
      (*m_retrieveCouplingDataFunction)(this, m_currentDateTime);
    }

    calculateHydraulicVariables();

    solve(m_timeStep);

    m_prevDateTime = m_currentDateTime;
    m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime = std::min(m_nextOutputTime + m_outputInterval / 86400.0 , m_endDateTime);
    }

    if(m_verbose)
    {
      printStatus();
    }
  }
}

void STSModel::prepareForNextTimeStep()
{

  m_minTemp = std::numeric_limits<double>::max();
  m_maxTemp = std::numeric_limits<double>::lowest();

  std::fill(m_maxSolute.begin(), m_maxSolute.end(), m_maxTemp);
  std::fill(m_minSolute.begin(), m_minSolute.end(), m_minTemp);

  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeHeatBalance(m_timeStep);
    m_totalHeatBalance += element->totalHeatBalance;
    m_totalRadiationHeatBalance += element->totalRadiationFluxesHeatBalance;
    m_totalEvaporationHeatBalance += element->totalEvaporativeHeatFluxesBalance;
    m_totalConvectiveHeatBalance += element->totalConvectiveHeatFluxesBalance;
    m_totalExternalHeatFluxBalance += element->totalExternalHeatFluxesBalance;

    element->prevTemperature.copy(element->temperature);

    m_minTemp = min(m_minTemp , element->temperature.value);
    m_maxTemp = max(m_maxTemp , element->temperature.value);

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->computeSoluteBalance(m_timeStep, j);
      m_totalSoluteMassBalance[j] += element->totalSoluteMassBalance[j];
      m_totalExternalSoluteFluxMassBalance[j] += element->totalExternalSoluteFluxesMassBalance[j];

      element->prevSoluteConcs[j].copy(element->soluteConcs[j]);

      m_minSolute[j] = min(m_minSolute[j] , (float)element->soluteConcs[j].value);
      m_maxSolute[j] = max(m_maxSolute[j] , (float)element->soluteConcs[j].value);
    }
  }

}

void STSModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_totalHeatBalance = 0.0;
  m_totalRadiationHeatBalance = 0.0;
  m_totalEvaporationHeatBalance = 0.0;
  m_totalConvectiveHeatBalance = 0.0;
  m_totalExternalHeatFluxBalance = 0.0;

  std::fill(m_totalSoluteMassBalance.begin(), m_totalSoluteMassBalance.end(), 0.0);
  std::fill(m_totalExternalSoluteFluxMassBalance.begin(), m_totalExternalSoluteFluxMassBalance.end(), 0.0);

  applyBoundaryConditions(m_currentDateTime);

  //Write initial output
  writeOutput();

  m_nextOutputTime += m_outputInterval / 86400.0;
}

void STSModel::applyBoundaryConditions(double dateTime)
{
  //reset external fluxes
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->externalHeatFluxes = 0.0;
    element->radiationFluxes = 0.0;

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->externalSoluteFluxes[j] = 0.0;
    }
  }

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->applyBoundaryConditions(dateTime);
  }
}

double STSModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;

  double maxCourantFactor = 0.0;// Δx / v (s^-1)

  if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
  {
    timeStep = m_minTimeStep;
    m_numCurrentInitFixedTimeSteps++;
  }
  else if(m_useAdaptiveTimeStep)
  {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size()  ; i++)
    {
      Element *element = m_elements[i];
      double dispersionFactor = element->computeDispersionFactor();


      if(dispersionFactor > maxCourantFactor)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        maxCourantFactor = dispersionFactor;
      }
    }

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;

  }

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(m_minTimeStep,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  return timeStep;
}

void STSModel::calculateHydraulicVariables()
{
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->calculateHydraulicVariables();

    if(m_useEvaporation)
    {
      element->computeDTDtEvaporation();

      if(m_useConvection)
      {
        element->computeDTDtConvection();
      }
    }
  }
}

void STSModel::solve(double timeStep)
{
  //Set initial input and output values to current values.
#ifdef USE_OPENMP
  //#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    m_solverCurrentValues[element->index] = element->temperature.value;
    m_solverOutputValues[element->index] = element->temperature.value;
  }

#ifdef USE_OPENMP
  //#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_solutes.size(); i++)
  {
    int startIndex = (i + 1) * m_elements.size();

    for(size_t j = 0 ; j < m_elements.size(); j++)
    {
      Element *element = m_elements[j];
      m_solverCurrentValues[startIndex + element->index] = element->soluteConcs[j].value;
      m_solverOutputValues[startIndex + element->index] = element->soluteConcs[j].value;
    }
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this;

  if(m_odeSolver->solve(m_solverCurrentValues.data(), m_solverSize, 0, timeStep,
                        m_solverOutputValues.data(), &STSModel::computeDYDt, &solverUserData))
  {
    m_currentDateTime = m_endDateTime;
    printf("CSH Solver failed \n");
  }
  else
  {

    //Apply computed values;
#ifdef USE_OPENMP
    //#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      double outputTemperature = m_solverOutputValues[element->index];
      element->temperature.value = outputTemperature;
    }

#ifdef USE_OPENMP
    //#pragma omp parallel for
#endif
    for(int i = 0 ; i < m_solutes.size(); i++)
    {
      int startIndex = (i + 1) * m_elements.size();

      for(size_t j = 0 ; j < m_elements.size(); j++)
      {
        Element *element = m_elements[j];
        element->soluteConcs[j].value = m_solverOutputValues[startIndex + element->index];
      }
    }
  }
}

void STSModel::computeDYDt(double t, double y[], double dydt[], void* userData)
{

  SolverUserData *solverUserData = (SolverUserData*) userData;
  STSModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
  //#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DTDt = element->computeDTDt(dt,y);
    dydt[element->index] = DTDt;
  }

#ifdef USE_OPENMP
  //#pragma omp parallel for
#endif
  for(int i = 0 ; i < modelInstance->m_solutes.size(); i++)
  {
    int startIndex = (i + 1) * modelInstance->m_elements.size();

    for(size_t j = 0 ; j < modelInstance-> m_elements.size(); j++)
    {
      Element *element = modelInstance->m_elements[j];
      double DSoluteDt = element->computeDSoluteDt(dt,&y[startIndex], i);
      dydt[startIndex + element->index] = DSoluteDt;
    }
  }
}
