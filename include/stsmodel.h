/*!
*  \file    stmproject.h
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU General Public License as published by the Free Software Foundation;
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

#ifndef STSMODEL_H
#define STSMODEL_H

#include "stscomponent_global.h"
#include "spatial/network.h"
#include "odesolver.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <QFileInfo>
#include <QTextStream>

class STSComponent;
struct Element;
struct ElementJunction;
class Edge;
class STSModel;

struct SolverUserData
{
    STSModel *model = nullptr;
    int variableIndex = -1;
};

class STSCOMPONENT_EXPORT STSModel : public QObject
{
    Q_OBJECT

    friend struct ElementJunction;
    friend struct Element;

  public:

    /*!
     * \brief STSModel - Constructor for the Computational engine for the Stream Temperature Model.
     */
    STSModel(STSComponent *component);

    /*!
     * \brief ~STSModel - Destructor for the Computational engine for the Stream Temperature Model.
     */
    ~STSModel();

    /*!
     * \brief minTimeStep - Minimum timestep for the model in seconds.
     * \return Returns the minimum timestep for the current model instance.
     */
    double minTimeStep() const;

    /*!
     * \brief setMinTimeStep - Sets the minimum timestep for the current model instance
     * \param timeStep - Minimum timestep in seconds.
     */
    void setMinTimeStep(double timeStep);

    /*!
     * \brief maxTimeStep
     * \return
     */
    double maxTimeStep() const;

    /*!
     * \brief setMaxTimeStep
     * \param timeStep
     */
    void setMaxTimeStep(double timeStep);

    /*!
     * \brief useAdaptiveTimeStep
     * \return
     */
    bool useAdaptiveTimeStep() const;

    /*!
     * \brief setUseAdaptiveTimeStep
     * \param use
     */
    void setUseAdaptiveTimeStep(bool use);

    /*!
     * \brief timeStepRelaxationFactor
     * \return
     */
    double timeStepRelaxationFactor() const;

    /*!
     * \brief setTimeStepRelaxationFactor
     * \param tStepRelaxFactor
     */
    void setTimeStepRelaxationFactor(double tStepRelaxFactor);

    /*!
     * \brief currentTimeStep
     * \return
     */
    double currentTimeStep() const;

    /*!
     * \brief startDateTime
     * \return
     */
    double startDateTime() const;

    /*!
     * \brief setStartDate
     * \param dateTime
     */
    void setStartDateTime(double dateTime);

    /*!
     * \brief endDateTime
     * \return
     */
    double endDateTime() const;

    /*!
     * \brief setEndDateTime
     * \param dateTime
     */
    void setEndDateTime(double dateTime);

    /*!
     * \brief outputInterval
     * \return
     */
    double outputInterval() const;

    /*!
     * \brief setOutputInterval
     * \param interval
     */
    void setOutputInterval(double interval);

    /*!
     * \brief currentDateTime
     * \return
     */
    double currentDateTime() const;

    /*!
     * \brief solver
     * \return
     */
    ODESolver *solver() const;


    /*!
     * \brief waterDensity
     * \return
     */
    double waterDensity() const;

    /*!
     * \brief setWaterDensity
     * \param value
     */
    void setWaterDensity(double value);

    /*!
     * \brief specificHeatCapacityWater
     * \return
     */
    double specificHeatCapacityWater() const;

    /*!
     * \brief setSpecificHeatCapacityWater
     * \param value
     */
    void setSpecificHeatCapacityWater(double value);

    /*!
     * \brief numSolutes
     * \return
     */
    int numSolutes() const;

    /*!
     * \brief setNumSolutes
     * \param numSolutes
     */
    void  setNumSolutes(int numSolutes);

    /*!
     * \brief setSoluteNuame
     * \param soluteIndex
     * \param soluteName
     */
    void setSoluteName(int soluteIndex, const std::string &soluteName);

    /*!
     * \brief solute
     * \param soluteIndex
     * \return
     */
    std::string solute(int soluteIndex) const;

    /*!
     * \brief numElementJunctions
     * \return
     */
    int numElementJunctions() const;

    /*!
     * \brief addControlVolumeNode
     * \param id
     * \param x
     * \param y
     * \param z
     * \return
     */
    ElementJunction *addElementJunction(const std::string &id, double x = 0, double y = 0, double z = 0);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(const std::string &id);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(int id);

    /*!
     * \brief getElementJunction
     * \param id
     * \return
     */
    ElementJunction *getElementJunction(const std::string &id) ;

    /*!
     * \brief getElementJunction
     * \param index
     * \return
     */
    ElementJunction *getElementJunction(int index) ;

    /*!
     * \brief numElements
     * \return
     */
    int numElements() const;

    /*!
     * \brief addElement
     * \param id
     * \param fromElement
     * \param toElement
     * \return
     */
    Element *addElement(const std::string &id, ElementJunction *upStream, ElementJunction *downStream);

    /*!
     * \brief removeElement
     * \param id
     */
    void deleteElement(const std::string &id);

    /*!
     * \brief removeElement
     * \param index
     */
    void deleteElement(int index);

    /*!
     * \brief getElement
     * \param id
     * \return
     */
    Element *getElement(const std::string &id);

    /*!
     * \brief getElement
     * \param index
     * \return
     */
    Element *getElement(int index);

    /*!
     * \brief discretizationFile
     * \return
     */
    QFileInfo discretizationFile() const;

    /*!
     * \brief setDiscretizationFile
     * \param discretizationFile
     */
    void setDiscretizationFile(const QFileInfo &discretizationFile);

    /*!
     * \brief hydrodynamicFile
     * \return
     */
    QFileInfo hydrodynamicFile() const;

    /*!
     * \brief setHydrodynamicFile
     * \param hydrodynamicFile
     */
    void setHydrodynamicFile(const QFileInfo &hydrodynamicFile);

    /*!
     * \brief outputCSVFile
     * \return
     */
    QFileInfo outputCSVFile() const;

    /*!
     * \brief setOutputCSVFile
     * \param outputFile
     */
    void setOutputCSVFile(const QFileInfo &outputFile);

    /*!
     * \brief initialize
     * \param errors
     * \return
     */
    bool initialize(std::list<std::string> &errors);

    /*!
     * \brief update
     */
    void update();

    /*!
     * \brief finalize
     * \param errors
     * \return
     */
    bool finalize(std::list<std::string> &errors);

  private:

    /*!
     * \brief initializeTimeVariables
     * \param errors
     * \return
     */
    bool initializeTimeVariables(std::list<std::string> &errors);

    /*!
     * \brief initializeElements
     * \param errors
     * \return
     */
    bool initializeElements(std::list<std::string> &errors);

    /*!
     * \brief initializeSolver
     * \param errors
     * \return
     */
    bool initializeSolver(std::list<std::string> & errors);

    /*!
     * \brief initializeInputFiles
     * \param errors
     * \return
     */
    bool initializeInputFiles(std::list<std::string> &errors);

    /*!
     * \brief intializeOutputFiles
     * \param errors
     * \return
     */
    bool initializeOutputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeCSVOutputFile
     * \param errors
     * \return
     */
    bool initializeCSVOutputFile(std::list<std::string> &errors);

    /*!
     * \brief initializeNetCDFOutputFile
     * \param errors
     * \return
     */
    bool initializeNetCDFOutputFile(std::list<std::string> &errors);

    /*!
     * \brief prepareForNextTimeStep
     */
    void prepareForNextTimeStep();

    /*!
     * \brief applyInitialConditions
     */
    void applyInitialConditions();

    /*!
     * \brief applyBoundaryConditions
     */
    void applyBoundaryConditions(double dateTime);

    /*!
     * \brief computeTimeStep
     * \return
     */
    double computeTimeStep();


    /*!
     * \brief solveHeat
     * \param timeStep
     */
    void solveHeatTransport(double timeStep);

    /*!
     * \brief solveSoluteContinuity
     * \param soluteIndex
     * \param timeStep
     */
    void solveSoluteTransport(int soluteIndex, double timeStep);

    /*!
     * \brief computeDTDt
     * \param model
     * \param variableIndex
     * \param t
     * \param y
     * \param dydt
     */
    static void computeDTDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief computeSoluteDYDt
     * \param t
     * \param y
     * \param dydt
     * \param userData
     */
    static void computeDSoluteDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief writeOutput
     */
    void writeOutput();

    /*!
     * \brief writeCSVOutput
     */
    void writeCSVOutput();

    /*!
     * \brief writeNetCDFOutput
     */
    void writeNetCDFOutput();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputFiles();

    /*!
     * \brief closeCSVOutputFile
     */
    void closeCSVOutputFile();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputNetCDFFile();

  private:

    //Solute variable names;
    std::vector<std::string> m_solutes; // Names of the solutes.

    //Time variables
    double m_timeStep, //seconds
    m_startDateTime, //MJD
    m_endDateTime, //MJD
    m_currentDateTime, //MJD
    m_maxTimeStep, //seconds
    m_minTimeStep, //seconds
    m_outputInterval, //seconds
    m_nextOutputTime,//MJD
    m_timeStepRelaxationFactor; //

    //Number of initial fixed timeSteps of the minimum timestep to use when using the adaptive time step;
    int m_numInitFixedTimeSteps,
    m_numCurrentInitFixedTimeSteps;

    bool m_useAdaptiveTimeStep,
    m_converged;


    //Element junctions
    std::vector<ElementJunction*> m_elementJunctions;
    std::unordered_map<std::string, ElementJunction*> m_elementJunctionsById; //added for fast lookup using identifiers instead of indexes.

    //1D Computational elements
    std::vector<Element*> m_elements;
    std::unordered_map<std::string, Element*> m_elementsById; //added for fast lookup using identifiers instead of indexes.

    //Solver Object
    ODESolver *m_solver = nullptr;

    //Global water properties
    double m_waterDensity, //kg/m^3
    m_cp;// 4187.0; // J/kg/C

    //File input and output
    QFileInfo m_discretizationFile, //Geometric discretization file specified using the SWMM version file input file format
    m_hydrodynamicFile, //File to specify hydrodynamic parameters, initial conditions, boundary conditions,
    m_outputCSVFileInfo; //Output CSV filepath

    QTextStream m_outputCSVStream; //Output CSV filestream

    //Parent component
    STSComponent *m_component;

};

#endif // STSMODEL_H