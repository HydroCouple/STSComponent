#include "stsmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "radiativefluxbc.h"
#include "hydraulicsbc.h"
#include "mainchannelbc.h"
#include "sourcebc.h"
#include "meteorologybc.h"
#include "mainchannelbc.h"
#include "temporal/timedata.h"
#include "temporal/timeseries.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencvar.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"

#include <QDir>
#include <QDate>
#include <cstdlib>
#include <errno.h>

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;

#endif

using namespace std;

bool STSModel::verbose() const
{
  return m_verbose;
}

void STSModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int STSModel::printFrequency() const
{
  return m_printFrequency;
}

void STSModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int STSModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void STSModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo STSModel::inputFile() const
{
  return m_inputFile;
}

void STSModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo STSModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void STSModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

QFileInfo STSModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void STSModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void STSModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("TimeStep (s): %f\tDateTime: %f\tTemp (°C) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalHeatBalance: %g (KJ)}", m_timeStep, m_currentDateTime,
           m_odeSolver->getIterations(), m_odeSolver->maxIterations(), m_minTemp, m_maxTemp, m_totalHeatBalance);

    printf("\n");

    m_currentPrintCount = 0;
  }
}

bool STSModel::initializeInputFiles(list<string> &errors)
{

  if (QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if (file.open(QIODevice::ReadOnly))
    {
      m_delimiters = QRegExp("(\\,|\\t|\\;|\\s+)");
      int currentFlag = -1;
      m_addedSoluteCount = 0;

      QTextStream streamReader(&file);
      int lineCount = 0;
      while (!streamReader.atEnd())
      {
        QString line = streamReader.readLine().trimmed();
        lineCount++;

        if (!line.isEmpty() && !line.isNull())
        {
          bool readSuccess = true;
          QString error = "";

          auto it = m_inputFileFlags.find(line.toStdString());

          if (it != m_inputFileFlags.cend())
          {
            currentFlag = it->second;
          }
          else if (!QStringRef::compare(QStringRef(&line, 0, 2), ";;"))
          {
            //commment do nothing
          }
          else
          {
            switch (currentFlag)
            {
              case 1:
                readSuccess = readInputFileOptionTag(line, error);
                break;
              case 2:
                readSuccess = readInputFileOutputTag(line, error);
                break;
              case 3:
                readSuccess = readInputFileSolutesTag(line, error);
                break;
              case 4:
                readSuccess = readInputFileElementJunctionsTag(line, error);
                break;
              case 5:
                readSuccess = readInputFileElementsTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileSourcesTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileHydraulicsTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileMeteorologyTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileMCBoundaryConditionTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileTimeSeriesTag(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
            break;
          }
        }
      }

      file.close();
    }
  }

  return true;
}

bool STSModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool STSModel::initializeCSVOutputFile(list<string> &errors)
{

  if (m_outputCSVFileInfo.isRelative())
  {
    m_outputCSVFileInfo = relativePathToAbsolute(m_outputCSVFileInfo);
  }

  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.absoluteDir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if (!file.isEmpty() && !file.isNull())
  {
    if (m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice(device);
    }

    if (m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Depth, Width, XSectionArea, "
                           "Temperature, TotalExternalHeatFluxesBalance, TotalRadFluxesHeatBalance, "
                           "TotalEvapFluxHeatBalance, TotalConvFluxHeatBalance, TotalHeatBalance, ExtRadiation, EvapRate, ConvRate, BackRadRate";

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        QString soluteName = QString::fromStdString(m_solutes[i]);
        m_outputCSVStream << ", " << soluteName << ", "
                          << "TotalAdvDispMassBalance_" + soluteName << ", "
                          << "TotalExternalFluxMassBalance_" + soluteName << ","
                          << "TotalMassBalance_" + soluteName ;
      }

      m_outputCSVStream << endl;
    }

    m_outputCSVStream.flush();

    return true;
  }

  return false;
}

bool STSModel::initializeNetCDFOutputFile(list<string> &errors)
{

#ifdef USE_NETCDF

  if (m_outputNetCDFFileInfo.isRelative())
  {
    m_outputNetCDFFileInfo = relativePathToAbsolute(m_outputNetCDFFileInfo);
  }

  if (m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() || m_outputNetCDFFileInfo.isDir())
  {
    return true;
  }
  else if (!m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() &&
           !m_outputNetCDFFileInfo.absoluteFilePath().isNull() &&
           !m_outputNetCDFFileInfo.absoluteDir().exists())
  {
    std::string message = "NetCDF output file directory does not exist: " + m_outputNetCDFFileInfo.absoluteFilePath().toStdString();
    errors.push_back(message);
    return false;
  }

  bool returnValue = false;

  closeOutputNetCDFFile();

  try
  {

    m_outNetCDFVariables.clear();

    m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);

    //time variable
    ThreadSafeNcDim timeDim = m_outputNetCDF->addDim("time");
    ThreadSafeNcVar timeVar = m_outputNetCDF->addVar("time", NcType::nc_DOUBLE, timeDim);
    timeVar.putAtt("time:long_name", "time");
    timeVar.putAtt("time:units", "days since 1858-11-17 0:0:0");
    timeVar.putAtt("time:calendar", "modified_julian");
    m_outNetCDFVariables["time"] = timeVar;

    //Add Solutes
    ThreadSafeNcDim solutesDim = m_outputNetCDF->addDim("solutes", m_solutes.size());


    //Add element junctions
    ThreadSafeNcDim junctionDim = m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers = m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("element_junction_id::long_name", "element junction identifier");

    ThreadSafeNcVar junctionX = m_outputNetCDF->addVar("x", NcType::nc_FLOAT, junctionDim);
    junctionX.putAtt("x:long_name", "junction x-coordinate");
    junctionX.putAtt("x:units", "m");

    ThreadSafeNcVar junctionY = m_outputNetCDF->addVar("y", NcType::nc_FLOAT, junctionDim);
    junctionY.putAtt("y:long_name", "junction y-coordinate");
    junctionY.putAtt("y:units", "m");

    ThreadSafeNcVar junctionZ = m_outputNetCDF->addVar("z", NcType::nc_FLOAT, junctionDim);
    junctionZ.putAtt("z:long_name", "junction z-coordinate");
    junctionZ.putAtt("z:units", "m");

    float *vertx = new float[m_elementJunctions.size()];
    float *verty = new float[m_elementJunctions.size()];
    float *vertz = new float[m_elementJunctions.size()];
    char **junctionIds = new char *[m_elementJunctions.size()];

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      ElementJunction *junction = m_elementJunctions[i];

      junctionIds[i] = new char[junction->id.size() + 1];
      std::strcpy(junctionIds[i], junction->id.c_str());

      vertx[i] = junction->x;
      verty[i] = junction->y;
      vertz[i] = junction->z;
    }

    junctionX.putVar(vertx);
    junctionY.putVar(verty);
    junctionZ.putVar(vertz);
    junctionIdentifiers.putVar(junctionIds);

    delete[] vertx;
    delete[] verty;
    delete[] vertz;

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      delete[] junctionIds[i];
    }

    delete[] junctionIds;

    //Add Elements
    ThreadSafeNcDim elementsDim = m_outputNetCDF->addDim("elements", m_elements.size());

    ThreadSafeNcVar elementIdentifiers = m_outputNetCDF->addVar("element_id", NcType::nc_STRING, elementsDim);
    elementIdentifiers.putAtt("element_id::long_name", "element identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;


    ThreadSafeNcVar elementFromJunction = m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("from_junction:long_name", "upstream junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction = m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("to_junction:long_name", "downstream junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;

    ThreadSafeNcVar element_x =  m_outputNetCDF->addVar("element_x", NcType::nc_FLOAT, elementsDim);
    element_x.putAtt("standard_name", "projection_x_coordinate");
    element_x.putAtt("long_name", "Element X Coordinate");
    element_x.putAtt("units", "m");
    m_outNetCDFVariables["element_x"] = element_x;

    ThreadSafeNcVar element_y =  m_outputNetCDF->addVar("element_y", NcType::nc_FLOAT, elementsDim);
    element_y.putAtt("standard_name", "projection_y_coordinate");
    element_y.putAtt("long_name", "Element Y Coordinate");
    element_y.putAtt("units", "m");
    m_outNetCDFVariables["element_y"] = element_y;

    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    float *elX = new float[m_elements.size()];
    float *elY = new float[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      std::strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->index;
      toJunctions[i] = element->downstreamJunction->index;

      elX[i] = element->x;
      elY[i] = element->y;
    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    element_x.putVar(elX);
    element_y.putVar(elY);

    delete[] fromJunctions;
    delete[] toJunctions;
    delete[] elX;
    delete[] elY;

    for (size_t i = 0; i < m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;


    //hydraulics variables
    ThreadSafeNcVar lengthVar =  m_outputNetCDF->addVar("length", "float",
                                                        std::vector<std::string>({"elements"}));
    lengthVar.putAtt("long_name", "Element Length");
    lengthVar.putAtt("units", "m");
    m_outNetCDFVariables["length"] = lengthVar;

    std::vector<float> lengths(m_elements.size());

    for(size_t i = 0; i < m_elements.size(); i++)
    {
      lengths[i] = static_cast<float>(m_elements[i]->length);
    }

    lengthVar.putVar(lengths.data());

    auto varOnOff = [this](const std::string& name) -> bool
    {
      return m_outNetCDFVariablesOnOff.find(name) != m_outNetCDFVariablesOnOff.end() ? m_outNetCDFVariablesOnOff[name] : true;
    };


    if((m_outNetCDFVariablesOnOff["depth"] = varOnOff("depth")))
    {
      ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      depthVar.putAtt("long_name", "Flow Depth");
      depthVar.putAtt("units", "m");
      m_outNetCDFVariables["depth"] = depthVar;
      m_outNetCDFVariablesIOFunctions["depth"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *depth = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          depth[i] = static_cast<float>(element->depth);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), depth);
        delete[] depth;
      };
    }

    if((m_outNetCDFVariablesOnOff["width"] = varOnOff("width")))
    {
      ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      widthVar.putAtt("long_name", "Top Width");
      widthVar.putAtt("units", "m");
      m_outNetCDFVariables["width"] = widthVar;
      m_outNetCDFVariablesIOFunctions["width"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *width = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          width[i] = static_cast<float>(element->width * element->beta);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), width);
        delete[] width;
      };
    }

    if((m_outNetCDFVariablesOnOff["xsection_area"] = varOnOff("xsection_area")))
    {
      ThreadSafeNcVar xsectAreaVar =  m_outputNetCDF->addVar("xsection_area", "float",
                                                             std::vector<std::string>({"time", "elements"}));
      xsectAreaVar.putAtt("long_name", "Flow Cross-Sectional Area");
      xsectAreaVar.putAtt("units", "m^2");
      m_outNetCDFVariables["xsection_area"] = xsectAreaVar;
      m_outNetCDFVariablesIOFunctions["xsection_area"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *xsection_area = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          xsection_area[i] = static_cast<float>(element->xSectionArea.value);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), xsection_area);
        delete[] xsection_area;
      };
    }

    if((m_outNetCDFVariablesOnOff["temperature"] = varOnOff("temperature")))
    {
      ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("temperature", "float",
                                                               std::vector<std::string>({"time", "elements"}));
      temperatureVar.putAtt("long_name", "Temperature");
      temperatureVar.putAtt("units", "°C");
      m_outNetCDFVariables["temperature"] = temperatureVar;
      m_outNetCDFVariablesIOFunctions["temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *temperature = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          temperature[i] = static_cast<float>(element->temperature.value);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), temperature);
        delete[] temperature;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_heat_balance"] = varOnOff("total_heat_balance")))
    {
      ThreadSafeNcVar totalHeatBalanceVar =  m_outputNetCDF->addVar("total_heat_balance", "float",
                                                                    std::vector<std::string>({"time"}));
      totalHeatBalanceVar.putAtt("long_name", "Total Heat Balance");
      totalHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_heat_balance"] = totalHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_evap_heat_balance"] = varOnOff("total_evap_heat_balance")))
    {
      ThreadSafeNcVar totalEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_evap_heat_balance", "float",
                                                                        std::vector<std::string>({"time"}));
      totalEvapHeatBalanceVar.putAtt("long_name", "Total Evaporation Heat Balance");
      totalEvapHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_evap_heat_balance"] = totalEvapHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_conv_heat_balance"] = varOnOff("total_conv_heat_balance")))
    {
      ThreadSafeNcVar totalConvHeatBalanceVar =  m_outputNetCDF->addVar("total_conv_heat_balance", "float",
                                                                        std::vector<std::string>({"time"}));
      totalConvHeatBalanceVar.putAtt("long_name", "Total Convection Heat Balance");
      totalConvHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_conv_heat_balance"] = totalConvHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["element_evap_heat_flux"] = varOnOff("element_evap_heat_flux")))
    {
      ThreadSafeNcVar elementEvapHeatFluxVar =  m_outputNetCDF->addVar("element_evap_heat_flux", "float",
                                                                       std::vector<std::string>({"time", "elements"}));

      elementEvapHeatFluxVar.putAtt("long_name", "Element Evaporation Heat Flux");
      elementEvapHeatFluxVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["element_evap_heat_flux"] = elementEvapHeatFluxVar;
      m_outNetCDFVariablesIOFunctions["element_evap_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *element_evap_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          element_evap_heat_flux[i] = static_cast<float>(element->evaporationHeatFlux);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), element_evap_heat_flux);
        delete[] element_evap_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["element_conv_heat_flux"] = varOnOff("element_conv_heat_flux")))
    {
      ThreadSafeNcVar elementConvHeatFluxVar =  m_outputNetCDF->addVar("element_conv_heat_flux", "float",
                                                                       std::vector<std::string>({"time", "elements"}));

      elementConvHeatFluxVar.putAtt("long_name", "Element Convective Heat Flux");
      elementConvHeatFluxVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["element_conv_heat_flux"] = elementConvHeatFluxVar;
      m_outNetCDFVariablesIOFunctions["element_conv_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *element_conv_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          element_conv_heat_flux[i] = static_cast<float>(element->convectionHeatFlux);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), element_conv_heat_flux);
        delete[] element_conv_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["element_mc_heat_flux"] = varOnOff("element_mc_heat_flux")))
    {
      ThreadSafeNcVar elementMCHeatFluxVar =  m_outputNetCDF->addVar("element_mc_heat_flux", "float",
                                                                       std::vector<std::string>({"time", "elements"}));

      elementMCHeatFluxVar.putAtt("long_name", "Element Main Channel Heat Exchange Flux");
      elementMCHeatFluxVar.putAtt("units", "J/s");
      m_outNetCDFVariables["element_mc_heat_flux"] = elementMCHeatFluxVar;
      m_outNetCDFVariablesIOFunctions["element_mc_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *element_conv_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          element_conv_heat_flux[i] = static_cast<float>(element->mcHeatExchangeFlux);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), element_conv_heat_flux);
        delete[] element_conv_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_radiation_flux_heat_balance"] = varOnOff("total_radiation_flux_heat_balance")))
    {
      ThreadSafeNcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_radiation_flux_heat_balance", "float",
                                                                                 std::vector<std::string>({"time"}));
      totalRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Radiation Flux Heat Balance");
      totalRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_radiation_flux_heat_balance"] = totalRadiationFluxHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_external_heat_flux_balance"] = varOnOff("total_external_heat_flux_balance")))
    {
      ThreadSafeNcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_external_heat_flux_balance", "float",
                                                                                std::vector<std::string>({"time"}));
      totalExternalHeatFluxBalanceVar.putAtt("long_name", "Total External Heat Flux Balance");
      totalExternalHeatFluxBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_external_heat_flux_balance"] = totalExternalHeatFluxBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_element_heat_balance"] = varOnOff("total_element_heat_balance")))
    {
      ThreadSafeNcVar totalElementHeatBalanceVar =  m_outputNetCDF->addVar("total_element_heat_balance", "float",
                                                                           std::vector<std::string>({"time", "elements"}));
      totalElementHeatBalanceVar.putAtt("long_name", "Total Element Heat Balance");
      totalElementHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_heat_balance"] = totalElementHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_heat_balance[i] = static_cast<float>(element->totalHeatBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_heat_balance);
        delete[] total_element_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_evap_heat_balance"] = varOnOff("total_element_evap_heat_balance")))
    {
      ThreadSafeNcVar totalElementEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_element_evap_heat_balance", "float",
                                                                               std::vector<std::string>({"time", "elements"}));
      totalElementEvapHeatBalanceVar.putAtt("long_name", "Total Element Evaporation Heat Balance");
      totalElementEvapHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_evap_heat_balance"] = totalElementEvapHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_evap_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_evap_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_evap_heat_balance[i] = static_cast<float>(element->totalEvaporativeHeatFluxesBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_evap_heat_balance);
        delete[] total_element_evap_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_conv_heat_balance"] = varOnOff("total_element_conv_heat_balance")))
    {
      ThreadSafeNcVar totalElementConvHeatBalanceVar =  m_outputNetCDF->addVar("total_element_conv_heat_balance", "float",
                                                                               std::vector<std::string>({"time", "elements"}));
      totalElementConvHeatBalanceVar.putAtt("long_name", "Total Element Convection Heat Balance");
      totalElementConvHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_conv_heat_balance"] = totalElementConvHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_conv_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_conv_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_conv_heat_balance[i] = static_cast<float>(element->totalConvectiveHeatFluxesBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_conv_heat_balance);
        delete[] total_element_conv_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_radiation_flux_heat_balance"] = varOnOff("total_element_radiation_flux_heat_balance")))
    {
      ThreadSafeNcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_element_radiation_flux_heat_balance", "float",
                                                                                        std::vector<std::string>({"time", "elements"}));
      totalElementRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Element Radiation Flux Heat Balance");
      totalElementRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_radiation_flux_heat_balance"] = totalElementRadiationFluxHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_radiation_flux_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_radiation_flux_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_radiation_flux_heat_balance[i] = static_cast<float>(element->totalRadiationFluxesHeatBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_radiation_flux_heat_balance);
        delete[] total_element_radiation_flux_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_external_heat_flux_balance"] = varOnOff("total_element_external_heat_flux_balance")))
    {
      ThreadSafeNcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_element_external_heat_flux_balance", "float",
                                                                                       std::vector<std::string>({"time", "elements"}));
      totalElementExternalHeatFluxBalanceVar.putAtt("long_name", "Total Element External Heat Flux Balance");
      totalElementExternalHeatFluxBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_external_heat_flux_balance"] = totalElementExternalHeatFluxBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_external_heat_flux_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_external_heat_flux_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_external_heat_flux_balance[i] = static_cast<float>(element->totalExternalHeatFluxesBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_external_heat_flux_balance);
        delete[] total_element_external_heat_flux_balance;
      };
    }

    if(m_solutes.size())
    {
      ThreadSafeNcVar solutes = m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
      solutes.putAtt("solute_names::long_name", "solute names");
      m_outNetCDFVariables["solutes"] = solutes;

      char **soluteNames = new char *[m_solutes.size()];

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        std::strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        delete[] soluteNames[i];
      }

      delete[] soluteNames;

      if((m_outNetCDFVariablesOnOff["solute_concentration"] = varOnOff("solute_concentration")))
      {
        ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "float",
                                                             std::vector<std::string>({"time", "solutes", "elements"}));
        solutesVar.putAtt("long_name", "Solute Concentration");
        solutesVar.putAtt("units", "kg/m^3");
        m_outNetCDFVariables["solute_concentration"] = solutesVar;
        m_outNetCDFVariablesIOFunctions["solute_concentration"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
        {
          if(elements.size())
          {
            int numSolutes = elements[0]->model->numSolutes();
            float *solute_concentration = new float[elements.size() * numSolutes];

            for (size_t i = 0; i < elements.size(); i++)
            {
              Element *element = elements[i];

              for (int j = 0; j < numSolutes; j++)
              {
                solute_concentration[i + j * elements.size()] = static_cast<float>(element->soluteConcs[j].value);
              }
            }
            variable.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)numSolutes, elements.size()}), solute_concentration);
            delete[] solute_concentration;
          }
        };
      }

      if((m_outNetCDFVariablesOnOff["total_solute_mass_balance"] = varOnOff("total_solute_mass_balance")))
      {
        ThreadSafeNcVar totalSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_solute_mass_balance", "float",
                                                                            std::vector<std::string>({"time", "solutes"}));
        totalSoluteMassBalanceVar.putAtt("long_name", "Total Solute Mass Balance");
        totalSoluteMassBalanceVar.putAtt("units", "kg");
        m_outNetCDFVariables["total_solute_mass_balance"] = totalSoluteMassBalanceVar;
      }

      if((m_outNetCDFVariablesOnOff["total_external_solute_flux_mass_balance"] = varOnOff("total_external_solute_flux_mass_balance")))
      {
        ThreadSafeNcVar totalExternalSoluteFluxMassBalanceVar = m_outputNetCDF->addVar("total_external_solute_flux_mass_balance", "float",
                                                                                       std::vector<std::string>({"time", "solutes"}));
        totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
        totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:units", "kg");
        m_outNetCDFVariables["total_external_solute_flux_mass_balance"] = totalExternalSoluteFluxMassBalanceVar;
      }

      if((m_outNetCDFVariablesOnOff["total_element_solute_mass_balance"] = varOnOff("total_element_solute_mass_balance")))
      {
        ThreadSafeNcVar totalElementSoluteMassBalanceVar = m_outputNetCDF->addVar("total_element_solute_mass_balance", "float",
                                                                                  std::vector<std::string>({"time", "solutes", "elements"}));
        totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:long_name", "total element solute mass balance");
        totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:units", "kg");
        m_outNetCDFVariables["total_element_solute_mass_balance"] = totalElementSoluteMassBalanceVar;
      }

      if((m_outNetCDFVariablesOnOff["total_element_external_solute_flux_mass_balance"] = varOnOff("total_element_external_solute_flux_mass_balance")))
      {
        ThreadSafeNcVar totalElementExternalSoluteFluxMassBalanceVar = m_outputNetCDF->addVar("total_element_external_solute_flux_mass_balance", "float",
                                                                                              std::vector<std::string>({"time", "solutes"}));
        totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
        totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:units", "kg");
        m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"] = totalElementExternalSoluteFluxMassBalanceVar;
      }
    }

    m_optionalOutputVariables.clear();
    m_optionalOutputVariables.reserve(m_outNetCDFVariablesOnOff.size());

    for (const auto& pair : m_outNetCDFVariablesOnOff)
    {
      if(pair.second && (m_outNetCDFVariablesIOFunctions.find(pair.first) != m_outNetCDFVariablesIOFunctions.end()))
        m_optionalOutputVariables.push_back(pair.first);
    }


    m_outputNetCDF->sync();

    return true;
  }
  catch (NcException &e)
  {
    std::string message = std::string(e.what());
    printf("%s\n", e.what());
    errors.push_back(message);
    return false;
  }

#endif

  return false;
}

bool STSModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  std::string optionsFlag = options[0].toStdString();
  auto it = m_optionsFlags.find(optionsFlag);

  if (it != m_optionsFlags.end())
  {
    int optionsIndex = it->second;

    switch (optionsIndex)
    {
      case 1:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            QString dateTimeString = options[1] + " " + options[2];
            if (SDKTemporal::DateTime::tryParse(dateTimeString, dateTime))
            {
              m_startDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
              m_currentDateTime = m_startDateTime;
            }
            else
            {
              foundError = true;
            }
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 2:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_endDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 3:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_outputInterval = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Report interval error";
            return false;
          }
        }
        break;
      case 4:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_maxTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Max time step error";
            return false;
          }
        }
        break;
      case 5:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_minTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Min time step error";
            return false;
          }
        }
        break;
      case 6:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numInitFixedTimeSteps = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of initial time step error";
            return false;
          }
        }
        break;
      case 7:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useAdaptiveTimeStep = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use adaptive time step error";
            return false;
          }
        }
        break;
      case 8:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_timeStepRelaxationFactor = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Time step relaxation factor error";
            return false;
          }
        }
        break;
      case 9:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int heatSolverMode = -1;

            if (it != m_solverTypeFlags.end())
              heatSolverMode = it->second;

            switch (heatSolverMode)
            {
              case 1:
                m_odeSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_odeSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_ADAMS);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::FUNCTIONAL);
                }
                break;
              case 4:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_BDF);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::NEWTON);
                  m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                }
                break;
              case 5:
                m_odeSolver->setSolverType(ODESolver::EULER);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver type error";
            return false;
          }
        }
        break;
      case 10:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver absolute tolerance error";
            return false;
          }
        }
        break;
      case 11:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver relative tolerance error";
            return false;
          }
        }
        break;
      case 12:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_linearSolverTypeFlags.find(code);

            int linearSolver = -1;

            if (it != m_linearSolverTypeFlags.end())
              linearSolver = it->second;

            switch (linearSolver)
            {
              case 1:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                break;
              case 2:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::FGMRES);
                break;
              case 3:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::Bi_CGStab);
                break;
              case 4:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::TFQMR);
                break;
              case 5:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::PCG);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Linear solver type error";
            return false;
          }
        }
        break;
      case 13:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_waterDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water density error";
            return false;
          }
        }
        break;
      case 14:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_cp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Specific heat capacity of water error";
            return false;
          }
        }
        break;
      case 15:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            int numSolutes = options[1].toInt(&parsed);

            if (parsed)
              setNumSolutes(numSolutes);

            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of solutes error";
            return false;
          }
        }
        break;
      case 16:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Verbose error";
            return false;
          }
        }
        break;
      case 17:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_flushToDiskFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Flush to disk frequency error";
            return false;
          }
        }
        break;
      case 18:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_printFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Print frequency error";
            return false;
          }
        }
        break;
      case 19:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useEvaporation = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use evaporation error";
            return false;
          }
        }
        break;
      case 20:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useConvection = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use convection error";
            return false;
          }
        }
        break;
      case 21:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_evapWindFuncCoeffA = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Evaporation coefficient error";
            return false;
          }
        }
        break;
      case 22:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_evapWindFuncCoeffB = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Evaporation coefficient error";
            return false;
          }
        }
        break;
      case 23:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_bowensCoeff = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Bowens coefficient error";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool STSModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  QString optionsFlag = options[0];

  if (options.size() == 2)
  {
    if (!QString::compare(optionsFlag, "csv", Qt::CaseInsensitive))
    {
      m_outputCSVFileInfo = QFileInfo(options[1]);
    }
    else if (!QString::compare(optionsFlag, "netcdf", Qt::CaseInsensitive))
    {
      m_outputNetCDFFileInfo = QFileInfo(options[1]);
    }
  }
  else
  {
    errorMessage = "Output file error";
    return false;
  }

  return true;
}

bool STSModel::readInputFileSolutesTag(const QString &line, QString &errorMessage)
{
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();
      m_addedSoluteCount++;
    }
    else
    {
      errorMessage = "Solute error count";
      return false;
    }
  }
  else
  {
    errorMessage = "Solute error";
    return false;
  }

  return true;
}

bool STSModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];

    bool workedX;
    bool workedY;
    bool workedZ;

    double x = columns[1].toDouble(&workedX);

    double y = columns[2].toDouble(&workedY);

    double z = columns[3].toDouble(&workedZ);

    if (workedX && workedY && workedZ)
    {
      addElementJunction(id.toStdString(), x, y, z);
    }
    else
    {
      errorMessage = "Junctions error";
      return false;
    }
  }
  else
  {
    errorMessage = "Junctions error";
    return false;
  }

  return true;
}

bool STSModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() > 10)
  {
    QString id = columns[0];
    QString fromId = columns[1];
    QString toId = columns[2];

    auto fromIt = m_elementJunctionsById.find(fromId.toStdString());
    auto toIt = m_elementJunctionsById.find(toId.toStdString());

    if (fromIt != m_elementJunctionsById.end() &&
        toIt != m_elementJunctionsById.end())
    {
      ElementJunction *ej1 = m_elementJunctionsById[fromId.toStdString()];
      ElementJunction *ej2 = m_elementJunctionsById[toId.toStdString()];

      Element::Bank bank = Element::Left;

      //      if(!QString::compare(columns[3],"Right", Qt::CaseInsensitive))
      //      {
      //        bank = Element::Right;
      //      }

      bool lengthOk;
      double length = columns[3].toDouble(&lengthOk);

      bool depthOk ;
      double depth = columns[4].toDouble(&depthOk);

      bool xsectionAreaOk ;
      double xsectionArea = columns[5].toDouble(&xsectionAreaOk);

      bool widthOk ;
      double width = columns[6].toDouble(&widthOk);

      bool betaOk ;
      double beta = columns[7].toDouble(&betaOk);

      bool tempOk ;
      double temp = columns[8].toDouble(&tempOk);

      bool mcTempOk ;
      double mcTemp = columns[9].toDouble(&mcTempOk);

      bool tempExchCoeffOk ;
      double tempExchCoeff = columns[10].toDouble(&tempExchCoeffOk);

      if (lengthOk && depthOk && xsectionAreaOk &&
          widthOk && betaOk && tempOk &&
          mcTempOk && tempExchCoeffOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->bank = bank;
        element->length = length;
        element->depth = depth;
        element->xSectionArea.value = xsectionArea;
        element->width = width;
        element->temperature.value = temp;
        element->mainChannelTemperature = mcTemp;
        element->temperatureExchangeCoefficient = tempExchCoeff;
        element->beta = beta;

        if (columns.size() > 10)
        {
          for (int i = 11; i < columns.size(); i = i+3)
          {
            bool soluteOk ;
            double solute = columns[i].toDouble(&soluteOk);

            bool mainChannelSoluteOk;
            double mainChannelSolute = columns[i+1].toDouble(&mainChannelSoluteOk);

            bool soluteExchangeCoeffOk;
            double soluteExchangeCoeff = columns[i+2].toDouble(&soluteExchangeCoeffOk);


            if (soluteOk && soluteExchangeCoeffOk &&
                mainChannelSoluteOk && i - 11 < (int)m_solutes.size())
            {
              element->soluteConcs[i - 11].value = solute;
              element->soluteExchangeCoefficients[i-11] = soluteExchangeCoeff;
              element->mainChannelSoluteConcs[i-11] = mainChannelSolute;
            }
            else
            {
              errorMessage = "Wrong initial solute parameters";
              return false;
            }
          }
        }
      }
      else
      {
        errorMessage = "";
        return false;
      }
    }
    else
    {
      errorMessage = "Wrong upstream or downstream junction";
      return false;
    }
  }

  return true;
}

bool STSModel::readInputFileSourcesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 7)
  {
    QString idFrom = columns[0];
    QString idTo = columns[2];

    bool okStart;
    double startFactor = columns[1].toDouble(&okStart);

    bool okEnd;
    double endFactor = columns[3].toDouble(&okEnd);

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if (okStart && okEnd && itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *elementFrom = itFrom->second;
      Element *elementTo = itTo->second;

      QString variable = columns[4].trimmed();
      QString type = columns[5];
      int variableType;

      int soluteIndex = -1;

      if (!QString::compare(variable, "HEAT", Qt::CaseInsensitive))
      {
        variableType = -1;
        soluteIndex = 0;
      }
      else
      {

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            variableType = i;
            soluteIndex = i;
            break;
          }
        }
      }

      if(soluteIndex > -1)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[6].toDouble(&valueOk);

          if (valueOk)
          {
            SourceBC *sourceBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            sourceBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = columns[6].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            SourceBC *sourceBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            sourceBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Source is is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Source is is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STSModel::readInputFileHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString idFrom = columns[0];
    QString idTo = columns[1];

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if ( itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      auto it = m_hydraulicVariableFlags.find(columns[2].toStdString());

      if (it != m_hydraulicVariableFlags.end())
      {
        int variableIndex = it->second;

        QString type = columns[3];

        if(type.toLower() == "VALUE")
        {
          bool valueOk;
          double value = columns[4].toDouble(&valueOk);

          if(valueOk)
          {
            HydraulicsBC *hydraulicsBC= new HydraulicsBC(itFrom->second, itTo->second, variableIndex, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            hydraulicsBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Value is invalid";
            return false;
          }
        }
        else if (type.toUpper() == "TIMESERIES")
        {
          std::string tsId = columns[4].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            HydraulicsBC *hydraulicsBC= new HydraulicsBC(itFrom->second, itTo->second, variableIndex, this);
            hydraulicsBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Timeseries was not found";
            return false;
          }
        }
        else
        {
          errorMessage = "Wrong timeseries or value type";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Could not find end/start cell";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

bool STSModel::readInputFileRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString idFrom = columns[0];
    QString idTo = columns[1];

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if ( itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {

      QString type = columns[2];

      if(type.toLower() == "VALUE")
      {
        bool valueOk;
        double value = columns[3].toDouble(&valueOk);

        if(valueOk)
        {
          RadiativeFluxBC *hydraulicsBC= new RadiativeFluxBC(itFrom->second, itTo->second, this);

          QUuid uid = QUuid::createUuid();
          QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

          ts->addRow(m_startDateTime, value);
          ts->addRow(m_endDateTime, value);

          m_timeSeries[ts->id().toStdString()] = ts;

          hydraulicsBC->setTimeSeries(ts);

          m_boundaryConditions.push_back(hydraulicsBC);
        }
        else
        {
          errorMessage = "Value is invalid";
          return false;
        }
      }
      else if (type.toUpper() == "TIMESERIES")
      {
        std::string tsId = columns[3].toStdString();
        auto tsIt = m_timeSeries.find(tsId);

        if (tsIt != m_timeSeries.end())
        {
          RadiativeFluxBC *hydraulicsBC= new RadiativeFluxBC(itFrom->second, itTo->second, this);
          hydraulicsBC->setTimeSeries(tsIt->second);
          m_boundaryConditions.push_back(hydraulicsBC);
        }
        else
        {
          errorMessage = "Timeseries was not found";
          return false;
        }
      }
      else
      {
        errorMessage = "Wrong timeseries or value type";
        return false;
      }

    }
    else
    {
      errorMessage = "Could not find end/start cell";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

bool STSModel::readInputFileMeteorologyTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_meteorologicalVariableFlags.find(variableType.toStdString());

      if (it != m_meteorologicalVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement,
                                                             variableIndex, this);
            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            meteorologyBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(meteorologyBC);
          }
          else
          {
            errorMessage = "Uniform meteorology value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if(tsIt != m_timeSeries.end())
          {
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement,
                                                             variableIndex, this);

            meteorologyBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(meteorologyBC);
          }
          else
          {
            errorMessage = "Specified uniform meteorology filepath does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform meteorology is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform meteorology boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STSModel::readInputFileMCBoundaryConditionTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString idFrom = columns[0];
    QString idTo = columns[1];

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *elementFrom = itFrom->second;
      Element *elementTo = itTo->second;

      QString variable = columns[2].trimmed();
      QString type = columns[3];
      int variableType;

      int soluteIndex = -1;

      if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
      {
        variableType = -1;
        soluteIndex = 0;
      }
      else
      {

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            variableType = i;
            soluteIndex = i;
            break;
          }
        }
      }

      if(soluteIndex > -1)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[4].toDouble(&valueOk);

          if (valueOk)
          {
            MainChannelBC *sourceBC = new MainChannelBC(elementFrom, elementTo, variableType, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            sourceBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = columns[6].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            MainChannelBC *sourceBC = new MainChannelBC(elementFrom, elementTo, variableType, this);
            sourceBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Source is is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Source is is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STSModel::readInputFileTimeSeriesTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    QFileInfo fileInfo(options[1].trimmed());

    if (fileInfo.isRelative())
      fileInfo = relativePathToAbsolute(fileInfo);

    if(QFile::exists(fileInfo.absoluteFilePath()))
    {
      QSharedPointer<TimeSeries> timeSeries(TimeSeries::createTimeSeries(options[0], fileInfo, this));

      if(!timeSeries.isNull())
      {
        m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
      }
      else
      {
        errorMessage = "Timeseries specified is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified filepath does not exist";
      return false;
    }
  }
  else
  {
    errorMessage = "TimeSeries must have two columns";
    return false;
  }

  return true;
}

void STSModel::writeOutput()
{
  m_currentflushToDiskCount++;

  if (m_currentflushToDiskCount >= m_flushToDiskFrequency)
  {
    m_flushToDisk = true;
    m_currentflushToDiskCount = 0;
  }
  else
  {
    m_flushToDisk = false;
  }

  writeCSVOutput();
  writeNetCDFOutput();
}

void STSModel::writeCSVOutput()
{
  if (m_outputCSVStream.device()->isOpen())
  {
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        << ", " << element->depth
                        << ", " << element->width
                        << ", " << element->xSectionArea.value
                        << ", " << element->temperature.value
                        << ", " << element->totalExternalHeatFluxesBalance
                        << ", " << element->totalRadiationFluxesHeatBalance
                        << ", " << element->totalEvaporativeHeatFluxesBalance
                        << ", " << element->totalConvectiveHeatFluxesBalance
                        << ", " << element->totalHeatBalance;


      m_outputCSVStream << ", " << element->radiationFluxes
                        << ", " << element->evaporationHeatFlux
                        << ", " << element->convectionHeatFlux;

      for (size_t j = 0; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value
                          << ", " << element->totalExternalSoluteFluxesMassBalance[j]
                             << ", " << element->totalSoluteMassBalance[j];
      }

      m_outputCSVStream << endl;
    }

    if (m_flushToDisk)
    {
      m_outputCSVStream.flush();
    }
  }
}

void STSModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF
  if (m_outputNetCDF)
  {

    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();

    //Set current dateTime
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);

    int nVars = static_cast<int>(m_optionalOutputVariables.size());

#ifdef USE_OPENMP
//#pragma omp parallel for
#endif
    for (int i = 0; i < nVars; i++)
    {
      std::string varName = m_optionalOutputVariables[static_cast<size_t>(i)];
      (m_outNetCDFVariablesIOFunctions[varName])(currentTime, m_outNetCDFVariables[varName], m_elements);
    }


    if(m_outNetCDFVariablesOnOff["total_heat_balance"])
      m_outNetCDFVariables["total_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_evap_heat_balance"])
      m_outNetCDFVariables["total_evap_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalEvaporationHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_conv_heat_balance"])
      m_outNetCDFVariables["total_conv_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalConvectiveHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_radiation_flux_heat_balance"])
      m_outNetCDFVariables["total_radiation_flux_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_external_heat_flux_balance"])
      m_outNetCDFVariables["total_external_heat_flux_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);


    if(numSolutes())
    {
      if(m_outNetCDFVariablesOnOff["total_solute_mass_balance"])
        m_outNetCDFVariables["total_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)numSolutes()}), m_totalSoluteMassBalance.data());

      if(m_outNetCDFVariablesOnOff["total_external_solute_flux_mass_balance"])
        m_outNetCDFVariables["total_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)numSolutes()}), m_totalExternalSoluteFluxMassBalance.data());
    }

    if (m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }
  }
#endif
}

void STSModel::closeOutputFiles()
{
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeCSVOutputFile();
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeOutputNetCDFFile();
    }
  }
}

void STSModel::closeCSVOutputFile()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void STSModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF
  if(m_outputNetCDF)
  {
    m_outputNetCDF->sync();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }
#endif
}

QFileInfo STSModel::relativePathToAbsolute(const QFileInfo &fileInfo)
{
  if (fileInfo.isRelative())
  {
    if (!m_inputFile.filePath().isEmpty() &&
        !m_inputFile.filePath().isNull() &&
        QFile::exists(m_inputFile.absoluteFilePath()))
    {
      QFileInfo absoluteFilePath = m_inputFile.absoluteDir().absoluteFilePath(fileInfo.filePath());

      if (absoluteFilePath.absoluteDir().exists())
      {
        return absoluteFilePath;
      }
    }
  }

  return fileInfo;
}

const unordered_map<string, int> STSModel::m_inputFileFlags({
                                                              {"[OPTIONS]", 1},
                                                              {"[OUTPUTS]", 2},
                                                              {"[SOLUTES]", 3},
                                                              {"[ELEMENTJUNCTIONS]", 4},
                                                              {"[ELEMENTS]", 5},
                                                              {"[SOURCES]", 6},
                                                              {"[HYDRAULICS]", 7},
                                                              {"[RADIATIVE_FLUXES]", 8},
                                                              {"[METEOROLOGY]", 9},
                                                              {"[MC_BOUNDARY_CONDITION]", 10},
                                                              {"[TIMESERIES]", 11},
                                                            });

const unordered_map<string, int> STSModel::m_optionsFlags({
                                                            {"START_DATETIME", 1},
                                                            {"END_DATETIME", 2},
                                                            {"REPORT_INTERVAL", 3},
                                                            {"MAX_TIME_STEP", 4},
                                                            {"MIN_TIME_STEP", 5},
                                                            {"NUM_INITIAL_FIXED_STEPS", 6},
                                                            {"USE_ADAPTIVE_TIME_STEP", 7},
                                                            {"TIME_STEP_RELAXATION_FACTOR", 8},
                                                            {"SOLVER", 9},
                                                            {"SOLVER_ABS_TOL", 10},
                                                            {"SOLVER_REL_TOL", 11},
                                                            {"LINEAR_SOLVER", 12},
                                                            {"WATER_DENSITY", 13},
                                                            {"WATER_SPECIFIC_HEAT_CAPACITY", 14},
                                                            {"NUM_SOLUTES", 15},
                                                            {"VERBOSE", 16},
                                                            {"FLUSH_TO_DISK_FREQ", 17},
                                                            {"PRINT_FREQ", 18},
                                                            {"EVAPORATION", 19},
                                                            {"CONVECTION", 20},
                                                            {"EVAP_WIND_FUNC_COEFF_A", 21},
                                                            {"EVAP_WIND_FUNC_COEFF_A", 22},
                                                            {"BOWENS_COEFF", 23},
                                                          });

const unordered_map<string, int> STSModel::m_solverTypeFlags({{"RK4", 1},
                                                              {"RKQS", 2},
                                                              {"ADAMS", 3},
                                                              {"BDF", 4}});

const unordered_map<string, int> STSModel::m_linearSolverTypeFlags({{"GMRES", 1},
                                                                    {"FGMRES", 2},
                                                                    {"Bi_CGStab", 3},
                                                                    {"TFQMR", 4},
                                                                    {"PCG", 5}});


const unordered_map<string, int> STSModel::m_hydraulicVariableFlags({{"DEPTH", 1},
                                                                     {"WIDTH", 2},
                                                                     {"XSECTION_AREA", 3}});

const unordered_map<string, int> STSModel::m_meteorologicalVariableFlags({{"RELATIVE_HUMIDITY", 1},
                                                                          {"AIR_TEMPERATURE", 2},
                                                                          {"WIND_SPEED", 3}});

const QRegExp STSModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
