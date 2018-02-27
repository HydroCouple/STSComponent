#include "stsmodel.h"
#include "element.h"
#include "elementjunction.h"

#include <QDir>

using namespace std;

QFileInfo STSModel::discretizationFile() const
{
  return m_discretizationFile;
}

void STSModel::setDiscretizationFile(const QFileInfo &discretizationFile)
{
  m_discretizationFile = discretizationFile;
}

QFileInfo STSModel::hydrodynamicFile() const
{
  return m_hydrodynamicFile;
}

void STSModel::setHydrodynamicFile(const QFileInfo &hydrodynamicFile)
{
  m_hydrodynamicFile = hydrodynamicFile;
}

QFileInfo STSModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void STSModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

bool STSModel::initializeInputFiles(list<string> &errors)
{
  return true;
}

bool STSModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool STSModel::initializeCSVOutputFile(list<string> &errors)
{
  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if(!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.dir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if(!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice( device);
    }

    if(m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Temperature";

      for(size_t i = 0; i < m_solutes.size(); i++)
      {
        m_outputCSVStream << ", " << QString::fromStdString(m_solutes[i]);
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
  return true;
}

void STSModel::writeOutput()
{
  writeCSVOutput();
  writeNetCDFOutput();
}

void STSModel::writeCSVOutput()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    for(size_t i = 0 ; i < m_elements.size() ; i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        <<  ", " << element->temperature.value;

      for(size_t j = 0 ; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value;
      }

      m_outputCSVStream << endl;
    }
  }
}

void STSModel::writeNetCDFOutput()
{

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
  if(m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void STSModel::closeOutputNetCDFFile()
{

}
