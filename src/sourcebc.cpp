/*!
*  \file    SourceBC.cpp
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
*  \todo Test transport on branching networks
*  \warning
*/

#include "sourcebc.h"
#include "element.h"
#include "stsmodel.h"
#include "temporal/timeseries.h"
#include "core/datacursor.h"

SourceBC::SourceBC(Element *startElement, double startElementLFactor,
                   Element *endElement, double endElementLFactor,
                   int variableIndex, STSModel *model)
  : m_model(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_startElementLFactor(startElementLFactor),
    m_endElementLFactor(endElementLFactor),
    m_variableIndex(variableIndex)
{
  m_dataCursor = new DataCursor();
}

SourceBC::~SourceBC()
{
  delete m_dataCursor;
}

void SourceBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);

  m_factors.clear();

  for(Element *element : m_profile)
  {
    if(element != m_startElement && element != m_endElement)
      m_factors[element] = element->length;
  }

  m_factors[m_startElement] = m_startElement->length * m_startElementLFactor;
  m_factors[m_endElement]   = m_endElement->length * m_endElementLFactor;
}

void SourceBC::prepare()
{

}

void SourceBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;


  if(m_timeSeries->interpolate(dateTime, 0, value))
  {
    switch (m_variableIndex)
    {
      case -1:
        for(Element *element : m_profile)
        {
          if(element == m_startElement)
          {
            element->externalHeatFluxes += value * element->length * m_endElementLFactor;

          }
          else if(element == m_endElement)
          {
            element->externalHeatFluxes += value * element->length * m_endElementLFactor;
          }
          else
          {
            element->externalHeatFluxes += value * element->length;
          }
        }
        break;
      default:
        for(Element *element : m_profile)
        {
          if(element == m_startElement)
          {
            element->externalSoluteFluxes[m_variableIndex] += value * element->length * m_endElementLFactor;

          }
          else if(element == m_endElement)
          {
            element->externalSoluteFluxes[m_variableIndex] += value * element->length * m_endElementLFactor;
          }
          else
          {
            element->externalSoluteFluxes[m_variableIndex] += value * element->length;
          }

          element->externalSoluteFluxes[m_variableIndex] += value * element->length;
        }
        break;
    }
  }

  if(m_timeSeries->numColumns() == static_cast<int>(m_profile.size()))
  {
    switch (m_variableIndex)
    {
      case -1:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, static_cast<int>(i), m_dataCursor, value))
          {
            Element *element = m_profile[i];
            element->externalHeatFluxes += value * m_factors[element];
          }
        }
        break;
      default:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, static_cast<int>(i), m_dataCursor, value))
          {
            Element *element = m_profile[i];
            double factor = m_factors[element];
            element->externalSoluteFluxes[m_variableIndex] += value * factor;
          }
        }
        break;
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      switch (m_variableIndex)
      {
        case -1:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];
            element->externalHeatFluxes += value * m_factors[element];
          }
          break;
        default:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];
            element->externalSoluteFluxes[m_variableIndex] += value * m_factors[element];
          }
          break;
      }
    }
  }
}

void SourceBC::clear()
{

}

Element *SourceBC::startElement() const
{
  return m_startElement;
}

void SourceBC::setStartElement(Element *element)
{
  m_startElement = element;
}

double SourceBC::startElementLFactor() const
{
  return m_startElementLFactor;
}

void SourceBC::setStartElementLFactor(double factor)
{
  m_startElementLFactor = factor;
}

Element *SourceBC::endElement() const
{
  return m_endElement;
}

void SourceBC::setEndElement(Element *element)
{
  m_endElement = element;
}

double SourceBC::endElementLFactor() const
{
  return m_endElementLFactor;
}

void SourceBC::setEndElementLFactor(double factor)
{
  m_endElementLFactor = factor;
}

QSharedPointer<TimeSeries> SourceBC::timeSeries() const
{
  return m_timeSeries;
}

void SourceBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
