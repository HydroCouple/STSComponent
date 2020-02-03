/*!
*  \file    MainChannelBC.cpp
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

#include "mainchannelbc.h"
#include "element.h"
#include "stsmodel.h"
#include "temporal/timeseries.h"
#include "core/datacursor.h"

MainChannelBC::MainChannelBC(Element *startElement,
                                         Element *endElement,
                                         int variableIndex,
                                         STSModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_variableIndex(variableIndex),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

MainChannelBC::~MainChannelBC()
{
 delete m_dataCursor;
}

void MainChannelBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void MainChannelBC::prepare()
{

}

void MainChannelBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;


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
            element->mainChannelTemperature = value;
          }
        }
        break;
      default:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, static_cast<int>(i), m_dataCursor, value))
          {
            Element *element = m_profile[i];
            element->mainChannelSoluteConcs[m_variableIndex] = value;
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
            element->mainChannelTemperature = value;
          }
          break;
        default:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];

            element->mainChannelSoluteConcs[m_variableIndex] = value;

          }
          break;
      }
    }
  }
}

void MainChannelBC::clear()
{

}

Element *MainChannelBC::startElement() const
{
  return m_startElement;
}

void MainChannelBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *MainChannelBC::endElement() const
{
  return m_endElement;
}

void MainChannelBC::setEndElement(Element *element)
{
  m_endElement = element;
}

QSharedPointer<TimeSeries> MainChannelBC::timeSeries() const
{
  return m_timeSeries;
}

void MainChannelBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
