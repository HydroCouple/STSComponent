/*!
*  \file    mcboundarycondition.h
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

#ifndef MCBOUNDARYCONDITION_H
#define MCBOUNDARYCONDITION_H

#include "stscomponent_global.h"
#include "iboundarycondition.h"

#include <QObject>
#include <QSharedPointer>

class TimeSeries;
class Element;
class STSModel;
class DataCursor;

class STSCOMPONENT_EXPORT MainChannelBC : public QObject,
    public virtual IBoundaryCondition
{
  public:

    MainChannelBC(Element *startElement,
                        Element *endElement,
                        int variableIndex,
                        STSModel *model);

    virtual ~MainChannelBC() override;

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:
    QSharedPointer<TimeSeries> m_timeSeries;
    std::vector<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_variableIndex;
    STSModel *m_model;
    DataCursor *m_dataCursor;
};


#endif // MCBOUNDARYCONDITION_H
