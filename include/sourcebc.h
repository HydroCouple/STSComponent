/*!
*  \file    nonpointsrctimeseriesbc.h
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

#ifndef NONPOINTSRCTIMESERIESBC_H
#define NONPOINTSRCTIMESERIESBC_H


#include "stscomponent_global.h"
#include "iboundarycondition.h"

#include <QObject>
#include <QSharedPointer>
#include <unordered_map>

class TimeSeries;
class Element;
class STSModel;
class DataCursor;

class STSCOMPONENT_EXPORT SourceBC :  public QObject,
    public virtual IBoundaryCondition
{
  public:

    SourceBC(Element *startElement, double startElementLFactor,
                            Element *endElement, double endElementLFactor,
                            int variableIndex, STSModel *model);

    virtual ~SourceBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    double startElementLFactor() const;

    void setStartElementLFactor(double factor);

    Element *endElement() const;

    void setEndElement(Element *element);

    double endElementLFactor() const;

    void setEndElementLFactor(double factor);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:
    std::vector<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    std::unordered_map<Element*, double> m_factors;
    double m_startElementLFactor, m_endElementLFactor;
    int m_variableIndex;
    QSharedPointer<TimeSeries> m_timeSeries;
    STSModel *m_model;
    DataCursor *m_dataCursor;
};


#endif // NONPOINTSRCTIMESERIESBC_H
