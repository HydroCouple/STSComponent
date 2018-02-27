/*!
 *  \file    stscomponent.cpp
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

#include "stdafx.h"
#include "stscomponent.h"


STSComponent::STSComponent(const QString &id, STSComponentInfo *modelComponentInfo)
  :AbstractModelComponent(id, modelComponentInfo)
{

}

STSComponent::~STSComponent()
{

}

QList<QString> STSComponent::validate()
{


  return QList<QString>();
}

void STSComponent::prepare()
{

}

void STSComponent::update(const QList<HydroCouple::IOutput *> &requiredOutputs)
{

}

void STSComponent::finish()
{

}

void STSComponent::intializeFailureCleanUp()
{

}

void STSComponent::createArguments()
{

}

bool STSComponent::initializeArguments(QString &message)
{

}

void STSComponent::createInputs()
{

}

void STSComponent::createOutputs()
{

}

void STSComponent::updateOutputValues(const QList<HydroCouple::IOutput *> &requiredOutputs)
{

}
