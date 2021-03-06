/*!
 *  \file    stmcomponentinfo.cpp
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

#include "stdafx.h"
#include "stscomponentinfo.h"
#include "stscomponent.h"
#include "spatial/geometryfactory.h"

using namespace HydroCouple;

STSComponentInfo::STSComponentInfo(QObject *parent)
  :AbstractModelComponentInfo(parent)
{
  GeometryFactory::registerGDAL();

  setId("Surface Transient Storage Temperature Model 1.0.0");
  setCaption("Surface Transient Storage Model Component");
  setIconFilePath("./../../resources/images/hydrocouplecomposer.png");
  setDescription("A one-dimensional channel heat and solute transport model.");
  setCategory("Hydrodyanmics/Heat");
  setCopyright("");
  setVendor("");
  setUrl("www.hydrocouple.org");
  setEmail("caleb.buahin@gmail.com");
  setVersion("1.0.0");

  QStringList documentation;
  documentation << "Several sources";
  setDocumentation(documentation);
}

STSComponentInfo::~STSComponentInfo()
{
}

IModelComponent *STSComponentInfo::createComponentInstance()
{
  QString id =  QUuid::createUuid().toString();
  STSComponent *component = new STSComponent(id, this);
  component->setDescription("STM Model Instance");
  return component;
}
