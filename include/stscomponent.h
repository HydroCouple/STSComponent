/*!
 *  \file    stscomponent.h
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

#ifndef STSCOMPONENT_H
#define STSCOMPONENT_H

#include "stscomponent_global.h"
#include "stscomponentinfo.h"
#include "temporal/abstracttimemodelcomponent.h"

class Dimension;
class Unit;
class Quantity;
class STSModel;
class HCGeometry;

class STSCOMPONENT_EXPORT STSComponent : public AbstractTimeModelComponent,
    public virtual HydroCouple::ICloneableModelComponent
{
    Q_OBJECT
    Q_INTERFACES(HydroCouple::ICloneableModelComponent)

  public:

    /*!
     * \brief STSComponent constructor
     * \param id Unique identifier for this component instance.
     * \param modelComponentInfo the parent ModelComponentInfo that generated this component instance.
     */
    STSComponent(const QString &id, STSComponentInfo* modelComponentInfo = nullptr);

    /*!
     * \brief ~STSComponent destructor
     */
    virtual ~STSComponent();

    /*!
     * \brief validate validates this component model instance
     * \return Returns a list of error messages.
     */
    QList<QString> validate() override;

    /*!
     * \brief prepare Prepares the model component instance.
     */
    void prepare() override;

    /*!
     * \brief update
     * \param requiredOutputs
     */
    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    /*!
     * \brief finish
     */
    void finish() override;

    /*!
     * \brief modelInstance
     * \return
     */
    STSModel *modelInstance() const;

    /*!
     * \brief parent
     * \return
     */
    HydroCouple::ICloneableModelComponent* parent() const override;

    /*!
     * \brief clone
     * \return
     */
    HydroCouple::ICloneableModelComponent* clone() override;

    /*!
     * \brief clones
     * \return
     */
    QList<HydroCouple::ICloneableModelComponent*> clones() const override;

  protected:

    /*!
     * \brief intializeFailureCleanUp
     */
    void initializeFailureCleanUp() override;

    /*!
     * \brief removeClone
     * \param component
     * \return
     */
    bool removeClone(STSComponent *component);

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    void createInputFileArguments();

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief initializeInputFilesArguments
     * \param message
     * \return
     */
    bool initializeInputFilesArguments(QString &message);

    /*!
     * \brief createGeometries
     */
    void createGeometries();

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    void createDepthInput();

    void createWidthInput();

    void createTemperatureInput();

    void createSoluteInput();

    void createHeatFluxInput();

    void createSoluteFluxInput();

    void createRadiativeFluxInput();

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    void createTemperatureOutput();

    void createSoluteConcOutput();

    void createDepthOutput();

    void createWidthOutput();

    void createWidthFractionOutput();

    void createHeatFluxOutput();

    void createSoluteFluxOutput();

    void createXSectionAreaOutput();

  private:

    IdBasedArgumentString *m_inputFilesArgument;

    Unit *m_radiationFluxUnit,
    *m_heatFluxUnit,
    *m_temperatureUnit,
    *m_soluteUnit,
    *m_soluteFluxUnit;

    Quantity *m_soluteConcQuantity,
             *m_soluteConcFluxQuantity;

    Dimension *m_timeDimension,
    *m_geometryDimension;

    STSModel *m_modelInstance;

    std::vector<QSharedPointer<HCGeometry>> m_elementGeometries;
    std::vector<QSharedPointer<HCGeometry>> m_elementJunctionGeometries;

    STSComponent *m_parent = nullptr;
    QList<HydroCouple::ICloneableModelComponent*> m_clones;
};

#endif //STSCOMPONENT_H
