#ifndef ELEMENTINPUT_H
#define ELEMENTINPUT_H

#include "stscomponent_global.h"
#include "spatiotemporal/timegeometryinput.h"
#include "spatiotemporal/timegeometrymultiinput.h"

#include <unordered_map>

class STSComponent;


class STSCOMPONENT_EXPORT ElementInput : public TimeGeometryInputDouble
{
    Q_OBJECT

  public:

    enum VariableType
    {
      Depth,
      Width,
      MCTemperature,
      MCSolute,
    };

    ElementInput(const QString &id,
                 Dimension *timeDimension,
                 Dimension *geometryDimension,
                 ValueDefinition *valueDefinition,
                 VariableType varType,
                 STSComponent *modelComponent);

    virtual ~ElementInput() override;

    /*!
     * \brief setProvider
     * \param provider
     */
    bool setProvider(HydroCouple::IOutput *provider) override;

    /*!
     * \brief canConsume
     * \param provider
     * \param message
     * \return
     */
    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    /*!
     * \brief retrieveValuesFromProvider
     */
    void retrieveValuesFromProvider() override;

    /*!
     * \brief applyData
     */
    void applyData() override;


    /*!
     * \brief variableType
     * \return
     */
    VariableType variableType() const;

    /*!
     * \brief setVariableType
     * \param variableType
     */
    void setVariableType(VariableType variableType);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    std::unordered_map<int,int> m_geometryMapping;
    std::unordered_map<int,double> m_geometryMappingOrientation;
    STSComponent *m_component;
    VariableType m_varType;
    int m_soluteIndex = -1;

};

class STSCOMPONENT_EXPORT ElementSourceInput : public  TimeGeometryMultiInputDouble
{
  public:

    enum SourceType
    {
      RadiativeFlux,
      HeatFlux,
      SoluteFlux
    };

    ElementSourceInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           SourceType srcType,
                           STSComponent *modelComponent);

    virtual ~ElementSourceInput() override;

    bool addProvider(HydroCouple::IOutput *provider) override;

    bool removeProvider(HydroCouple::IOutput *provider) override;

    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    void retrieveValuesFromProvider() override;

    void applyData() override;

    SourceType sourceType() const;

    void setSourceType(SourceType srcType);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    STSComponent *m_component;
    SourceType m_srcType;
    std::unordered_map<HydroCouple::IOutput*, std::unordered_map<int,int>> m_geometryMapping;
    int m_soluteIndex;
};


#endif // ELEMENTINPUT_H
