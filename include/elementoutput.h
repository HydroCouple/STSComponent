#ifndef ELEMENTOUTPUT_H
#define ELEMENTOUTPUT_H

#include "stscomponent_global.h"
#include "spatiotemporal/timegeometryoutput.h"

class STSComponent;

class STSCOMPONENT_EXPORT ElementOutput: public TimeGeometryOutputDouble
{

    Q_OBJECT

  public:

    enum VariableType
    {
      Temperature,
      SoluteConc,
      Depth,
      Width,
      WidthFraction,
      XSectionArea,
      MCHeatFlux,
      MCSoluteFlux
    };

    ElementOutput(const QString &id,
                  Dimension *timeDimension,
                  Dimension *geometryDimension,
                  ValueDefinition *valueDefinition,
                  VariableType variableType,
                  STSComponent *modelComponent);


    virtual ~ElementOutput();

    void updateValues(HydroCouple::IInput *querySpecifier) override;

    void updateValues() override;

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    STSComponent *m_component;
    int m_variableType;
    int m_soluteIndex;
};

#endif // ELEMENTOUTPUT_H
