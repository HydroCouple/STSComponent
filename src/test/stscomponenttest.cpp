#include "stdafx.h"
#include "test/stscomponenttest.h"
#include "odesolver.h"
#include "stsmodel.h"
#include "elementjunction.h"
#include "element.h"
#include "variable.h"


void STSComponentTest::kuparuk_river1()
{
  QBENCHMARK_ONCE
  {

    //    //Error messages
    //    std::list<std::string> errors;

    //    //Create model instance
    //    STSModel *model = new STSModel(nullptr);

    //    //Read input filel
    //    model->setInputFile(QFileInfo("../../examples/kuparuk_river1/kuparuk_river1.inp"));

    //    //initialize model
    //    if(model->initialize(errors))
    //    {

    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }
    //    else
    //    {
    //      for(std::string error : errors)
    //      {
    //        printf("%s\n", error.c_str());
    //      }
    //    }

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void STSComponentTest::kuparuk_river2()
{
  QBENCHMARK_ONCE
  {

    //Error messages
    std::list<std::string> errors;

    //Create model instance
    STSModel *model = new STSModel(nullptr);

    //Read input filel
    model->setInputFile(QFileInfo("../../examples/kuparuk_river2/kuparuk_river2.inp"));

    //initialize model
    if(model->initialize(errors))
    {

      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }
    else
    {
      for(std::string error : errors)
      {
        printf("%s\n", error.c_str());
      }
    }


    //finalize model
    model->finalize(errors);

    delete model;
  }
}
