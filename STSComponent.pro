#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2018
#License GNU Lesser General Public License (see <http: //www.gnu.org/licenses/> for details).
#The STSComponent is a stream surface transient storage zone temperature and solute transport model.

TEMPLATE = lib
VERSION = 1.0.0
TARGET = STSComponent
QT -= gui
QT += testlib


DEFINES += STSCOMPONENT_LIBRARY
DEFINES += USE_OPENMP
DEFINES += USE_MPI
DEFINES += USE_CVODE
DEFINES += USE_NETCDF
#DEFINES += USE_CVODE_OPENMP

#Compile as library or executable
contains(DEFINES,STSCOMPONENT_LIBRARY){
  TEMPLATE = lib
  message("Compiling HTSComponent as library")
} else {
  TEMPLATE = app
  CONFIG-=app_bundle
  message("Compiling HTSComponent as application")
}

CONFIG += c++11

linux{
CONFIG += debug_and_release
}

#Added for faster compilation
*msvc* { # visual studio spec filter
      QMAKE_CXXFLAGS += /MP /O2
}

PRECOMPILED_HEADER = ./include/stdafx.h

INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include \
               ./../ODESolver/include


HEADERS += ./include/stdafx.h\
           ./include/stscomponent_global.h \
           ./include/stscomponent.h \
           ./include/stscomponentinfo.h \
           ./include/stsmodel.h \
           ./include/iboundarycondition.h \
           ./include/elementjunction.h \
           ./include/test/stscomponenttest.h \
           ./include/variable.h \
           ./include/element.h \
           ./include/iboundarycondition.h \
           ./include/radiativefluxbc.h \
           ./include/hydraulicsbc.h \
           ./include/pointsrctimeseriesbc.h \
           ./include/sourcebc.h \
           ./include/meteorologybc.h \
           ./include/mainchannelbc.h \
           ./include/elementinput.h \
           include/elementoutput.h

SOURCES +=./src/stdafx.cpp \
          ./src/stscomponent.cpp \
          ./src/stscomponentinfo.cpp \
          ./src/elementinput.cpp\
          ./src/main.cpp \
          ./src/element.cpp \
          ./src/stsmodel.cpp \
          ./src/elementjunction.cpp \
          ./src/test/stscomponenttest.cpp \
          ./src/stsmodelio.cpp \
          ./src/stscompute.cpp \
          ./src/radiativefluxbc.cpp \
          ./src/hydraulicsbc.cpp \
          ./src/sourcebc.cpp \
          ./src/meteorologybc.cpp \
          ./src/mainchannelbc.cpp \
          src/elementoutput.cpp

macx{

    INCLUDEPATH += /usr/local \
                   /usr/local/include

    contains(DEFINES, USE_CVODE){

    contains(DEFINES, USE_CVODE){
    message("CVODE enabled")
    LIBS += -L/usr/local/lib -lsundials_cvode
     }

    contains(DEFINES, USE_NETCDF){

    message("NetCDF enabled")
    LIBS += -L/usr/local/lib -lnetcdf-cxx4
    }

    LIBS += -L/usr/local/lib -lnetcdf-cxx4

    contains(DEFINES,USE_OPENMP){

        QMAKE_CC = /usr/local/opt/llvm/bin/clang
        QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
        QMAKE_LINK = /usr/local/opt/llvm/bin/clang++

        QMAKE_CFLAGS+= -fopenmp
        QMAKE_LFLAGS+= -fopenmp
        QMAKE_CXXFLAGS+= -fopenmp

        INCLUDEPATH += /usr/local/opt/llvm/lib/clang/5.0.0/include
        LIBS += -L /usr/local/opt/llvm/lib -lomp

      message("OpenMP enabled")
     } else {
      message("OpenMP disabled")
     }

    contains(DEFINES,USE_MPI){

        QM
        QMAKE_CXX = /usr/local/bin/mpicxx
        QMAKE_LINK = /usr/local/bin/mpicxx

        QMAKE_CFLAGS += $$system(/usr/local/bin/mpicc --showme:compile)
        QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
        QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)

        LIBS += -L/usr/local/lib -lmpi

        message("MPI enabled")
     } else {

      message("MPI disabled")
     }
}

linux{

INCLUDEPATH += /usr/include \
               ../gdal/include

    contains(DEFINES,UTAH_CHPC){

         INCLUDEPATH += /uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.4.1/include \
                        ../netcdf-cxx4-4.3.0/installdir/include \
                        ../hypre/build/include


         LIBS += -L/uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/lib -lhdf5 \
                 -L/uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/lib -lnetcdf_c++4 \
                 -L../hypre/build/lib -lHYPRE

        contains(DEFINES,USE_CVODE){

            message("CVODE enabled")

            INCLUDEPATH += ../sundials-3.1.1/instdir/include
            LIBS += -L../sundials-3.1.1/instdir/lib -lsundials_cvode
        }

         message("Compiling on CHPC")
     }

    contains(DEFINES,USE_OPENMP){

    QMAKE_CFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    QMAKE_CXXFLAGS += -fopenmp

    LIBS += -L/usr/lib/x86_64-linux-gnu -lgomp

      message("OpenMP enabled")
     } else {

      message("OpenMP disabled")
     }
}

win32{

    #Windows vspkg package manager installation path
    VCPKGDIR = C:/vcpkg/installed/x64-windows

    INCLUDEPATH += $${VCPKGDIR}/include \
                   $${VCPKGDIR}/include/gdal

    CONFIG(debug, debug|release) {
    LIBS += -L$${VCPKGDIR}/debug/lib -lgdald
        } else {
    LIBS += -L$${VCPKGDIR}/lib -lgdal
    }


    contains(DEFINES, USE_CVODE){
    message("CVODE enabled")
    CONFIG(debug, debug|release) {
        message("CVODE debug")

        LIBS += -L$${VCPKGDIR}/debug/lib -lsundials_cvode

        } else {

        LIBS += -L$${VCPKGDIR}/lib -lsundials_cvode

        }
    }

    contains(DEFINES, USE_NETCDF){
    message("NetCDF enabled")
    CONFIG(release, debug|release) {
        LIBS += -L$${VCPKGDIR}/lib -lnetcdf \
                -L$${VCPKGDIR}/lib -lnetcdf-cxx4
        } else {
        LIBS += -L$${VCPKGDIR}/debug/lib -lnetcdf \
                -L$${VCPKGDIR}/debug/lib -lnetcdf-cxx4
        }
    }
    
    contains(DEFINES,USE_OPENMP){

        QMAKE_CFLAGS += -openmp
        QMAKE_LFLAGS += -openmp
        QMAKE_CXXFLAGS += -openmp
        QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
        QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

        message("OpenMP enabled")
     } else {

      message("OpenMP disabled")
     }

    contains(DEFINES,USE_MPI){
       message("MPI enabled")

        CONFIG(debug, debug|release) {
            LIBS += -L$${VCPKGDIR}/debug/lib -lmsmpi
          } else {
            LIBS += -L$${VCPKGDIR}/lib -lmsmpi
        }

    } else {
      message("MPI disabled")
    }

    QMAKE_CXXFLAGS += /MP
    QMAKE_LFLAGS += /incremental /debug:fastlink
}

CONFIG(debug, debug|release) {

    win32 {
       QMAKE_CXXFLAGS += /MDd /O2
    }

    macx {
       QMAKE_CXXFLAGS += -O3
    }

    linux {
       QMAKE_CXXFLAGS += -O3
    }

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{

    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/ &&"
    QMAKE_POST_LINK += "cp -a ./../ODESolver/build/debug/*ODESolver.* ./build/debug/";

    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.1.0.0 \
            -L./../ODESolver/build/debug -lODESolver

  }

   linux{

    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.so.1.0.0
     }

   win32{

    QMAKE_POST_LINK += "copy ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK1
     }
}

CONFIG(release, debug|release) {


   win32 {
    QMAKE_CXXFLAGS += /MD
   }

    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib/macx -lHydroCoupleSDK \
            -L./../ODESolver/lib/macx -lODESolver
   }

   linux{
    LIBS += -L./../HydroCoupleSDK/lib/linux -lHydroCoupleSDK \
            -L./../ODESolver/lib/linux -lODESolver
   }

   win32{
    LIBS += -L./../HydroCoupleSDK/lib/win32 -lHydroCoupleSDK1 \
            -L./../ODESolver/lib/win32 -lODESolver1
   }


     contains(DEFINES,STSCOMPONENT_LIBRARY){
         #MacOS
         macx{
             DESTDIR = lib/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./lib/macx/ &&"
             QMAKE_POST_LINK += "cp -a ./../ODESolver/lib/macx/*ODESolver.* ./lib/macx/";
          }

         #Linux
         linux{
             DESTDIR = lib/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./lib/linux/ &&"
             QMAKE_POST_LINK += "cp -a ./../ODESolver/lib/linux/*ODESolver.* ./lib/linux/";
          }

         #Windows
         win32{
             DESTDIR = lib/win32
             QMAKE_POST_LINK += "copy /B .\..\HydroCoupleSDK\lib\win32\HydroCoupleSDK* .\lib\win32 &&"
             QMAKE_POST_LINK += "copy /B .\..\ODESolver\lib\win32\ODESolver* .\lib\win32"
          }
     } else {
         #MacOS
         macx{
             DESTDIR = bin/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./bin/macx/ &&"
             QMAKE_POST_LINK += "cp -a ./../ODESolver/lib/macx/*ODESolver.* ./bin/macx/";
          }

         #Linux
         linux{
             DESTDIR = bin/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./bin/linux/ &&"
             QMAKE_POST_LINK += "cp -a ./../ODESolver/lib/linux/*ODESolver.* ./bin/linux/";
          }

         #Windows
         win32{
             DESTDIR = bin/win32
             QMAKE_POST_LINK += "copy /B .\..\HydroCoupleSDK\lib\win32\HydroCoupleSDK* .\bin\win32 &&"
             QMAKE_POST_LINK += "copy /B .\..\ODESolver\lib\win32\ODESolver* .\bin\win32"
          }
     }
  }
}
