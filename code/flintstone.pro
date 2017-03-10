TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
SOURCES += main.cpp \
    basis.cpp \
    gc.cpp \
    Coulomb_Functions.cpp

HEADERS += \
    basis.h \
    maths.h \
    gc.h \
    Coulomb_Functions.hpp

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
LIBS += -larmadillo -llapack -lblas
