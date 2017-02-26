TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

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
