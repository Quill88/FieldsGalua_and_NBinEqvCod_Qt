TEMPLATE = app
TARGET = FieldsGaluaQt

QT += core
CONFIG += debug console c++11
CONFIG -= app_bundle
DEFINES += _CONSOLE
TEMPLATE = app
INCLUDEPATH += $$PWD \
	$$PWD/gmp/include \
    $$PWD/debug
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
LIBS += -L$$PWD/gmp/lib/ -lgmpxx -lgmp
DEPENDPATH += $$PWD/gmp/include
MOC_DIR += $$PWD/GeneratedFiles/debug
OBJECTS_DIR += $$PWD/debug
UI_DIR += $$PWD/GeneratedFiles
RCC_DIR += $$PWD/GeneratedFiles

HEADERS += $$PWD/galuafield.h \
    $$PWD/galuarow.h \
    $$PWD/koder.h \
    $$PWD/nBinEqvCod.h \
    $$PWD/nBinEqvVec.h \
    $$PWD/smart_fact.h
SOURCES += $$PWD/galuafield.cpp \
    $$PWD/galuarow.cpp \
    $$PWD/koder.cpp \
    $$PWD/main.cpp \
    $$PWD/nBinEqvCod.cpp \
    $$PWD/nBinEqvVec.cpp
