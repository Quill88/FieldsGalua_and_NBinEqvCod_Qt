TEMPLATE = app
TARGET = FieldsGaluaQt
DESTDIR = ./debug
QT += core
CONFIG += debug console
DEFINES += _CONSOLE
INCLUDEPATH += . \
    ./debug \
    $(QTDIR)/mkspecs/win32-msvc2013
DEPENDPATH += .
MOC_DIR += ./GeneratedFiles/debug
OBJECTS_DIR += debug
UI_DIR += ./GeneratedFiles
RCC_DIR += ./GeneratedFiles
HEADERS += ./galuafield.h \
    ./galuarow.h \
    ./koder.h \
    ./nBinEqvCod.h \
    ./nBinEqvVec.h
SOURCES += ./galuafield.cpp \
    ./galuarow.cpp \
    ./koder.cpp \
    ./main.cpp \
    ./nBinEqvCod.cpp \
    ./nBinEqvVec.cpp
