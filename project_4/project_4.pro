TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


LIBS += -llapack -larmadillo -lblas
SOURCES += main.cpp \
    lib.cpp

HEADERS += \
    lib.h
