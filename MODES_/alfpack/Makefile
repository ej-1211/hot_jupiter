

# This Makefile will build ALFPACK on most unix systems.  To use it,
# name this file "Makefile", and put it in the same directory as the 
# ALFPACK source code, and simply type "make".
#
#
# Set the following path to the directory where you want the ALFPACK
# binary to be installed.
#
FC = gfortran
ALFPACK	= .

LIBDIR	= $(ALFPACK)

OBJECTS = \
	 alfk.o dalfk.o dlfma.o dlfma1.o dlfpt.o drecn.o dtrih.o lfim.o\
	 lfim1.o lfin.o lfin1.o lfma.o lfma1.o lfmb.o lfmb1.o lfna.o\
	 lfna1.o lfnb.o lfnb1.o lfnc.o lfnc1.o lfp.o lfp1.o lfpt.o qp.o\
	 recm.o trih.o trivec.o dlfnc.o dlfnc1.o

FSOURCES = \
	 alfk.f dalfk.f dlfma.f dlfma1.f dlfpt.f drecn.f dtrih.f lfim.f\
	 lfim1.f lfin.f lfin1.f lfma.f lfma1.f lfmb.f lfmb1.f lfna.f\
	 lfna1.f lfnb.f lfnb1.f lfnc.f lfnc1.f lfp.f lfp1.f lfpt.f qp.f\
	 recm.f trih.f trivec.f dlfnc.f dlfnc1.f

OTHR_SOURCES = Makefile README

FLOAT_OPT =

DEBUG_OPT = -O2

FFLAGS	= $(FLOAT_OPT) $(DEBUG_OPT)

CFLAGS	= $(FLOAT_OPT) $(DEBUG_OPT)

TARGET	= $(LIBDIR)/libalfpack.a

all:	$(TARGET)

$(TARGET): $(OBJECTS)
	ar -vru $(TARGET) $(OBJECTS)
	ranlib $(TARGET)

#
# When ALFPACK has been successfully installed, you will want to clean
# up the ALFPACK directory.  To do this, type: "make clean".
#
clean:
	rm -f $(OBJECTS) core
