CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) 
LIBS = -L$(CFITSIO) -lcfitsio -lm
GLIBS = 
GLIBS += 
OBJECTS = subtractImages.o 
HEADERS = globalConstants.h

ALL : subtractImages.exe
	@echo "Listo!"

subtractImages.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o subtractImages.exe $(LIBS) $(GLIBS) $(CFLAGS)

subtractImages.o : subtractImages.cc $(HEADERS)
	$(CPP) -c subtractImages.cc -o subtractImages.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
