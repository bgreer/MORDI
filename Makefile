


MKL = -L/central/intel/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
FFTW = -I/shared/fftw/include -L/shared/fftw/lib -lfftw3
TCLAP = -I/home8/begr7169/SOFTWARE/tclap-1.2.1/include

LIBS = ${FFTW} ${TCLAP} -lm $(MKL) -lrt
FLAG = -O2 -ip -ipo -openmp -Weffc++

CC = mpic++

EXECUTABLE = mcd

OBJECTS = main.o

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(FLAG) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)

$(OBJECTS) : *.cpp *.h
	$(CC) -c $(FLAG) *.cpp $(LIBS)

clean : 
	rm -f $(OBJECTS) $(EXECUTABLE)
