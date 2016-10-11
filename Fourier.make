CC = gcc


##CFLAGS = -c -O2 -I/home/$(USER)/local/include/ -I/home/$(USER)/local2/include/ -I/usr/include/  
##CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/include/ -I/usr/include/
##LFLAGS = -lm -L/home/$(USER)/local/lib -L/home/$(USER)/local2/lib -Wl,-R /home/$(USER)/local/lib 

CFLAGS = -c -O2 -I/home/$(USER)/local/include/ -I/usr/include/  

CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/include/ -I/usr/include/

LFLAGS = -lm -L/home/$(USER)/local/lib -Wl,-R /home/$(USER)/local/lib 


Fourier:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lgsl -lgslcblas -lfftw3 -lnttw -lfrtw -lm -o  $@.x		

debug: 
	echo Estoy compilando para debugear $@.c
	$(CC) $(CFLAGSDEBUG) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lgsl -lgslcblas -lfftw3 -lnttw -lfrtw -o $@.x


clean:
	rm -rf *.~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
