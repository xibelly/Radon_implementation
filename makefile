CC = gcc

CFLAGS = -c -O0  -I/home/$(USER)/local/include/ -I/usr/include/

#-I/home/$(USER)/Madagascar/rsf/include

CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/include/ -I/usr/include/

LFLAGS = -L/home/$(USER)/local/lib -L/home/$(USER)/Madagascar/rsf/lib

#-Wl,-R /home/$(USER)/local/lib /home/$(USER)/Madagascar/rsf/lib


myradon2:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o


radon_rays_v2:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) myradon2.o $@.o $(LFLAGS) -lm -lfftw3 -o  $@.x		

debug: 
	echo Estoy compilando para debugear $@.c
	$(CC) $(CFLAGSDEBUG) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lfftw3 -o $@.x


clean:
	rm -rf *.~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
