CC = gcc

CFLAGS = -c -O1 -I/home/$(USER)/local/src/RSFSRC/include -I/home/$(USER)/local/include/ -I/usr/include/

CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/src/RSFSRC/include -I/home/$(USER)/local/include/ -I/usr/include/

LFLAGS = -L/home/$(USER)/local/src/RSFSRC/lib  -L/home/$(USER)/local/lib -L/home/$(USER)/Madagascar/rsf/lib -Wl,-R /home/$(USER)/local/lib 

#ctoeplitz_reg:
#	echo Estoy compilando $@.c
#	$(CC) $(CFLAGS) $@.c -o $@.o


#myradon2:
#	echo Estoy compilando $@.c
#	$(CC) $(CFLAGS) $@.c -o $@.o
#	$(CC) $@.o $(LFLAGS) -lm -lrsf -o  $@.x

radon_rays:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lfftw3 -o  $@.x		

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
