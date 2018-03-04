######################################################################
#   .____    .___  ________  _________._________________________     #
#   |    |   |   |/  _____/ /   _____/|   \_   _____/\__    ___/     #
#   |    |   |   /   \  ___ \_____  \ |   ||    __)    |    |        #
#   |    |___|   \    \_\  \/        \|   ||     \     |    |        #
#   |_______ \___|\______  /_______  /|___|\___  /     |____|        #
#           \/           \/        \/          \/                    #
#                                                                    #
######################################################################
CC  = g++
CFLAGS	=-c -O3 -Wall

## Similar path to be used on Ubuntu

INCLUDE =-I/home/ambrish/bin/openbabel/include/openbabel-2.0
LDFLAGS =-rdynamic /home/ambrish/bin/openbabel/lib/libopenbabel.so -Wl,-rpath,/home/ambrish/bin/openbabel/lib/

## Similar path To be used on MacOS
#INCLUDE =-I/Users/roya/bin/include/openbabel-2.0
#LDFLAGS =-lopenbabel -L/Users/roya/bin/lib

LIBS	=-lm
EXE = LIGSIFT

OBJ = main.o sm.o pharm.o coor.o alignment.o rmsd.o lap.o comm.o eigen.o


all: $(EXE)

LIGSIFT: $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o LIGSIFT  $(LDFLAGS)
	@rm -f *.o

main.o: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) main.cpp 

sm.o: Smallmolecule.cpp
	$(CC) $(CFLAGS) $(INCLUDE) Smallmolecule.cpp -o sm.o

alignment.o: Alignment.cpp 
	$(CC) $(CFLAGS) Alignment.cpp -o alignment.o

pharm.o: Pharmacophores.cpp
	$(CC) $(CFLAGS) Pharmacophores.cpp -o pharm.o

rmsd.o: RMSD.cpp
	$(CC) $(CFLAGS) RMSD.cpp -o rmsd.o

coor.o: Coordinate.cpp
	$(CC) $(CFLAGS) Coordinate.cpp -o coor.o

lap.o: LAP.cpp
	$(CC) $(CFLAGS) LAP.cpp -o lap.o

comm.o: Common.cpp
	$(CC) $(CFLAGS) Common.cpp -o comm.o

eigen.o: Eigen.cpp
	$(CC) $(CFLAGS) Eigen.cpp -o eigen.o

clean:
	@rm -f *.o  LIGSIFT
