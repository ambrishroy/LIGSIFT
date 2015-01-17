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

OBJ = main.o sm.o pharm.o coor.o alignment.o rmsd.o lap.o comm.o eigen.o
CFLAGS = -c -I/usr/local/include/openbabel-2.0
LDFLAGS = -lopenbabel -L/usr/local/lib
EXE = LIGSIFT

all: $(EXE)

LIGSIFT: $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -O3 -o LIGSIFT
	@rm -f *.o

main.o: main.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -O3 main.cpp

sm.o: Smallmolecule.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c Smallmolecule.cpp -O3 -o sm.o

alignment.o: Alignment.cpp 
	$(CC) -c Alignment.cpp -O3 -o alignment.o

pharm.o: Pharmacophores.cpp
	$(CC) -c Pharmacophores.cpp -O3 -o pharm.o

rmsd.o: RMSD.cpp
	$(CC) -c RMSD.cpp -o rmsd.o

coor.o: Coordinate.cpp
	$(CC) -c Coordinate.cpp -O3 -o coor.o

lap.o: LAP.cpp
	$(CC) -c LAP.cpp -O3 -o lap.o

comm.o: Common.cpp
	$(CC) -c Common.cpp -O3 -o comm.o

eigen.o: Eigen.cpp
	$(CC) -c Eigen.cpp -O3 -o eigen.o

clean:
	@rm -f *.o  LIGSIFT
