#Compiladores
CC=g++
CXX=g++

GFT_DIR=./gft
DISF_LIB = ./disf
# FLAGS= -Wall -O3 -lpthread -msse
# #-march=native 

# LINKS= -lz -lm -fopenmp

FLAGS= -Wall -g

LINKS= -lpthread -lz -lm 

#Bibliotecas
GFTLIB  = -L$(GFT_DIR)/lib -lgft
GFTFLAGS  = -I$(GFT_DIR)/include


#Rules
all: run

libgft:
	$(MAKE) -C $(GFT_DIR)

streamISF: streamISF.cpp intersection.cpp libgft
	make -f MakefileDisf.make
	$(CXX) $(FLAGS) $(GFTFLAGS) \
		streamISF.cpp intersection.cpp $(GFTLIB) -o streamISF $(LINKS) -I $(DISF_LIB)/include -I $(DISF_LIB)/externals -L $(DISF_LIB)/lib -ldisf -g

run: clean streamISF
#	./streamISF ./datasets/soccer output/soccer 5000 100 150 1 0 1 50
#	./streamISF ./datasets/girl output/girl 5000 100 100 1 0 1 50
	./streamISF ./dataset/handwaving1 output/handwaving1 3000 30 100 1 0 1 50
#	./streamISF ./dataset/walking5 output/walking5 3000 30 100 1 0 1 50

clean:
	make -f MakefileDisf.make clean
	$(RM) *~ *.o streamISF