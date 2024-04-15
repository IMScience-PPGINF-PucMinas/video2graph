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
	# make -f MakefileDisf.make
	$(CXX) $(FLAGS) $(GFTFLAGS) \
		streamISF.cpp intersection.cpp $(GFTLIB) -o streamISF $(LINKS) -I $(DISF_LIB)/include -I $(DISF_LIB)/externals -L $(DISF_LIB)/lib -ldisf -g

run: streamISF
	./streamISF ./dataset/handwaving1 3000 30 100 1 0 1 50 1
	./streamISF ./dataset/handwaving2 3000 30 100 1 0 1 50 1
	./streamISF ./dataset/handwaving3 3000 30 100 1 0 1 50 1
	./streamISF ./dataset/handwaving4 3000 30 100 1 0 1 50 1
	./streamISF ./dataset/handwaving5 3000 30 100 1 0 1 50 1
	./streamISF ./dataset/walking1 3000 30 100 1 0 1 50 2
	./streamISF ./dataset/walking2 3000 30 100 1 0 1 50 2
	./streamISF ./dataset/walking3 3000 30 100 1 0 1 50 2
	./streamISF ./dataset/walking4 3000 30 100 1 0 1 50 2
	./streamISF ./dataset/walking5 3000 30 100 1 0 1 50 2

clean:
	# make -f MakefileDisf.make clean
	$(RM) *~ *.o streamISF