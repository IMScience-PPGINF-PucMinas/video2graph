DEMO_DIR = ./disf

BIN_DIR = $(DEMO_DIR)/bin
STB_DIR = $(DEMO_DIR)/externals
SRC_DIR = $(DEMO_DIR)/src
LIB_DIR = $(DEMO_DIR)/lib
INCLUDE_DIR = $(DEMO_DIR)/include
OBJ_DIR = $(DEMO_DIR)/obj

CC = gcc
CFLAGS = -g -pg -Wall -fPIC -std=gnu11 -pedantic -Wno-unused-result -O3 
LIBS = -lm

HEADER_INC = -I $(STB_DIR) -I $(INCLUDE_DIR)
LIB_INC = -L $(LIB_DIR) -ldisf

all: obj
	$(eval ALL_OBJS := $(wildcard $(OBJ_DIR)/*.o))
	ar csr $(LIB_DIR)/libdisf.a $(ALL_OBJS)

obj: \
$(OBJ_DIR)/Utils.o \
$(OBJ_DIR)/IntList.o \
$(OBJ_DIR)/Color.o \
$(OBJ_DIR)/PrioQueue.o \
$(OBJ_DIR)/Image.o \
$(OBJ_DIR)/DISF.o 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INCLUDE_DIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@ $(HEADER_INC) $(LIBS)

$@.c: $@.c
	$(CC) $(CFLAGS) $@.c -o $(BIN_DIR)/$@ $(HEADER_INC) $(LIB_INC) $(LIBS)

clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(LIB_DIR)/*.a
	rm -rf $(BIN_DIR)/*
