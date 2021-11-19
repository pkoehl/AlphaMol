# ================= Libraries and Includes ====================================

INOUT_INCLUDE_DIR = ./InOut/include
STRUCT_INCLUDE_DIR = ./DataStruct/include
DELCX_INCLUDE_DIR = ./Delcx/include
ALFCX_INCLUDE_DIR = ./Alphacx/include
VOLUMES_INCLUDE_DIR = ./Volumes/include

INCLUDE_DIRS = -I$(INC_DIR) -I$(INOUT_INCLUDE_DIR) -I$(STRUCT_INCLUDE_DIR)\
		-I$(DELCX_INCLUDE_DIR) -I$(ALFCX_INCLUDE_DIR) -I$(VOLUMES_INCLUDE_DIR)
GMP_LIB_DIR =
LIB_DIRS = 

GMP_LIBS = -lgmp
LIBS = $(GMP_LIBS) -lstdc++

# ================= Project Directories ====================================

INC_DIR = ./project/include
SRC_DIR = ./project/src
OBJ_DIR = ./project/src
BIN_DIR = ./bin

# ================= Project Name ===========================================

EXT=
NAME=AlphaMol
NAMEOBJ=$(OBJ_DIR)/$(NAME).o
NAMEBIN=$(BIN_DIR)/$(NAME)$(EXT)

# ================= Compilers and Flags ====================================

CC 		:= gcc
CFLAGS 		:= -c -O3 -ansi -Wall -Werror -pedantic 
CPP 		:= g++
CPPFLAGS 	:= -c -O3 -Wdeprecated-declarations -Wuninitialized -ansi -Werror -Wunused -std=c++11
FC 		:= gfortran
FFLAGS 		:= -c -O3 -fcray-pointer

LOAD_LIB_PATH =

LD_FLAGS = -O3

# ================= Pattern rules ==========================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIRS) $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $< -o $@

# ================= Compile source code ====================================

OBJECTS = \
$(NAMEOBJ)

# ================= Generate Executable ====================================

$(NAMEBIN) : $(OBJECTS)
	$(CPP) -o $(NAMEBIN) $(LD_FLAGS) $(OBJECTS) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

all: $(OBJECTS)
	$(CPP) -o $(NAMEBIN) $(LD_FLAGS) $(OBJECTS) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

clean:
	touch $(OBJ_DIR)/junk.o; rm -f $(OBJ_DIR)/*.o $(NAMEBIN)

$(OBJECTS):
