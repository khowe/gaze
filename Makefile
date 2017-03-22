
SRC = ./src
INC = ./include
OBJ = ./
BIN = ./


INCPATH = -I$(INC)
LIB =  -lm -lexpat 

OBJS =	$(OBJ)/util.o \
	$(OBJ)/structure.o \
	$(OBJ)/str_parse.o \
	$(OBJ)/options.o \
	$(OBJ)/engine.o \
	$(OBJ)/info.o \
	$(OBJ)/output.o \
	$(OBJ)/g_features.o \
	$(OBJ)/gff.o \
	$(OBJ)/g_engine.o \
	$(OBJ)/sequence.o \
	$(OBJ)/gaze.o

CC = gcc 
CFLAGS =  -c -O2 -Wall 


gaze : $(BIN)/gaze

all: $(BIN)/gaze

$(BIN)/gaze : $(OBJS)
	$(CC) -o $@ $(OBJS) $(LIB)

$(OBJ)/util.o : $(SRC)/util.c $(INC)/util.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/util.o $(SRC)/util.c

$(OBJ)/structure.o : $(SRC)/structure.c $(INC)/structure.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/structure.o $(SRC)/structure.c

$(OBJ)/str_parse.o : $(SRC)/str_parse.c $(INC)/str_parse.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/str_parse.o $(SRC)/str_parse.c

$(OBJ)/options.o : $(SRC)/options.c $(INC)/options.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/options.o $(SRC)/options.c

$(OBJ)/g_features.o : $(SRC)/g_features.c $(INC)/g_features.h 
	$(CC) $(CFLAGS) $(INCPATH) -c -o $(OBJ)/g_features.o $(SRC)/g_features.c

$(OBJ)/engine.o : $(SRC)/engine.c $(INC)/engine.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/engine.o $(SRC)/engine.c

$(OBJ)/info.o : $(SRC)/info.c $(INC)/info.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/info.o $(SRC)/info.c

$(OBJ)/output.o : $(SRC)/output.c $(INC)/output.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/output.o $(SRC)/output.c

$(OBJ)/g_engine.o : $(SRC)/g_engine.c $(INC)/g_engine.h
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/g_engine.o $(SRC)/g_engine.c

$(OBJ)/sequence.o : $(SRC)/sequence.c $(INC)/sequence.h
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/sequence.o $(SRC)/sequence.c

$(OBJ)/gff.o : $(SRC)/gff.c $(INC)/gff.h
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/gff.o $(SRC)/gff.c

$(OBJ)/gaze.o : $(SRC)/gaze.c
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/gaze.o $(SRC)/gaze.c

# clean up

clean :
	rm -f $(OBJ)/*.o $(BIN)/gaze

