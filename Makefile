
# If you have installed glib and/or expat locally,
# you will have to amend the following lines

# put the location of glib.h here. 
# If centrally installed, comment this line out.
GLIBINC = -I/usr/local/include
# put the location of the glibconfig.h here. 
# If centrally installed, comment this line out.
GLIBCONFIGINC = -I/usr/local/lib/glib/include 
# Put the location of expat.h file here. 
# If centrally installed, comment this line out.
EXPATINC = -I$(HOME)/libraries/osf/include

# put the location of the glib and expat libraries themselves here.
# If centrally installed, comment these lines out
EXPATLIB = -L$(HOME)/libraries/osf/lib
GLIBLIB = 

###################################################
##### shouldn't need to change anything below here 
###################################################

OBJ = ./obj
SRC = ./src
BIN = ./bin
INC = ./include

INCPATH = -I$(INC) $(GLIBINC) $(GLIBCONFIGINC) $(EXPATINC)

LIB = $(EXPATLIB) $(GLIBLIB) -lm -lglib -lexpat 

OBJS = 	$(OBJ)/structure.o \
	$(OBJ)/str_parse.o \
	$(OBJ)/options.o \
	$(OBJ)/engine.o \
	$(OBJ)/info.o \
	$(OBJ)/output.o \
	$(OBJ)/features.o \
	$(OBJ)/g_engine.o  \
	$(OBJ)/gaze.o

#TRACE_LEV = -DTRACE=1

#CC = gcc 
#CFLAGS =  -c -O2 -Wall 

CC = cc
CFLAGS = -c -O2 -fullwarn


gaze : $(BIN)/gaze

all: $(BIN)/gaze

$(BIN)/gaze : $(OBJS)
	$(CC) -o $@ $(OBJS) $(LIB)

$(OBJ)/structure.o : $(SRC)/structure.c $(INC)/structure.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/structure.o $(SRC)/structure.c

$(OBJ)/str_parse.o : $(SRC)/str_parse.c $(INC)/str_parse.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/str_parse.o $(SRC)/str_parse.c

$(OBJ)/options.o : $(SRC)/options.c $(INC)/options.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/options.o $(SRC)/options.c

$(OBJ)/features.o : $(SRC)/features.c $(INC)/features.h 
	$(CC) $(CFLAGS) $(INCPATH) -c -o $(OBJ)/features.o $(SRC)/features.c

$(OBJ)/engine.o : $(SRC)/engine.c $(INC)/engine.h
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/engine.o $(SRC)/engine.c

$(OBJ)/info.o : $(SRC)/info.c $(INC)/info.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/info.o $(SRC)/info.c

$(OBJ)/output.o : $(SRC)/output.c $(INC)/output.h 
	$(CC) $(CFLAGS) $(INCPATH) -o $(OBJ)/output.o $(SRC)/output.c

$(OBJ)/g_engine.o : $(SRC)/g_engine.c $(INC)/g_engine.h
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/g_engine.o $(SRC)/g_engine.c

$(OBJ)/gaze.o : $(SRC)/gaze.c
	$(CC) $(CFLAGS) $(TRACE_LEV) $(INCPATH) -o $(OBJ)/gaze.o $(SRC)/gaze.c

# clean up

clean :
	rm -f $(OBJ)/*.o

