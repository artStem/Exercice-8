CC = g++
CFLAGS = -std=c++0x -Wall -g
EXEC_NAME = Exercice8
INCLUDES =
LIBS =
OBJ_FILES = Exercice8_Majeux_Stemplowski.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES) *.out*

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

