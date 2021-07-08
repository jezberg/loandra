CC = g++
CFLAGS = -std=c++11 -O2 -Wall -Wextra -Wshadow -g -Wfatal-errors
OBJFILES = preprocessor.o inputreader.o outputreader.o preprocessedinstance.o trace.o utility.o probleminstance.o timer.o clause.o log.o AMSLEX.o touchedlist.o preprocessorinterface.o cardinalityconstraint.o

all: preprocessor

preprocessor: main.o $(OBJFILES)
	$(CC) $(CFLAGS) main.o $(OBJFILES) -o preprocessor

.PHONY: lib
lib: preprocessorinterface.o $(OBJFILES)
	mkdir -p lib
	ar rvs lib/libpreprocessor.a preprocessorinterface.o $(OBJFILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -MMD -c $< -o $@

-include *.d

clean: 
	rm -f lib/*.a
	rm -f *.o *.d
