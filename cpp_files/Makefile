#To compile the c++ codeS, do "make" in a terminal
#To compile landscape, do "make landscape.x"
#To compile mcmc_output, do "make mcmc_output.x"

CC        = g++
SOURCES = libraries/library_droso.cpp
OBJECTS = library_droso.o
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.x )

all: $(OBJECTS) $(OBJS)

%.x: %.cpp
	$(CC) $< $(SOURCES) -o $@

# compilation of the libraries
$(OBJECTS): $(SOURCES)
	$(CC)  -fpic -c $< -o $@

clean:
	rm -r *.x *.o *.out
