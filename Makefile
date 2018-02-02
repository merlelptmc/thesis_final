#To compile the c++ codeS, do "make" in a terminal
#To compile landscape, do "make landscape.x"
#To compile mcmc_output, do "make mcmc_output.x"

CC        = g++ -std=c++0x
SRC_1     = mcmc_output.cpp
OUT_CPP_1 = mcmc_output.x
SRC_2	  = mfsa_output.cpp
OUT_CPP_2 = mfsa_output.x
SRC_3     = landscape.cpp
OUT_CPP_3 = landscape.x

SOURCES = library_droso.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(OUT_CPP_1) $(OUT_CPP_2) $(OUT_CPP_3)

$(OUT_CPP_1): $(OBJECTS) $(SRC_1) $(SOURCES)
	$(CC) $(SRC_1) $(OBJECTS) -o $(OUT_CPP_1)

$(OUT_CPP_2): $(OBJECTS) $(SRC_2) $(SOURCES)
	$(CC) $(SRC_2) $(OBJECTS) -o $(OUT_CPP_2)

$(OUT_CPP_3): $(OBJECTS) $(SRC_3) $(SOURCES)
	$(CC) $(SRC_1) $(OBJECTS) -o $(OUT_CPP_3)
	
# compilation of the libraries
$%.o : $%.cpp
	$(CC)  -fpic -c    $< -o $@
	
clean:
	rm -r *.x *.o *.out
