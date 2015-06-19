CXX=g++ -std=c++0x # -cxxlib-icc-abr

CXXFLAGS=  -march=native -O3 -Wall  -Wextra -Wshadow
#CXXFLAGS=  -g  -pg -Wall  -Wextra -Wshadow 

OUT_CODE=pop
SRC= BP_recursion_zeroT.cpp main.cpp field.cpp
LIB= BP_recursion.h field.h data.h MersenneTwister.h

OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")

all:	$(OUT_CODE)

$(OUT_CODE):$(OBJ) 
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ -lstdc++ -lm


$(OBJ): $(SRC1)
	$(CXX) $(CXXFLAGS) -o $@ -c $(@:.o=.cpp)

tar: 
	tar -cjvf population_dynamics_$(DATE).tz2 $(SRC) $(LIB) Makefile
