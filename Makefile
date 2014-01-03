COMPILER = g++ 
COMPILER_FLAGS = -c -g -O3 -march=native $(WARNINGS)
WARNINGS = -Wall -Werror -Wextra -pedantic 
LINKER = g++ 
LINKER_FLAGS = -g

TAGGER = tagger 

all : $(TAGGER)

$(TAGGER) : pos-tagger.o backward.o forward.o viterbi.o log_functions.o
	$(LINKER) -g log_functions.o viterbi.o forward.o backward.o pos-tagger.cpp -o $@

pos-tagger.o : pos-tagger.cpp backward.o forward.o viterbi.o log_functions.o
	$(COMPILER) $(COMPILER_FLAGS) $<

%.o : %.cpp %.h
	$(COMPILER) $(COMPILER_FLAGS) $< -o $@

clean:
	-rm -rf *.o $(TAGGER)
