CPPS=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
TARGETS=$(patsubst %.cpp,%,${CPPS})
TARGETSOPTIM=$(patsubst %.cpp,%-optim,${CPPS})
FLAGS=
INCS=-fopenmp
LIBS=-lnewntl -lfplll -lmpfr -lgmp -lm

all: ${TARGETS}

optim: ${TARGETSOPTIM}

%: %.cpp ${HEADERS} 
	g++ -g3 -O0 -Wall ${FLAGS} ${INCS} -o $@ $< ${LIBS}

%-optim: %.cpp ${HEADERS}
	g++ -Wall -O9 -DNDEBUG ${FLAGS} ${INCS} -o $@ $< ${LIBS}


clean:
	rm -f *.o ${TARGETS} ${TARGETSOPTIM}
