CC = g++
CFLAGS =-Wall -O2 -g

all: blockclust
blockclust: src/BlockClust.o src/BlockGroup.o src/Block.o src/Read.o src/BlockEdge.o src/GspanCreator.o src/FeatureComputer.o src/BlockGroupAnnotator.o
	${CC} ${CFLAGS} -o $@ $^

clean:
	rm src/*.o

cleanest: clean
	rm blockclust
