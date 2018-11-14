CXX=g++
CXXFLAGS+=-Wall -O2 -g

all: blockclust
blockclust: src/BlockClust.o src/BlockGroup.o src/Block.o src/Read.o src/BlockEdge.o src/GspanCreator.o src/FeatureComputer.o src/BlockGroupAnnotator.o
	${CXX} ${CXXFLAGS} -o $@ $^

clean:
	rm blockclust
	rm src/*.o
