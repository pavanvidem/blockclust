CC = g++
#CFLAGS =-c -Wall

all: BlockClust
BlockClust: BlockClust.o BlockGroup.o Block.o Read.o BlockEdge.o GspanCreator.o FeatureComputer.o BlockClust.o BlockGroupAnnotator.o 
BlockGroupAnnotator: BlockGroupAnnotator.o
#	${CC} ${CFLAGS} -o $@ $^
	$(CC) -o $@ $^

clean:
	rm *.o

cleanest: clean
	rm BlockClust
