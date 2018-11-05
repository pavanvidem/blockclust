/*
 * Block.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "Block.h"

Block::Block(){
	blockId = 0;
	blockStart = 0;
	blockEnd = 0;
	blockExpression = 0;
}

Block::Block(int blkId, unsigned long blkStart, unsigned long blkEnd, double blkExpression){
	blockId = blkId;
	blockStart = blkStart;
	blockEnd = blkEnd;
	blockExpression = blkExpression;
}

Block::~Block(){
}

unsigned long Block::getStart(){ return blockStart; }
unsigned long Block::getEnd(){ return blockEnd; }
double Block::getLength(){ return blockEnd - blockStart + 1; }
double Block::getExpression(){ return blockExpression; }
double Block::getArea(){ return getLength()*getExpression();}
vector<Read> Block::getReads(){ return reads; }
void Block::setReads(vector<Read> blockReads){ reads = blockReads; }

////////////////////////////////////////////

double Block::readLengthEntropy(){
	map<double, unsigned int> readLengthMap;
	float numberOfReads=0;
	vector<Read> readsList = getReads();
	vector<Read>::iterator readIterator;
	for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
		double readLength = readIterator->getLength();
		readLengthMap[readLength]++;
		numberOfReads++;
	}
	double readLengthEntropy = 0;
	for( map<double, unsigned int>::iterator ii=readLengthMap.begin(); ii!=readLengthMap.end(); ++ii){
		unsigned int lengthFrequency = ii->second;
		double lengthProbability = lengthFrequency/numberOfReads;
		readLengthEntropy += lengthProbability*(log2(lengthProbability));
	}
	readLengthEntropy *= -1;
	readLengthMap.clear();
	return readLengthEntropy;
}

////////////////////////////////////////////

double Block::readExpressionEntropy(){
	map<double, unsigned int> readExpressionMap;
	float numberOfReads=0;
	vector<Read> readsList = getReads();
	vector<Read>::iterator readIterator;
	for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
		double readExpression = readIterator->getExpression();
		double rouded = floor(readExpression * 100 + 0.5) / 100;
		readExpressionMap[rouded]++;
		numberOfReads++;
	}
	double readExpressionEntropy = 0;
	for( map<double, unsigned int>::iterator ii=readExpressionMap.begin(); ii!=readExpressionMap.end(); ++ii){
		unsigned int expressionFrequency = ii->second;
		double expressionProbability = expressionFrequency/numberOfReads;
		readExpressionEntropy += expressionProbability*(log2(expressionProbability));
	}
	readExpressionEntropy *= -1;
	readExpressionMap.clear();
	return readExpressionEntropy;
}

////////////////////////////////////////////

double Block::minReadLength(){
	vector<Read> readsList = getReads();
	vector<Read>::iterator readIterator;
	double minReadLength = 0;
	for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
		if(minReadLength == 0 || readIterator->getLength() < minReadLength)
			minReadLength = readIterator->getLength();
		
	}
	readsList.clear();

	return minReadLength;
}

////////////////////////////////////////////

double Block::noOfMultiMapped(){
	double multiMapped = 0;
	vector<Read> readsList = getReads();
	vector<Read>::iterator readIterator;
	for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
		double noOfMappings = readIterator->getNumberOfMappings();
		double readCount = readIterator->getReadCount();
		if(noOfMappings > 1)
			multiMapped ++;	
	}
	readsList.clear();

	return multiMapped;
}

////////////////////////////////////////////

double Block::medianReadExpression(){
	vector<Read> readsList = getReads();
	vector<Read>::iterator readIterator;
	vector<double> readExpressions;
	for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
		double readExpression = readIterator->getExpression();
		readExpressions.push_back(readExpression);
	}
	readsList.clear();

	typedef vector<double>::size_type vec_sz;
	vec_sz size = readExpressions.size();
	if (size == 0){
		cerr << "median of an empty vector" << endl;
		exit(0);
	}
	sort(readExpressions.begin(), readExpressions.end());
	vec_sz mid = size/2;
	return size % 2 == 0 ? (readExpressions[mid] + readExpressions[mid-1]) / 2 : readExpressions[mid];
}



