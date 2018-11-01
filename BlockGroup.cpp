/*
 * BlockGroup.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */


#include "BlockGroup.h"

BlockGroup::BlockGroup(){
	blockGroupId = 0;
	blockGroupStart = 0;
	blockGroupEnd = 0;
	blockGroupExpression = 0;
	blockGroupReadCount = 0;
	blockGroupClass = "";
}


BlockGroup::~BlockGroup(){
}

unsigned int BlockGroup::getId(){
	return blockGroupId;
}
void BlockGroup::setId(unsigned int id){
	blockGroupId = id;
}
unsigned long BlockGroup::getStart(){
	return blockGroupStart;
}
void BlockGroup::setStart(unsigned long start){
	blockGroupStart = start;
}
unsigned long BlockGroup::getEnd(){
	return blockGroupEnd;
}
void BlockGroup::setEnd(unsigned long end){
	blockGroupEnd = end;
}
double BlockGroup::getLength(){
	return blockGroupEnd - blockGroupStart +1;
}
double BlockGroup::getExpression(){
	return blockGroupExpression;
}
void BlockGroup::setExpression(double expr){
	blockGroupExpression = expr;
}
unsigned int BlockGroup::getReadCount(){
	return blockGroupReadCount;
}
void BlockGroup::setReadCount(unsigned int count){
	blockGroupReadCount = count;
}
string BlockGroup::getClass(){
	return blockGroupClass;
}
void BlockGroup::setClass(string classLabel){
	blockGroupClass = classLabel;
}
vector<Block> BlockGroup::getBlocks(){
	return blocks;
}
void BlockGroup::setBlocks(vector<Block> blks){
	blocks = blks;
}
double BlockGroup::getNumberOfBlocks(){
	return blocks.size();
}
vector<BlockEdge> BlockGroup::getBlockEdges(){
	return blockEdges;
}
void BlockGroup::setBlockEdges(vector<BlockEdge> blkedgs){
	blockEdges = blkedgs;
}

////////////////////////////////////////////

double BlockGroup::startPositionEntropy(){
	vector<Block> blocksList = getBlocks();
	vector<Block>::iterator blockIterator;
	map<double, unsigned long> readStartPositonMap;
	float numberOfReads=0;
	for (blockIterator = blocksList.begin(); blockIterator != blocksList.end(); ++blockIterator ){
		vector<Read> readsList = blockIterator->getReads();
		vector<Read>::iterator readIterator;
		for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
			unsigned int readStart = readIterator->getStart();
			readStartPositonMap[readStart]++;
			numberOfReads++;
		}
		readsList.clear();
	}
	blocksList.clear();

	return computeEntropy(readStartPositonMap, numberOfReads);
}

////////////////////////////////////////////

double BlockGroup::endPositionEntropy(){
	vector<Block> blocksList = getBlocks();
	vector<Block>::iterator blockIterator;
	map<double, unsigned long> readEndPositonMap;
	float numberOfReads=0;
	for (blockIterator = blocksList.begin(); blockIterator != blocksList.end(); ++blockIterator ){
		vector<Read> readsList = blockIterator->getReads();
		vector<Read>::iterator readIterator;
		for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
			unsigned int readEnd = readIterator->getEnd();
			readEndPositonMap[readEnd]++;
			numberOfReads++;
		}
		readsList.clear();
	}
	blocksList.clear();

	return computeEntropy(readEndPositonMap, numberOfReads);
}

////////////////////////////////////////////

double BlockGroup::readLengthEntropy(){
	vector<Block> blocksList = getBlocks();
	vector<Block>::iterator blockIterator;
	map<double, unsigned long> readLengthMap;
	float numberOfReads=0;
	for ( blockIterator = blocksList.begin(); blockIterator != blocksList.end(); ++blockIterator ){
		vector<Read> readsList = blockIterator->getReads();
		vector<Read>::iterator readIterator;
		for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
			unsigned int readLength = readIterator->getLength();
			readLengthMap[readLength]++;
			numberOfReads++;
		}
		readsList.clear();
	}
	blocksList.clear();

	return computeEntropy(readLengthMap, numberOfReads);
}

////////////////////////////////////////////

double BlockGroup::medianReadExpression(){
	vector<Block> blocksList = getBlocks();
	vector<Block>::iterator blockIterator;
	vector<double> readExpressions;
	for(blockIterator = blocksList.begin(); blockIterator != blocksList.end(); ++blockIterator ){
		vector<Read> readsList = blockIterator->getReads();
		vector<Read>::iterator readIterator;
		for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
			double readExpression = readIterator->getExpression();
			readExpressions.push_back(readExpression);
		}
		readsList.clear();
	}
	blocksList.clear();

	typedef vector<double>::size_type vec_sz;
	vec_sz size = readExpressions.size();
	if (size == 0){
		cerr << "empty read expressions" << endl;
		exit(0);
	}
	sort(readExpressions.begin(), readExpressions.end());
	vec_sz mid = size/2;
	return size % 2 == 0 ? (readExpressions[mid] + readExpressions[mid-1]) / 2 : readExpressions[mid];
}


////////////////////////////////////////////

double BlockGroup::q1(){
	vector<Block> blocksList = getBlocks();
	vector<Block>::iterator blockIterator;
	vector<double> readExpressions;
	for(blockIterator = blocksList.begin(); blockIterator != blocksList.end(); ++blockIterator ){
		vector<Read> readsList = blockIterator->getReads();
		vector<Read>::iterator readIterator;
		for ( readIterator = readsList.begin(); readIterator != readsList.end(); ++readIterator ){
			double readExpression = readIterator->getExpression();
			readExpressions.push_back(readExpression);
		}
		readsList.clear();
	}
	blocksList.clear();

	typedef vector<double>::size_type vec_sz;
	vec_sz size = readExpressions.size();
	if (size == 0){
		cerr << "empty read expressions" << endl;
		exit(0);
	}
	sort(readExpressions.begin(), readExpressions.end());
	vec_sz mid = size/2;


	double quant_one;
	if (size % 2 == 0){
		if ((mid-1)/2 == 0)
			quant_one = (readExpressions[(mid)/2] + readExpressions[((mid)/2)-1])/2;
		else
			quant_one = readExpressions[(mid-1)/2];
	}
	else
		quant_one = readExpressions[mid/2];

	return quant_one;
}


////////////////////////////////////////////

double BlockGroup::computeEntropy(map<double, unsigned long> readInfoMap, float numberOfReads){
	double entropy = 0;
	for( map<double,unsigned long>::iterator ii=readInfoMap.begin(); ii!=readInfoMap.end(); ++ii){
		unsigned int infoFrequency = ii->second;
		double infoProbability = infoFrequency/numberOfReads;
		entropy += infoProbability*(log2(infoProbability));
	}
	entropy *= -1;

	readInfoMap.clear();

	return entropy;
}
