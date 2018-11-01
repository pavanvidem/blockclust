/*
 * Read.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "Read.h"

Read::Read(){
	readId  = "";
	readStart = 0;
	readEnd = 0;
	readExpression = 0;
	numberOfMappings = 0;
	readCount = 0;
}

Read::Read(string rdId, unsigned long rdStart, unsigned long rdEnd, double rdExpression, double noOfMappings, double rdCount){
	readId = rdId;
	readStart = rdStart;
	readEnd = rdEnd;
	readExpression = rdExpression;
	numberOfMappings = noOfMappings;
	readCount = rdCount;
}

Read::~Read(){
	// TODO Auto-generated destructor stub
}

string Read::getId(){ return readId; }
unsigned long Read::getStart(){ return readStart; }
unsigned long Read::getEnd(){ return readEnd; }
unsigned short Read::getLength(){ return readEnd - readStart +1; }
double Read::getExpression(){ return readExpression; }
double Read::getNumberOfMappings(){ return numberOfMappings; }
double Read::getReadCount(){ return readCount; }
double Read::getArea(){ return getLength()*getExpression(); }


