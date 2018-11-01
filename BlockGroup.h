/*
 * BlockGroup.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef BLOCKGROUP_H_
#define BLOCKGROUP_H_
#include <string>
#include <stdlib.h>
#include <sstream>
#include <cstdio>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include "Block.h"
#include "BlockEdge.h"
using namespace std;

class BlockGroup{
	public:
		BlockGroup();
		virtual ~BlockGroup();
		unsigned int getId();
		void setId(unsigned int);
		unsigned long getStart();
		void setStart(unsigned long);
		unsigned long getEnd();
		void setEnd(unsigned long);
		double getLength();
		double getExpression();
		void setExpression(double);
		unsigned int getReadCount();
		void setReadCount(unsigned int);
		string getClass();
		void setClass(string);
		vector<Block> getBlocks();
		void setBlocks(vector<Block>);
		vector<BlockEdge> getBlockEdges();
		double getNumberOfBlocks();
		
		void setBlockEdges(vector<BlockEdge>);
		double startPositionEntropy();
		double endPositionEntropy();
		double readLengthEntropy();
		double medianReadExpression();
		double q1();
		double computeEntropy(map<double, unsigned long>, float);

	private:
		unsigned int blockGroupId;
		unsigned long blockGroupStart;
		unsigned long blockGroupEnd;
		double blockGroupExpression;
		unsigned int blockGroupReadCount;
		string blockGroupClass;
		vector<Block> blocks;
		vector<BlockEdge> blockEdges;
};

#endif /* BLOCKGROUP_H_ */
