/*
 * Block.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef BLOCK_H_
#define BLOCK_H_
#include <iostream>
#include <vector>
#include <stdio.h>
#include "Read.h"
#include <map>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <numeric>
using namespace std;

class Block{
	public:
		Block();
		Block(int blkId, unsigned long blkStart, unsigned long blkEnd, double blkExpression);
		virtual ~Block();
		unsigned long getStart();
		unsigned long getEnd();
		double getLength();
		double getExpression();
		double getArea();
		vector<Read> getReads();
		void setReads(vector<Read>);

		double readLengthEntropy();
		double readExpressionEntropy();
		double minReadLength();
		double noOfMultiMapped();
		double medianReadExpression();

	private:
    	int blockId;
    	unsigned long blockStart;
    	unsigned long blockEnd;
    	double blockExpression;
    	vector<Read> reads;
};

#endif /* BLOCK_H_ */
