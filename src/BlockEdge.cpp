/*
 * BlockEdge.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "BlockEdge.h"

BlockEdge::BlockEdge(){
}

BlockEdge::BlockEdge(Block lBlock, Block rBlock){
	leftBlock = lBlock;
	rightBlock = rBlock;
}

BlockEdge::~BlockEdge(){
}

Block BlockEdge::getLeftBlock(){ return leftBlock; }
Block BlockEdge::getRightBlock(){ return rightBlock; }

/////////////////////////////////////////

double BlockEdge::blockContiguity(){
	unsigned long leftBlockStart = getLeftBlock().getStart();
	unsigned long rightBlockStart = getRightBlock().getStart();
	unsigned long leftBlockEnd = getLeftBlock().getEnd();
	unsigned long rightBlockEnd = getRightBlock().getEnd();


	double edgeLength = 0;
	long overlapLength = 0;
	double percentageOverlap = 0;

	if(leftBlockStart <= rightBlockStart){
		if(rightBlockStart < leftBlockEnd){
			if(leftBlockEnd < rightBlockEnd){
				edgeLength = rightBlockEnd-leftBlockStart+1;
				overlapLength = leftBlockEnd-rightBlockStart+1;
			}
			else{
				edgeLength = leftBlockEnd -leftBlockStart+1;
				overlapLength = rightBlockEnd-rightBlockStart+1;
			}
		}
		else{
			edgeLength = rightBlockEnd-leftBlockStart+1;
			overlapLength = leftBlockEnd-rightBlockStart+1;
		}
	}
	else{
		if(leftBlockStart < rightBlockEnd){
			if(rightBlockEnd < leftBlockEnd){
				edgeLength = leftBlockEnd-rightBlockStart+1;
				overlapLength = rightBlockEnd-leftBlockStart+1;
			}
			else{
				edgeLength = rightBlockEnd -rightBlockStart+1;
				overlapLength = leftBlockEnd-leftBlockStart+1;
			}
		}
		else{
			edgeLength = leftBlockEnd-rightBlockStart+1;
			overlapLength = rightBlockEnd-leftBlockStart+1;
		}
	}
	percentageOverlap = (overlapLength/edgeLength)*100;
	return percentageOverlap;
}

/////////////////////////////////////////

double BlockEdge::medianExpressionDiff(){
	double leftBlockMedianExpression = getLeftBlock().medianReadExpression();
	double rightBlockMedianExpression = getRightBlock().medianReadExpression();
	double medianExpressionDiff =  leftBlockMedianExpression > rightBlockMedianExpression ? leftBlockMedianExpression - rightBlockMedianExpression : rightBlockMedianExpression - leftBlockMedianExpression;
	return medianExpressionDiff;
}
