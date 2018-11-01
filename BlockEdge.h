/*
 * BlockEdge.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef BLOCKEDGE_H_
#define BLOCKEDGE_H_
#include "Block.h"

class BlockEdge{
	public:
		BlockEdge();
		BlockEdge(Block, Block);
		virtual ~BlockEdge();
		Block getLeftBlock();
		Block getRightBlock();

		double blockContiguity();
		double medianExpressionDiff();

	private:
		Block leftBlock;
		Block rightBlock;
};

#endif /* BLOCKEDGE_H_ */
