/*
 * BlockClust.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef BLOCKCLUST_H_
#define BLOCKCLUST_H_
#include <string>
#include <sstream>
#include "algorithm"
#include <iostream>
#include <time.h>
#include "BlockGroup.h"
#include "FeatureComputer.h"
#include "GspanCreator.h"
#include "BlockGroupAnnotator.h"
using namespace std;

class BlockClust{
public:
	BlockClust();
	virtual ~BlockClust();
	void aucroc(const char* similarityMatrix, map <int, string> blockGroupClassMap, ofstream &ROC, string config);
	vector<string> split(const string& delim, const string& str);
private:
};

#endif /* BLOCKCLUST_H_ */
