/*
 * FeatureComputer.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef FEATURECOMPUTER_H_
#define FEATURECOMPUTER_H_
#include "BlockClust.h"
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <string.h>
#include "BlockGroup.h"
#include "Block.h"
#include "BlockEdge.h"

using namespace std;

class FeatureComputer{
	public:
		FeatureComputer();
		virtual ~FeatureComputer();
		void parseBlockBusterOutput(string bboFile);
		void buildBlockEdges();
		map<string, map<int, string> > parseBinBoundaries(string configFile);
		vector<string> getBinBoundaries(string feat);
		void writeFeatures(string outDir);
		string discretize(vector<string>, double);
		string join( const string,  const vector<string>);
		void init(string bboFile, string configFile, string outDir, unsigned short nrOfBins);
		void init(string bboFile, string outDir, unsigned short nrOfBins);
		void computeBinBoundaries(unsigned short nrOfBins, string outDir);
		void sort(vector<double> &v);
		vector<double> splice(vector<double> &v, unsigned int offset, unsigned int length);
		void writeBounds(unsigned short nrOfBins, unsigned short totalBins, double previousBinMax, vector<double> v, ofstream &BB);
		map<int, string> getBlockGroupClassMap();
		vector<BlockGroup> getBlockGroupList();
		map<string, map<int, string> > getBinBoundariesMap();
		vector<string> split(const string& delim, const string& str);
		void writeToFile(string path, string content);

		typedef double (BlockGroup::*blockGroupFunction)();
		typedef map<string, blockGroupFunction> blockGroupFunctionMap;
		blockGroupFunctionMap bg_featFunctionMap;

		typedef double (Block::*blockFunction)();
		typedef map<string, blockFunction> blockFunctionMap;
		blockFunctionMap b_featFunctionMap;

		typedef double (BlockEdge::*blockEdgeFunction)();
		typedef map<string, blockEdgeFunction> blockEdgeFunctionMap;
		blockEdgeFunctionMap be_featFunctionMap;

	private:
		map<string, map<int, string> > binBoundariesMap;
		vector<BlockGroup> blockGroupsList;
		map <int, string> blockGroupClassMap;
};

#endif /* FEATURECOMPUTER_H_ */
