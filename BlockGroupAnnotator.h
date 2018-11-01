/*
 * BlockGroupAnnotator.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef BLOCKGROUPANNOTATOR_H_
#define BLOCKGROUPANNOTATOR_H_
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

class BlockGroupAnnotator{
	public:
		BlockGroupAnnotator();
		virtual ~BlockGroupAnnotator();
		void parseAnnotationsBED(string acceptAnnotationsBEDFile, string rejectAnnotationsBEDFile);
		string annotateBlockGroups(string acceptAnnotationsBEDFile, string rejectAnnotationsBEDFile, string bboFile, string outputDir);
		void parseBlockbusterOut(string bboFile);
		void writeAnnotatedBlockGroups(string outFile);
		float overlapAvg(string overlapRatios);
		vector<string> split(const string& delim, const string& str);
	private:
		map<string, string> acceptAnnotationsMap;
		map<string, string> rejectAnnotationsMap;
		map<string, vector<string> > blockGroupReadsMap;
		map<string, map<string,string> > blockGroupAnnotationsMap;
		vector<string> rejectBlockGroupsList;
};

#endif /* BLOCKGROUPANNOTATOR_H_ */
