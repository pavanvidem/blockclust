/*
 * GspanCreator.h
 *
 *  Created on: Nov 29, 2012
 *      Author: videmp
 */

#ifndef GSPANCREATOR_H_
#define GSPANCREATOR_H_
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
#include <vector>
#include <map>

using namespace std;

class GspanCreator{
public:
	GspanCreator();
	virtual ~GspanCreator();
	void createGspanFile(string configuration, string featuresFile, string outputDir);
	vector<string> split(const string& delim, const string& str);
private:
};

#endif /* GSPANCREATOR_H_ */

