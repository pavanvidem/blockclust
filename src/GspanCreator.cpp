/*
 * GspanCreator.cpp
 *
 *  Created on: Nov 29, 2012
 *      Author: videmp
 */

#include "GspanCreator.h"

GspanCreator::GspanCreator(){
}

GspanCreator::~GspanCreator(){
}

void GspanCreator::createGspanFile(string configuration, string featuresFile, string outputDir){
    char blockGroupFeatures[1000] = "";
    char blockFeatures[1000]= "";
    char blockEdgeFeatures[1000]= "";

    sscanf(configuration.c_str(), "%s %s %s", blockGroupFeatures, blockFeatures, blockEdgeFeatures);
    int count = 0;
    for (int i = 0; i < 1000; i++){
    	 if (blockEdgeFeatures[i] == ',')
    		 count++;
		 if(blockFeatures[i] == ',')
			 count++;
    }
	unsigned short seqenceDegree = count + 2;
    vector<string> bgf = split(",", blockGroupFeatures);
    vector<string> bf = split(",", blockFeatures);
    vector<string> bef = split(",", blockEdgeFeatures);

    ifstream FF;
    FF.open(featuresFile.c_str());
	if(!FF.is_open()){
	    cerr << "Error opening file '" << featuresFile << "'!!!" << endl;
	    exit(1);
	}
    ofstream OUT;
    OUT.open((outputDir + "/discretized.gspan").c_str());
	if(!OUT.is_open()){
	    cerr << "Error opening file '" << outputDir << "/discretized.gspan'!!!" << endl;
	    exit(1);
	}
	string LINE;
	getline(FF, LINE);

	map <unsigned short, string> featureNamesMap;
	vector<string> featureNames = split("\t", LINE);
	for(unsigned int i=0; i<featureNames.size(); i++){
		featureNamesMap[i] = featureNames[i];
	}

	map<unsigned short, string>::iterator it;
	while(!FF.eof()){
		getline(FF, LINE);
		if(!LINE.empty()){
			map <string, vector<string> > featureValuesMap;
			vector<string> featureValues = split("\t", LINE);
			for(unsigned int i=1; i<featureValues.size(); i++){
				featureValuesMap[featureNamesMap.find(i)->second] = split(",", featureValues[i]);
			}
			for(unsigned int i=0; i<bgf.size(); i++){
				OUT << featureValuesMap.find(bgf[i])->second[0];
				for(unsigned short j=1; j<seqenceDegree; j++){
					OUT << "0";
				}
			}

			OUT << "|";
			size_t numberOfBlocks = featureValuesMap.find(bf[0])->second.size();
			for(unsigned int i=0; i<numberOfBlocks; i++){
				for(unsigned int j=0; j<bf.size(); ++j){
					OUT << featureValuesMap.find(bf[j])->second[i];
				}
				for(unsigned int j=0; j<bef.size(); j++){
					if(i < featureValuesMap.find(bef[j])->second.size()){
						OUT << featureValuesMap.find(bef[j])->second[i];
					}
					else
						OUT << "N";
				}
			}
			OUT << endl;
			featureValuesMap.clear();
		}
	}
	FF.close();
	OUT.close();
	featureNamesMap.clear();
	bgf.clear(); bf.clear(); bef.clear();
    return;
}

/////////////////////////////////////////////

vector<string> GspanCreator::split(const string& delim, const string& str){
    size_t start_pos = 0;
    size_t match_pos;
    size_t substr_length;
    vector<string> result;
    while((match_pos = str.find(delim, start_pos)) != string::npos){
    	substr_length = match_pos - start_pos;
        if (substr_length > 0){
            result.push_back(str.substr(start_pos, substr_length));
        }
        start_pos = match_pos + delim.length();
    }
    substr_length = str.length() - start_pos;
    if (substr_length > 0){
        result.push_back(str.substr(start_pos, substr_length));
    }
    return result;
}
