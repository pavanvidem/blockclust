/*
 * FeatureComputer.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "FeatureComputer.h"


FeatureComputer::FeatureComputer(){
    bg_featFunctionMap["SPE"] = &BlockGroup::startPositionEntropy;
    bg_featFunctionMap["EPE"] = &BlockGroup::endPositionEntropy;
    bg_featFunctionMap["RLE"] = &BlockGroup::readLengthEntropy;
    bg_featFunctionMap["GMDX"] = &BlockGroup::medianReadExpression;
    bg_featFunctionMap["Q1"] = &BlockGroup::q1;

    b_featFunctionMap["BRLE"] = &Block::readLengthEntropy;
    b_featFunctionMap["BREE"] = &Block::readExpressionEntropy;
    b_featFunctionMap["BL"] = &Block::getLength;
    b_featFunctionMap["MINRL"] = &Block::minReadLength;
    b_featFunctionMap["BMM"] = &Block::noOfMultiMapped;

    be_featFunctionMap["BC"] = &BlockEdge::blockContiguity;
    be_featFunctionMap["MDXD"] = &BlockEdge::medianExpressionDiff;
}

FeatureComputer::~FeatureComputer(){
}

/////////////////////////////////////////////

map<int, string> FeatureComputer::getBlockGroupClassMap(){
	return blockGroupClassMap;
}
vector<BlockGroup> FeatureComputer::getBlockGroupList(){
	return blockGroupsList;
}
map<string, map<int, string> > FeatureComputer::getBinBoundariesMap(){
	return binBoundariesMap;
}

/////////////////////////////////////////////

void FeatureComputer::parseBlockBusterOutput(string bboFile){
	string LINE;
	ifstream BBO;
	vector<Block> blocksList;
	vector<Read> readsList;
	BlockGroup blockGroup;
	stringstream bgClass;
	unsigned long blockStart, blockEnd;
	unsigned int currentBlockId = 0;
	double blockExpression;
	BBO.open (bboFile.c_str());
	if(!BBO.is_open()){
	    cout << "Error opening file '" << bboFile << "'!!!" << endl;
	    exit(1);
	}
	int count=1;
	char blockGroupChrom[50], blockGroupStrand[5], readChrom[50], readStrand[5], readId[100], annotation[200]="", annotationOverlap[200] = "", blockGroupClass[200] = "";
	while(!BBO.eof()){
		getline(BBO, LINE);
		if(LINE == "")
			continue;
		unsigned int  blockGroupId,  readCount, blockCount, blockId;
		unsigned long  blockGroupStart, blockGroupEnd, readStart, readEnd;
		double blockGroupExpression, readExpression;
		if(LINE.substr(0,1) == ">"){
			if(currentBlockId > 0){
				Block block(currentBlockId, blockStart, blockEnd, blockExpression);
				block.setReads(readsList);
				blocksList.push_back(block);
				blockGroup.setBlocks(blocksList);
				blockGroupsList.push_back(blockGroup);
			}
			sscanf(LINE.c_str(), ">cluster_%d %s %lu %lu %s %lf %d %d %s %s %s", &blockGroupId, blockGroupChrom, &blockGroupStart, &blockGroupEnd, blockGroupStrand, &blockGroupExpression, &readCount, &blockCount, annotation, blockGroupClass, annotationOverlap);
			blockGroup.setId(blockGroupId);
			blockGroup.setStart(blockGroupStart);
			blockGroup.setEnd(blockGroupEnd);
			blockGroup.setExpression(blockGroupExpression);
			blockGroup.setReadCount(readCount);
			blockGroup.setClass(blockGroupClass);
			blocksList.clear();
			readsList.clear();
			blockStart = 4294967295;
			blockEnd = 0;
			currentBlockId = 1;
			blockExpression = 0;
			bgClass << blockGroupClass << ":" << count;
			blockGroupClassMap[count] = bgClass.str();
			bgClass.str("");
			count++;
		}
		else{
			vector<string> fields = split("\t", LINE);
			sscanf(fields[0].c_str(),"%s", readChrom);
			sscanf(fields[1].c_str(),"%lu", &readStart);
			sscanf(fields[2].c_str(),"%lu", &readEnd);
			sscanf(fields[3].c_str(),"%s", readId);
			sscanf(fields[4].c_str(),"%lf", &readExpression);
			sscanf(fields[5].c_str(),"%s", readStrand);
			sscanf(fields[6].c_str(),"%d", &blockId);

			vector<string> readIdFields = split("|", readId);

			double noOfMappings = atof(readIdFields.at(readIdFields.size()-1).c_str());
			double readCount = atof(readIdFields.at(readIdFields.size()-2).c_str());
			if(blockId == currentBlockId){
				Read read(readId, readStart, readEnd, readExpression, noOfMappings, readCount);
				readsList.push_back(read);
			}
			else{
				Block block(currentBlockId, blockStart, blockEnd, blockExpression);
				block.setReads(readsList);
				blocksList.push_back(block);
				blockGroup.setBlocks(blocksList);
				readsList.clear();
				blockStart = 4294967295;
				blockEnd = 0;
				blockExpression = 0;
				Read read(readId, readStart, readEnd, readExpression, noOfMappings, readCount);
				readsList.push_back(read);
			}
			if(readStart < blockStart){
				blockStart = readStart;
			}
			if(blockEnd < readEnd){
				blockEnd = readEnd;
			}
			currentBlockId = blockId;
			blockExpression += readExpression;
		}
	}
	// create last block of final block group
	Block block(currentBlockId, blockStart, blockEnd, blockExpression);
	block.setReads(readsList);
	blocksList.push_back(block);
	blockGroup.setBlocks(blocksList);
	blockGroupsList.push_back(blockGroup);
	BBO.close();
	readsList.clear();
	blocksList.clear();
	return;
}


/////////////////////////////////////////////

/* Test mode */
void FeatureComputer::init(string bboFile, string configFile, string outDir, unsigned short nrOfBins){
	cout << "parsing Blockbuster output file: " << bboFile << endl;
	parseBlockBusterOutput(bboFile);
	cout << "Finished\n";
	buildBlockEdges();
	cout << "Built edges\n";
    cout << "parsing config file...\n";
	binBoundariesMap = parseBinBoundaries(configFile);
	cout << "Done \n";
    cout << "Writing features...\n";
	writeFeatures(outDir);
	cout << "Done \n";
    return;
}

/////////////////////////////////////////////

/* Train mode */
void FeatureComputer::init(string bboFile, string outDir, unsigned short nrOfBins){
	stringstream bb;
	bb << outDir << "/bin_boundaries." << nrOfBins;
	cout << "parsing  Blockbuster output file: " << bboFile << endl;
	parseBlockBusterOutput(bboFile);
	cout << "Finished\n";
	buildBlockEdges();
	cout << "Built edges\n";
    cout << "computing bin boundaries...\n";
    computeBinBoundaries(nrOfBins, bb.str().c_str());
    cout << "Done \n";
    cout << "parsing bin boundaries file...\n";
    binBoundariesMap = parseBinBoundaries(bb.str().c_str());
	cout << "Done \n";
	cout << "Writing features...\n";
    writeFeatures(outDir);
    cout << "Done \n";
    return;
}

/////////////////////////////////////////////

void FeatureComputer::computeBinBoundaries(unsigned short nrOfBins, string bbFile){
	ofstream BB;
	BB.open(bbFile.c_str());
	if(!BB.is_open()){
	    cerr << "Error opening file '" << bbFile << "'!!!" << endl;
	    exit(1);
	}

	BB << "Bin_boundaries:" << endl;

	for(blockGroupFunctionMap::iterator gg=bg_featFunctionMap.begin(); gg!=bg_featFunctionMap.end(); ++gg){
		vector<double> featureVales;
		for(unsigned int i=0; i<blockGroupsList.size(); i++){
			double featValue = (blockGroupsList[i].*(gg->second))();
			featureVales.push_back(featValue);
		}
		sort(featureVales);
		BB << ">" << gg->first << endl;
		writeBounds(nrOfBins, nrOfBins, -INFINITY, featureVales, BB);
	}

	for(blockFunctionMap::iterator bb=b_featFunctionMap.begin(); bb!=b_featFunctionMap.end(); ++bb){
		vector<double> featureVales;
		for(unsigned int g=0; g<blockGroupsList.size(); g++){
			vector<Block> blocks = blockGroupsList[g].getBlocks();
			for(unsigned int b=0; b<blocks.size(); b++){
				double featValue = (blocks[b].*(bb->second))();
				featureVales.push_back(featValue);
			}
			sort(featureVales);
		}
		BB << ">" << bb->first << endl;
		writeBounds(nrOfBins, nrOfBins, -INFINITY, featureVales, BB);
	}

	for(blockEdgeFunctionMap::iterator ee=be_featFunctionMap.begin(); ee!=be_featFunctionMap.end(); ++ee){
		vector<double> featureVales;
		for(unsigned int g=0; g<blockGroupsList.size(); g++){
			vector<BlockEdge> blockEdges = blockGroupsList[g].getBlockEdges();
			for(unsigned int e=0; e<blockEdges.size(); e++){
				double featValue = (blockEdges[e].*(ee->second))();
				featureVales.push_back(featValue);
			}
			sort(featureVales);
		}
		BB << ">" << ee->first << endl;
		writeBounds(nrOfBins, nrOfBins, -INFINITY, featureVales, BB);
	}
}

/////////////////////////////////////////////

void FeatureComputer::buildBlockEdges(){
	for(unsigned int i=0; i<blockGroupsList.size(); i++){
		vector<Block> blocks = blockGroupsList[i].getBlocks();
		vector<BlockEdge> blockEdges;
		for (unsigned int j=0; j<blocks.size()-1; j++){
				BlockEdge blockEdge(blocks[j], blocks[j+1]);
				blockEdges.push_back(blockEdge);
		}
		blockGroupsList[i].setBlockEdges(blockEdges);
	}
}

/////////////////////////////////////////////

map<string, map<int, string> > FeatureComputer::parseBinBoundaries(string configFile){
	map<string, map<int, string> > binBoundariesMap;
	string LINE;
	ifstream BB;
	BB.open (configFile.c_str());
	if(!BB.is_open()){
	    cerr << "Error opening file '" << configFile << "'!!!" << endl;
	    exit(1);
	}
	string currentFeature;
	bool readBoundaries=false;
	while(!BB.eof()){
		getline(BB, LINE);
		if(readBoundaries){
			if(LINE.substr(0,1) == ">")
				currentFeature = LINE.erase (0,1);
			else if(!LINE.empty()){
				int bin; char start[120], end[120];
				sscanf(LINE.c_str(), "%d\t%s\t%s", &bin, start, end);
				binBoundariesMap[currentFeature][bin] = string(start) + "\t" + string(end);
			}
		}
		if(LINE == "Bin_boundaries:"){
			readBoundaries = true;
		}
	}
	BB.close();
	return binBoundariesMap;
}

/////////////////////////////////////////////

vector<string> FeatureComputer::getBinBoundaries(string feat){
	vector<string> boundaries;
	for(unsigned int i=1; i<=binBoundariesMap[feat].size(); i++){
		boundaries.push_back(binBoundariesMap[feat][i]);
	}
	return boundaries;
}

/////////////////////////////////////////////

string FeatureComputer::discretize(vector<string> binBoundaries, double non_discretized){
	int discretized;
	int binCount = 1;
	for (unsigned int i=0; i< binBoundaries.size(); i++){
		double start_f, end_f;
		if(string::npos != binBoundaries[i].find("-inf\t")){
			sscanf(binBoundaries[i].c_str(), "-inf %lf", &end_f);
			start_f = -INFINITY;
			if(non_discretized >= start_f && non_discretized <= end_f){
				discretized = binCount;
			}
		}
		else if(string::npos != binBoundaries[i].find("\tinf")){
			sscanf(binBoundaries[i].c_str(), "%lf inf", &start_f);
			end_f = INFINITY;
			if(non_discretized > start_f && non_discretized <= end_f){
				discretized = binCount;
			}
		}
		else{
			sscanf(binBoundaries[i].c_str(), "%lf %lf", &start_f, &end_f);
			if(non_discretized > start_f && non_discretized <= end_f){
				discretized = binCount;
			}
		}
		binCount++;
	}
	stringstream d;
	d << discretized;
	return d.str();
}

/////////////////////////////////////////////

void FeatureComputer::writeBounds(unsigned short nrOfBins, unsigned short totalBins, double previousBinMax, vector<double> v, ofstream &BB){
    float totalInstances = v.size();
    float sizeOfBin = totalInstances/nrOfBins;
	sizeOfBin = round(sizeOfBin);
	vector<double> currentBin = splice(v, 0, (unsigned int)sizeOfBin);
	vector<double> leftBin = v;

	double currentBinMax = currentBin[currentBin.size()-1];
	int moved = 0;
	for (unsigned int i=0; i<leftBin.size(); i++){
		if(leftBin[i] == currentBinMax){
			currentBin.push_back(leftBin[i]);
			moved++;
		}
	}
	vector<double> newLeftBin = splice(leftBin, moved, leftBin.size()-moved);
	unsigned int binCount = totalBins - nrOfBins + 1;
	if(nrOfBins == 2){
		if(totalBins == 2){
			double n_inf = -INFINITY;
			BB << binCount << "\t" << n_inf << "\t" << currentBinMax << endl;
		}
		else{
			BB << binCount << "\t" << previousBinMax << "\t" << currentBinMax << endl;
			previousBinMax = currentBinMax;
		}
		double p_inf = INFINITY;
		binCount++;
		if(previousBinMax != -INFINITY){
			BB << binCount <<"\t" << previousBinMax << "\t" << p_inf << endl;
		}
	}
	else if(totalBins == nrOfBins){
		double n_inf = -INFINITY;
		BB << binCount << "\t" << n_inf << "\t"<< currentBinMax << endl;
		nrOfBins--;
		previousBinMax = currentBinMax;
		writeBounds(nrOfBins, totalBins, previousBinMax, newLeftBin, BB);
	}
	else{
		BB << binCount << "\t" << previousBinMax << "\t" << currentBinMax << endl;
		nrOfBins--;
		previousBinMax = currentBinMax;
		writeBounds(nrOfBins, totalBins, previousBinMax, newLeftBin, BB);
	}
}


/////////////////////////////////////////////

void FeatureComputer::writeFeatures(string outDir){
    stringstream bgAnnotationsFileContent;
    stringstream ndfFileContent;
    stringstream dfFileContent;

    // Header line
	ndfFileContent << "blockgroup";
	dfFileContent << "blockgroup";
	for(blockGroupFunctionMap::iterator gg=bg_featFunctionMap.begin(); gg!=bg_featFunctionMap.end(); ++gg){
		ndfFileContent << "\t" << gg->first;
		dfFileContent << "\t" << gg->first;
	}
	for(blockFunctionMap::iterator bb=b_featFunctionMap.begin(); bb!=b_featFunctionMap.end(); ++bb){
		ndfFileContent << "\t" << bb->first;
		dfFileContent << "\t" << bb->first;
	}
	for(blockEdgeFunctionMap::iterator ee=be_featFunctionMap.begin(); ee!=be_featFunctionMap.end(); ++ee){
		ndfFileContent << "\t" << ee->first;
		dfFileContent << "\t" << ee->first;
	}
	ndfFileContent << "\n";
	dfFileContent << "\n";

	for(unsigned int g=0; g<blockGroupsList.size(); g++){
		vector<Block> blocks = blockGroupsList[g].getBlocks();
		vector<BlockEdge> blockEdges = blockGroupsList[g].getBlockEdges();

		bgAnnotationsFileContent << "blockgroup_" << blockGroupsList[g].getId() << "\t" << blockGroupsList[g].getClass() << endl;
		ndfFileContent	<< "blockgroup_" << blockGroupsList[g].getId();
		dfFileContent << "blockgroup_" << blockGroupsList[g].getId();

    	// whole graph (blockgroup) features //
		for(blockGroupFunctionMap::iterator gg=bg_featFunctionMap.begin(); gg!=bg_featFunctionMap.end(); ++gg){
			vector<string> featureBinBoudaries = getBinBoundaries(gg->first);
			double featValue = (blockGroupsList[g].*(gg->second))();
			string d_featValue = discretize(featureBinBoudaries, featValue);
			ndfFileContent << "\t" << featValue;
			dfFileContent << "\t" << d_featValue;
		}

		// graph vertex (block) features //
		for(blockFunctionMap::iterator bb=b_featFunctionMap.begin(); bb!=b_featFunctionMap.end(); ++bb){
			ndfFileContent << "\t";
			dfFileContent << "\t";
			stringstream s; vector<string> blockDiscretizedFeat, blockFeat;
			for(unsigned int b=0; b<blocks.size(); b++){
				vector<string> featureBinBoudaries = getBinBoundaries(bb->first);
				double featValue = (blocks[b].*(bb->second))();
				string d_featValue = discretize(featureBinBoudaries, featValue);
				s << featValue;
				blockFeat.push_back(s.str());
				s.str("");
				blockDiscretizedFeat.push_back(d_featValue);
			}
			string d_feat = join (",", blockDiscretizedFeat);
			string nd_feat = join (",", blockFeat);
			ndfFileContent << nd_feat;
			dfFileContent << d_feat;
		}

		// graph edges (between blocks) features //
		for(blockEdgeFunctionMap::iterator ee=be_featFunctionMap.begin(); ee!=be_featFunctionMap.end(); ++ee){
			ndfFileContent << "\t";
			dfFileContent << "\t";
			stringstream s; vector<string> blockDiscretizedFeat, blockFeat;
			for(unsigned int e=0; e<blockEdges.size(); e++){
				vector<string> featureBinBoudaries = getBinBoundaries(ee->first);
				double featValue = (blockEdges[e].*(ee->second))();
				string d_featValue = discretize(featureBinBoudaries, featValue);
				s << featValue;
				blockFeat.push_back(s.str());
				s.str("");
				blockDiscretizedFeat.push_back(d_featValue);
			}
			string d_feat = join (",", blockDiscretizedFeat);
			string nd_feat = join (",", blockFeat);
			ndfFileContent << nd_feat;
			dfFileContent << d_feat;
		}
		ndfFileContent << endl;
		dfFileContent << endl;
	}
	writeToFile(outDir + "/blockgroup_annotations.txt", bgAnnotationsFileContent.str());
	writeToFile(outDir + "/non_discretized.feat", ndfFileContent.str());
	writeToFile(outDir + "/discretized.feat", dfFileContent.str());
	return;
}

/////////////////////////////////////////////

void FeatureComputer::writeToFile(string path, string content){
	ofstream FH;
	FH.open (path.c_str());
	if(!FH.is_open()){
	    cerr << "Error opening file '" << path << "'!!!" << endl;
	    exit(1);
	}
	FH << content;
	FH.flush(); FH.close(); FH.clear();
	return;
}

/////////////////////////////////////////////

vector<string> FeatureComputer::split(const string& delim, const string& str){
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

/////////////////////////////////////////////

vector<double> FeatureComputer::splice(vector<double> &v, unsigned int offset, unsigned int length){
	vector<double> slicedVector;
	for(unsigned int i=offset; i<offset+length; i++){
		slicedVector.push_back(v[i]);
	}
	v.erase(v.begin()+offset, v.begin()+offset+length);
	return slicedVector;
}

/////////////////////////////////////////////

string FeatureComputer::join(const string delim, const vector<string> v){
    stringstream s;
    for(unsigned int i=0; i<v.size(); i++){
    	if(i != 0)
    		s << delim;
    	s << v[i];
    }
    return s.str();
}

/////////////////////////////////////////////

void FeatureComputer::sort(vector<double> &v){
     int i, j;
     double temp;
     int vSize = v.size();
     for (i=0; i< vSize; i++){  // element to be compared
          for(j = (i+1); j < vSize; j++){   // rest of the elements
                if (v[i] > v[j]){          // descending order
					temp= v[i];          // swap
					v[i] = v[j];
					v[j] = temp;
               }
          }
     }
     return;
}
