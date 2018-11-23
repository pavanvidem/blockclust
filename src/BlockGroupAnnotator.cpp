/*
 * BlockGroupAnnotator.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "BlockGroupAnnotator.h"
#include <getopt.h>
BlockGroupAnnotator::BlockGroupAnnotator(){
}

BlockGroupAnnotator::~BlockGroupAnnotator(){
}


/////////////////////////////////////////////

string BlockGroupAnnotator::annotateBlockGroups(string acceptAnnotationsBEDFile, string rejectAnnotationsBEDFile, string bboFile, string outputDir){
	parseAnnotationsBED(acceptAnnotationsBEDFile, rejectAnnotationsBEDFile);
	parseBlockbusterOut(bboFile.c_str());
	stringstream s;
	s << outputDir << "/" << "annotated.bbo";
	string annotatedBBOFile = s.str();
	writeAnnotatedBlockGroups(annotatedBBOFile.c_str());
	s.str("");
	s << outputDir << "/" << "annotated.bbo";
	annotatedBBOFile = s.str().c_str();
	return annotatedBBOFile;
}


void BlockGroupAnnotator::parseAnnotationsBED(string acceptAnnotationsBEDFile, string rejectAnnotationsBEDFile){
	ifstream ANNOT;
	// Accepted annotations - ncRNAs
	ANNOT.open(acceptAnnotationsBEDFile.c_str());
	if(!ANNOT.is_open()){
	    cerr << "Error opening file '" << acceptAnnotationsBEDFile << "'!!!" << endl;
	    exit(1);
	}
	string LINE;
	while(!ANNOT.eof()){
		getline(ANNOT, LINE);
		char chrom[50], annotation[100], score[50], strand[5];
		unsigned long start, end;
		sscanf(LINE.c_str(), "%s\t%lu\t%lu\t%s\t%s\t%s", chrom, &start, &end, annotation, score, strand);
		stringstream s;
		s << chrom << "\t" << start<< "\t" << end << "\t" << strand;
		acceptAnnotationsMap[s.str()] = annotation;
	}
	ANNOT.close();
	ANNOT.clear();

	// Reject annotations - all other gene annoatations
	ANNOT.open(rejectAnnotationsBEDFile.c_str());
	if(!ANNOT.is_open()){
	    cerr << "Error opening file '" << rejectAnnotationsBEDFile << "'!!!" << endl;
	    exit(1);
	}

	while(!ANNOT.eof()){
		getline(ANNOT, LINE);
		char chrom[50], annotation[100], score[50], strand[5];
		unsigned long start, end;
		sscanf(LINE.c_str(), "%s\t%lu\t%lu\t%s\t%s\t%s", chrom, &start, &end, annotation, score, strand);
		stringstream s;
		s << chrom << "\t" << start<< "\t" << end << "\t" << strand;
		rejectAnnotationsMap[s.str()] = annotation;
	}
	ANNOT.close();
	ANNOT.clear();
}

void BlockGroupAnnotator::parseBlockbusterOut(string bboFile){
	ifstream BBO;
	BBO.open(bboFile.c_str());
	if(!BBO.is_open()){
	    cerr << "Error opening file '" << bboFile << "'!!!" << endl;
	    exit(1);
	}
	string LINE;
	string currentBlockGroup="";
	string blockGroup;
	vector<string> reads;
	map<string, string> annotationBlockGroupsMap;
	unsigned int unknownCount = 0;
	bool firstBlock = false;
	char blockGroupChrom[50], blockGroupStrand[5];
	unsigned int  blockGroupId,  readCount, blockCount;
	unsigned long  blockGroupStart, blockGroupEnd;
	double blockGroupExpression;
	short blockGroupLength = 0;
	while(!BBO.eof()){
		getline(BBO, LINE);
		if(LINE == "")
			continue;
		if(LINE.substr(0,1) == ">"){
			vector<string> fields = split("\t", LINE);
			sscanf(fields[0].c_str(),">cluster_%d", &blockGroupId);
			sscanf(fields[1].c_str(),"%s", blockGroupChrom);
			sscanf(fields[2].c_str(),"%lu", &blockGroupStart);
			sscanf(fields[3].c_str(),"%lu", &blockGroupEnd);
			sscanf(fields[4].c_str(),"%s", blockGroupStrand);
			sscanf(fields[5].c_str(),"%lf", &blockGroupExpression);
			sscanf(fields[6].c_str(),"%d", &readCount);
			sscanf(fields[7].c_str(),"%d", &blockCount);
			blockGroupLength = blockGroupEnd - blockGroupStart + 1;
			// ignore block groups containing only 1 block
			if(blockCount < 2 || blockGroupLength > 200 ||  blockGroupLength < 50) //|| blockGroupExpression < 50
				continue;
			if(firstBlock){
				blockGroupReadsMap[currentBlockGroup] = reads;
				reads.clear();
			}
			firstBlock =true;
			blockGroup = fields[0];
			for (int i=1;i<=7;i++)
				blockGroup.append("\t").append(fields[i]);
			map<string, string>::iterator mapIterator;
			int numberOfOverlaps = 0;
			for(mapIterator = acceptAnnotationsMap.begin(); mapIterator != acceptAnnotationsMap.end(); ++mapIterator){
				char annotationChrom[50], annoationStrand[5];
				unsigned long annotationStart, annotationEnd;
				sscanf(mapIterator->first.c_str(), "%s\t%lu\t%lu\t%s", annotationChrom, &annotationStart, &annotationEnd, annoationStrand);
				int annotationLength = annotationEnd - annotationStart + 1;
				long overlapStart = 0;
				long overlapEnd = 0;

				if(strcmp(annotationChrom,blockGroupChrom))
					continue;

				if(annotationStart <= blockGroupStart && annotationEnd <= blockGroupEnd && blockGroupStart < annotationEnd){
					overlapStart = blockGroupStart;
					overlapEnd = annotationEnd;
				}
				if(annotationStart > blockGroupStart && annotationEnd > blockGroupEnd && annotationStart < blockGroupEnd){
					overlapStart = annotationStart;
					overlapEnd = blockGroupEnd;
				}
				if(annotationStart <= blockGroupStart && annotationEnd >= blockGroupEnd){
					overlapStart = blockGroupStart;
					overlapEnd = blockGroupEnd;
				}
				if(annotationStart >= blockGroupStart && annotationEnd <= blockGroupEnd){
					overlapStart = annotationStart;
					overlapEnd = annotationEnd;
				}

				stringstream currentAnnotation;
				currentAnnotation << annotationChrom<<"\t"<<annotationStart<<"\t"<<annotationEnd<<"\t"<<annoationStrand;

				if(overlapStart == 0 || overlapEnd == 0)
					continue;

				numberOfOverlaps++;
				float overlapLength = overlapEnd - overlapStart + 1;
				float annotationOverlapRatio = overlapLength/annotationLength;
				float blockGroupOverlapRatio = overlapLength/blockGroupLength;
				stringstream currentOverlapRatios;
				currentOverlapRatios << blockGroupOverlapRatio<<":"<<annotationOverlapRatio;
				if(annotationOverlapRatio < 0.7 || blockGroupOverlapRatio < 0.7)
					continue;

				if(blockGroupAnnotationsMap.count(blockGroup) > 0){
					string previousOverlap = blockGroupAnnotationsMap.find(blockGroup)->second.begin()->second;
					float currentOverlapRatiosAvg = overlapAvg(currentOverlapRatios.str());
					float existingOverlapRatiosAvg = overlapAvg(previousOverlap);
					if(currentOverlapRatiosAvg > existingOverlapRatiosAvg){
						blockGroupAnnotationsMap.erase(blockGroup);
						blockGroupAnnotationsMap[blockGroup][currentAnnotation.str()] = currentOverlapRatios.str();
					}
				}
				else
					blockGroupAnnotationsMap[blockGroup][currentAnnotation.str()] = currentOverlapRatios.str();

				if(annotationBlockGroupsMap.count(currentAnnotation.str()) > 0){
					string existingBlockGroup = annotationBlockGroupsMap.find(currentAnnotation.str())->second;

					if(blockGroupAnnotationsMap.find(existingBlockGroup)->second.count(currentAnnotation.str()) >0 ){
						string existingOverlapRatios = blockGroupAnnotationsMap.find(existingBlockGroup)->second.find(currentAnnotation.str())->second;
						float currentOverlapRatiosAvg = overlapAvg(currentOverlapRatios.str());
						float existingOverlapRatiosAvg = overlapAvg(existingOverlapRatios);
						if(currentOverlapRatiosAvg > existingOverlapRatiosAvg){
							annotationBlockGroupsMap[currentAnnotation.str()] = blockGroup;
							rejectBlockGroupsList.push_back(existingBlockGroup);
						}
						else{
							rejectBlockGroupsList.push_back(blockGroup);
						}
					}
				}
				else
					annotationBlockGroupsMap[currentAnnotation.str()] = blockGroup;
			}
			if(numberOfOverlaps == 0){
				for(mapIterator = rejectAnnotationsMap.begin(); mapIterator != rejectAnnotationsMap.end(); ++mapIterator){
					char annotationChrom[50], annoationStrand[5];
					unsigned long annotationStart, annotationEnd;
					sscanf(mapIterator->first.c_str(), "%s\t%lu\t%lu\t%s", annotationChrom, &annotationStart, &annotationEnd, annoationStrand);
					if(!strcmp(annotationChrom,blockGroupChrom)){
						if((annotationStart < blockGroupStart && blockGroupStart < annotationEnd) || (blockGroupStart < annotationStart && annotationStart < blockGroupEnd)){
							numberOfOverlaps++;
							continue;
						}
					}
				}
			}
			if(numberOfOverlaps == 0){
				unknownCount++;
				stringstream s1;
				s1 << "unknown_" << unknownCount;
				stringstream s2;
				s2 << "unknown_"<< unknownCount << ":unknown";
				blockGroupAnnotationsMap[blockGroup][s1.str()] = s2.str();
			}
			currentBlockGroup = blockGroup;
		}
		else{
			if(blockCount < 2 || blockGroupLength > 200 ||  blockGroupLength < 50) //|| blockGroupExpression < 50
				continue;
			reads.push_back(LINE);
		}
	}
	blockGroupReadsMap[currentBlockGroup] = reads;
	BBO.close();
	BBO.clear();
	reads.clear();
	annotationBlockGroupsMap.clear();
	return;
}

float BlockGroupAnnotator::overlapAvg(string overlapRatios){
	float blockGroupOverlapRatio, annotationOverlapRatio;
	sscanf(overlapRatios.c_str(), "%f:%f", &blockGroupOverlapRatio, &annotationOverlapRatio);
	return (blockGroupOverlapRatio+annotationOverlapRatio)/2;
}

void BlockGroupAnnotator::writeAnnotatedBlockGroups(string outFile){
	ofstream OUT;
	OUT.open(outFile.c_str());
	if(!OUT.is_open()){
	    cerr << "Error opening file '" << outFile << "'!!!" << endl;
	    exit(1);
	}
	vector<string> reads;
	map<string, map<string,string> >::iterator it1;
	map<string,string>::iterator it2;
	for(it1 = blockGroupAnnotationsMap.begin(); it1 != blockGroupAnnotationsMap.end(); ++it1){
		if (std::find(rejectBlockGroupsList.begin(), rejectBlockGroupsList.end(), it1->first) != rejectBlockGroupsList.end()){
			continue;
		}
		OUT << it1->first;
		for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
			string annotation = "";
			if (std::string::npos != it2->first.find("unknown"))
				annotation = it2->second;
			else
				annotation = acceptAnnotationsMap.find(it2->first)->second;
			OUT << "\t" << annotation;
			size_t f1 = annotation.find(":");
			size_t f2 = annotation.find(":",f1+1,1);

		  	string family = annotation.substr(f1+1,f2-f1-1);
//			while(f != string::npos){
//				annotation.erase(0, f+1);
//				f = annotation.find(":");
//			}
			OUT << "\t" << family << "\t" << it2->second << endl;
		}
		reads = blockGroupReadsMap.find(it1->first)->second;
		for(unsigned int r=0; r<reads.size(); r++)
			OUT << reads[r] << endl;
		reads.clear();
	}
	OUT.close();
	return;
}

/////////////////////////////////////////////

vector<string> BlockGroupAnnotator::split(const string& delim, const string& str){
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

