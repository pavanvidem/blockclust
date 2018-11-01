/*
 * BlockClust.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "BlockClust.h"
#include <iostream>
using namespace std;

BlockClust::BlockClust(){
}

BlockClust::~BlockClust(){
}

int main (int argc, char **argv){
	string testFile, configFile, outDir, rejectBED, acceptBED;
    unsigned short bitSize = 0;
    bool reportAUC=false;
    int index;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "c:o:a:b:t:q:r:")) != -1)
        switch (c){
        case 'a':
        	acceptBED = optarg;
            break;
        case 'r':
        	rejectBED = optarg;
            break;
        case 'o':
        	outDir = optarg;
            break;
        case 'b':
        	sscanf(optarg, "%hu", &bitSize);
            break;
        case 't':
        	testFile = optarg;
            break;
        case 'q':
            reportAUC = true;
            break;  
        case 'c':
            configFile = optarg;
            break;
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
            return 1;
		default:
            break;
	}
    for (index = optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);

    if(testFile.empty() || configFile.empty() || acceptBED.empty() || rejectBED.empty() || outDir.empty()){
    	fprintf (stderr, "options '-t', '-c', '-a', '-r' and '-o' are mandatory in clustering mode\n");
    	return 1;
    }

    printf ("Blockbuster output	: %s\nConfiguration		: %s\nAccept annotations	: %s\nReject annotations	: %s\nOutput directory	: %s\n\n",
    		testFile.c_str(), configFile.c_str(), acceptBED.c_str(), rejectBED.c_str(), outDir.c_str());

    unsigned short nrOfBins = 0, radius = 0;
    char blockGroupFeatures[20] = "";
    char blockFeatures[50]= "";
    char blockEdgeFeatures[50]= "";
	string LINE;
	ifstream BB;
	BB.open(configFile.c_str());
	if(!BB.is_open()){
	    cerr << "Error opening file '" << configFile << "'!!!" << endl;
	    exit(1);
	}
	while(!BB.eof()){
		getline(BB, LINE);
		sscanf(LINE.c_str(), "Discretization_level: %hu", &nrOfBins);
		sscanf(LINE.c_str(), "Radius: %hu", &radius);
		sscanf(LINE.c_str(), "Feature_combination: %s %s %s", blockGroupFeatures, blockFeatures, blockEdgeFeatures);
	}
	BB.close();
	if (strcmp(blockGroupFeatures, "") == 0 || strcmp(blockFeatures, "") == 0 || strcmp(blockEdgeFeatures, "") == 0 || nrOfBins == 0 || radius == 0){
		cerr << "Please check the configuration file!\n";
		exit(0);
	}

    int count = 0;
    for (int i = 0; i < 50; i++){
    	 if (blockEdgeFeatures[i] == ',')
    		 count++;
		 if(blockFeatures[i] == ',')
			 count++;
    }
    unsigned short sequenceDegree = count + 2;
    unsigned short distance = 2*radius+1;
    BlockGroupAnnotator ba;
    string annotatedBBOFile = ba.annotateBlockGroups(acceptBED, rejectBED, testFile, outDir);
    FeatureComputer fc;
    fc.init(annotatedBBOFile,configFile,outDir,nrOfBins);
    GspanCreator gc;
    stringstream s;
    s << blockGroupFeatures << " " << blockFeatures << " " << blockEdgeFeatures;
    string configuration = s.str();
	gc.createGspanFile(configuration, outDir + "/discretized.feat", outDir);

	s.str("");
	distance = 2*radius+1;
	cout << "Radius: " << radius << "\n" << "Distance: " << distance << "\n" << "Bit size: "<< bitSize << endl;
	s << "EDeN -i "<< outDir <<"/discretized.gspan -f SEQUENCE -M "<<sequenceDegree<< " -b "<< bitSize << " -a MATRIX -r "<<radius<<" -d "<< distance << " -g DIRECTED -y "<<outDir<<" >> "<<outDir<<"/EDeN.log";
	system(s.str().c_str()); //EDeN call
	s.str("");
    if(reportAUC){
	    s << outDir << "/roc_measures_n_" << nrOfBins << "_r_" << radius << ".out";
        ofstream ROC;
	    ROC.open(s.str().c_str());
	    if(!ROC.is_open()){
	        cerr << "Error opening file '" << s.str() << "'!!!" << endl;
	        exit(1);
	    }
	    s.str("");    
	    s << outDir<< "/matrix";
	    BlockClust bc;
	    bc.aucroc(s.str().c_str(), fc.getBlockGroupClassMap(), ROC, configuration);
	    ROC.close();
	}
    return 0;
}


/////////////////////////////////////////////

void BlockClust::aucroc(const char* similarityMatrix, map<int, string> blockGroupClassMap, ofstream &ROC, string config){
	ifstream SIM;
	SIM.open (similarityMatrix);
	if(!SIM.is_open()){
	    cerr << "Error opening file '" << similarityMatrix << "'!!!" << endl;
	    exit(1);
	}
	int count=1;

	map <string, vector<float> > classSpecificAUCs;
	string LINE;
	ofstream TEMP;
	while(!SIM.eof()){
		getline(SIM, LINE);
		if(LINE == "")
			continue;
		vector<string> fields = split(" ", LINE);
		string rowClass = blockGroupClassMap.find(count)->second;
		stringstream scount;
		scount << ":" << count;
		size_t f = rowClass.find(scount.str());
		while(f != string::npos){
			rowClass.replace(f, rowClass.length(), "");
			f = rowClass.find(scount.str());
		}
	    map<float, map<int, string> > similaritiesMap;
		for(unsigned int i=0; i<fields.size(); i++){
			similaritiesMap[atof(fields[i].c_str())][i] = blockGroupClassMap.find(i+1)->second;
		}
		scount.str("");
		map<float, map<int, string> >::iterator it1;
		map<int, string>::iterator it2;
		string similarities = "";
		for(it1 = similaritiesMap.begin(); it1 != similaritiesMap.end(); ++it1){
			int cnt = 1;
			stringstream sim;
			sim << it1->first;
			for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
				string fieldClass = it2->second;
				size_t f = fieldClass.find(":");
				while(f != string::npos){
					fieldClass.replace(f, fieldClass.length(), "");
					f = fieldClass.find(":");
				}
				if(fieldClass == rowClass){
					similarities.append("1").append(" ").append(sim.str()).append("\n");
				}
				else{
					similarities.append("0").append(" ").append(sim.str()).append("\n");
				}
			}
		    cnt++;
			sim.clear();
		}
		stringstream temp;
		temp << similarityMatrix <<"."<< count << ".rocin";
		TEMP.open(temp.str().c_str());
		if(!TEMP.is_open()){
		    cerr << "Error opening file '" << temp.str().c_str() << "'!!!" << endl;
		    exit(1);
		}
		TEMP << similarities;
		TEMP.close();
		stringstream cmd;
		cmd << "perf -ROC < " << temp.str().c_str();
	    FILE* pipe = popen(cmd.str().c_str(), "r");
	    if (!pipe) fprintf (stderr, "ERROR");
	    char buffer[128];
	    std::string result = "";
	    while(!feof(pipe)){
	    	if(fgets(buffer, 128, pipe) != NULL)
	    		result += buffer;
	    }
	    pclose(pipe);
	    cmd.str("");
	    cmd << "rm " << temp.str().c_str();
	    system(cmd.str().c_str());
	    cmd.str("");
	    temp.str("");
	    float aucroc;
	    sscanf(result.c_str(), "ROC\t%f\n", &aucroc);
	    classSpecificAUCs[rowClass].push_back(aucroc);
		count++;
	}
	SIM.close();

	float allClassAUCsum = 0;
	map<string, vector<float> >::iterator it;
	unsigned int knownBGCount = 0;
	unsigned int unknownBGCount = 0;

	ROC << "ncRNA class" << "\t" << "Nr. of ncRNAs" << "\tAUC ROC" << endl;
	cout << "ncRNA class" << "\t" << "Nr. of ncRNAs" << "\tAUC ROC" << endl;
	
	for(it = classSpecificAUCs.begin(); it != classSpecificAUCs.end(); ++it){
	    float classAUCSum = 0;
		string annotation = it->first;
        if(annotation != "unknown"){
		    for (unsigned int i=0; i<it->second.size(); i++){
			    classAUCSum += it->second[i];
			    allClassAUCsum += it->second[i];
			    knownBGCount++;
		    }
		    float classSpecificAvgAUC = classAUCSum/it->second.size();
		    ROC << annotation << "\t\t" << it->second.size() << "\t\t" << classSpecificAvgAUC << endl;
		    cout << annotation << "\t\t" << it->second.size() << "\t\t" << classSpecificAvgAUC << endl;
		}
		else{
		    for (unsigned int i=0; i<it->second.size(); i++){
    		    unknownBGCount++;
		    }
		}
	}
	
	float allClassAvgAUC = allClassAUCsum/knownBGCount;

	ROC << "Average\t\t" << knownBGCount << "\t\t" << allClassAvgAUC << endl;
	cout << "Average\t\t" << knownBGCount << "\t\t" << allClassAvgAUC << endl;
    cout << endl;

//    cout << "unknown" << "\t" << unknownBGCount << endl;
	classSpecificAUCs.clear();
	return;
}

/////////////////////////////////////////////

vector<string> BlockClust::split(const string& delim, const string& str){
	std::size_t start_pos = 0;
	std::size_t match_pos;
	std::size_t substr_length;
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
