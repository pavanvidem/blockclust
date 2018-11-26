/*
 * BlockClust.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#include "BlockClust.h"
#include <iostream>
#include <getopt.h>
using namespace std;

BlockClust::BlockClust(){
}

BlockClust::~BlockClust(){
}

void print_usage() {
    printf("Efficient clustering and classification of non-coding RNAs from short read RNA-seq profiles\n");
    printf("-------------------------------------------------------------------------------------------\n");
    std::cout <<
        "Usage: blockclust \n"
            "       -i, --in       [blockbuster output]\n"
            "       -a, --accept   [accept annotations]\n"
            "       -r, --reject   [reject annotations]\n"
            "       -c, --config   [config file]\n"
            "       -o, --out      [output dir]\n"
            "       --help     Show help\n";
    exit(1);
}

int main (int argc, char* argv[]){
    int opt = 0;
    string bboFile, configFile, acceptBED, rejectBED, outDir;
    int index = 0;

    static struct option long_options[] = {
        {"in",      required_argument,  0,  'i' },
        {"accept",  required_argument,  0,  'a' },
        {"reject",  required_argument,  0,  'r' },
        {"config",  required_argument,  0,  'c' },
        {"out",     required_argument,  0,  'o' },
        {"help",    no_argument,        0,  'h' },
        {0,         0,                  0,  0   }
    };

    while ((opt = getopt_long(argc, argv, "a:c:h:i:o:r:",
                   long_options, &index )) != -1) {
        switch (opt) {
            case 'a':
                acceptBED = optarg;
                break;
            case 'r':
                rejectBED = optarg;
                break;
            case 'o':
                outDir = optarg;
                break;
            case 'i':
                bboFile = optarg;
                break;
            case 'c':
                configFile = optarg;
                break;
            case '?':
            case 'h':
            default:
                print_usage();
                break;
        }
    }

    if (argc == 1)
        print_usage();

    for (index = optind; index < argc; index++){
        printf ("Non-option argument %s\n", argv[index]);
        exit(1);
    }

     if (bboFile.empty() || configFile.empty() || acceptBED.empty() || rejectBED.empty() || outDir.empty()){
    	cerr << "options '-i', '-c', '-a', '-r' and '-o' are mandatory in ANALYSIS mode" << endl;
        cerr << "Use blockclust --help for more information" << endl;
    	exit(1);
    }

    printf ("Blockbuster output	: %s\nConfiguration		: %s\nAccept annotations	: %s\nReject annotations	: %s\nOutput directory	: %s\n\n",
    		bboFile.c_str(), configFile.c_str(), acceptBED.c_str(), rejectBED.c_str(), outDir.c_str());

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
    string annotatedBBOFile = ba.annotateBlockGroups(acceptBED, rejectBED, bboFile, outDir);
    FeatureComputer fc;
    fc.init(annotatedBBOFile,configFile,outDir,nrOfBins);
    GspanCreator gc;
    stringstream s;
    s << blockGroupFeatures << " " << blockFeatures << " " << blockEdgeFeatures;
    string configuration = s.str();
	gc.createGspanFile(configuration, outDir + "/discretized.feat", outDir);

	s.str("");
	distance = 2*radius+1;
	cout << "Radius: " << radius << "\n" << "Distance: " << distance << endl;
	s << "EDeN -i "<< outDir <<"/discretized.gspan -f SEQUENCE -M "<<sequenceDegree<< " -b 15 -a MATRIX -r "<<radius<<" -d "<< distance << " -g DIRECTED -y "<<outDir<<" >> "<<outDir<<"/EDeN.log";
	system(s.str().c_str()); //EDeN call
	s.str("");
    return 0;
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
