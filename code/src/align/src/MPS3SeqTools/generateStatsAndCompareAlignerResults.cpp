// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//incorporate 7 individual programs into one
// for stats of alignments, junctions and compare them with ground truth.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
//#include <omp.h>
#include "../general/index_info.h"
#include "../general/read_block_test.h"
#include "../general/bwtmap_info.h"
#include "../general/DoubleAnchorScore.h"
#include "../general/sbndm.h"
#include "../general/splice_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/alignInferJunctionHash_info_vec.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 17)
	{
		cout << "Executable inputIndex outputFolder threadNum dataType('R'or'S') SJ_offset SJ_min_size SJ_max_size " << endl;
		cout << "SJsupportNum_max readNum(not pairNum!) readLength_max " << endl;
		cout << "inputGroundTruth/Annotation_name inputGroundTruth/Annotation_SJ_file " << endl;
		cout << "AlignerName_1 SAM_1 AlignerName_2 SAM_2 (... AlignerName_N, SAM_N)" << endl;
		exit(1)
	}

	string outputFolder = argv[2];
	outputFolder += "/";
	string log_file_str = outputFolder + "log";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	ofstream log_ofs(log_file_str.c_str());

	bool realData_simulatedData_bool;
	string realData_simulatedData_Str = argv[4];
	if(realData_simulatedData_Str == "R")
		realData_simulatedData_bool = true;
	else if(realData_simulatedData_bool == "S")
		realData_simulatedData_bool = false;
	else
	{
		cout << "Please specify the dataType: R for RealData, S for Simulated Data" << endl;
		exit(1);
	}

	string threadNumStr = argv[3];
	string offsetStr = argv[5];
	string SJ_min_size_str = argv[6];
	string SJ_max_size_str = argv[7];
	string SJsupportNum_max_str = argv[8];
	string readNum_max_str = argv[9];
	string readLength_max_str = argv[10];
	int threadNum = atoi(threadNumStr.c_str());
	int offset = atoi(offsetStr.c_str());
	int SJ_min_size = atoi(SJ_min_size_str.c_str());
	int SJ_max_size = atoi(SJ_max_size_str.c_str());
	int SJsupportNum_max = atoi(SJsupportNum_max_str.c_str());
	int readLength_max = atoi(readLength_max_str.c_str());
	int readNum_max = atoi(readNum_max_str.c_str());
	string groundTruthName = argv[11];
	string groundTruthSJfileStr = argv[12];
	for(int tmp = 13; tmp < argc; tmp++)
	{
		
	}
}
