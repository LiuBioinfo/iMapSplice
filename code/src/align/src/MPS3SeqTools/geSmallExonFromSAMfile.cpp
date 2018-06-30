// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/smallExon_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <OutputSmallExon> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}

	string inputIndexStr = argv[1]; inputIndexStr += "/";
	string parameter_file = inputIndexStr += "_parameter";
	ifstream parameter_file_ifs(parameter_file.c_str());

	string outputSmallExonstr = argv[2];

	vector<string> inputSAMfileVec;
	for(int tmp = 3; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	cout << "start to initaite alignInferJunctionHashInfo " << endl;

	SmallExonHash_Info* smallExonHashInfo = new SmallExonHash_Info();
	//cout << "output initial first..." << endl;
	int chromNum = indexInfo->returnChromNum();
	smallExonHashInfo->initiateSmallExonHashInfo(chromNum);
	//cout << "output initial ..." << endl;
	cout << "start to insert SmallExon into SmallExonHash" << endl;
	// insert SmallExon into SmallExonHash
	smallExonHashInfo->insertSmallExonFromAlignmentFileVec(inputSAMfileVec, indexInfo);
	cout << "start to output smallExonHashInfo" << endl;
	// output SmallExonHash info
	smallExonHashInfo->outputSmallExonHashInfo(indexInfo, outputSmallExonstr);

	delete smallExonHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_file_ifs.close();
	return 0;
}
