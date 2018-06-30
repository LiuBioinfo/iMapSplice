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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "./../general/decoyAlignHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexInfoPath inputSimulatedTransReadAlignResults outputFilePrefix" << endl;
		exit(1);
	}
	string inputSimulatedTransReadAlignResults = argv[2];
	string outputFilePrefix = argv[3];
	string outputFile_decoySNPmapPos = outputFilePrefix + "_decoy_pos.txt";
	string outputFile_decoyAlignment = outputFilePrefix + "_decoy_alignment.txt";

	cout << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	//parameter_ifs.close();
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	parameter_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;	

	DecoyAlignHash_Info tmpDecoyAlignHashInfo;
	cout << "start to do initiate_decoyAlignInfoIndexMapVec_decoyAlignAreaMapVec" << endl;
	tmpDecoyAlignHashInfo.initiate_decoyAlignInfoIndexMapVec_decoyAlignAreaMapVec(chromNum);
	cout << "start to do generateDecoyAlignHash" << endl;
	tmpDecoyAlignHashInfo.generateDecoyAlignHash_fromSimulatedTransReadAlignResults(
		inputSimulatedTransReadAlignResults, indexInfo);
	cout << "start to do output decoySNPcorrespondingMapPos" << endl;
	tmpDecoyAlignHashInfo.output_decoySNPcorrespondingMapPos(outputFile_decoySNPmapPos, indexInfo);
	cout << "start to do output decoyAlignment" << endl;
	tmpDecoyAlignHashInfo.output_decoyAlignment(outputFile_decoyAlignment, indexInfo);
	cout << "all jobs done ......" << endl;
	delete indexInfo;
	return 0;
}