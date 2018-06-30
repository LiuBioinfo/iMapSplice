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

#include "../../../general/index_info.h"
#include "general/starChimericJunc_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputFusionChimercOutputPath outputFolderPath" << endl;
		exit(1);
	}
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;

	string outputFolderPath = argv[3];
	cout << "creating results folder ...." << endl;
	outputFolderPath += "/";
	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());	

	string inputFusionChimercOutputPath = argv[2];
	StarChimericJuncHash_Info* starChimericJuncHashInfo 
		= new StarChimericJuncHash_Info();
	starChimericJuncHashInfo->initiate(chromNum);	
	starChimericJuncHashInfo->generateFromStarChimericJuncOutputFile(
		inputFusionChimercOutputPath, indexInfo);

	string output_fusionJunc_includingEncompassedJunc 
		= outputFolderPath + "junc_includingEncompassed.txt";
	string output_fusionJunc_onlySpannedJunc
		= outputFolderPath + "junc_onlySpanned.txt";
	starChimericJuncHashInfo->output_includingEncompassed(
		output_fusionJunc_includingEncompassedJunc, indexInfo);
	starChimericJuncHashInfo->output_onlySpanned(
		output_fusionJunc_onlySpannedJunc, indexInfo);

	starChimericJuncHashInfo->memoryFree();
	delete starChimericJuncHashInfo;
	delete indexInfo;
	return 0;
}