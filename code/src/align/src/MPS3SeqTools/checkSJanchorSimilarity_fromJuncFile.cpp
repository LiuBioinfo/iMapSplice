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
#include "../general/alignInferJunctionHash_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <InputSJ> <OutputSJ_anchorSimilarity>" << endl;
		exit(1);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	int maxReadBaseNumInPathStructure = 30;
	string outputSJanchorSimilarityStr = argv[3];
	string outputSJanchorSimilarityStr_stats = outputSJanchorSimilarityStr + ".stats";
	string inputSJfileStr = argv[2];

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);

	cout << "start to insert SJ into SJmap" << endl;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_onlyChrNamePos(inputSJfileStr, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to getAlter spliceSites" << endl;
	int alterSpliceSiteOffset = 0;
	alignInferJunctionHashInfo->compareExtensionAnchorSimilarity(indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to output SJ map" << endl;	
	// output SJmap
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_onlyExtensionAnchorSimilarity(
		indexInfo, outputSJanchorSimilarityStr);

	alignInferJunctionHashInfo->outputExtensionSpliceSiteAnchorSimilarity(
		outputSJanchorSimilarityStr_stats);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	delete alignInferJunctionHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}