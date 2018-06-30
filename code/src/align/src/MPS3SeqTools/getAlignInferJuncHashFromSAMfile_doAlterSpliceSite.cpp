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

clock_t loadIndex_begin, loadIndex_end, loadIndex_cost = 0,
	getAlignInferJuncHash_begin, getAlignInferJuncHash_end, getAlignInferJuncHash_cost = 0,
	getSJhashInfo_begin, getSJhashInfo_end, getSJhashInfo_cost = 0,
	doExtension_begin, doExtension_end, doExtension_cost = 0,
	checkAlterSpliceSite_begin, checkAlterSpliceSite_end, checkAlterSpliceSite_cost = 0,
	outputSJ_begin, outputSJ_end, outputSJ_cost = 0,
	refineAlignment_begin, refineAlignment_end, refineAlignment_cost = 0;

time_t nowtime;
struct tm *local;


int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <OutputSJ> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	int maxReadBaseNumInPathStructure = 30;
	string outputSJstr = argv[2];

	vector<string> inputSAMfileVec;
	for(int tmp = 3; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}

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
	//cout << "output initial first..." << endl;
	//alignInferJunctionHashInfo->outputAlignInferInfoHashInfo(indexInfo, outputSJstr);
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	//cout << "output initial ..." << endl;
	//alignInferJunctionHashInfo->outputAlignInferInfoHashInfo(indexInfo, outputSJstr);
	cout << "start to insert SJ into SJmap" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfo->insertJuncFromAlignmentFileVec(
		inputSAMfileVec, indexInfo, maxReadBaseNumInPathStructure);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to getAlter spliceSites" << endl;
	int alterSpliceSiteOffset = 0;
	alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity(
		SJhashInfo, alterSpliceSiteOffset, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	cout << "start to output SJ map" << endl;	
	// output SJmap
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
		indexInfo, outputSJstr);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... current time ......" << endl << endl; 

	delete alignInferJunctionHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}
