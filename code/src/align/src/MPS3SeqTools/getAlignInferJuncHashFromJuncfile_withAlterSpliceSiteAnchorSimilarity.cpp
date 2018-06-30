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

int main(int argc, char**argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <OutputSJ> <inputJunc_1>" << endl;
		exit(1);
	}
	//string alterSpliceSiteOffsetStr = argv[4];
	//int alterSpliceSiteOffset = atoi(alterSpliceSiteOffsetStr.c_str());

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

	string outputSJstr = argv[2];
	//ofstream sj_ofs(outputSJstr.c_str());

	//vector<string> inputJUNCfileVec;
	//for(int tmp = 3; tmp < argc; tmp++)
	//{
	//	string tmpJUNCstr = argv[tmp];
	//	inputJUNCfileVec.push_back(tmpJUNCstr);
	//}
	string inputJuncFile = argv[3];

	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfo->insertJuncFromJuncFile_withAlterSpliceSiteAnchorSimilarity(
		inputJuncFile, indexInfo);
	
	// cout << "start to initaite SJhashInfo" << endl;
	// SJhash_Info* SJhashInfo = new SJhash_Info();
	// SJhashInfo->initiateAreaAndStringHash(chromNum);
	// cout << "start to convert 2 SJhashInfo" << endl;
	// alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);

	// cout << "start to getAlter spliceSites" << endl;

	// alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity(
	// 	SJhashInfo, alterSpliceSiteOffset, indexInfo);
	// cout << "start to output SJ map" << endl;
	// output SJmap
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, outputSJstr);

	delete alignInferJunctionHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}
