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
		cout << "Executable <InputIndexInfo> <inputSAM> <OutputFolder> " << endl;
		exit(1);
	}
	int maxReadBaseNumInPathStructure = 30;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	

	string outputSJpath_highConfidence = outputFolderStr + "highConfidence.alignInferJuncHash";
	string outputSJpath_lowConfidence = outputFolderStr + "lowConfidence.alignInferJuncHash";
	string outputSJpath_lowConfidence_kept = outputFolderStr + "lowConfidence_kept.alignInferJuncHash";
	string outputSJpath_lowConfidence_fitlerOut = outputFolderStr + "lowConfidence_filterOut.alignInferJuncHash";


	vector<string> inputSAMfileVec;
	string tmpSAMstr = argv[2];
	inputSAMfileVec.push_back(tmpSAMstr);

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

	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);
	cout << "start to getAlter spliceSites" << endl;
	int alterSpliceSiteOffset = 0;
	alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity(
		SJhashInfo, alterSpliceSiteOffset, indexInfo);

	cout << "start to output SJ map" << endl;	
	// output SJmap
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, outputSJpath_all);
	cout << "start to output SJ filtering results ..." << endl;
	alignInferJunctionHashInfo->outputSJ_extended(outputSJpath_extended, indexInfo);
	alignInferJunctionHashInfo->outputSJ_nonExtended_nonAlterSpliceSite(outputSJpath_nonExtended_nonAlterSpliceSite, indexInfo);
	alignInferJunctionHashInfo->outputSJ_nonExtended_withAlterSpliceSite(outputSJpath_nonExtended_withAlterSpliceSite, indexInfo);
	//alignInferJunctionHashInfo->outputSJ_nonExtended_withAlterSpliceSite_valid(outputSJpath_nonExtended_withAlterSpliceSite_valid, indexInfo);
	//alignInferJunctionHashInfo->outputSJ_nonExtended_withAlterSpliceSite_invalid(outputSJpath_nonExtended_withAlterSpliceSite_invalid, indexInfo);
	alignInferJunctionHashInfo->outputSJ_kept(outputSJpath_kept, indexInfo);
	alignInferJunctionHashInfo->outputSJ_filterOut(outputSJpath_filterOut, indexInfo);
	alignInferJunctionHashInfo->outputSJ_nonExtended_multiAnchorExactlyTheSame(outputSJpath_nonExtended_multiAnchorExactlyTheSame, indexInfo);
	
	delete alignInferJunctionHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}
