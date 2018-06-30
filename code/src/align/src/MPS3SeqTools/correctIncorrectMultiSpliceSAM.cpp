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
//#include <hash_map>
#include <map>
#include <set>

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
#include "../general/alignmentToJunc_supportNum.h"
#include "../general/fixSingleAnchorNWDP_info.h"
#include "../general/refineAlignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable indexFolderPath inputSAMfilePath outputFolderPath" << endl;
		exit(1);
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

	string inputSJalignInferHash = argv[2];
	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo
		= new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	alignInferJunctionHashInfo->insertJuncFromJuncFile(
		inputSJalignInferHash, indexInfo);
	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);
	cout << "start to getAlter spliceSites" << endl;
	alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity(SJhashInfo, alterSpliceSiteOffset, indexInfo);
	
	cout << "start to creat output folder ......" << endl;
	string inputSAMfilePath = argv[3];
	string outputFolderPath = argv[4];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";	

	string outputFinalSamPath = outputFolder + "final.sam";
	string outputExtensionAtSpliceSite = outputFolder + "extensionAtSpliceSite.sam";
	string outputExtensionAtSpliceSite_ori = outputFolder + "extensionAtSpliceSite_ori.sam";
	string outputAlterSpliceSite = outputFolder + "alterSpliceSite.sam";
	string outputAlterSpliceSite_ori = outputFolder + "alterSpliceSite_ori.sam"

	ofstream finalSam_ofs(outputFinalSamPath.c_str());
	ofstream outputExtensionAtSpliceSite_ofs(outputExtensionAtSpliceSite.c_str());
	ofstream outputExtensionAtSpliceSite_ori_ofs(outputExtensionAtSpliceSite_ori.c_str());
	ofstream outputAlterSpliceSite_ofs(outputAlterSpliceSite.c_str());
	ofstream outputAlterSpliceSite_ori_ofs(outputAlterSpliceSite_ori.c_str());


	finalSam_ofs.close();
	ofstream outputExtensionAtSpliceSite_ofs.close();
	ofstream outputExtensionAtSpliceSite_ori_ofs.close();
	ofstream outputAlterSpliceSite_ofs.close();
	ofstream outputAlterSpliceSite_ori_ofs.close();
	delete SJhashInfo;
	delete alignInferJunctionHashInfo;
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}