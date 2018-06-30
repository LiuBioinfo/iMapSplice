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
	//string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	//char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	//cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	//indexInfo->readGenome(chrom);
	//cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	//cout << "start to load every chromosome" << endl;
	//indexInfo->initiate();	
	//cout << "start to initiate chrNameIndexArray" << endl;
	//indexInfo->initiateChrNameIndexArray(1000);
	//cout << "finish loading chromosomes" << endl;

	// string inputSJalignInferHash = argv[2];
	// cout << "start to initaite alignInferJunctionHashInfo " << endl;
	// AlignInferJunctionHash_Info* alignInferJunctionHashInfo
	// 	= new AlignInferJunctionHash_Info();
	// int chromNum = indexInfo->returnChromNum();
	// alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	// cout << "start to insert SJ into SJmap" << endl;
	// alignInferJunctionHashInfo->insertJuncFromJuncFile(
	// 	inputSJalignInferHash, indexInfo);
	
	cout << "start to creat output folder ......" << endl;
	string inputSAMfilePath = argv[2];
	ifstream sam_ifs(inputSAMfilePath.c_str());

	string outputFolder = argv[3];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";	

	string compatibleSamPath = outputFolder + "compatiblePair.sam";
	string incompatibleSamPath = outputFolder + "incompatiblePair.sam";
	ofstream compatibleSam_ofs(compatibleSamPath.c_str());
	ofstream incompatibleSam_ofs(incompatibleSamPath.c_str());

	RefineAlignment_Info refineAlignmentInfo;
	refineAlignmentInfo.detectIncompatiblePairedAlignment(
		sam_ifs,
		compatibleSam_ofs,
		incompatibleSam_ofs, indexInfo);

	compatibleSam_ofs.close();
	incompatibleSam_ofs.close();
	//delete alignInferJunctionHashInfo;
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}