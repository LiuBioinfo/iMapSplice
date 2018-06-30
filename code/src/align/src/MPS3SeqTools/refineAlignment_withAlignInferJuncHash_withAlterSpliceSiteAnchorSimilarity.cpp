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
#include "../general/alignInferJunctionHash_info.h"
#include "../general/fixSingleAnchorNWDP_info.h"
#include "../general/fixDoubleAnchorMatch_info.h"
#include "../general/refineAlignment_info.h"
//#include "shortAnchorSJalignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexInfo> <inputGroundTruthSJ> <inputSAM> <outputFolder>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputJuncFile = argv[2];
	string inputSAMpath = argv[3];
	string outputFolder = argv[4];

	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";
	string output_final_sam_file = outputFolder + "output_final.sam";
	string output_corrected_sam_file_original = outputFolder + "corrected.original.sam";
	string output_corrected_sam_file = outputFolder + "corrected.sam";
	string output_log_file = outputFolder + "process.log";
	ofstream log_ofs(output_log_file.c_str());

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
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfo->insertJuncFromJuncFile_withAlterSpliceSiteAnchorSimilarity(
		inputJuncFile, indexInfo);
	cout << "end of get alignmentToJunc from SJfile" << endl;


	cout << "start to correct original alignment file " << endl;
	RefineAlignment_Info refineAlignmentInfo;
	refineAlignmentInfo.refineAlignment(
		inputSAMpath, output_final_sam_file,
		output_corrected_sam_file, output_corrected_sam_file_original,
		alignInferJunctionHashInfo, indexInfo);

	delete alignInferJunctionHashInfo;
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}