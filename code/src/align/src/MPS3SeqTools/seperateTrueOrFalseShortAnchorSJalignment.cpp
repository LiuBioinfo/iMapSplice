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
#include "shortAnchorSJalignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexInfo> <inputGroundTruthSJ> <outputFolder> <inputSAM> " << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string outputFolder = argv[3];
	string inputSAMpath = argv[4];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";
	string output_correctShortAnchorSJ_sam_file = outputFolder + "correctShortAnchorSJ.sam";
	string output_incorrectShortAnchorSJ_sam_file = outputFolder + "incorrectShortAnchorSJ.sam";
	string output_correctShortAnchorSJ_seqExtension_file = outputFolder + "correctShortAnchorSJ_seqExtension.txt";
	string output_incorrectShortAnchorSJ_seqExtension_file = outputFolder + "incorrectShortAnchorSJ_seqExtension.txt";

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

	AlignmentToJunc_supportNum_Info* align2juncInfo = new AlignmentToJunc_supportNum_Info();
	int chromNum = indexInfo->returnChromNum();
	align2juncInfo->initiateAlignmentToJuncInfo(chromNum);	
	string inputSJpath = argv[2];
	string invalidSJpath = outputFolder + "invalidSJ";
	cout << "getAlignmentToJuncInfoFromSJfile file 1 ..." << endl;
	int min_intron_size_int = 50;
	int max_intron_size_int = 300000;
	int validSJnum, invalidSJnum;
	align2juncInfo->getAlignmentToJuncInfoFromSJfile(
		inputSJpath, indexInfo, 
		min_intron_size_int, max_intron_size_int, 
		validSJnum, invalidSJnum, invalidSJpath);
	ShortAnchorSJalignment_Info shortAnchorSJalignmentInfo;
	shortAnchorSJalignmentInfo.extractShortAnchorSJ_seqExtension(
		inputSAMpath, indexInfo, align2juncInfo,
		output_correctShortAnchorSJ_sam_file,
		output_incorrectShortAnchorSJ_sam_file,
		output_correctShortAnchorSJ_seqExtension_file,
		output_incorrectShortAnchorSJ_seqExtension_file);
	shortAnchorSJalignmentInfo;
	delete align2juncInfo;
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}