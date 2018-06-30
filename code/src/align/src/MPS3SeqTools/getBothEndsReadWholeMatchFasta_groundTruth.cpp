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
#include "../general/splice_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/fixSingleAnchorNWDP_info.h"
#include "../general/fixDoubleAnchorMatch_info.h"
#include "../general/refineAlignment_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <inputSAM> <OutputFastaFilePrefix> " << endl;
		exit(1);
	}
	string inputSAMpath = argv[2];
	ifstream sam_ifs(inputSAMpath.c_str());
	string outputPrefix = argv[3];
	string outputFa_1 = outputPrefix + "_1.fa";
	string outputFa_2 = outputPrefix + "_2.fa";
	ofstream fa_1_ofs(outputFa_1.c_str());
	ofstream fa_2_ofs(outputFa_2.c_str());

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

	RefineAlignment_Info refineAlignmentInfo;
	refineAlignmentInfo.outputWholeMatchPairReadInFasta_groundTruth(
		sam_ifs, fa_1_ofs, fa_2_ofs, indexInfo);

	delete indexInfo;
	sam_ifs.close();
	fa_1_ofs.close();
	fa_2_ofs.close();
	parameter_ifs.close();
	return 0;
}
