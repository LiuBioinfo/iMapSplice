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
#include "../general/bioEventAlignInferHash_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <OutputBioEventAlignInferHashFile> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}

	int maxReadBaseNumInPathStructure = 30;

	string inputIndexStr = argv[1]; inputIndexStr += "/";
	string parameter_file = inputIndexStr += "_parameter";
	ifstream parameter_file_ifs(parameter_file.c_str());

	string outputBioEventAlignInferHashFileStr = argv[2];

	vector<string> inputSAMfileVec;
	for(int tmp = 3; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	cout << "start to initaite bioEventAlignInferHashInfo " << endl;

	BioEventAlignInferHash_Info* bioEventAlignInferHashInfo = new BioEventAlignInferHash_Info();
	int chromNum = indexInfo->returnChromNum();
	bioEventAlignInferHashInfo->initiateBioEventAlignInferHashInfo(chromNum);
	cout << "start to insert bioEventAlignInfer into bioEventAlignInferHashInfo" << endl;
	// insert bioEventAlignInfer into bioEventAlignInferHashInfo
	bioEventAlignInferHashInfo->insertBioEventAlignInferInfoFromAlignmentFileVec(
		inputSAMfileVec, indexInfo, maxReadBaseNumInPathStructure);
	cout << "start to output bioEventAlignInferHashInfo" << endl;
	
	// output bioEventAlignInferHashInfo
	bioEventAlignInferHashInfo->outputBioEventAlignInferHashInfo(
		indexInfo, outputBioEventAlignInferHashFileStr);

	delete bioEventAlignInferHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_file_ifs.close();
	return 0;
}
