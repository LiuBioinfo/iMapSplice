// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: SAM files
// output: SJ  
// (chr_name pos_donerEnd pos_acceptorStart SJ_name 
//	support_num strand max_seg_left max_seg_right
//	SJ_type flankString)

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

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <OutputSJ> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}

	string inputIndexStr = argv[1]; inputIndexStr += "/";
	string parameter_file = inputIndexStr += "_parameter";
	ifstream parameter_file_ifs(parameter_file.c_str());

	string outputSJstr = argv[2];
	//ofstream sj_ofs(outputSJstr.c_str());

	vector<string> inputSAMfileVec;
	for(int tmp = 3; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	AlignmentToJunc_supportNum_Info* align2juncInfo = new AlignmentToJunc_supportNum_Info();
	int chromNum = indexInfo->returnChromNum();
	align2juncInfo->initiateAlignmentToJuncInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	// insert SJ into SJmap
	align2juncInfo->insertJuncFromAlignmentFileVec(inputSAMfileVec, indexInfo);
	cout << "start to output SJ map" << endl;
	// output SJmap
	align2juncInfo->outputSJmapVec(outputSJstr, indexInfo, true);

	delete align2juncInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_file_ifs.close();
	return 0;
}