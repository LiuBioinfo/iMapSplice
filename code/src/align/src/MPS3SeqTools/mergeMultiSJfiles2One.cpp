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

#include "general/read_block_test.h"
#include "general/index_info.h"
#include "general/splice_info.h"
#include "general/alignmentToJunc_supportNum.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 6)
	{
		//cout << "Executable <InputIndexInfo> <SJ_file_1> <SJ_file_2> <offset> <outputCompareFile> <min_intron_size> <max_intron_size>" << endl; // SJ_file_1 is the ground truth
		cout << "Executable <InputIndexInfo> <min_intron_size> <max_intron_size> <outputSJfile> <inputSJ_1> (<inputSJ_2> ...)" << endl;
		exit(1);
	}
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string min_intron_size_str = argv[2];
	string max_intron_size_str = argv[3];
	int min_intron_size = atoi(min_intron_size_str.c_str());
	int max_intron_size = atoi(max_intron_size_str.c_str());
	string mergedSJfile = argv[4];
	//ofstream mergedSJfile_ofs(mergedSJfile.c_str());
	AlignmentToJunc_supportNum_Info* mergedJuncInfo = new AlignmentToJunc_supportNum_Info();
	mergedJuncInfo->initiateAlignmentToJuncInfo(chromNum);
	//int juncInfo_merged_validSJnum, juncInfo_merged_invalidSJnum;
	string mergedSJfile_report = mergedSJfile + ".report.log";
	ofstream mergedSJfile_report_ofs(mergedSJfile_report.c_str());
	for(int tmp = 5; tmp < argc; tmp++)
	{
		string tmpSJfile = argv[tmp];
		int juncInfo_tmp_validSJnum, juncInfo_tmp_inValidSJnum;
		string invalidJuncFile_tmp = mergedSJfile + ".invalidSJ." + int_to_str(tmp-4);
		mergedJuncInfo->getAlignmentToJuncInfoFromSJfile(tmpSJfile, indexInfo,
			min_intron_size, max_intron_size, juncInfo_tmp_validSJnum, juncInfo_tmp_inValidSJnum,
			invalidJuncFile_tmp);
		mergedSJfile_report_ofs << "compareSJ: " << tmp-4 << endl 
			<< "validSJ: " << juncInfo_tmp_validSJnum << endl
			<< "invalidSJ: " << juncInfo_tmp_inValidSJnum << endl;
		cout << "getAlignmentToJuncInfoFrom " << tmpSJfile << endl;
		mergedSJfile_report_ofs << "getAlignmentToJuncInfoFrom " << tmpSJfile << endl;
	}

	mergedJuncInfo->outputSJmapVec(mergedSJfile, indexInfo, true);

	delete mergedJuncInfo;
	mergedSJfile_report_ofs.close();
	delete indexInfo;
	//mergedSJfile_ofs.close();
	parameter_ifs.close();
}	