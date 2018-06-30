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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable <InputIndexFolderPath> <GroundTruthJuncFile> <ToCompareJuncFile> ";
		cout << " <SJsizeMin> <SJmaxMin> <outputFolder>" << endl;
		exit(1);
	}
	string outputFolderStr = argv[6];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "GroundTruthJuncFile: " << argv[2] << endl;
	log_ofs << "ToCompareJuncFile: " << argv[3] << endl;
	log_ofs << "SJsizeMin: " << argv[4] << endl;
	log_ofs << "SJmaxMin: " << argv[5] << endl;
	log_ofs << "outputFolder: " << argv[6] << endl;

	cout << "defining output files ......" << endl;
	log_ofs << endl << "defining output files ......" << endl;
	string outputCompareFileStr_total = outputFolderStr + "total.compare";
	ofstream compare_ofs(outputCompareFileStr_total.c_str());
	string foundSJ_file = outputFolderStr + "foundJunc.junc";
	string unfoundSJ_file = outputFolderStr + "unfoundSJ.junc";
	ofstream found_ofs(foundSJ_file.c_str());
	ofstream unfound_ofs(unfoundSJ_file.c_str());

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string min_intron_size_str = argv[4];
	int min_intron_size_int = atoi(min_intron_size_str.c_str());
	log_ofs << "min_intron_size_int: " << min_intron_size_int << endl;
	string max_intron_size_str = argv[5];
	int max_intron_size_int = atoi(max_intron_size_str.c_str());	
	log_ofs << "max_intron_size_int: " << max_intron_size_int << endl;
	cout << "generating 2 sam2alignInferJuncHash" << endl;
	log_ofs << "generating 2 sam2alignInferJuncHash" << endl;

	string juncHash_file_1_str = argv[2];
	cout << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	log_ofs << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_1 = new AlignInferJunctionHash_Info();
	juncHash_1->initiateAlignInferJunctionInfo(chromNum);
	juncHash_1->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_1_str, indexInfo);

	int junc_total_num = 0;
	int junc_valid_num = 0;
	int junc_found_num = 0;
	int junc_unfound_num = 0;
	string juncHash_file_2_str = argv[3];
	ifstream junc_2_ofs(juncHash_file_2_str.c_str());
	while(!junc_2_ofs.eof())
	{
		string tmpJuncStr;
		getline(junc_2_ofs, tmpJuncStr);
		if(tmpJuncStr == "")
			break;
		junc_total_num ++;
		int tabLoc_1 = tmpJuncStr.find("\t");
		int tabLoc_2 = tmpJuncStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpJuncStr.find("\t", tabLoc_2 + 1);
		string tmpJunc_chrName = tmpJuncStr.substr(0, tabLoc_1);
		string tmpJunc_startPosStr = tmpJuncStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpJunc_endPosStr = tmpJuncStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_1 - 1);
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		int tmpJunc_distance = tmpJunc_endPos - tmpJunc_startPos - 1;
		if((tmpJunc_distance >= min_intron_size_int)&&(tmpJunc_distance <= max_intron_size_int))
			junc_valid_num ++;
		int tmpIndex = juncHash_1->searchAndReturnAlignInferInfoVecIndex(tmpJunc_chrNameInt,
			tmpJunc_startPos, tmpJunc_endPos);
		if(tmpIndex < 0)
		{	
			junc_unfound_num ++;
			unfound_ofs << tmpJuncStr << endl;
		}
		else
		{
			junc_found_num ++;
			found_ofs << tmpJuncStr << endl;
		}
	}
	compare_ofs << "Total_junc_num: " << junc_total_num << endl;
	compare_ofs << "Valid_junc_num: " << junc_valid_num << endl;
	compare_ofs << "Found_junc_num: " << junc_found_num << endl;
	compare_ofs << "Unfound_junc_num: " << junc_unfound_num << endl;
	junc_2_ofs.close();
	return 0;
}	