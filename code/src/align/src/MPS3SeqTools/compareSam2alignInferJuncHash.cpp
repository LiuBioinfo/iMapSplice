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

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable <IndexInput> <juncFile_1> <juncFile_2> <offset> <SJsizeMin> <SJsizeMax> <outputFolder> " << endl;
		exit(1);
	}
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string offsetStr = argv[4];
	int offset = atoi(offsetStr.c_str());
	string min_intron_size_str = argv[5];
	string max_intron_size_str = argv[6];
	int min_intron_size_int = atoi(min_intron_size_str.c_str());
	int max_intron_size_int = atoi(max_intron_size_str.c_str());	
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[7];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputCompareFileStr = outputFolderStr + "juncHash.compare";
	ofstream compare_ofs(outputCompareFileStr.c_str());
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());
	string foundSJ_inJuncFile_1 = outputFolderStr + "foundJunc_inJuncFile_1.junc";
	string unfoundSJ_inJuncFile_1 = outputFolderStr + "unfoundSJ_inJuncFile_1.junc";
	string foundSJ_inJuncFile_2 = outputFolderStr + "foundJunc_inJuncFile_2.junc";
	string unfoundSJ_inJuncFile_2 = outputFolderStr + "unfoundSJ_inJuncFile_2.junc";

	string juncHash_file_1_str = argv[2];
	string juncHash_file_2_str = argv[3];
	cout << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_1 = new AlignInferJunctionHash_Info();
	AlignInferJunctionHash_Info* juncHash_2 = new AlignInferJunctionHash_Info();
	juncHash_1->initiateAlignInferJunctionInfo(chromNum);
	juncHash_2->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to read 2 juncfiles ...." << endl;
	juncHash_1->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_1_str, indexInfo);
	juncHash_2->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_2_str, indexInfo);
	cout << "start to compare junctions" << endl;
	int juncNum_total_1 = 0;
	int juncNum_total_2 = 0;
	int juncNum_validIntronSize_1 = 0;
	int juncNum_validIntronSize_2 = 0;
	int juncNum_invalidIntronSize_1 = 0;
	int juncNum_invalidIntronSize_2 = 0;
	int foundJuncNum_in1stJuncHash_withinOffset = 0;
	int foundJuncNum_in2ndJuncHash_withinOffset = 0;
	int unfoundJuncNum_in1stJuncHash_withinOffset = 0;
	int unfoundJuncNum_in2ndJuncHash_withinOffset = 0;	

	juncHash_1->countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise(
		juncHash_2, offset, min_intron_size_int, max_intron_size_int,
		juncNum_total_1, juncNum_validIntronSize_1, juncNum_invalidIntronSize_1,
		foundJuncNum_in1stJuncHash_withinOffset, unfoundJuncNum_in1stJuncHash_withinOffset,
		foundSJ_inJuncFile_1, unfoundSJ_inJuncFile_1, indexInfo);

	juncHash_2->countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise(
		juncHash_1, offset, min_intron_size_int, max_intron_size_int,
		juncNum_total_2, juncNum_validIntronSize_2, juncNum_invalidIntronSize_2,
		foundJuncNum_in2ndJuncHash_withinOffset, unfoundJuncNum_in2ndJuncHash_withinOffset,
		foundSJ_inJuncFile_2, unfoundSJ_inJuncFile_2, indexInfo);

	compare_ofs << "juncNum_total_1: " << juncNum_total_1 << endl;
	compare_ofs << "juncNum_validIntronSize_1: " << juncNum_validIntronSize_1 << endl;
	compare_ofs << "juncNum_invalidIntronSize_1: " << juncNum_invalidIntronSize_1 << endl;
	compare_ofs << endl;
	compare_ofs << "juncNum_total_2: " << juncNum_total_2 << endl;
	compare_ofs << "juncNum_validIntronSize_2: " << juncNum_validIntronSize_2 << endl;
	compare_ofs << "juncNum_invalidIntronSize_2: " << juncNum_invalidIntronSize_2 << endl;
	compare_ofs << endl;

	compare_ofs << "foundJuncNum_in1stJuncHash_withinOffset: " 
		<< foundJuncNum_in1stJuncHash_withinOffset << endl;
	compare_ofs << "unffoundJuncNum_in1stJuncHash_withinOffset: " 
		<< unfoundJuncNum_in1stJuncHash_withinOffset << endl;
	compare_ofs << endl;
	compare_ofs << "foundJuncNum_in2ndJuncHash_withinOffset: " 
		<< foundJuncNum_in2ndJuncHash_withinOffset << endl;
	compare_ofs << "unfoundJuncNum_in2ndJuncHash_withinOffset: " 
		<< unfoundJuncNum_in2ndJuncHash_withinOffset << endl;
	compare_ofs << endl;

	double sensitivity 
		= ((double)(foundJuncNum_in1stJuncHash_withinOffset)/(double)(juncNum_validIntronSize_1)) * 100;
	double specificity
		= ((double)(foundJuncNum_in2ndJuncHash_withinOffset)/(double)(juncNum_validIntronSize_2)) * 100;

	compare_ofs << "sensitivity: " << sensitivity << "\t" 
		<< foundJuncNum_in1stJuncHash_withinOffset << " / " << juncNum_validIntronSize_1 << endl;
	compare_ofs << "specificity: " << specificity << "\t"
		<< foundJuncNum_in2ndJuncHash_withinOffset << " / " << juncNum_validIntronSize_2 << endl;

	delete juncHash_1;
	delete juncHash_2;

	return 0;
}