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
	if(argc != 8)
	{
		cout << "Executable <InputIndexFolderPath> <SJfromAnnotationFile> <ToCompareJuncFile> ";
		cout << "<offset> <SJsizeMin> <SJmaxMin> <outputFolder>" << endl;
		exit(1);
	}
	string outputFolderStr = argv[7];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "SJfromAnnotationFile: " << argv[2] << endl;
	log_ofs << "ToCompareJuncFile: " << argv[3] << endl;
	log_ofs << "offset: " << argv[4] << endl;
	log_ofs << "SJsizeMin: " << argv[5] << endl;
	log_ofs << "SJmaxMin: " << argv[6] << endl;
	log_ofs << "outputFolder: " << argv[7] << endl;

	cout << "defining output files ......" << endl;
	log_ofs << endl << "defining output files ......" << endl;
	string outputCompareFileStr_total = outputFolderStr + "total.compare";
	ofstream compare_ofs(outputCompareFileStr_total.c_str());
	string foundSJ_inJuncFile_1 = outputFolderStr + "recalled_annotatedSJ.junc";
	string unfoundSJ_inJuncFile_1 = outputFolderStr + "missed_annotatedSJ.junc";
	string foundSJ_inJuncFile_2 = outputFolderStr + "correct.junc";
	string unfoundSJ_inJuncFile_2 = outputFolderStr + "incorrect.junc";

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	cout << "assigning parameters ......" << endl;
	log_ofs << "assigning parameters ......" << endl;
	string offsetStr = argv[4];
	int offset = atoi(offsetStr.c_str());
	log_ofs << "offset: " << offset << endl;
	string min_intron_size_str = argv[5];
	int min_intron_size_int = atoi(min_intron_size_str.c_str());
	log_ofs << "min_intron_size_int: " << min_intron_size_int << endl;
	string max_intron_size_str = argv[6];
	int max_intron_size_int = atoi(max_intron_size_str.c_str());	
	log_ofs << "max_intron_size_int: " << max_intron_size_int << endl;
	cout << "generating 2 sam2alignInferJuncHash" << endl;
	log_ofs << "generating 2 sam2alignInferJuncHash" << endl;
	string juncHash_file_1_str = argv[2];
	string juncHash_file_2_str = argv[3];
	cout << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	log_ofs << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_1 = new AlignInferJunctionHash_Info();
	AlignInferJunctionHash_Info* juncHash_2 = new AlignInferJunctionHash_Info();
	juncHash_1->initiateAlignInferJunctionInfo(chromNum);
	juncHash_2->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to read 2 juncfiles ...." << endl;
	log_ofs << "start to read 2 juncfiles ...." << endl;
	juncHash_1->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_1_str, indexInfo);
	juncHash_2->insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(juncHash_file_2_str, indexInfo);
	cout << "start to compare junctions" << endl;
	log_ofs << "start to compare junctions" << endl;
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
	juncHash_1->countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise_outputChrNamePosOnly(
		juncHash_2, offset, min_intron_size_int, max_intron_size_int,
		juncNum_total_1, juncNum_validIntronSize_1, juncNum_invalidIntronSize_1,
		foundJuncNum_in1stJuncHash_withinOffset, unfoundJuncNum_in1stJuncHash_withinOffset,
		foundSJ_inJuncFile_1, unfoundSJ_inJuncFile_1, indexInfo);
	juncHash_2->countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise_outputSupNumAnchorSize(
		juncHash_1, offset, min_intron_size_int, max_intron_size_int,
		juncNum_total_2, juncNum_validIntronSize_2, juncNum_invalidIntronSize_2,
		foundJuncNum_in2ndJuncHash_withinOffset, unfoundJuncNum_in2ndJuncHash_withinOffset,
		foundSJ_inJuncFile_2, unfoundSJ_inJuncFile_2, indexInfo);
	cout << "start to output junction comparing results " << endl;
	log_ofs << "start to output junction comparing results " << endl;
	compare_ofs << "juncNum_total_1: " << juncNum_total_1 << endl;
	compare_ofs << "juncNum_validIntronSize_1: " << juncNum_validIntronSize_1 << endl;
	compare_ofs << "juncNum_invalidIntronSize_1: " << juncNum_invalidIntronSize_1 << endl;
	compare_ofs << endl;
	compare_ofs << "juncNum_total_2: " << juncNum_total_2 << endl;
	compare_ofs << "juncNum_validIntronSize_2: " << juncNum_validIntronSize_2 << endl;
	compare_ofs << "juncNum_invalidIntronSize_2: " << juncNum_invalidIntronSize_2 << endl;
	compare_ofs << endl;
	compare_ofs << "foundJuncNum_in1stJuncHash_withinOffset: " << foundJuncNum_in1stJuncHash_withinOffset << endl;
	compare_ofs << "unffoundJuncNum_in1stJuncHash_withinOffset: " << unfoundJuncNum_in1stJuncHash_withinOffset << endl;
	compare_ofs << endl;
	compare_ofs << "foundJuncNum_in2ndJuncHash_withinOffset: " << foundJuncNum_in2ndJuncHash_withinOffset << endl;
	compare_ofs << "unfoundJuncNum_in2ndJuncHash_withinOffset: " << unfoundJuncNum_in2ndJuncHash_withinOffset << endl;
	compare_ofs << endl;
	double sensitivity = ((double)(foundJuncNum_in1stJuncHash_withinOffset)/(double)(juncNum_validIntronSize_1)) * 100;
	double specificity = ((double)(foundJuncNum_in2ndJuncHash_withinOffset)/(double)(juncNum_validIntronSize_2)) * 100;
	compare_ofs << "sensitivity: " << sensitivity << " " << foundJuncNum_in1stJuncHash_withinOffset << "/" << juncNum_validIntronSize_1 << endl;
	compare_ofs << "specificity: " << specificity << " " << foundJuncNum_in2ndJuncHash_withinOffset << "/" << juncNum_validIntronSize_2 << endl;		
	cout << "All jobs done ! " << endl;
	log_ofs << "All jobs done ! " << endl;	
	delete juncHash_1;
	delete juncHash_2;
	delete indexInfo;
	compare_ofs.close();
	parameter_ifs.close();
	log_ofs.close();
}	