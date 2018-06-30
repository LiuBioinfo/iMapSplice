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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"

using namespace std;

int getSupNumFromJuncStr(string& tmpJuncStr)
{
	int tabLoc_1 = tmpJuncStr.find("\t");
	int tabLoc_2 = tmpJuncStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpJuncStr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpJuncStr.find("\t", tabLoc_3 + 1);
	string supNumStr = tmpJuncStr.substr(tabLoc_4 + 1, tabLoc_4 - tabLoc_3 - 1);
	return atoi(supNumStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 10)
	{
		cout << "Executable <InputIndexFolderPath> <GroundTruthJuncFile> <ToCompareJuncFile> ";
		cout << "<offset> <SJsizeMin> <SJmaxMin> <SupNumThresMax> <outputFolder> <toolName>" << endl;
		exit(1);
	}
	string toolName = argv[9];
	string outputFolderStr = argv[8];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	
	cout << "defining output files ......" << endl;
	log_ofs << "defining output files ......" << endl;
	string outputCompareFileStr_total = outputFolderStr + "total.compare";
	ofstream compare_ofs(outputCompareFileStr_total.c_str());
	string outputCompareFileStr_varySupNumThres_TPR_FPR = outputFolderStr + "varySupNumThres_TPR_FPR.compare";
	ofstream compare_varySupNum_TPR_FPR_ofs(outputCompareFileStr_varySupNumThres_TPR_FPR.c_str());
	string outputCompareFileStr_varySupNumThres_trueNum_falseNum = outputFolderStr + "varySupNumThres_trueNum_falseNum.compare";
	ofstream compare_varySupNum_trueNum_falseNum_ofs(outputCompareFileStr_varySupNumThres_trueNum_falseNum.c_str());
	string outputCompareFileStr_varySupNumThres_raw = outputFolderStr + "varySupNumThres_raw.compare";
	ofstream compare_varySupNum_raw_ofs(outputCompareFileStr_varySupNumThres_raw.c_str());

	string foundSJ_inJuncFile_1 = outputFolderStr + "foundJunc_inJuncFile_1.junc";
	string unfoundSJ_inJuncFile_1 = outputFolderStr + "unfoundSJ_inJuncFile_1.junc";
	string foundSJ_inJuncFile_2 = outputFolderStr + "foundJunc_inJuncFile_2.junc";
	string unfoundSJ_inJuncFile_2 = outputFolderStr + "unfoundSJ_inJuncFile_2.junc";

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
	string supNumThresMax_str = argv[7];
	int sup_num_thres_max = atoi(supNumThresMax_str.c_str());
	log_ofs << "sup_num_thres_max: " << sup_num_thres_max << endl;

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
	juncHash_1->insertJuncFromJuncFile_chrNamePos_supportNum(juncHash_file_1_str, indexInfo);
	juncHash_2->insertJuncFromJuncFile_chrNamePos_supportNum(juncHash_file_2_str, indexInfo);
	
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

	cout << "start to compare junctions with varying supporting num" << endl;
	log_ofs << "start to compare junctions with varying supporting num" << endl;
	vector<int> foundJuncNumVec_in1stJuncHash_withinOffset_varySupNum;
	for(int tmp = 0; tmp < sup_num_thres_max; tmp++)
		foundJuncNumVec_in1stJuncHash_withinOffset_varySupNum.push_back(0);
	juncHash_1->compareWithAnotherAlignInferJuncHash_juncWise_varySupNumMaxInTheOtherJuncHash(
		juncHash_2, offset, min_intron_size_int, max_intron_size_int, sup_num_thres_max,
		foundJuncNumVec_in1stJuncHash_withinOffset_varySupNum, indexInfo);
	
	cout << "start to generate foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum " << endl;
	log_ofs << "start to generate foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum " << endl;
	vector<int> foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum;
	for(int tmp = 0; tmp < sup_num_thres_max; tmp++)
		foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum.push_back(0);
	ifstream foundSJ_inJuncFile_2_ifs(foundSJ_inJuncFile_2.c_str());
	while(!foundSJ_inJuncFile_2_ifs.eof())
	{
		string tmpStr;
		getline(foundSJ_inJuncFile_2_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpSupNum = getSupNumFromJuncStr(tmpStr);
		if(tmpSupNum > sup_num_thres_max)
			tmpSupNum = sup_num_thres_max;
		for(int tmp = 1; tmp <= tmpSupNum; tmp++)
			foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum[tmp-1] ++;
	}
	foundSJ_inJuncFile_2_ifs.close();
	
	cout << "start to  generate unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum" << endl;
	log_ofs << "start to generate unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum" << endl;
	vector<int> unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum;
	for(int tmp = 0; tmp < sup_num_thres_max; tmp++)
		unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum.push_back(0);
	ifstream unfoundSJ_inJuncFile_2_ifs(unfoundSJ_inJuncFile_2.c_str());
	while(!unfoundSJ_inJuncFile_2_ifs.eof())
	{
		string tmpStr;
		getline(unfoundSJ_inJuncFile_2_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpSupNum = getSupNumFromJuncStr(tmpStr);
		if(tmpSupNum > sup_num_thres_max)
			tmpSupNum = sup_num_thres_max;
		for(int tmp = 1; tmp <= tmpSupNum; tmp++)
			unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum[tmp-1] ++;		
	}
	unfoundSJ_inJuncFile_2_ifs.close();

	cout << "start to output junction comparing results with varying supporting num" << endl;
	log_ofs << "start to output junction comparing results with varying supporting num" << endl;	
	for(int tmp = 0; tmp < sup_num_thres_max; tmp++)
	{
		//cout << "tmp: " << tmp << endl;
		int tmp_sup_num_thres = tmp + 1;
		int tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum = foundJuncNumVec_in1stJuncHash_withinOffset_varySupNum[tmp];
		int tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum = foundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum[tmp];
		int tmp_unfoundJuncNum_in2ndJuncHash_withinOffset_varySupNum = unfoundJuncNumVec_in2ndJuncHash_withinOffset_varySupNum[tmp];
		int tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum 
			= tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum + tmp_unfoundJuncNum_in2ndJuncHash_withinOffset_varySupNum;
		//cout << "tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum: " << endl << tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum << endl;
		//cout << "tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum: " << endl << tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum << endl;
		//cout << "tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum: " << endl << tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum << endl;
		double tmp_sensitivity_varySupNum 
			= ((double)tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum / (double)juncNum_validIntronSize_1) * 100;
		double tmp_specificity_varySupNum 
			= ((double)tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum / 
				(double)tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum) * 100;
		compare_varySupNum_raw_ofs << "SupNum[" << tmp_sup_num_thres << "]:" << endl;
		compare_varySupNum_raw_ofs << "sensitivity: " << tmp_sensitivity_varySupNum << " "
			<< tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum << "/" << juncNum_validIntronSize_1 << endl;
		compare_varySupNum_raw_ofs << "specificity: " << tmp_specificity_varySupNum << " "
			<< tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum << "/" << tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum << endl;
		compare_varySupNum_TPR_FPR_ofs << tmp_sup_num_thres << "\t" << tmp_sensitivity_varySupNum << "\t" 
			<< 100 - tmp_specificity_varySupNum << "\t" << toolName << endl;
		compare_varySupNum_trueNum_falseNum_ofs << tmp_sup_num_thres << "\t" << tmp_foundJuncNum_in1stJuncHash_withinOffset_varySupNum << "\t"
			<< tmp_totalJuncNum_in2ndJuncHash_withinOffset_varySupNum - tmp_foundJuncNum_in2ndJuncHash_withinOffset_varySupNum << "\t" 
			<< toolName << endl;
	}

	delete juncHash_1;
	delete juncHash_2;
	delete indexInfo;
	compare_ofs.close();
	compare_varySupNum_trueNum_falseNum_ofs.close();
	compare_varySupNum_TPR_FPR_ofs.close();
	compare_varySupNum_raw_ofs.close();
	parameter_ifs.close();
	log_ofs.close();
}