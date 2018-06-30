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

#include "../constantDefinitions.h"
#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../phase2/spliceJunctionHash_info.h"
#include "../general/insertAnnotatedSJ2matchThroughEnd_info.h"
//#include "../general/fixSingleAnchorNWDP_info.h"
//#include "../general/fixDoubleAnchorMatch_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexInfo> <inputAnnotatedSJ> ";
		cout << "<inputSAM> <outputFolder>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputAnnotatedSJpath = argv[2];
	string inputSAMpath = argv[3];
	string outputFolderStr = argv[4];

	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_file = outputFolderStr + "log";
	ofstream log_ofs(log_file.c_str());

	string annotatedSJpath = outputFolderStr + "annotated.alignInferJuncHash";
	string headerSection_file = outputFolderStr + "headerSection";
	string finalSAM_file = outputFolderStr + "final.sam";
	string refinedSAM_file = outputFolderStr + "refined.sam";
	string refinedSAM_ori_file = outputFolderStr + "refined_ori.sam";
	string keptSAM_file = outputFolderStr + "kept.sam";

	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;
	int chromNum = indexInfo->returnChromNum();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePosOnly(
		inputAnnotatedSJpath, indexInfo);	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	log_ofs << "start to convert 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to output alignInferJuncHash with chrNamePos only ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to output alignInferJuncHash with chrNamePos only ..." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_chrNamePosOnly(
		indexInfo, annotatedSJpath);	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to refine alignment with annotated SJs ......" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to refine alignment with annotated SJs ......" << endl;

	ifstream sam_ifs(inputSAMpath.c_str());
	ofstream headerSection_ofs(headerSection_file.c_str());
	ofstream finalSAM_ofs(finalSAM_file.c_str());
	ofstream refinedSAM_ofs(refinedSAM_file.c_str());
	ofstream refinedSAM_ori_ofs(refinedSAM_ori_file.c_str());
	ofstream keptSAM_ofs(keptSAM_file.c_str());
	while(!sam_ifs.eof())
	{
		string tmpLineStr;
		getline(sam_ifs, tmpLineStr);
		if((sam_ifs.eof())||(tmpLineStr == ""))
			break;
		//cout << "tmpLineStr: " << tmpLineStr << endl;
		if(tmpLineStr.at(0) == '@')
		{
			headerSection_ofs << tmpLineStr << endl;
			continue;
		}
		// initiate, if reads unmapped, then output directly
		InsertAnnotatedSJ2matchThroughEnd_Info* tmpInsertAnnotatedSJ2samInfo
			= new InsertAnnotatedSJ2matchThroughEnd_Info();
		bool initiateMappedRead_bool 
			= tmpInsertAnnotatedSJ2samInfo->initiateWithSamStr(
				tmpLineStr, indexInfo);
		if(!initiateMappedRead_bool)
		{
			finalSAM_ofs << tmpLineStr << endl;
			keptSAM_ofs << tmpLineStr << endl;
			delete tmpInsertAnnotatedSJ2samInfo;
			continue;
		}

		// try to insert annotated SJs 2 alignment
		bool insertAnnotatedSJ2alignment_bool 
			= tmpInsertAnnotatedSJ2samInfo->insertAnnotatedSJ2matchThroughEnd(
				alignInferJunctionHashInfo, SJhashInfo,
				indexInfo);
		//cout << "insertAnnotatedSJ2alignment_bool: " << insertAnnotatedSJ2alignment_bool << endl;
		if(!insertAnnotatedSJ2alignment_bool)
		{
			finalSAM_ofs << tmpLineStr << endl;
			keptSAM_ofs << tmpLineStr << endl;
		}
		else
		{
			// bool tmpRefinedAlignmentCompleteOrNot 
			// 	= tmpInsertAnnotatedSJ2samInfo->completeOrNot();
			// if(tmpRefinedAlignmentCompleteOrNot )	
			// {
				string tmpRefinedSamStr 
					= tmpInsertAnnotatedSJ2samInfo->returnRefinedSamStr(indexInfo);
				finalSAM_ofs << tmpRefinedSamStr << endl;
				refinedSAM_ofs << tmpRefinedSamStr << endl;
				refinedSAM_ori_ofs << tmpLineStr << endl;
			// }
			// else
			// {
			// 	finalSAM_ofs << tmpLineStr << endl;
			// 	keptSAM_ofs << tmpLineStr << endl;				
			// }
		}
		delete tmpInsertAnnotatedSJ2samInfo;
	}

	parameter_ifs.close();
	delete indexInfo;
	finalSAM_ofs.close();
	refinedSAM_ofs.close();
	refinedSAM_ori_ofs.close();
	keptSAM_ofs.close();
	return 0;
}	