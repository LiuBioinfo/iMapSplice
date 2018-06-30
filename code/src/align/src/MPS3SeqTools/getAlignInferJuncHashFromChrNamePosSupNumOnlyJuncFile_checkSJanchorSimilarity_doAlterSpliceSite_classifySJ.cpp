// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: alignerInferJunc raw file,  
// format: chrName chrStartPos chrEndPos JuncName supportNum
// Note: can specify the mode of classification,
// forRemapping: keep low support and noncanonical ones
// finalClassification: filter out low support and noncanonical ones

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

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexInfo> <inputAlignInferJuncHashFile> <OutputFolderPath> <Mode (forRemapping/final)>" << endl;
		exit(1);
	}
	int lowSupportNum_threshold = 2;
	int offsetInt_alterSpliceSite = 5;
	string SJ_true_sign_str = "T";
	string SJ_false_sign_str = "F";
	bool mode_forRemapping_or_finalClassification_bool = false;
	string modeStr = argv[4];
	if(modeStr == "forRemapping")
		mode_forRemapping_or_finalClassification_bool = true;
	else if(modeStr == "final")
		mode_forRemapping_or_finalClassification_bool = false;
	else
	{
		cout << "invalid mode: " << modeStr << endl;
		cout << "please choose forRemapping or final" << endl; 
		exit(1);
	}

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
	int chromNum = indexInfo->returnChromNum();
	cout << "finish loading chromosomes" << endl;

	string inputAlignInferJuncHashPath = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	

	string outputPath_log = outputFolderStr + "log.txt";
	ofstream log_ofs(outputPath_log.c_str());
	string outputSJpath_original = outputFolderStr + "original.alignerInferJuncHash";
	string outputSJpath_all_classified = outputFolderStr + "all_classified.alignInferJuncHash";
	string outputSJpath_all_withSeqSimilarity = outputFolderStr + "all_withSeqSimilarity.alignerInferJuncHash";
	string outputSJpath_kept = outputFolderStr + "kept.alignInferJuncHash";
	string outputSJpath_filterOut = outputFolderStr + "filterOut.alignInferJuncHash";
	string outputSJpath_filterOut_canBeExtendedDirectly = outputFolderStr + "filterOut_canBeExtendedDirectly.alignInferJuncHash";
	string outputSJpath_filterOut_canMove2alterSpliceSite = outputFolderStr + "filterOut_canMove2alterSpliceSite.alignerInferJuncHash";
	string outputSJpath_filterOut_lowSupNumNonCanonical = outputFolderStr + "filterOut_lowSupNumNonCanonical.alignerInferJuncHash";
	//ofstream originalSJ_ofs(outputSJpath_original.c_str());
	// ofstream allClassifiedSJ_ofs(outputSJpath_all_classified.c_str());
	// ofstream keptSJ_ofs(outputSJpath_kept.c_str());
	// ofstream filterOutSJ_ofs(outputSJpath_filterOut.c_str());
	// ofstream filterOutSJ_canBeExtendedDirectly_ofs(outputSJpath_filterOut_canBeExtendedDirectly.c_str());
	// ofstream filterOutSJ_canMove2alterSpliceSite_ofs(outputSJpath_filterOut_canMove2alterSpliceSite.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to insert SJs from JuncFile into alignInferJuncHash" << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to insert SJs from JuncFile into alignInferJuncHash" << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePos_supportNum(
		inputAlignInferJuncHashPath, indexInfo);


	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to output alignInferJuncHash without anchor sequence similarity ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to output alignInferJuncHash without anchor sequence similarity ..." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo(
		indexInfo, outputSJpath_original);


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
	cout << "start to get extension and alter splice sites possibility only for low support SJs" << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to get extension and alter splice sites possibility only for low support SJs" << endl;
	int alterSpliceSiteOffset = 0;
	int defaultAnchorSize = 30;
	alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity_allSJ_withDefaultAnchorSize(
		SJhashInfo, alterSpliceSiteOffset, defaultAnchorSize, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to output alignInferJuncHash after checking extension and alterSpliceSite anchor sequence similarity ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to output alignInferJuncHash after checking extension and alterSpliceSite anchor sequence similarity ..." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
		indexInfo, outputSJpath_all_withSeqSimilarity);

	if(mode_forRemapping_or_finalClassification_bool) // classifcation for remapping, keep low support and noncanonical ones
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "..."; 
		cout << "start to output alignInferJuncHash after classifying for remapping" << endl;
		log_ofs << endl << "[" << asctime(local) << "..."; 
		log_ofs << "start to output alignInferJuncHash after classifying for remapping" << endl;
		// for remapping, so would not throw out low support noncanonical SJs
		alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_classified_forRemapping( 
			indexInfo, outputSJpath_all_classified, outputSJpath_kept, 
			outputSJpath_filterOut, outputSJpath_filterOut_canBeExtendedDirectly, 
			outputSJpath_filterOut_canMove2alterSpliceSite, log_ofs);
	}
	else
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "..."; 
		cout << "start to output alignInferJuncHash after the final classifying " << endl;
		log_ofs << endl << "[" << asctime(local) << "..."; 
		log_ofs << "start to output alignInferJuncHash after the final classifying (including )" << endl;
		// for finalClassification, so would throw out low support noncanonical SJs
		alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_classified_finalClassification( 
			indexInfo, outputSJpath_all_classified, outputSJpath_kept, 
			outputSJpath_filterOut, outputSJpath_filterOut_canBeExtendedDirectly, 
			outputSJpath_filterOut_canMove2alterSpliceSite,
			outputSJpath_filterOut_lowSupNumNonCanonical,
			lowSupportNum_threshold, log_ofs);		
	}

	log_ofs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}