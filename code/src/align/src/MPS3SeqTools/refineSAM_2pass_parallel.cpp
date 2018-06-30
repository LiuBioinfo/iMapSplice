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

time_t nowtime;
struct tm *local;

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/alignInferJunctionHash_info_vec.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <threadNum> <OutputSJ> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}
	//int maxReadBaseNumInPathStructure = 30;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string threadNumStr = argv[2];
	int threadNum_int = atoi(threadNumStr.c_str());

	int alignInferJuncHashInfoVecSize = threadNum_int;

	vector<string> inputSAMfileVec;
	for(int tmp = 4; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initiate chrNameIndexArray" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initiate chrNameIndexArray" << endl;	
	indexInfo->initiateChrNameIndexArray(1000);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite merged alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite merged alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_merged 
		= new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo_merged->initiateAlignInferJunctionHashInfo(chromNum);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite alignInferJunctionHashInfo vec " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite alignInferJunctionHashInfo vec " << endl;	
	AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec 
		= new AlignInferJunctionHash_Info_Vec();
	alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(alignInferJuncHashInfoVecSize, chromNum);
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to insert SJ into SJmap from SAM file" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to insert SJ into SJmap from SAM file" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(
		inputSAMfileVec, indexInfo, alignInferJuncHashInfoVecSize, log_ofs);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to merge alignInferJuncHashInfo in vec" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to merge alignInferJuncHashInfo in vec" << endl;
	alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(
		alignInferJunctionHashInfo_merged, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initiate SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initiate SJhashInfo" << endl;	
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo_merged->convert2SJhashInfo(SJ, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to output SJ map" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to output SJ map" << endl;	
	// output SJmap
	alignInferJunctionHashInfo_merged->outputAlignInferInfoHashInfo_chrNamePos_supportNum(
		indexInfo, outputSJstr);


	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "end of running getAlignInferJuncFromSAMfile_chrNamePos" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "end of running getAlignInferJuncFromSAMfile_chrNamePos" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to check flankstring case" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to check flankstring case" << endl;
	// generate_flankStringCase_parallel



	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to getAnchorSimilarity_extension" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to getAnchorSimilarity_extension" << endl;
	// getAnchorSimilarity_extension_parallel

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to getAnchorSimilarity_alterSpliceSite" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to getAnchorSimilarity_alterSpliceSite" << endl;
	// getAnchorSimilarity_alterSpliceSite_parallel?


	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to clip alignments with alignInferJuncHash" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to clip alignments with alignInferJuncHash" << endl;
	// start to clip alignments ... //


	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "end of program" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "end of program" << endl;
	alignInferJunctionHashInfoVec->freeMemory();
	delete alignInferJunctionHashInfoVec;
	delete alignInferJunctionHashInfo_merged;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}