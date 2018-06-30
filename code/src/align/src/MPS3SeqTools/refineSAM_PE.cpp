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

#include "../general/otherFunc.h"
#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/fixSingleAnchorNWDP_info.h"
#include "../general/fixDoubleAnchorMatch_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <inputSAM> <outputFolder>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputSAMpath = argv[2];
	string outputFolderStr = argv[3];

	int maxReadBaseNumInPathStructure = 30;
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_file = outputFolderStr + "log";
	string headerSection_file = outputFolderStr + "headerSection";
	string alignInferJuncHashFile_file_beforeSeqSimi = outputFolderStr + "alignInferJuncHash.junc.beforeSeqSimi";
	string alignInferJuncHashFile_file_afterSeqSimi = outputFolderStr + "alignInferJuncHash.junc.afterSeqSimi";
	string final_sam_file = outputFolderStr + "final.sam";
	string multi_sam_file_ori = outputFolderStr + "multi_ori.sam";
	string unmap_sam_file_ori = outputFolderStr + "unmap_ori.sam";
	string unique_final_sam_file = outputFolderStr + "unique_final.sam";
	string unique_corrected_sam_file_original = outputFolderStr + "unique_corrected.original.sam";
	string unique_corrected_sam_file = outputFolderStr + "unique_corrected.sam";
	string multi_final_sam_file = outputFolderStr + "multi_final.sam";
	string multi_corrected_sam_file_original = outputFolderStr + "multi_corrected.original.sam";
	string multi_corrected_sam_file = outputFolderStr + "multi_corrected.sam";
	ofstream log_ofs(log_file.c_str());
	ofstream unique_final_sam_ofs(unique_final_sam_file.c_str());
	ofstream unique_corrected_sam_original_ofs(unique_corrected_sam_file_original.c_str());
	ofstream unique_corrected_sam_ofs(unique_corrected_sam_file.c_str());

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
	cout << "finish loading chromosomes" << endl;
	log_ofs << "finish loading chromosomes" << endl;

	vector<string> inputSAMfileVec;
	inputSAMfileVec.push_back(inputSAMpath);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to initiate refineSAMvecInfo" << endl;
	log_ofs << "start to initiate refineSAMvecInfo" << endl;
	RefineSAM_PE_Vec_Info* refineSAMvecInfo_PE = new RefineSAM_PE_Vec_Info();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to insert SJ of alignments into alignInferJuncHash" << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to insert SJ of alignments into alignInferJuncHash" << endl;	
	// insert SJ into SJmap, filter can be extended unique alignments
	alignInferJunctionHashInfo->insertJuncFromAlignmentFileVec_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_PE(
		inputSAMfileVec, indexInfo, maxReadBaseNumInPathStructure,
		multi_sam_file_ori, unmap_sam_file_ori, unique_final_sam_ofs,
		headerSection_file,
		//unique_corrected_sam_original_ofs, unique_corrected_sam_ofs,
		refineSAMvecInfo_PE);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to output alignInferJuncHash without anchor sequence similarity ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to output alignInferJuncHash without anchor sequence similarity ..." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo(
		indexInfo, alignInferJuncHashFile_file_beforeSeqSimi);

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
	alignInferJunctionHashInfo->getAlterSpliceSites_compareAnchorSimilarity_onlyLowSupportSJ(
		SJhashInfo, alterSpliceSiteOffset, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to output alignInferJuncHash after checking extension and alterSpliceSite anchor sequence similarity ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to output alignInferJuncHash after checking extension and alterSpliceSite anchor sequence similarity ..." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
		indexInfo, alignInferJuncHashFile_file_afterSeqSimi);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to refine and output low support and confidence reads ...." << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to refine and output low support and confidence reads ...." << endl;
	alignInferJunctionHashInfo->correct_output_SAM_withLowSupportLowConfidenceSJ_PE(refineSAMvecInfo_PE, 
		unique_final_sam_ofs, unique_corrected_sam_original_ofs, unique_corrected_sam_ofs, indexInfo);

	unique_final_sam_ofs.close();
	unique_corrected_sam_original_ofs.close();
	unique_corrected_sam_ofs.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to correct and refine mulit-alignments" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to correct and refine mulit-alignments" << endl;
	alignInferJunctionHashInfo->refineMultiAlignment_PE(multi_sam_file_ori, multi_final_sam_file,
		multi_corrected_sam_file_original, multi_corrected_sam_file, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to generate final SAM file ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to generate final SAM file ..." << endl;
	string cat_cmd = "cat " + headerSection_file + " " + unique_final_sam_file 
		+ " " + multi_final_sam_file + " " + unmap_sam_file_ori + " > " + final_sam_file;
	system(cat_cmd.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "end of generating final SAM file ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "end of generating final SAM file ..." << endl;

	log_ofs.close();
	delete indexInfo;
	refineSAMvecInfo_PE->freeMemory();
	delete refineSAMvecInfo_PE;
	delete SJhashInfo;
	delete alignInferJunctionHashInfo;
	return 0;
}