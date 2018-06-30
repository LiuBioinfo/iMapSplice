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
	//string final_sam_file = outputFolderStr + "final.sam";
	string multi_sam_file_ori = outputFolderStr + "multi_ori.sam";
	string unmap_sam_file_ori = outputFolderStr + "unmap_ori.sam";
	string unique_final_sam_file = outputFolderStr + "unique_final.sam";
	//string unique_corrected_sam_file_original = outputFolderStr + "unique_corrected.original.sam";
	//string unique_corrected_sam_file = outputFolderStr + "unique_corrected.sam";
	//string multi_final_sam_file = outputFolderStr + "multi_final.sam";
	//string multi_corrected_sam_file_original = outputFolderStr + "multi_corrected.original.sam";
	//string multi_corrected_sam_file = outputFolderStr + "multi_corrected.sam";
	ofstream log_ofs(log_file.c_str());
	ofstream unique_final_sam_ofs(unique_final_sam_file.c_str());
	//ofstream unique_corrected_sam_original_ofs(unique_corrected_sam_file_original.c_str());
	//ofstream unique_corrected_sam_ofs(unique_corrected_sam_file.c_str());

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
	RefineSAM_SE_Vec_Info* refineSAMvecInfo_SE = new RefineSAM_SE_Vec_Info();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "start to parseSAM" << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "start to parseSAM" << endl;	
	// insert SJ into SJmap, filter can be extended unique alignments
	alignInferJunctionHashInfo->parseSAMfileVec_SE(
		inputSAMfileVec, indexInfo, maxReadBaseNumInPathStructure,
		multi_sam_file_ori, unmap_sam_file_ori, unique_final_sam_ofs,
		headerSection_file,
		//unique_corrected_sam_original_ofs, unique_corrected_sam_ofs,
		refineSAMvecInfo_SE);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "..."; 
	cout << "end of parseSAM" << endl;
	log_ofs << endl << "[" << asctime(local) << "..."; 
	log_ofs << "end of parseSAM" << endl;	

	parameter_ifs.close();
	unique_final_sam_ofs.close();
	log_ofs.close();
	delete indexInfo;
	refineSAMvecInfo_SE->freeMemory();
	delete refineSAMvecInfo_SE;
	//delete SJhashInfo;
	delete alignInferJunctionHashInfo;
	return 0;
}