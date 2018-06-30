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
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>


#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../phase1/arrayQueue.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"
#include "../../phase2/spliceJunctionHash_info.h"
#include "../../phase2/unmapEnd_info.h"
#include "../../phase2/unfixedHead.h"
#include "../../phase2/unfixedTail.h"
#include "../../phase2/incompleteLongHead.h"
#include "../../phase2/incompleteLongTail.h"
#include "../../phase2/sam2junc.h"
#include "../../fixHeadTail.h"
#include "../../phase2/fixOneEndUnmapped.h"
#include "../../fixPhase1.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"
#include "../../general/alignmentToJunc.h"
#include "./candidateFusionRegionDetection/fusionRegion_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable GlobalIndex LocalIndex inputAlignmentInfo_unpaired ";
		cout << "inputAlignmentInfo_incomplete outputFolder thread_num fasta_or_fastq" << endl;
		exit(1);
	}
	string globalIndexStr = argv[1];
	string localIndexStr = argv[2];
	string inputAlignInfo_unpaired = argv[3]; // candidate encompassing reads
	string inputAlignInfo_incompletePaired = argv[4]; // candidate spanning reads

	string outputFolderStr = argv[5];
	string threadNumStr = argv[6];
	string fasta_or_fastq_str = argv[7];

	cout << "alignAll_fixFusion starts ..." << endl;

	int normalRecordNum_fixFusion = 1;//500000;
	int threads_num = atoi(threadNumStr.c_str());
	bool SE_or_PE_bool = false;
	bool fasta_or_fastq_bool;
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta")||(fasta_or_fastq_str == "FASTA"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq")||(fasta_or_fastq_str == "FASTQ"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "input file format invalid" << endl;
		exit(1);
	}

	cout << "creat folders and files ......" << endl;	
	string mkdirOutputCommand = "mkdir -p " + outputFolderStr;
	system(mkdirOutputCommand.c_str());
	string outputFilePrefix = outputFolderStr + "/";

	string log_file = outputFilePrefix + "process.log";
	ofstream log_ofs(log_file.c_str());

	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);

	log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;


	string candidateFusionRegionDetection_folder = outputFilePrefix + "candidateFusionRegionDetection";
	string mkdir_candidateFusionRegionDetection_folder = "mkdir -p " + candidateFusionRegionDetection_folder;
	system(mkdir_candidateFusionRegionDetection_folder.c_str());
	candidateFusionRegionDetection_folder += "/";

	string candidateFusionRegionDetection_regionInfo_folder
		= candidateFusionRegionDetection_folder + "candidateRegionInfo";
	string mkdir_candidateFusionRegionDetection_regionInfo_folder
		= "mkdir -p " + candidateFusionRegionDetection_regionInfo_folder;
	system(mkdir_candidateFusionRegionDetection_regionInfo_folder.c_str());

	string candidateFusionRegionDetection_peAlignInfo_complete_unpair_folder
		= candidateFusionRegionDetection_folder + "complete_unpair";
	string mkdir_candidateFusionRegionDetection_peAlignInfo_complete_unpair_folder
		= "mkdir -p " + candidateFusionRegionDetection_peAlignInfo_complete_unpair_folder;
	system(mkdir_candidateFusionRegionDetection_peAlignInfo_complete_unpair_folder.c_str());

	string candidateFusionRegionDetection_peAlignInfo_incomplete_unpair_folder
		= candidateFusionRegionDetection_folder + "incomplete_unpair";
	string mkdir_candidateFusionRegionDetection_peAlignInfo_incomplete_unpair_folder
		= "mkdir -p " + candidateFusionRegionDetection_peAlignInfo_incomplete_unpair_folder;
	system(mkdir_candidateFusionRegionDetection_peAlignInfo_incomplete_unpair_folder.c_str());

	string candidateFusionRegionDetection_peAlignInfo_incomplete_pair_folder
		= candidateFusionRegionDetection_folder + "incomplete_pair";
	string mkdir_candidateFusionRegionDetection_peAlignInfo_incomplete_pair_folder
		= "mkdir -p " + candidateFusionRegionDetection_peAlignInfo_incomplete_pair_folder;
	system(mkdir_candidateFusionRegionDetection_peAlignInfo_incomplete_pair_folder.c_str());

	string candidateFusionRegionDetection_regionInfo_file 
		= candidateFusionRegionDetection_regionInfo_folder + "/candidateFusionRegion.txt";
	ofstream candiFusionRegion_ofs(candidateFusionRegionDetection_regionInfo_file.c_str());

	log_ofs << "start to initiate genomicRegionVecInfo " << endl;

	GenomicRegionVec_Info* genomicRegionVecInfo = new GenomicRegionVec_Info();
	genomicRegionVecInfo->initiate_withIndexInfo(indexInfo);
	FusionRegionVec_Info* fusionRegionVecInfo = new FusionRegionVec_Info();

	ifstream unpair_alignInfo_ifs(inputAlignInfo_unpaired.c_str());

	int normalRecordNum = normalRecordNum_fixFusion; //1000000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;// = normalRecordNum;

	log_ofs << "start to proess unpaired alignments " << endl;
	for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
	{
		//cout << "tmpTurn: " << tmpTurn << endl;
		vector<string> line1StrVec(normalRecordNum);
		vector<string> line2StrVec(normalRecordNum);
		vector<string> line3StrVec(normalRecordNum);
		vector<string> line4StrVec(normalRecordNum);
		vector<string> line5StrVec(normalRecordNum);
		vector<string> line6StrVec(normalRecordNum);
		vector<string> line7StrVec(normalRecordNum);
		vector<string> line8StrVec(normalRecordNum);
		vector<string> line9StrVec(normalRecordNum);
		vector<string> line10StrVec(normalRecordNum);

		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		realRecordNum = normalRecordNum;
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
			string line1, line2, line3, line4, line5, line6, line7, 
				line8, line9, line10, line11;	
			if(unpair_alignInfo_ifs.eof())
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
			}				
			getline(unpair_alignInfo_ifs, line11);
			if(unpair_alignInfo_ifs.eof())
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
			}		
			getline(unpair_alignInfo_ifs, line1);
			getline(unpair_alignInfo_ifs, line2);
			getline(unpair_alignInfo_ifs, line3);
			getline(unpair_alignInfo_ifs, line4);
			getline(unpair_alignInfo_ifs, line5);
			getline(unpair_alignInfo_ifs, line6);
			getline(unpair_alignInfo_ifs, line7);
			getline(unpair_alignInfo_ifs, line8);
			getline(unpair_alignInfo_ifs, line9);
			getline(unpair_alignInfo_ifs, line10);
			//getline(inputUnfixedHeadTailRecord_ifs, line11);

			line1StrVec[recordNumTmp] = line1;
			line2StrVec[recordNumTmp] = line2;
			line3StrVec[recordNumTmp] = line3;
			line4StrVec[recordNumTmp] = line4;
			line5StrVec[recordNumTmp] = line5;
			line6StrVec[recordNumTmp] = line6;
			line7StrVec[recordNumTmp] = line7;
			line8StrVec[recordNumTmp] = line8;
			line9StrVec[recordNumTmp] = line9;
			line10StrVec[recordNumTmp] = line10;
		}
		omp_set_num_threads(threads_num);
		//omp_set_num_threads(1);
		#pragma omp parallel for			
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			////////////  parse long head reads record after 1-mapping process  ////////
			int threadNO = omp_get_thread_num();
			PE_Read_Info peReadInfo;// = new PE_Read_Info();
			PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
			peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(
				line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
				line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
				line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
				line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, 
				indexInfo, fasta_or_fastq_bool, SE_or_PE_bool);
			//cout << "line2StrVec[tmpOpenMP]: " << line2StrVec[tmpOpenMP] << endl;
			//cout << endl << "***********************************************************" << endl;
			//cout << "readName: " << endl << peReadInfo.returnReadName_1() << endl;	
			fusionRegionVecInfo->checkAndGenerateCandidateFusionRegionInfo_fromUnpairCompleteReadPair(
				peAlignInfo, genomicRegionVecInfo, indexInfo);

			peAlignInfo->memoryFree();
			delete peAlignInfo;
		}
	}

	int validFusionRegionPair_supportNum_min = 10;
	fusionRegionVecInfo->output_fusionRegionVecInfo(candiFusionRegion_ofs, genomicRegionVecInfo, indexInfo,
		validFusionRegionPair_supportNum_min);

	fusionRegionVecInfo->memoryFree();
	delete fusionRegionVecInfo;
	genomicRegionVecInfo->memoryFree();
	delete genomicRegionVecInfo;	
	delete indexInfo;
	candiFusionRegion_ofs.close();
	return 0;
}