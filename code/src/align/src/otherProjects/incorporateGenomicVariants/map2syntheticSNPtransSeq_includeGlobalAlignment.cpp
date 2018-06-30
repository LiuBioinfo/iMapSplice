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
#include "../../general/gene_count_vec.h"
#include "../../general/transcript2geneMap_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "general/syntheticSNPtransSeq_info.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info.h"
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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndex inputGlobalIndex inputRead_1 inputRead_2 outputFolder threadNum" << endl;
		exit(1);
	}
	bool fasta_or_fastq_bool = true;
	bool InputAsFastq = false;
	//bool SE_or_PE_bool = true;
	bool checkQualSeqForReadSegSeq = false;
	bool Do_cirRNA = false;
	bool annotation_provided_bool = false;
	bool Do_annotation_only_bool = false;
	bool Do_extendHeadTail_phase1 = true;
	Annotation_Info* annotationInfo = new Annotation_Info();

	string threadNumStr = argv[6];
	int threads_num = atoi(threadNumStr.c_str());

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());	

	cout << "start to load indexes" << endl;
	string indexStr = argv[1];
	indexStr.append("/");
	string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);	
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	cout << "SyntheticSNPtransSeq index files loaded" << endl;
	
	string indexStr_global = argv[2];
	indexStr_global.append("/");
	string preIndexArrayPreStr_global = indexStr_global;
	preIndexArrayPreStr_global.append("/");
	string preIndexMapLengthArrayStr_global = preIndexArrayPreStr_global; preIndexMapLengthArrayStr_global.append("_MapLength"); ifstream preIndexMapLengthArray_ifs_global(preIndexMapLengthArrayStr_global.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr_global = preIndexArrayPreStr_global; preIndexIntervalStartArrayStr_global.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs_global(preIndexIntervalStartArrayStr_global.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr_global = preIndexArrayPreStr_global; preIndexIntervalEndArrayStr_global.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs_global(preIndexIntervalEndArrayStr_global.c_str(), ios::binary);
	int* preIndexMapLengthArray_global; preIndexMapLengthArray_global = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs_global.read((char*)preIndexMapLengthArray_global, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray_global; preIndexIntervalStartArray_global = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs_global.read((char*)preIndexIntervalStartArray_global, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray_global; preIndexIntervalEndArray_global = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs_global.read((char*)preIndexIntervalEndArray_global, PreIndexSize * sizeof(int));

	string chrom_bit_file_global = indexStr_global; chrom_bit_file_global.append("_chrom"); ifstream chrom_bit_file_ifs_global(chrom_bit_file_global.c_str(),ios::binary);
	string parameter_file_global = indexStr_global; parameter_file_global.append("_parameter"); ifstream parameter_file_ifs_global(parameter_file_global.c_str(),ios::binary);
	Index_Info* indexInfo_global = new Index_Info(parameter_file_ifs_global, log_ofs);
	char *chrom_global; chrom_global = (char*)malloc((indexInfo_global->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs_global.read((char*)chrom_global, (indexInfo_global->returnIndexSize()) * sizeof(char)); 
	indexInfo_global->readGenome(chrom_global);
	indexInfo_global->initiate();	
	indexInfo_global->initiateChrNameIndexArray(1000);
	string SA_file_global = indexStr_global; SA_file_global.append("_SA"); 
	string lcpCompress_file_global = indexStr_global; lcpCompress_file_global.append("_lcpCompress"); 
	string childTab_file_global = indexStr_global; childTab_file_global.append("_childTab"); 
	string verifyChild_file_global = indexStr_global; verifyChild_file_global.append("_detChild"); 	
    unsigned int *sa_global; sa_global = (unsigned int*)malloc((indexInfo_global->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab_global; childTab_global = (unsigned int*)malloc((indexInfo_global->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress_global; lcpCompress_global = (BYTE*)malloc((indexInfo_global->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild_global; verifyChild_global = (BYTE*)malloc((indexInfo_global->returnIndexSize()) * sizeof(BYTE)); 	
	ifstream SA_file_ifs_global(SA_file_global.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs_global(lcpCompress_file_global.c_str(),ios::binary);
	ifstream childTab_file_ifs_global(childTab_file_global.c_str(),ios::binary);
	ifstream verifyChild_file_ifs_global(verifyChild_file_global.c_str(),ios::binary);
	SA_file_ifs_global.read((char*)sa_global, (indexInfo_global->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress_file_ifs_global.read((char*)lcpCompress_global, (indexInfo_global->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs_global.read((char*)childTab_global, (indexInfo_global->returnIndexSize()) * sizeof(unsigned int));
	verifyChild_file_ifs_global.read((char*)verifyChild_global, (indexInfo_global->returnIndexSize()) * sizeof(BYTE));
	cout << "Global index files loaded" << endl;


	string InputReadFile = argv[3];
	string InputReadFile_PE = argv[4];
	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;
	int normalRecordNum = 1;//normalRecordNum_1stMapping; //1000000;//1500000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	int readPairNum = 0;
	int readTotalNum = 0;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);

	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}

	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... map2syntheticSNPtransSeq process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... map2syntheticSNPtransSeq process starts ......" << endl << endl; 

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();
	for(tmpTurn = 0; ; tmpTurn++) // used to control turns to process
	{
		if(EndOfRecord)
			break;		
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		getline(inputRead_ifs, line1); // readName_1
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		getline(inputRead_ifs, line2); // readSeq_1
    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;    		
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    			readQualSeq1Vec[recordNumTmp] = line4;
    		}
    		getline(inputRead_PE_ifs, line1_PE); // readName_2
    		readName2Vec[recordNumTmp] = line1_PE.substr(1);
    		getline(inputRead_PE_ifs, line2_PE);
    		line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
    		readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
    		if(InputAsFastq)
    		{
    			getline(inputRead_PE_ifs, line3_PE);
    			getline(inputRead_PE_ifs, line4_PE);
    			readQualSeq2Vec[recordNumTmp] = line4_PE;
    		}
		}
		readTotalNum += realRecordNum;
		log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		log_ofs << "start to fix reads, turn: " << tmpTurn + 1 << endl;

		omp_set_num_threads(threads_num);
		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();

			///////////////////////////     map to synthetic SNP trans seq      ////////////////////////////////// 
			string nullStr = "*";
			PE_Read_Info readInfo_end1, readInfo_end2; //= new PE_Read_Info();
    		FixPhase1Info fixPhase1Info_end1, fixPhase1Info_end2;
    		PE_Read_Alignment_Info peAlignInfo_syntheticSNPtransSeq_end1, peAlignInfo_syntheticSNPtransSeq_end2;
    		PE_Read_Alignment_Info peAlignInfo_syntheticSNPtransSeq_end1_converted2wholeGenome, 
    			peAlignInfo_syntheticSNPtransSeq_end2_converted2wholeGenome;
			/////////////////////////////   end 1  ////////////////////////////////////////////////
			readInfo_end1.initiateReadInfo(readName1Vec[tmpOpenMP], nullStr, readSeq1Vec[tmpOpenMP], 
				nullStr, readQualSeq1Vec[tmpOpenMP], nullStr, fasta_or_fastq_bool, true);
			fixPhase1Info_end1.fixPhase1_segInfo_map2syntheticSNPtransSeq(
				sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, readInfo_end1);
			fixPhase1Info_end1.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_end1, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
			fixPhase1Info_end1.fixPhase1_gapInfo(readInfo_end1, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
			peAlignInfo_syntheticSNPtransSeq_end1.initiatePeAlignInfo(
				fixPhase1Info_end1.pathInfo_Nor1, fixPhase1Info_end1.pathInfo_Rcm1, 
				fixPhase1Info_end1.pathInfo_Nor2, fixPhase1Info_end1.pathInfo_Rcm2, indexInfo, true);
			peAlignInfo_syntheticSNPtransSeq_end1.chooseBestAlignment_final_SE();
			// peAlignInfo_syntheticSNPtransSeq_end1_converted2wholeGenome.convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
			// 	peAlignInfo_syntheticSNPtransSeq_end1, indexInfo, indexInfo_global, readInfo_end1, true);
			/////////////////////////////   end 2  ////////////////////////////////////////////////
			readInfo_end2.initiateReadInfo(readName2Vec[tmpOpenMP], nullStr, readSeq2Vec[tmpOpenMP], 
				nullStr, readQualSeq2Vec[tmpOpenMP], nullStr, fasta_or_fastq_bool, true);
			fixPhase1Info_end2.fixPhase1_segInfo_map2syntheticSNPtransSeq(
				sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, readInfo_end2);				
			fixPhase1Info_end2.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_end2, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);			
			fixPhase1Info_end2.fixPhase1_gapInfo(readInfo_end2, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
			peAlignInfo_syntheticSNPtransSeq_end2.initiatePeAlignInfo(
				fixPhase1Info_end2.pathInfo_Nor1, fixPhase1Info_end2.pathInfo_Rcm1, 
				fixPhase1Info_end2.pathInfo_Nor2, fixPhase1Info_end2.pathInfo_Rcm2, indexInfo, true);			
			peAlignInfo_syntheticSNPtransSeq_end2.chooseBestAlignment_final_SE();		
			// peAlignInfo_syntheticSNPtransSeq_end2_converted2wholeGenome.convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
			// 	peAlignInfo_syntheticSNPtransSeq_end2, indexInfo, indexInfo_global, readInfo_end2, true);
			cout << endl << "*********************************" << endl;
			cout << "###############################" << endl; fixPhase1Info_end1.coutFixedPathInfo_SE(readInfo_end1, indexInfo);
			//cout << "## after choosing best alignment for read_1: ##" << endl; peAlignInfo_syntheticSNPtransSeq_end1.cout_seAlignVec_final();
			cout << "###############################" << endl; fixPhase1Info_end2.coutFixedPathInfo_SE(readInfo_end2, indexInfo);
			//cout << "## after choosing best alignment for read_2: ##" << endl; peAlignInfo_syntheticSNPtransSeq_end2.cout_seAlignVec_final();
			fixPhase1Info_end1.memoryFree();
			fixPhase1Info_end2.memoryFree();

			////////////////////////////////   map to whole genome   //////////////////////////////////
			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP], 
				readSeq2Vec[tmpOpenMP], readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], fasta_or_fastq_bool, false);
			//cout << endl << "read_name_1: " << readName1Vec[tmpOpenMP] << endl;
			//cout << "read_name_2: " << readName2Vec[tmpOpenMP] << endl;
    		FixPhase1Info fixPhase1Info;
			fixPhase1Info.fixPhase1_segInfo(sa_global, lcpCompress_global, childTab_global, chrom_global, verifyChild_global, indexInfo_global, 
				preIndexMapLengthArray_global, preIndexIntervalStartArray_global, preIndexIntervalEndArray_global, readInfo, 
				repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, false);
			// cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info.coutSegInfo_Nor1(indexInfo_global);
			// cout << "segInfo_Rcm1: " << endl; fixPhase1Info.coutSegInfo_Rcm1(indexInfo_global);
			// cout << "segInfo_Nor2: " << endl; fixPhase1Info.coutSegInfo_Nor2(indexInfo_global);
			// cout << "segInfo_Rcm2: " << endl; fixPhase1Info.coutSegInfo_Rcm2(indexInfo_global);	
			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo_global, readInfo, annotation_provided_bool, Do_annotation_only_bool, 
				annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, false);
			// cout << "getPossiPathInfo_Nor1: " << endl; fixPhase1Info.coutPossiPathStr_Nor1();
			// cout << "getPossiPathInfo_Rcm1: " << endl; fixPhase1Info.coutPossiPathStr_Rcm1();
			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo_global, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, false);
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo(fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo_global, false);
			cout << "###############################" << endl; cout << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl;
			peAlignInfo.includeSyntheticSNPtransSeqPeAlignInfo(peAlignInfo_syntheticSNPtransSeq_end1, 
				peAlignInfo_syntheticSNPtransSeq_end2, indexInfo, indexInfo_global, readInfo);
			cout << "## After including SyntheticSNPtransSeqPeAlignInfo ##" << endl; cout << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl;
			peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(), readInfo.returnReadSeqLength_2());
			cout << "after filtering and choosing best alignments: " << endl; peAlignInfo.cout_finalAlignPair();		

			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();			
		} // read file end		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads, turn: " << tmpTurn + 1 << endl;
	}
	cout << endl << "finish gene count on total reads" << endl;
	cout << "readTotalNum: " << readTotalNum << endl;
	log_ofs << "readTotalNum: " << readTotalNum << endl;	

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	log_ofs.close();
	return 0;
}