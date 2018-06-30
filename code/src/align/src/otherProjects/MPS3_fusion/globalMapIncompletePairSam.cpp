// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//Note: this program is used to do global mapping for headSoftclipped sequence in
// uniquePairedAlignment to detect fusion
// Note: no header section in input SAM file
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
#include "general/incompleteUniquePairedAlignment2detectFusion_info.h"
#include "general/refineFusionBreakPoint_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexInfoPath ";
		cout << "inputIncompleteUniquePairedAlignmentPath ";
		cout << "outputFolder threads_num fasta_or_fastq min_fusion_distance" << endl;
		exit(1);
	}
	int anchorSeqSize2output = 15;
	string positiveStrand = "+";
	string negativeStrand = "-";

	int offsetForAdjustingFusionBreakPoint_max = 5;
	string minFusionDistanceStr = argv[6];
	int min_fusion_distance = atoi(minFusionDistanceStr.c_str());
	int normalRecordNum_1stMapping = 200000;
	bool fasta_or_fastq_bool = true;
	string fasta_or_fastq_str = argv[5];
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "Please set the correct format, fasta or fastq" << endl;
		exit(1);
	}
	bool Do_cirRNA = true;
	Do_cirRNA = false;
	bool annotation_provided_bool = false;
	bool Do_annotation_only_bool = false;
	bool Do_extendHeadTail_phase1 = true;
	bool checkQualSeqForReadSegSeq = false;
	//Do_extendHeadTail_phase1 = false;
	string threadNumStr = argv[4];
	int threads_num = atoi(threadNumStr.c_str());
	omp_set_num_threads(threads_num);

	int readTotalNum = 0;

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	string outputDirStr = outputFolderStr + "/";
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
   	string settingsLogStr = outputDirStr + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	settings_log_ofs << "CommandLine:" << endl;
   	for(int tmp = 1; tmp < argc; tmp++)
   		settings_log_ofs << "\t" << argv[tmp];
   	settings_log_ofs << endl;
   	string progressLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());

   	string runtimeLogStr = outputDirStr + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());

   	string oriSAM_outputPath = outputDirStr + "/kept.sam";
   	string nonFusionSAM_outputPath = outputDirStr + "/nonFusion.sam";
   	string fusionSAM_outputPath = outputDirStr + "/fusion.sam";
   	string fusionSAM_outputPath_withBreakPoint = outputDirStr + "fusionReadsWithBreakPoint_raw.txt";
    //string fusionSAM_outputPath_withBreakPoint_adjusted = outputDirStr + "fusionReadsWithBreakPoint_adjusted.txt";
	string fusionSAM_outputPath_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering
		= outputDirStr + "fusionReadsWithBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering.txt";
	string fusionSAM_outputPath_withBreakPoint_adjusted_withFusionAnchorSeq
		= outputDirStr + "fusionReadsWithBreakPoint_adjusted_withFusionAnchorSeq.txt";
	string fusionSAM_outputPath_withBreakPoint_adjusted_withMapRegion
		= outputDirStr + "fusionReadsWithBreakPoint_adjusted_withMapRange.txt";

   	ofstream oriSAM_ofs(oriSAM_outputPath.c_str());
   	ofstream nonFusionSAM_ofs(nonFusionSAM_outputPath.c_str());
   	ofstream fusionSAM_ofs(fusionSAM_outputPath.c_str());
   	ofstream fusionSAM_withBreakPoint_ofs(fusionSAM_outputPath_withBreakPoint.c_str());
   	//ofstream fusionSAM_withBreakPoint_adjusted_ofs(fusionSAM_outputPath_withBreakPoint_adjusted.c_str());
   	ofstream fusionSAM_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering_ofs(
   		fusionSAM_outputPath_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering.c_str());
   	ofstream fusionSAM_withBreakPoint_adjusted_withFusionAnchorSeq_ofs(
   		fusionSAM_outputPath_withBreakPoint_adjusted_withFusionAnchorSeq.c_str());
   	ofstream fusionSAM_withBreakPoint_adjusted_withMapRegion_ofs(
   		fusionSAM_outputPath_withBreakPoint_adjusted_withMapRegion.c_str());

	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) 
		<< "... MPS fusion starts to do globalMapping for outerSoftClipUniquePairedAlignment ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) 
		<< "... MPS fusion starts to do globalMapping for outerSoftClipUniquePairedAlignment ......" << endl << endl; 


	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 

    string indexStr = argv[1];
    string preIndexArrayPreStr = indexStr; 	
    preIndexArrayPreStr.append("/");
    indexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; 
	preIndexMapLengthArrayStr.append("_MapLength"); 
	ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; 
	preIndexIntervalStartArrayStr.append("_IntervalStart"); 
	ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; 
	preIndexIntervalEndArrayStr.append("_IntervalEnd"); 
	ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; 
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); 
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; 
	preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; 
	preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	log_ofs << "finish loading preIndex ..." << endl;
	string chrom_bit_file = indexStr; 
	chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; 
	parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	log_ofs << "index: " << indexStr << endl;
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
 	
	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 	

	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;

	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;
	log_ofs << "start to load SA" << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "start to load childTab " << endl;
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load detChild" << endl;
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;
	runtime_log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;	

	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	Annotation_Info* annotationInfo = new Annotation_Info();
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////	
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///


	string inputIncompleteUniquePairedAlignmentPath = argv[2];
	ifstream incompleteUniquePairedAlignment_ifs(inputIncompleteUniquePairedAlignmentPath.c_str());
	int normalRecordNum = normalRecordNum_1stMapping;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	vector<string> inputSamStr1Vec(normalRecordNum);
	vector<string> inputSamStr2Vec(normalRecordNum);
	vector<string> outputOriSamStrVec(normalRecordNum);
	vector<string> outputNonfusionSamStrVec(normalRecordNum);
	vector<string> outputFusionSamStrVec(normalRecordNum);
	vector<string> outputFusionBreakPointStrVec(normalRecordNum);
	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) 
		<< "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) 
		<< "... 1st mapping process starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) 
		<< "... 1st mapping process starts ......" << endl << endl; 

	string line1, line2;

	for(tmpTurn = 0; 
		//tmpTurn <= 300     //used to control # of rounds to process
		; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		//cout << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) 
			<< "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if(incompleteUniquePairedAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
    		}
    		getline(incompleteUniquePairedAlignment_ifs, line1); // readName_1
    		if(incompleteUniquePairedAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		inputSamStr1Vec[recordNumTmp] = line1;
    		getline(incompleteUniquePairedAlignment_ifs, line2);		
    		inputSamStr2Vec[recordNumTmp] = line2;
		}

		readTotalNum += realRecordNum;

		runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;

		cout << endl << "[" << asctime(local) 
			<< "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		cout << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << "realRecordNum: " << realRecordNum << endl;
		cout << "threads_num: " << threads_num << endl;
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			//cout << "start to get threadNO " << endl;
			int threadNO = omp_get_thread_num();
			bool SE_or_PE_bool = true;
			//cout << "threadNO: " << threadNO << endl;
			IncompleteUniquePairedAlignment2detectFusion_Info incompleteUniquePairedAlignmentInfo;
			bool incompleteUniquePaired_bool
				= incompleteUniquePairedAlignmentInfo.initiateWith2samStr_globalMapOuterUnfixedEnd2detectFusionBreakPoint(
					inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
			RefineFusionBreakPoint_Info tmpRefineInfo;
			bool leftReadHeadUnfixed_bool = false;
			PE_Read_Info readInfo_unfixedLeftReadHead;
			FixPhase1Info fixPhase1Info_unfixedLeftReadHead;
			PE_Read_Alignment_Info peAlignInfo_unfixedLeftReadHead;

			bool rightReadTailUnfixed_bool = false;
			PE_Read_Info readInfo_unfixedRightReadTail;
			FixPhase1Info fixPhase1Info_unfixedRightReadTail;
			PE_Read_Alignment_Info peAlignInfo_unfixedRightReadTail;

			if(incompleteUniquePaired_bool)
			{
				leftReadHeadUnfixed_bool = incompleteUniquePairedAlignmentInfo.leftReadHeadUnfixedOrNot();
				rightReadTailUnfixed_bool = incompleteUniquePairedAlignmentInfo.rightReadTailUnfixedOrNot();
				cout << "leftReadHeadUnfixed_bool: " << leftReadHeadUnfixed_bool << endl;
				cout << "rightReadTailUnfixed_bool: " << rightReadTailUnfixed_bool << endl;				
				if(leftReadHeadUnfixed_bool)
				{
					string tmpReadName_1 = incompleteUniquePairedAlignmentInfo.returnReadName_1();
					string tmpUnfixedLeftReadHead_readSeq = incompleteUniquePairedAlignmentInfo.returnUnfixedLeftReadHead_readSeq();
					string tmpUnfixedLeftReadHead_qualSeq = incompleteUniquePairedAlignmentInfo.returnUnfixedLeftReadHead_qualSeq(fasta_or_fastq_bool);
					string tmpStr;
					readInfo_unfixedLeftReadHead.initiateReadInfo(tmpReadName_1, tmpStr, tmpUnfixedLeftReadHead_readSeq, tmpStr, 
						tmpUnfixedLeftReadHead_qualSeq, tmpStr, fasta_or_fastq_bool, true);
					fixPhase1Info_unfixedLeftReadHead.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
						preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray, 
						readInfo_unfixedLeftReadHead, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
					fixPhase1Info_unfixedLeftReadHead.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedLeftReadHead, 
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
					fixPhase1Info_unfixedLeftReadHead.fixPhase1_gapInfo(readInfo_unfixedLeftReadHead, indexInfo, Do_cirRNA,
						Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
					cout << "output map info: " << endl; fixPhase1Info_unfixedLeftReadHead.coutDebugInfo(readInfo_unfixedLeftReadHead, indexInfo, SE_or_PE_bool);
					peAlignInfo_unfixedLeftReadHead.initiatePeAlignInfo(fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor1,
						fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm1, fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor2,
						fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm2, indexInfo, true);
					cout << endl << endl << "peAlignInfo: \n" << peAlignInfo_unfixedLeftReadHead.returnPeAlignInfoStr() << endl << endl;
					peAlignInfo_unfixedLeftReadHead.chooseBestAlignment_final_SE();
				}
				if(rightReadTailUnfixed_bool)
				{
					string tmpReadName_2 = incompleteUniquePairedAlignmentInfo.returnReadName_2();
					string tmpUnfixedRightReadTail_readSeq = incompleteUniquePairedAlignmentInfo.returnUnfixedRightReadTail_readSeq();
					string tmpUnfixedRightReadTail_qualSeq = incompleteUniquePairedAlignmentInfo.returnUnfixedRightReadTail_qualSeq(fasta_or_fastq_bool);
					string tmpStr;
					readInfo_unfixedRightReadTail.initiateReadInfo(tmpReadName_2, tmpStr, tmpUnfixedRightReadTail_readSeq, tmpStr, 
						tmpUnfixedRightReadTail_qualSeq, tmpStr, fasta_or_fastq_bool, true);
					fixPhase1Info_unfixedRightReadTail.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
						preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray, 
						readInfo_unfixedRightReadTail, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
					fixPhase1Info_unfixedRightReadTail.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedRightReadTail, 
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
					fixPhase1Info_unfixedRightReadTail.fixPhase1_gapInfo(readInfo_unfixedRightReadTail, indexInfo, Do_cirRNA,
						Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
					cout << "output map info: " << endl; fixPhase1Info_unfixedRightReadTail.coutDebugInfo(readInfo_unfixedLeftReadHead, indexInfo, SE_or_PE_bool);
					peAlignInfo_unfixedRightReadTail.initiatePeAlignInfo(fixPhase1Info_unfixedRightReadTail.pathInfo_Nor1,
						fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm1, fixPhase1Info_unfixedRightReadTail.pathInfo_Nor2,
						fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm2, indexInfo, true);
					cout << endl << endl << "peAlignInfo: \n" << peAlignInfo_unfixedRightReadTail.returnPeAlignInfoStr() << endl << endl;
					peAlignInfo_unfixedRightReadTail.chooseBestAlignment_final_SE();
				}
				cout << "start to refine fusionBreakPoint" << endl;
				tmpRefineInfo.refineFusionBreakPoint(InsertedBaseNumMaxBetweenFusedTranscript,
					min_fusion_distance, incompleteUniquePairedAlignmentInfo, 
					leftReadHeadUnfixed_bool, readInfo_unfixedLeftReadHead, peAlignInfo_unfixedLeftReadHead,
					rightReadTailUnfixed_bool, readInfo_unfixedRightReadTail, peAlignInfo_unfixedRightReadTail, indexInfo);
				cout << endl << "end of refining fusion break point" << endl;
			}

			outputOriSamStrVec[tmpOpenMP] = "";			
			outputNonfusionSamStrVec[tmpOpenMP] = "";
			outputFusionSamStrVec[tmpOpenMP] = "";
			outputFusionBreakPointStrVec[tmpOpenMP] = "";			

			if(incompleteUniquePaired_bool)
			{
				bool leftReadHeadUnfixed_fixed_bool = tmpRefineInfo.return_leftReadHeadUnfixed_fixed_bool();
				bool rightReadTailUnfixed_fixed_bool = tmpRefineInfo.return_rightReadTailUnfixed_fixed_bool();
				if(leftReadHeadUnfixed_fixed_bool || rightReadTailUnfixed_fixed_bool)
				{
					outputFusionSamStrVec[tmpOpenMP] = tmpRefineInfo.returnFusionSamStr();
					outputFusionBreakPointStrVec[tmpOpenMP] = tmpRefineInfo.returnFusionBreakPointStr();
				}
				else
					outputNonfusionSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
			}
			else
				outputOriSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];

			fixPhase1Info_unfixedLeftReadHead.memoryFree();
			fixPhase1Info_unfixedRightReadTail.memoryFree();
			peAlignInfo_unfixedLeftReadHead.memoryFree();
			peAlignInfo_unfixedRightReadTail.memoryFree();	
		} // read file end
		
		#ifdef CAL_TIME
		align_end = clock();
		align_cost = align_cost + align_end - align_begin;
		output_begin = clock();
		#endif

		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << "start to output ... turn: " << tmpTurn+1 << endl;
	
		// for(int tmp = 0; tmp < realRecordNum; tmp++)
		// {	
		// 	if(outputOriSamStrVec[tmp] != "")
		// 		oriSAM_ofs << outputOriSamStrVec[tmp] << endl;		
		// 	if(outputNonfusionSamStrVec[tmp] != "")
		// 		nonFusionSAM_ofs << outputNonfusionSamStrVec[tmp] << endl;
		// 	if(outputFusionSamStrVec[tmp] != "")
		// 		fusionSAM_ofs << outputFusionSamStrVec[tmp] << endl;
		// 	if(outputFusionSamStrVec_withBreakPoint[tmp] != "")
		// 		fusionSAM_withBreakPoint_ofs << outputFusionSamStrVec_withBreakPoint[tmp] << endl;
		// 	//if(outputFusionSamStrVec_withBreakPoint_adjusted[tmp] != "")
		// 	//	fusionSAM_withBreakPoint_adjusted_ofs << outputFusionSamStrVec_withBreakPoint_adjusted[tmp] << endl;
		// 	if(outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmp] != "")
		// 		fusionSAM_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering_ofs
		// 			<< outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmp] << endl;
		// 	if(outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmp] != "")
		// 		fusionSAM_withBreakPoint_adjusted_withFusionAnchorSeq_ofs 
		// 			<< outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmp] << endl;
		// 	if(outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmp] != "")
		// 		fusionSAM_withBreakPoint_adjusted_withMapRegion_ofs
		// 			<< outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmp] << endl;
		// }
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
		cout << endl << "[" << asctime(local) 
			<< "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}
	settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	runtime_log_ofs << endl << "[" << asctime(local) 
		<< "finish running globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion ! " << endl << endl;	
	cout << endl << "[" << asctime(local) 
		<< "finish running globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion ! " << endl << endl;	

	incompleteUniquePairedAlignment_ifs.close();
	delete indexInfo;
	free(preIndexMapLengthArray);
	free(preIndexIntervalStartArray);
	free(preIndexIntervalEndArray);
	free(sa); 
	free(childTab);
	free(lcpCompress);
	free(verifyChild);
	oriSAM_ofs.close();
	nonFusionSAM_ofs.close();
	fusionSAM_ofs.close();
	fusionSAM_withBreakPoint_ofs.close();
	//fusionSAM_withBreakPoint_adjusted_ofs.close();
	fusionSAM_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering_ofs.close();
	fusionSAM_withBreakPoint_adjusted_withFusionAnchorSeq_ofs.close();
	fusionSAM_withBreakPoint_adjusted_withMapRegion_ofs.close();

	SA_file_ifs.close();
	lcpCompress_file_ifs.close();
	childTab_file_ifs.close();
	verifyChild_file_ifs.close();
	preIndexMapLengthArray_ifs.close();
	preIndexIntervalStartArray_ifs.close();
	preIndexIntervalEndArray_ifs.close();
	parameter_file_ifs.close();
	chrom_bit_file_ifs.close();
	stats_ofs.close();
	runtime_log_ofs.close();
	settings_log_ofs.close();
    log_ofs.close();
	return 0;
}