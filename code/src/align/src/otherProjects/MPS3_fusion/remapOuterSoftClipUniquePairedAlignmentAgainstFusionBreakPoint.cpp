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
//#include "../../general/otherFunc2.h"
#include "general/fusionBreakPointHash_info.h"
#include "general/incompleteUniquePairedAlignment2detectFusion_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexInfoPath inputFusionBreakPointPath ";
		cout << " inputOuterSoftClipPairedAlignmentPath ";
		cout << " outputFolder threads_num fasta_or_fastq " << endl;
		exit(1);
	}

	string threads_num_str = argv[5];
	int threads_num = atoi(threads_num_str.c_str());

	int normalRecordNum_1stMapping = 200000;
	bool fasta_or_fastq_bool = true;
	string fasta_or_fastq_str = argv[6];
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "Please set the correct format, fasta or fastq" << endl;
		exit(1);
	}	

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	string outputDirStr = outputFolderStr + "/";
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
   	string logStr = outputDirStr + "/log.txt";
   	ofstream log_ofs(logStr.c_str());
	log_ofs << "Command Line:";
	for(int tmp = 0; tmp < argc; tmp++)
	{
		log_ofs << "\t" << argv[tmp];
	}
	log_ofs << endl << "************************" << endl;

	string globalIndexStr = argv[1];
	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	log_ofs << "index: " << indexStr << endl;
	int chromNum = indexInfo->returnChromNum();
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


   	string oriSAM_outputPath = outputDirStr + "/kept.sam";
   	string nonFusionSAM_outputPath = outputDirStr + "/nonFusion.sam";
   	string fusionSAM_outputPath = outputDirStr + "/fusion.sam";
   	string fusionReadsWithBreakPoint_outputPath = outputDirStr + "/fusion_breakPointInfo.txt";
   	string updatedFusionBreakPointHash_outputPath = outputDirStr + "/fusionBreakPointHashInfo_updated.txt";
   	//string fusionSAM_outputPath_withBreakPoint_raw = outputDirStr + "fusionReadsWithBreakPoint_raw.txt";
   	//string fusionSAM_outputPath_withBreakPoint_adjusted = outputDirStr + "fusionReadsWithBreakPoint_adjusted.txt";
   	ofstream oriSAM_ofs(oriSAM_outputPath.c_str());
   	ofstream nonFusionSAM_ofs(nonFusionSAM_outputPath.c_str());
   	ofstream fusionSAM_ofs(fusionSAM_outputPath.c_str());
   	ofstream fusionReadsWithBreakPoint_ofs(fusionReadsWithBreakPoint_outputPath.c_str());
   	//ofstream updatedFusionBreakPointHash_ofs(updatedFusionBreakPointHash_outputPath.c_str());
   	//ofstream fusionSAM_withBreakPoint_adjusted_ofs(fusionSAM_outputPath_withBreakPoint_adjusted.c_str());

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
	cout << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash ......" << endl << endl; 

	FusionBreakPointHash_Info* fusionBreakPointHashInfo_merged = new FusionBreakPointHash_Info();
	fusionBreakPointHashInfo_merged->initiateWithChromNum(chromNum);
	string inputFusionBreakPointPath = argv[2];
	fusionBreakPointHashInfo_merged->generateFusionBreakPointHashInfo_fromFuionJuncFile(
		inputFusionBreakPointPath, indexInfo);

	cout << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHashInfoVec ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHashInfoVec ......" << endl << endl;	

	vector<FusionBreakPointHash_Info*> fusionBreakPointHashInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo 
			= new FusionBreakPointHash_Info();
		tmpFusionBreakPointHashInfo->initiateWithChromNum(chromNum);
		//tmpFusionBreakPointHashInfo->copyFrom(fusionBreakPointHashInfo_merged); // toImplement
		tmpFusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromFuionJuncFile(
			inputFusionBreakPointPath, indexInfo);
		tmpFusionBreakPointHashInfo->clearSupportNum();
		tmpFusionBreakPointHashInfo->clearSupportNum_encompassing();
		fusionBreakPointHashInfoVec.push_back(tmpFusionBreakPointHashInfo);
	}

	string inputIncompleteUniquePairedAlignmentPath = argv[3];
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
	vector<string> outputFusionBreakPointInfoStrVec(normalRecordNum);
	//vector<string> outputFusionSamStrVec_withBreakPoint_raw(normalRecordNum);
	//vector<string> outputFusionSamStrVec_withBreakPoint_adjusted(normalRecordNum);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) 
		<< "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) 
		<< "... 1st mapping process starts ......" << endl << endl; 

	int readTotalNum = 0;
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

		cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		

		cout << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		cout << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << "realRecordNum: " << realRecordNum << endl;
		cout << "threads_num: " << threads_num << endl;
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			//cout << "start to get threadNO " << endl;
			int threadNO = omp_get_thread_num();
			IncompleteUniquePairedAlignment2detectFusion_Info incompleteUniquePairedAlignmentInfo;
			//cout << endl << endl << "*********************************************" << endl;
			//cout << "start to initiate incompleteUniquePairedAlignmentInfo" << endl;
			bool incompleteUniquePaired_bool
				= incompleteUniquePairedAlignmentInfo.initiateWith2samStr_remapOuterUnfixedEndAgainstFusionBreakPoint(
					inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
			//cout << "inputSAM_1: " << endl << inputSamStr1Vec[tmpOpenMP] << endl;
			//cout << "inputSAM_2: " << endl << inputSamStr2Vec[tmpOpenMP] << endl;
			//cout << "incompleteUniquePaired_bool: " << incompleteUniquePaired_bool << endl;
			bool leftReadHeadUnfixed_bool = false;
			bool rightReadTailUnfixed_bool = false;
			bool leftReadHeadUnfixed_fixed_bool = false;
			bool rightReadTailUnfixed_fixed_bool = false;
			//bool fusionRemappingFixed_bool = false;
			RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1
				= new RemapAgainstFusionBreakPoint_Info();
			RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2
				= new RemapAgainstFusionBreakPoint_Info();
			RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1
				= new RemapAgainstFusionBreakPoint_Info();
			RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2
				= new RemapAgainstFusionBreakPoint_Info();				
			if(incompleteUniquePaired_bool)
			{
				leftReadHeadUnfixed_bool 
					= incompleteUniquePairedAlignmentInfo.leftReadHeadUnfixedOrNot();
				rightReadTailUnfixed_bool
					= incompleteUniquePairedAlignmentInfo.rightReadTailUnfixed_bool();
				//cout << "leftReadHeadUnfixed_bool: " << leftReadHeadUnfixed_bool << endl;
				//cout << "rightReadTailUnfixed_bool: " << rightReadTailUnfixed_bool << endl;
				if(leftReadHeadUnfixed_bool)
				{	
					incompleteUniquePairedAlignmentInfo.remapUnfixedHeadAgainstFusionBreakPoint_leftRead(
						tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1,
						tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2,
						fusionBreakPointHashInfo_merged, indexInfo);
					leftReadHeadUnfixed_fixed_bool 				
						= (tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1->returnResultsSize()
							+ tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2->returnResultsSize() == 1);
					//cout << "leftReadHeadUnfixed_fixed_bool: " << leftReadHeadUnfixed_fixed_bool << endl;
				}
				if(rightReadTailUnfixed_bool)
				{
					incompleteUniquePairedAlignmentInfo.remapUnfixedTailAgainstFusionBreakPoint_rightRead(
						tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1,
						tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2,
						fusionBreakPointHashInfo_merged, indexInfo);
					rightReadTailUnfixed_fixed_bool
						= (tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1->returnResultsSize()
							+ tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2->returnResultsSize() == 1);
					//cout << "rightReadTailUnfixed_fixed_bool: " << rightReadTailUnfixed_fixed_bool << endl; 
				}
			}
			//cout << "start to get results str " << endl;
			outputOriSamStrVec[tmpOpenMP] = "";			
			outputNonfusionSamStrVec[tmpOpenMP] = "";
			outputFusionSamStrVec[tmpOpenMP] = "";
			outputFusionBreakPointInfoStrVec[tmpOpenMP] = "";
			bool fusionRemappingFixed_bool = (leftReadHeadUnfixed_fixed_bool||rightReadTailUnfixed_fixed_bool);
			if(incompleteUniquePaired_bool)
			{
				if(fusionRemappingFixed_bool)
				{
					string tmpFusionSamStr
						= "fixME:currentlyNull";
						// = incompleteUniquePairedAlignmentInfo.returnFixedFusionSamStr(
						// 	leftReadHeadUnfixed_fixed_bool, 
						// 	tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1,
						// 	tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2, 
						// 	rightReadTailUnfixed_fixed_bool, 
						// 	tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1,
						// 	tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2,
						// 	indexInfo);
					outputFusionSamStrVec[tmpOpenMP] = tmpFusionSamStr;
					//cout << "start to get tmpFusionBreakPointInfoStr. .." << endl;
					string tmpFusionBreakPointInfoStr
						= incompleteUniquePairedAlignmentInfo.returnFusionBreakPointInfoStr_supportNumIncrement(
							leftReadHeadUnfixed_fixed_bool, 
							tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1,
							tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2, 
							rightReadTailUnfixed_fixed_bool, 
							tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1,
							tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2,
							fusionBreakPointHashInfoVec[threadNO], indexInfo);
					//cout << "end of getting tmpFusionBreakPointInfoStr ..." << endl;
					outputFusionBreakPointInfoStrVec[tmpOpenMP] = tmpFusionBreakPointInfoStr;
					//cout << "end of inserting into outputFusionBreakPointInfoStrVec" << endl;
					// string tmpFusionBreakPointInfoStr_adjusted
					// 	= incompleteUniquePairedAlignmentInfo.returnFusionBreakPointInfoStr_adjusted(
					// 		leftReadHeadUnfixed_fixed_bool, 
					// 		tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1,
					// 		tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2, 
					// 		rightReadTailUnfixed_fixed_bool, 
					// 		tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1,
					// 		tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2,
					// 		indexInfo);
					// outputFusionSamStrVec_withBreakPoint_raw[tmpOpenMP] = tmpFusionBreakPointInfoStr_raw;
					// outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP] = tmpFusionBreakPointInfoStr_adjusted;
				}
				else
				{
					outputNonfusionSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
				}
			}
			else
			{
				outputOriSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
			}

			delete tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene1;
			delete tmpRemapInfo_unfixedHead_leftRead_oriGeneAsGene2;
			delete tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene1;
			delete tmpRemapInfo_unfixedTail_rightRead_oriGeneAsGene2;
		}	

		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << "start to output ... turn: " << tmpTurn+1 << endl;

		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{	
			if(outputOriSamStrVec[tmp] != "")
				oriSAM_ofs << outputOriSamStrVec[tmp] << endl;		
			if(outputNonfusionSamStrVec[tmp] != "")
				nonFusionSAM_ofs << outputNonfusionSamStrVec[tmp] << endl;
			if(outputFusionSamStrVec[tmp] != "")
				fusionSAM_ofs << outputFusionSamStrVec[tmp] << endl;
			if(outputFusionBreakPointInfoStrVec[tmp] != "")
				fusionReadsWithBreakPoint_ofs << outputFusionBreakPointInfoStrVec[tmp] << endl;
			// if(outputFusionSamStrVec_withBreakPoint_adjusted[tmp] != "")
			// 	fusionSAM_withBreakPoint_ofs << outputFusionSamStrVec_withBreakPoint_adjusted[tmp] << endl;
		}
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}

	cout << "start to do support number sum up" << endl;
	// sum up support numbers
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		fusionBreakPointHashInfo_merged->addSupportNumFromAnotherFusionBreakPointHashInfo(
			fusionBreakPointHashInfoVec[tmp]);
	cout << "start to do anchor size update" << endl;
	// update anchor size
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		fusionBreakPointHashInfo_merged->updateAnchorSizeFromAnotherFusionBreakPointHashInfo(
			fusionBreakPointHashInfoVec[tmp]);	
	cout << "start to do XM update" << endl;
	// update XM
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		fusionBreakPointHashInfo_merged->updateXMfromAnotherFusionBreakPointHashInfo(
			fusionBreakPointHashInfoVec[tmp]);	
	cout << "start to output fusionBreakPointHash" << endl;
	// output
	fusionBreakPointHashInfo_merged->outputFusionBreakPointHashInfoStr(
		updatedFusionBreakPointHash_outputPath, indexInfo);

	cout << "start to do memory free !" << endl;

	fusionBreakPointHashInfo_merged->memoryFree();
	delete fusionBreakPointHashInfo_merged;
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
	{
		fusionBreakPointHashInfoVec[tmp]->memoryFree();
		delete fusionBreakPointHashInfoVec[tmp];
	}	
   	log_ofs.close();
   	delete indexInfo;
	cout << "end of running remapOutSoftClipUniquePairedAlignmentAgainstFusionBreakPoint !" << endl;
	return 0;
}

