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
#include "general/uniqueUnpairedAlignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputIndexInfoPath inputFusionBreakPointPath ";
		cout << " inputUniqueUnpairedReads ";
		cout << " outputFolder threads_num fasta_or_fastq strandedOrNotStranded" << endl;
		exit(1);
	}
	time_t nowtime;
	struct tm *local;

	bool strandedOrNotStranded_bool;
	string strandedOrNotStrandedStr = argv[7];
	if((strandedOrNotStrandedStr == "stranded")||(strandedOrNotStrandedStr == "Stranded"))
		strandedOrNotStranded_bool = true;
	else if((strandedOrNotStrandedStr == "nonstranded")||(strandedOrNotStrandedStr == "nonStranded")
		||(strandedOrNotStrandedStr == "Nonstranded")||(strandedOrNotStrandedStr == "NonStranded"))
		strandedOrNotStranded_bool = false;
	else
	{
		cout << "error in command, strandedOrNotStrandedStr: " << strandedOrNotStrandedStr << endl;
		cout << "it should be: stranded or nonStranded" << endl;
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
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	int chromNum = indexInfo->returnChromNum();

   	string oriSAM_outputPath = outputDirStr + "/kept.sam";
   	string nonFusionSAM_outputPath = outputDirStr + "/nonFusion.sam";
   	string fusionEncompassingSAM_outputPath = outputDirStr + "/fusionEncompassing.sam";
   	string fusionEncompassingReads_breakPoint_outputPath = outputDirStr + "/fusionEncompassing_breakPoint.txt";
   	string updatedFusionBreakPointHash_outputPath = outputDirStr + "/fusionBreakPointHashInfo_updated.txt";
   	//string fusionSAM_outputPath_withBreakPoint_adjusted = outputDirStr + "fusionReadsWithBreakPoint_adjusted.txt";
   	ofstream oriSAM_ofs(oriSAM_outputPath.c_str());
   	ofstream nonFusionSAM_ofs(nonFusionSAM_outputPath.c_str());
   	ofstream fusionEncompassingSAM_ofs(fusionEncompassingSAM_outputPath.c_str());
   	ofstream fusionEncompassingReadsWithreakPoint_ofs(fusionEncompassingReads_breakPoint_outputPath.c_str());
   	//ofstream fusionSAM_withBreakPoint_adjusted_ofs(fusionSAM_outputPath_withBreakPoint_adjusted.c_str());

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to generate FusionBreakPointHash_merged ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash_merged ......" << endl << endl; 

	FusionBreakPointHash_Info* fusionBreakPointHashInfo_merged
		= new FusionBreakPointHash_Info();
	fusionBreakPointHashInfo_merged->initiateWithChromNum(chromNum);
	string inputFusionBreakPointPath = argv[2];
	if(strandedOrNotStranded_bool)
		fusionBreakPointHashInfo_merged->generateFusionBreakPointHashInfo_fromFuionJuncFile(
			inputFusionBreakPointPath, indexInfo);
	else
		fusionBreakPointHashInfo_merged->generateFusionBreakPointHashInfo_fromNonStrandedFuionJuncFile(
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
		if(strandedOrNotStranded_bool)
			tmpFusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromFuionJuncFile(
					inputFusionBreakPointPath, indexInfo);
		else
			tmpFusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromNonStrandedFuionJuncFile(
					inputFusionBreakPointPath, indexInfo);			
		tmpFusionBreakPointHashInfo->clearSupportNum();
		tmpFusionBreakPointHashInfo->clearSupportNum_encompassing();
		fusionBreakPointHashInfoVec.push_back(tmpFusionBreakPointHashInfo);
	}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to count encompassing unique unpaired reads ......" << endl << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to count encompassing unique unpaired reads ......" << endl << endl;

	string inputUniqueUnpairedAlignmentPath = argv[3];
	ifstream uniqueUnpairedAlignments_ifs(inputUniqueUnpairedAlignmentPath.c_str());
	int normalRecordNum = normalRecordNum_1stMapping;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	vector<string> inputSamStr1Vec(normalRecordNum);
	vector<string> inputSamStr2Vec(normalRecordNum);
	vector<string> outputOriSamStrVec(normalRecordNum);
	vector<string> outputNonfusionSamStrVec(normalRecordNum);
	vector<string> outputFusionEncompassingSamStrVec(normalRecordNum);
	vector<string> outputFusionEncompassingReadsWithBreakPointStrVec(normalRecordNum);
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
    		if(uniqueUnpairedAlignments_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
    		}
    		getline(uniqueUnpairedAlignments_ifs, line1); // readName_1
    		if(uniqueUnpairedAlignments_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
    		}
    		inputSamStr1Vec[recordNumTmp] = line1;
    		getline(uniqueUnpairedAlignments_ifs, line2);		
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

			UniqueUnpairedAlignment_Info tmpUniqueUnpairedAlignmentInfo;
			//cout << endl << endl << "************************************************" << endl;
			//cout << "start to initiateWith2samStr ...." << endl;
			bool uniqueUnpaired_bool
				= tmpUniqueUnpairedAlignmentInfo.initiateWith2samStr(
					inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
			//cout << "uniqueUnpaired_bool: " << uniqueUnpaired_bool << endl;	
			int fusionBreakPointEncompassed_found_index = -1;
			if(uniqueUnpaired_bool)
			{
				//cout << "readName: " << tmpUniqueUnpairedAlignmentInfo.returnReadName_withoutSignForEnd() << endl;
				//cout << "start to searchFusionBreakPointEncompassedAndAddSupport ......" << endl;
				fusionBreakPointEncompassed_found_index
					= tmpUniqueUnpairedAlignmentInfo.searchFusionBreakPointEncompassedAndAddSupport(
						fusionBreakPointHashInfo_merged);
				//cout << "fusionBreakPointEncompassed_found_index: " << fusionBreakPointEncompassed_found_index << endl;
				//cout << "start to do supportNumIncrementWithIndex_encompassing ...." << endl;
				if(fusionBreakPointEncompassed_found_index >= 0)
					fusionBreakPointHashInfoVec[threadNO]->supportNumIncrementWithIndex_encompassing(
						fusionBreakPointEncompassed_found_index);
			}
			//cout << "start to generate results strVec....." << endl;
			outputOriSamStrVec[tmpOpenMP] = "";			
			outputNonfusionSamStrVec[tmpOpenMP] = "";
			outputFusionEncompassingSamStrVec[tmpOpenMP] = "";
			outputFusionEncompassingReadsWithBreakPointStrVec[tmpOpenMP] = "";
			if(uniqueUnpaired_bool)
			{
				if(fusionBreakPointEncompassed_found_index >= 0)
				{
					string tmpFusionEncompassingSamStr 
						= inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
					outputFusionEncompassingSamStrVec[tmpOpenMP] = tmpFusionEncompassingSamStr;
					
					string tmpReadName = tmpUniqueUnpairedAlignmentInfo.returnReadName_withoutSignForEnd();
					string tmpChrNameStr_gene1 = fusionBreakPointHashInfoVec[threadNO]->returnChrNameStrWithIndex_1(
						fusionBreakPointEncompassed_found_index, indexInfo);
					string tmpChrNameStr_gene2 = fusionBreakPointHashInfoVec[threadNO]->returnChrNameStrWithIndex_2(
						fusionBreakPointEncompassed_found_index, indexInfo);
					int tmpBreakPointPos_gene1 = fusionBreakPointHashInfoVec[threadNO]->returnBreakPointWithIndex_1(
						fusionBreakPointEncompassed_found_index);
					int tmpBreakPointPos_gene2 = fusionBreakPointHashInfoVec[threadNO]->returnBreakPointWithIndex_2(
						fusionBreakPointEncompassed_found_index);
					string tmpStrand_gene1 = fusionBreakPointHashInfoVec[threadNO]->returnStrandWithIndex_1(
						fusionBreakPointEncompassed_found_index);
					string tmpStrand_gene2 = fusionBreakPointHashInfoVec[threadNO]->returnStrandWithIndex_2(
						fusionBreakPointEncompassed_found_index);
					string tmpFusionEncompassingReadsWithBreakPointInfoStr
						= tmpReadName + "\t" + tmpChrNameStr_gene1 + "\t" + tmpChrNameStr_gene2 + "\t"
							+ int_to_str(tmpBreakPointPos_gene1) + "\t" + int_to_str(tmpBreakPointPos_gene2) + "\t"
							+ tmpStrand_gene1 + "\t" + tmpStrand_gene2;
					// fusionBreakPointHashInfoVec[threadNO]->addSupportNum_withChrNameStrBreakPointPos_encompassing(
					// 	1, tmpChrNameStr_gene1, tmpChrNameStr_gene2, tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, indexInfo);
					//fusionBreakPointHashInfoVec[threadNO]->addSupportNum_withIndexInBreakPointInfoVec_encompassing(1, fusionBreakPointEncompassed_found_index);		
					outputFusionEncompassingReadsWithBreakPointStrVec[tmpOpenMP]
						= tmpFusionEncompassingReadsWithBreakPointInfoStr;
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
			if(outputFusionEncompassingSamStrVec[tmp] != "")
				fusionEncompassingSAM_ofs << outputFusionEncompassingSamStrVec[tmp] << endl;
			if(outputFusionEncompassingReadsWithBreakPointStrVec[tmp] != "")
				fusionEncompassingReadsWithreakPoint_ofs << outputFusionEncompassingReadsWithBreakPointStrVec[tmp] << endl;
		}
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}

	// sum up encompassing support numbers and output
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		fusionBreakPointHashInfo_merged->addSupportNumFromAnotherFusionBreakPointHashInfo_encompassing(
			fusionBreakPointHashInfoVec[tmp]);

	fusionBreakPointHashInfo_merged->outputFusionBreakPointHashInfoStr(
		updatedFusionBreakPointHash_outputPath, indexInfo);


	uniqueUnpairedAlignments_ifs.close();
	fusionBreakPointHashInfo_merged->memoryFree();
	delete fusionBreakPointHashInfo_merged;
	for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
	{
		fusionBreakPointHashInfoVec[tmp]->memoryFree();
		delete fusionBreakPointHashInfoVec[tmp];
	}
   	log_ofs.close();
   	delete indexInfo;
	return 0;
}

