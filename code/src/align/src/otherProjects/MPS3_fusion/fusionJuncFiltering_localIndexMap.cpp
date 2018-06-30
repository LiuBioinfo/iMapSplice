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
#include "general/incompleteUniquePairedAlignment2detectFusion_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputGlobalIndexPath inputLocalIndexPath ";
		//cout << "leftHead_or_rightTail";
		cout << "inputAdjustedBreakPointFileWithSAM inputJuncFile ";
		cout << "outputFolder threads_num fasta_or_fastq " << endl;
		exit(1);
	}
	int normalRecordNum = 1000000;

	string globalIndexPath = argv[1];
	string localIndexPath = argv[2];
	string inputAdjustedBreakPointFileWithSAM = argv[3];
	// /////////////////////////////////////////////////////////////////////////////////
	// /////////////////    	specify leftHead or rightTail      //////////////////////
	// /////////////////////////////////////////////////////////////////////////////////	
	// string leftHead_or_rightTail_str = argv[3];
	// bool leftHead_or_rightTail_bool = false;
	// if(leftHead_or_rightTail_str == "leftHead")
	// 	leftHead_or_rightTail_bool = true;
	// else if(leftHead_or_rightTail_str == "rightTail")
	// 	leftHead_or_rightTail_bool = false;
	// else
	// {
	// 	cout << "Please set the correct info, leftHead or rightTail" << endl;
	// 	exit(1);
	// }	
	/////////////////////////////////////////////////////////////////////////////////
	////////////////////////    	creating folders      ///////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
   	string settingsLogStr = outputFolderStr + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	settings_log_ofs << "CommandLine:" << endl;
   	for(int tmp = 1; tmp < argc; tmp++)
   		settings_log_ofs << "\t" << argv[tmp] << endl;
   	settings_log_ofs << endl;
   	string adjustedBreakPoint_afterLocalIndexMapFiltering_outputFilePath 
   		= outputFolderStr + "/adjustedBreakPointVec_afterLocalIndexMapFiltering.txt";
   	ofstream adjustedBreakPoint_afterLocalIndexMapFiltering_ofs(
   		adjustedBreakPoint_afterLocalIndexMapFiltering_outputFilePath.c_str());
   	string peAlignSAM_fusionOtherEnd_outputFilePath
   		= outputFolderStr + "/fusionOtherEnd.sam";
   	ofstream peAlignSAM_fusionOtherEnd_ofs(peAlignSAM_fusionOtherEnd_outputFilePath.c_str());
	/////////////////////////////////////////////////////////////////////////////////
	//////////////////    	check threads and fasta or fastq      ///////////////////
	/////////////////////////////////////////////////////////////////////////////////
	string threadNumStr = argv[6];
	int threads_num = atoi(threadNumStr.c_str());
	omp_set_num_threads(threads_num);
	bool fasta_or_fastq_bool = true;
	string fasta_or_fastq_str = argv[7];
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "Please set the correct format, fasta or fastq" << endl;
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////
	/////////////////////    	load global index      ///////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	string indexStr = globalIndexPath + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	/////////////////////////////////////////////////////////////////////////////
	////////////////////// set some initial parameters //////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	string juncfile_alignInferHash = outputFolderStr + "/juncFile_alignInferJuncHash.txt";
	bool spliceJunctionHashExists;
	int chromNum = indexInfo->returnChromNum();
	string inputJuncFile = argv[4];
	Annotation_Info* annotationInfo = new Annotation_Info();
	bool annotation_provided_bool = false;
	bool Do_annotation_only_bool = false;
	bool checkQualSeqForReadSegSeq = false;
	bool checkQualSeqForShortAnchorSeqToTargetMap = false;

	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();	
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_onlyChrNamePos(inputJuncFile, indexInfo);
	cout << "start to do convert2SJhashInfo" << endl;
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());		
	alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
	int junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
	if(junctionNum_in_alignInferJuncHashInfo == 0)
		spliceJunctionHashExists = false;
	settings_log_ofs << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
	cout << "start to output alignInferJunctionHashInfo ...." << endl;
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_chrNamePos_supportNum(indexInfo, juncfile_alignInferHash);
	delete alignInferJunctionHashInfo;

	bool Do_cirRNA = true;
	Do_cirRNA = false;
	bool Do_extendHeadTail_fixOneEndUnmapped = true;
	//Do_extendHeadTail_fixOneEndUnmapped = false;
	bool Do_extendHeadTail_fixHeadTail = true;
	//Do_extendHeadTail_fixHeadTail = false;	
	bool Do_fixHeadTail_remapping = true;
	//Do_fixHeadTail_remapping = false;
	bool Do_fixHeadTail_greedyMapping = true;
	//Do_fixHeadTail_greedyMapping = false;
	bool Do_fixHeadTail_remappingAndTargetMapping = true;
	//Do_fixHeadTail_remappingAndTargetMapping = false;
	bool Do_fixHeadTail_remappingAgain = true;
	Do_fixHeadTail_remappingAgain = false;
	bool Do_fixHeadTail_extend2end_finalStep = true;
	//Do_fixHeadTail_extend2end_finalStep = false;
	bool Do_fixHeadTail_extend2end_fixIndel = true;
	//Do_fixHeadTail_extend2end_fixIndel = false;




	/////////////////////////////////////////////////////////////////////////////
	/////////////////////    	load local index      ///////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;
	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;
	string secondLevelIndexStr = localIndexPath + "/";
	settings_log_ofs << "start to load second-level index ..." << endl;
	int secondLevelIndexNO = 0;
	for(int tmpChrNO = 0; tmpChrNO < indexInfo->returnChromNum(); tmpChrNO ++)
	{
		for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
		{
			char tmpFileNumChar[4];
			sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
			string tmpFileNumStr = tmpFileNumChar;			
			string inputIndexFileStr = secondLevelIndexStr + "/" + indexInfo->returnChrNameStr(tmpChrNO) + "/" + tmpFileNumStr + "/";

			string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
			ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
			string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
			ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);
			string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
			ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
			string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
			ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
			string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
			ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					
			
			int sizeOfIndex = indexInfo->returnSecondLevelIndexNormalSize() + 1;
			char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
			for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
			secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
			if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);

			bool No_ATGC_Bool = true;
			for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
			{
				char ch = tmpSecondLevelChrom[tmpMallocSpace];
				if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
				{
					No_ATGC_Bool = false;
					break;
				}
			}				
			if(No_ATGC_Bool)
				indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
			secondLevelChrom.push_back(tmpSecondLevelChrom);
			unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
			secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
			secondLevelSa.push_back(tmpSecondLevelSa);
			BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
			secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
			secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);					
			unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
			secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
			secondLevelChildTab.push_back(tmpSecondLevelChildTab);
			BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
			secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
			secondLevelDetChild.push_back(tmpSecondLevelDetChild);
			secondLevelChrom_file_ifs.close();
			secondLevelSA_file_ifs.close();
			secondLevelLcpCompress_file_ifs.close();
			secondLevelChildTab_file_ifs.close();
			secondLevelDetChild_file_ifs.close();							
			secondLevelIndexNO ++;
		}
		settings_log_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
	}
	settings_log_ofs << "finish loading ALL 2nd-level index !" << endl;
	settings_log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
	/////////////////////////////////////////////////////////////////////////////
	/////////////  filtering based on possibiligy of local index map ////////////
	/////////////////////////////////////////////////////////////////////////////
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;// = normalRecordNum;
	ifstream inputAdjustedBreakPointFileWithSAM_ifs(inputAdjustedBreakPointFileWithSAM.c_str());
	for(tmpTurn = 0; ; tmpTurn ++)
	{
		vector<string> line1StrVec(normalRecordNum);
		vector<string> line2StrVec(normalRecordNum);
		vector<string> line3StrVec(normalRecordNum);
		vector<string> line4StrVec(normalRecordNum);	
		vector<string> adjustedBreakPointVec_afterLocalIndexMapFiltering(normalRecordNum);
		vector<string> PeAlignSamVec_fusionOtherEnd(normalRecordNum);
		
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		settings_log_ofs << endl << "[" << asctime(local) << "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
		cout << endl << "[" << asctime(local) << "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;		
		realRecordNum = normalRecordNum;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
			string line1, line2, line3, line4;
			if(inputAdjustedBreakPointFileWithSAM_ifs.eof())
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;				
			}
			getline(inputAdjustedBreakPointFileWithSAM_ifs, line1);
			if(inputAdjustedBreakPointFileWithSAM_ifs.eof())
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;				
			}
			getline(inputAdjustedBreakPointFileWithSAM_ifs, line2);
			getline(inputAdjustedBreakPointFileWithSAM_ifs, line3);
			getline(inputAdjustedBreakPointFileWithSAM_ifs, line4);
			line1StrVec[recordNumTmp] = line1;
			line2StrVec[recordNumTmp] = line2;
			line3StrVec[recordNumTmp] = line3;
			line4StrVec[recordNumTmp] = line4;
		}
		
		settings_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		settings_log_ofs << endl << "[" << asctime(local) << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
		settings_log_ofs << endl << "[" << asctime(local) << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;					
		cout << endl << "[" << asctime(local) << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
		cout << endl << "[" << asctime(local) << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;				

		omp_set_num_threads(threads_num);
		#pragma omp parallel for			
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			int tmp_multiMapSeg_maxLength = 0;
			bool SE_or_PE_bool = false;
			string tmpAdjustedBreakPointInfoStr = line1StrVec[tmpOpenMP];
			string tmpOriPeSamStr_1 = line2StrVec[tmpOpenMP];
			string tmpOriPeSamStr_2 = line3StrVec[tmpOpenMP];
			string tmpDetectedFusionOtherEndSamStr = line4StrVec[tmpOpenMP];
			//cout << "tmpAdjustedBreakPointInfoStr: " << endl << tmpAdjustedBreakPointInfoStr << endl;
			//cout << "tmpOriPeSamStr_1: " << endl << tmpOriPeSamStr_1 << endl;
			//cout << "tmpOriPeSamStr_2: " << endl << tmpOriPeSamStr_2 << endl;
			//cout << "tmpDetectedFusionOtherEndSamStr: " << endl << tmpDetectedFusionOtherEndSamStr << endl;
			//AdjustedBreakPoint_Info* tmpAdjustedBreakPointInfo = new AdjustedBreakPoint_Info();
			//tmpAdjustedBreakPointInfo->initiate(tmpAdjustedBreakPointInfoStr, indexInfo);
			//cout << "start to generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_forLocalIndexMap2filterFusion ...." << endl;
			PE_Read_Info peReadInfo;
			PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
			peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_forLocalIndexMap2filterFusion(
				tmpOriPeSamStr_1, tmpOriPeSamStr_2, tmpDetectedFusionOtherEndSamStr, //leftHead_or_rightTail_bool,
				peReadInfo, indexInfo, fasta_or_fastq_bool, tmp_multiMapSeg_maxLength);
			cout << endl << "*********************************************" << endl;
			cout << "tmpReadName_1: " << peReadInfo.returnReadName_1() << endl;
			cout << "before Fixing peAlignInfo: " << endl;
			cout << peAlignInfo->returnPeAlignInfoStr() << endl;
			//cout << "start to do fixOneEndUnmappedInfo " << endl;
			FixOneEndUnmappedInfo* fixOneEndUnmappedInfo = new FixOneEndUnmappedInfo();
			fixOneEndUnmappedInfo->fixOneEndUnmapped(peReadInfo, peAlignInfo, secondLevelChrom, secondLevelSa,
					secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
					Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, 
					annotationInfo, MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq);			

			FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
			//cout << "start to do fixHeadTail_areaAndStringHash_new_remappingOnly " << endl;
			if(Do_fixHeadTail_remapping)
			{
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa, secondLevelLcpCompress,
					secondLevelChildTab, secondLevelDetChild, spliceJunctionHashExists,
					indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);
				// fixHeadTailInfo->fixHeadTail_shortAnchorRemappingOnly_withAlignInfer(
				// 	peReadInfo, peAlignInfo, SJ, spliceJunctionHashExists,
				// 	indexInfo, alignInferJunctionHashInfo);
			}
			//cout << "start to do fixHeadTail_areaAndStringHash_new_greedyMappingOnly " << endl;		
			if(Do_fixHeadTail_greedyMapping)
			{
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, 
					secondLevelChildTab, secondLevelDetChild, spliceJunctionHashExists, indexInfo, 
					Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool);
			}
			//cout << "start to do fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping" << endl;
			if(Do_fixHeadTail_remappingAndTargetMapping)
			{
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, 
					secondLevelChildTab, secondLevelDetChild, spliceJunctionHashExists, indexInfo, 
					Do_extendHeadTail_fixHeadTail, checkQualSeqForShortAnchorSeqToTargetMap, SE_or_PE_bool);	
			}
			if(Do_fixHeadTail_remappingAgain)
			{
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa, secondLevelLcpCompress,
					secondLevelChildTab, secondLevelDetChild, spliceJunctionHashExists,
					indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);	
			}
			if(Do_fixHeadTail_extend2end_finalStep)
				fixHeadTailInfo->fixHeadTail_extend2end_finalStepForAligner(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);
			if(Do_fixHeadTail_extend2end_fixIndel)
				fixHeadTailInfo->fixHeadTail_extend2end_fixIndel(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);	
			peAlignInfo->removeDuplicateMismatch(SE_or_PE_bool);
			peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();

			PeAlignSamVec_fusionOtherEnd[tmpOpenMP] = "";

			bool pairExistsBool = peAlignInfo->finalPairExistsBool();
			cout << "after Fixing peAlignInfo: " << endl;
			cout << peAlignInfo->returnPeAlignInfoStr() << endl;
			if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
			{
				// bool align_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
				// if(align_lowScore_bool)
				// 	PeAlignSamVec_fusionOtherEnd[tmpOpenMP] = peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);
				// else
					PeAlignSamVec_fusionOtherEnd[tmpOpenMP] = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool, tmp_multiMapSeg_maxLength);	
			}
			else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				PeAlignSamVec_fusionOtherEnd[tmpOpenMP] = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo, fasta_or_fastq_bool);

			//delete tmpAdjustedBreakPointInfo;
			delete fixOneEndUnmappedInfo;
			delete fixHeadTailInfo;
		}

		nowtime = time(NULL);
		local = localtime(&nowtime);
		settings_log_ofs << endl << "[" << asctime(local) << "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
		settings_log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		cout << endl << "[" << asctime(local) << "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			if(adjustedBreakPointVec_afterLocalIndexMapFiltering[tmp] != "")
				adjustedBreakPoint_afterLocalIndexMapFiltering_ofs << adjustedBreakPointVec_afterLocalIndexMapFiltering[tmp] << endl;
			if(PeAlignSamVec_fusionOtherEnd[tmp] != "")
				peAlignSAM_fusionOtherEnd_ofs << PeAlignSamVec_fusionOtherEnd[tmp] << endl;
		}

		nowtime = time(NULL);
		local = localtime(&nowtime);
		settings_log_ofs << endl << "[" << asctime(local) << "finish output, turn: " << tmpTurn+1 << endl << endl;
		cout << endl << "[" << asctime(local) << "finish output, turn: " << tmpTurn+1 << endl << endl;
	}

	cout << "end of fusionJuncFiltering_localIndexMap" << endl;
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl;  
	settings_log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl;  

	peAlignSAM_fusionOtherEnd_ofs.close();
	adjustedBreakPoint_afterLocalIndexMapFiltering_ofs.close();
	inputAdjustedBreakPointFileWithSAM_ifs.close();

	delete annotationInfo;
	delete indexInfo;
	return 0;
}