// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef MAPSPLICE_FUSION_INFO_H
#define MAPSPLICE_FUSION_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "fusionBreakPointHash_info.h"
#include "uniqueUnpairedAlignment_info.h"
#include "remapAgainstFusionBreakPoint_info.h"
#include "incompleteUniquePairedAlignment2detectFusion_info.h"

using namespace std;

typedef map<int, int> FusionJuncEndPos2fusionInfoVecIndexMap;
typedef map<int, FusionJuncEndPos2fusionInfoVecIndexMap > FusionJuncPosPair2fusionInfoVecIndexMap;

class MapSplice_Fusion_Info
{
private:
	// initialMapResultsFolder
	string input_folder_path;
	// outputFolder
	string output_folder_path;

	// global map outerSoftClip paired alignment
	// input
	string input_incompletePairedAlignmentPath_globalMapOuterSoftClipUniquePairedSAM;
	// output
	string output_folder_globalMapOuterSoftClipUniquePairedSAM;
	string output_oriSAM_globalMapOuterSoftClipUniquePairedSAM;
	string output_nonFusionSAM_globalMapOuterSoftClipUniquePairedSAM;
	string output_fusionSAM_globalMapOuterSoftClipUniquePairedSAM;
	string output_breakPoint_globalMapOuterSoftClipUniquePairedSAM;
	string output_breakPoint_adjusted_globalMapOuterSoftClipUniquePairedSAM;

	// generate fusion junc from adjusted fusion break point results file
	// input
	string input_breakPoint_adjusted_generateFusionJuncFromAdjustedBreakPointFile;
	// output
	string output_rawFusionJunc_generateFusionJuncFromAdjustedBreakPointFile;
	//string output_filteredFusionJunc_generateFusionJuncFromAdjustedBreakPointFile;

	// remap outerSoftClipp paired alignment
	// input
	string input_fusionJunc_remapOuterSoftClipUniquePairedSAM;
	string input_incompletePairedAlignmentPath_remapOuterSoftClipUniquePairedSAM;
	// output
	string output_folder_remapOuterSoftClipUniquePairedSAM;
	string output_oriSAM_remapOuterSoftClipUniquePairedSAM;
	string output_nonFusionSAM_remapOuterSoftClipUniquePairedSAM;
	string output_fusionSAM_remapOuterSoftClipUniquePairedSAM;
	string output_breakPoint_remapOuterSoftClipUniquePairedSAM;
	string output_updatedFusionJunc_remapOuterSoftClipUniquePairedSAM;

	// count encompassing fusion reads
	// input
	string input_fusionJunc_countEncompassingReads;
	string input_unpairedSAM_countEncompassingReads;
	// output
	string output_folder_countEncompassingReads;
	string output_oriSAM_countEncompassingReads;
	string output_nonFusionSAM_countEncompassingReads;
	string output_fusionSAM_countEncompassingReads;
	string output_fusionReadsEncompassingBreakPoint_countEncompassingReads;
	string output_updatedFusionJunc_countEncompassingReads;


	int normalRecordNum_globalMap;
	int normalRecordNum_remap;
	int normalRecordNum_countEncompassingReads;
public:
	MapSplice_Fusion_Info()
	{
		normalRecordNum_globalMap = 200000;
		normalRecordNum_remap = 200000;
		normalRecordNum_countEncompassingReads = 200000;
	}

	void initiateInputFilePathFromInitialMapResults(
		string& initialMapResultsFolderPath)
	{
		input_folder_path = initialMapResultsFolderPath;
	}

	void initiateInputOutputFilePath_globalMapPairedAlignment()
	{
		// input
		input_incompletePairedAlignmentPath_globalMapOuterSoftClipUniquePairedSAM
			= input_folder_path 
				+ "/phase2_output/fixHeadTail_incomplete_pair.sam"; 		
		// output
		output_folder_globalMapOuterSoftClipUniquePairedSAM 
			= output_folder_path + "/globalMap";
		string output_folder_path_globalMap 
			= output_folder_globalMapOuterSoftClipUniquePairedSAM;
		cout << "creating results folder for global map...." << endl;
		string mkdir_globalMap = "mkdir -p " + output_folder_path_globalMap;
		system(mkdir_globalMap.c_str());				
		output_oriSAM_globalMapOuterSoftClipUniquePairedSAM
			= output_folder_path_globalMap + "/kept.sam";
		output_nonFusionSAM_globalMapOuterSoftClipUniquePairedSAM
			= output_folder_path_globalMap + "/nonFusion.sam";
		output_fusionSAM_globalMapOuterSoftClipUniquePairedSAM
			= output_folder_path_globalMap + "/fusion.sam";
		output_breakPoint_globalMapOuterSoftClipUniquePairedSAM
			= output_folder_path_globalMap + "/breakPoint_raw.txt";
		output_breakPoint_adjusted_globalMapOuterSoftClipUniquePairedSAM
			= output_folder_path_globalMap + "/breakPoint_adjusted.txt";
	}

	bool adjustFusionBreakPoint_searchForCanonicalOnesWithOffset(
		bool breakPointInEnd1OrEnd2_bool,
		bool leftGene_NorOrRcm_bool, bool rightGene_NorOrRcm_bool,// here, left & right indicates the location in read seq.
		string& leftChrNameStr, string& rightChrNameStr,
		int leftBreakPointPos, int rightBreakPointPos,
		int leftAnchorLength, int rightAnchorLength,
		string& canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1,
		string& canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2,
		int& canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1,
		int& canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2, 
		string& canonicalAndStrandIdenfiedStrand_1,
		string& canonicalAndStrandIdenfiedStrand_2,
		int& canonicalAndStrandIdenfiedAnchorLength_1, 
		int& canonicalAndStrandIdenfiedAnchorLength_2,	
		Index_Info* indexInfo, int offset_max)
	{
		int tmp_offset_max = offset_max;
		int shorterAnchorLength = leftAnchorLength;
		if(rightAnchorLength < leftAnchorLength)
			shorterAnchorLength = rightAnchorLength;
		if(tmp_offset_max > shorterAnchorLength - 1)
			tmp_offset_max = shorterAnchorLength - 1;

		//return false;
		int leftChrNameInt = indexInfo->convertStringToInt(leftChrNameStr);
		int rightChrNameInt = indexInfo->convertStringToInt(rightChrNameStr);
		// string rawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
		// 	leftChrNameInt, rightChrNameInt, leftBreakPointPos, rightBreakPointPos,
		// 	leftGene_NorOrRcm_bool, rightGene_NorOrRcm_bool, breakPointInEnd1OrEnd2_bool);
		vector<int> candiOffsetVec;

		candiOffsetVec.push_back(0);
		for(int tmp = 0; tmp < tmp_offset_max; tmp++)
		{	
			candiOffsetVec.push_back(tmp+1);
			candiOffsetVec.push_back(-1-tmp);
		}
		if(leftGene_NorOrRcm_bool && rightGene_NorOrRcm_bool)
		{
			//if((rawFusionJuncFlankString == "GTAG")||(rawFusionJuncFlankString == "CTAC"))
			//	return true;//false;
			for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetVec.size(); tmpOffsetIndex ++)
			{
				int tmpOffset = candiOffsetVec[tmpOffsetIndex];
				string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt,
					leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
					true, true, breakPointInEnd1OrEnd2_bool);
				if(tmpRawFusionJuncFlankString == "GTAG")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
					canonicalAndStrandIdenfiedStrand_1 = "+";
					canonicalAndStrandIdenfiedStrand_2 = "+";
					canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset;
					return true;
				}
				else if(tmpRawFusionJuncFlankString == "CTAC")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
					canonicalAndStrandIdenfiedStrand_1 = "-";
					canonicalAndStrandIdenfiedStrand_2 = "-";
					canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset;
					return true;
				}
				else
				{}
			}
			return false;	
		}
		else if((!leftGene_NorOrRcm_bool) && (!rightGene_NorOrRcm_bool))
		{
			//if((rawFusionJuncFlankString == "AGGT")||(rawFusionJuncFlankString == "ACCT"))
			//	return false;
			for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetVec.size(); tmpOffsetIndex ++)
			{
				int tmpOffset = candiOffsetVec[tmpOffsetIndex];
				string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt,
					leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
					false, false, breakPointInEnd1OrEnd2_bool);
				if(tmpRawFusionJuncFlankString == "GTAG")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
					canonicalAndStrandIdenfiedStrand_1 = "+";
					canonicalAndStrandIdenfiedStrand_2 = "+";
					canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset;
					return true;
				}
				else if(tmpRawFusionJuncFlankString == "CTAC")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
					canonicalAndStrandIdenfiedStrand_2 = "-";
					canonicalAndStrandIdenfiedStrand_1 = "-";
					canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset;				
					return true;
				}
				else
				{}
			}
			return false;	
		}
		else if((!leftGene_NorOrRcm_bool) && rightGene_NorOrRcm_bool) // Rev For
		{
			if(!breakPointInEnd1OrEnd2_bool)
			{
				for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetVec.size(); tmpOffsetIndex ++)
				{
					int tmpOffset = candiOffsetVec[tmpOffsetIndex];
					string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
						leftChrNameInt, rightChrNameInt,
						leftBreakPointPos + tmpOffset, rightBreakPointPos - tmpOffset,
						false, true, (breakPointInEnd1OrEnd2_bool));
					if(tmpRawFusionJuncFlankString == "GTCT")
					{
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos - tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
						canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset;
						canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset;
						canonicalAndStrandIdenfiedStrand_2 = "-";
						canonicalAndStrandIdenfiedStrand_1 = "+";
						return true;
					}
					else if(tmpRawFusionJuncFlankString == "CTGT")
					{
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos - tmpOffset;					
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
						canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset;
						canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset;
						canonicalAndStrandIdenfiedStrand_2 = "-";
						canonicalAndStrandIdenfiedStrand_1 = "+";
						return true;				
					}
					else
					{}
				}		
			}
			else
			{
				for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetVec.size(); tmpOffsetIndex ++)
				{
					int tmpOffset = candiOffsetVec[tmpOffsetIndex];
					string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
						leftChrNameInt, rightChrNameInt,
						leftBreakPointPos - tmpOffset, rightBreakPointPos + tmpOffset,
						false, true, (breakPointInEnd1OrEnd2_bool));
					if(tmpRawFusionJuncFlankString == "AGAC")
					{
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos - tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
						canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset;
						canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset;
						canonicalAndStrandIdenfiedStrand_2 = "+";
						canonicalAndStrandIdenfiedStrand_1 = "-";
						return true;
					}
					else if(tmpRawFusionJuncFlankString == "ACAG")
					{
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos - tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
						canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
						canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset;
						canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset;
						canonicalAndStrandIdenfiedStrand_2 = "+";
						canonicalAndStrandIdenfiedStrand_1 = "-";
						return true;			
					}
					else
					{}
				}
			}
			return false;	
		}
		else // REV FOR
		{
			//if((rawFusionJuncFlankString == "ACAG")||(rawFusionJuncFlankString == "AGAC"))
			//	return false;
			/*for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetVec.size(); tmpOffsetIndex ++)
			{
				int tmpOffset = candiOffsetVec[tmpOffsetIndex];
				string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt,
					leftBreakPointPos + tmpOffset, rightBreakPointPos - tmpOffset,
					false, true);
				if((tmpRawFusionJuncFlankString == "ACAG"))
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos - tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength - tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedStrand_1 = "-";
					canonicalAndStrandIdenfiedStrand_2 = "+";	
					return true;
				}
				else if((tmpRawFusionJuncFlankString == "AGAC"))
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos - tmpOffset;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength - tmpOffset;
					canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength + tmpOffset;
					canonicalAndStrandIdenfiedStrand_1 = "-";
					canonicalAndStrandIdenfiedStrand_2 = "+";	
					return true;				
				}
			}*/
			cout << "error in adjustFusionBreakPoint_searchForCanonicalOnesWithOffset from ";
			cout << "globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion.cpp" << endl;
			return false;	
		}
	}

	void globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion(
		Index_Info* indexInfo,
		unsigned int* sa, BYTE* lcpCompress,
		unsigned int* childTab, char* chrom,
		BYTE* verifyChild, int* preIndexMapLengthArray, 
		unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray,		
		//RepeatRegion_Info* repeatRegionInfo,
		int offsetForAdjustingFusionBreakPoint_max, int threads_num, 
		bool fasta_or_fastq_bool, int min_fusion_distance, ofstream& log_ofs)
	{
		time_t nowtime;
		struct tm *local;

	   	ofstream oriSAM_ofs(
	   		output_oriSAM_globalMapOuterSoftClipUniquePairedSAM.c_str());
	   	ofstream nonFusionSAM_ofs(
	   		output_nonFusionSAM_globalMapOuterSoftClipUniquePairedSAM.c_str());
	   	ofstream fusionSAM_ofs(
	   		output_fusionSAM_globalMapOuterSoftClipUniquePairedSAM.c_str());
	   	ofstream breakPoint_ofs(
	   		output_breakPoint_globalMapOuterSoftClipUniquePairedSAM.c_str());
	   	ofstream breakPoint_adjusted_ofs(
	   		output_breakPoint_adjusted_globalMapOuterSoftClipUniquePairedSAM.c_str());

		//int offsetForAdjustingFusionBreakPoint_max = 5;
		//int min_fusion_distance = atoi(min_fusion_distance_str.c_str());
		int normalRecordNum_1stMapping = normalRecordNum_globalMap;
		int readTotalNum = 0;

		bool Do_cirRNA = true;
		Do_cirRNA = false;
		bool annotation_provided_bool = false;
		bool Do_annotation_only_bool = false;
		bool Do_extendHeadTail_phase1 = true;
		bool checkQualSeqForReadSegSeq = false;
		//////////////   start to load annotation  ///////////////////////
		Annotation_Info* annotationInfo = new Annotation_Info();
		//////////////   end of loading annotation ///////////////////////		

		ifstream incompleteUniquePairedAlignment_ifs(
			input_incompletePairedAlignmentPath_globalMapOuterSoftClipUniquePairedSAM.c_str());
		int normalRecordNum = normalRecordNum_globalMap;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;

		vector<string> inputSamStr1Vec(normalRecordNum);
		vector<string> inputSamStr2Vec(normalRecordNum);
		vector<string> outputOriSamStrVec(normalRecordNum);
		vector<string> outputNonfusionSamStrVec(normalRecordNum);
		vector<string> outputFusionSamStrVec(normalRecordNum);
		vector<string> outputFusionSamStrVec_withBreakPoint(normalRecordNum);
		vector<string> outputFusionSamStrVec_withBreakPoint_adjusted(normalRecordNum);
		vector< RepeatRegion_Info* > repeatRegionInfoVec;
		for(int tmp = 0; tmp < threads_num; tmp++)
		{
			RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
			repeatRegionInfoVec.push_back(repeatRegionInfo);
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 

		string line1, line2;
		for(tmpTurn = 0; 
			//tmpTurn <= 300     //used to control # of rounds to process
			; tmpTurn++)
		{		
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
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
				int threadNO = omp_get_thread_num();
				bool SE_or_PE_bool = true;
				//cout << "threadNO: " << threadNO << endl;
				IncompleteUniquePairedAlignment2detectFusion_Info incompleteUniquePairedAlignmentInfo;
				//cout << "start to initiate incompleteUniquePairedAlignmentInfo" << endl;
				bool incompleteUniquePaired_bool
					= incompleteUniquePairedAlignmentInfo.initiateWith2samStr_globalMapOuterUnfixedEnd2detectFusionBreakPoint(
						inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
				//cout << "incompleteUniquePaired_bool: " << incompleteUniquePaired_bool << endl;
				// if(incompleteUniquePaired_bool)
				// {	
				// 	cout << "readName_1: " << incompleteUniquePairedAlignmentInfo.returnReadName_1() << endl;
				// 	cout << "readName_2: " << incompleteUniquePairedAlignmentInfo.returnReadName_2() << endl;
				// }
				bool leftReadHeadUnfixed_bool = false;
				PE_Read_Info readInfo_unfixedLeftReadHead;
				FixPhase1Info fixPhase1Info_unfixedLeftReadHead;
				PE_Read_Alignment_Info peAlignInfo_unfixedLeftReadHead;
				bool leftReadHeadUnfixed_fixed_bool = false;
				bool leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool = false;
				string leftReadHeadUnfixed_fixed_SAMstr;

				bool rightReadTailUnfixed_bool = false;
				PE_Read_Info readInfo_unfixedRightReadTail;
				FixPhase1Info fixPhase1Info_unfixedRightReadTail;
				PE_Read_Alignment_Info peAlignInfo_unfixedRightReadTail;
				bool rightReadTailUnfixed_fixed_bool = false;
				bool rightReadTailUnfixed_fixed_Nor_or_Rcm_bool = false;
				string rightReadTailUnfixed_fixed_SAMstr;

				if(incompleteUniquePaired_bool)
				{
					leftReadHeadUnfixed_bool 
						= incompleteUniquePairedAlignmentInfo.leftReadHeadUnfixedOrNot();
					rightReadTailUnfixed_bool
						= incompleteUniquePairedAlignmentInfo.rightReadTailUnfixedOrNot();
					//cout << "leftReadHeadUnfixed_bool: " << leftReadHeadUnfixed_bool << endl;
					//cout << "rightReadTailUnfixed_bool: " << rightReadTailUnfixed_bool << endl;

					if(leftReadHeadUnfixed_bool)
					{
						//cout << "start to do leftReadHeadUnfixed_bool" << endl;
						string tmpReadName_1 
							= incompleteUniquePairedAlignmentInfo.returnReadName_1();
						//cout << "tmpReadName_1: " << tmpReadName_1 << endl;
						string tmpUnfixedLeftReadHead_readSeq
							= incompleteUniquePairedAlignmentInfo.returnUnfixedLeftReadHead_readSeq();
						//cout << "tmpUnfixedLeftReadHead_readSeq: " << tmpUnfixedLeftReadHead_readSeq << endl;
						string tmpUnfixedLeftReadHead_qualSeq
							= incompleteUniquePairedAlignmentInfo.returnUnfixedLeftReadHead_qualSeq(fasta_or_fastq_bool);
						//cout << "tmpUnfixedLeftReadHead_qualSeq: " << tmpUnfixedLeftReadHead_qualSeq << endl;
						string tmpStr;
						//cout << "start to initiateReadInfo " << endl;
						readInfo_unfixedLeftReadHead.initiateReadInfo(tmpReadName_1, tmpStr, tmpUnfixedLeftReadHead_readSeq, tmpStr, 
							tmpUnfixedLeftReadHead_qualSeq, tmpStr, fasta_or_fastq_bool, true);
						//cout << "start to do fixPhase1_segInfo" << endl;
						fixPhase1Info_unfixedLeftReadHead.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
							preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray, 
							readInfo_unfixedLeftReadHead, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
						//cout << "start to do fixPhase1_pathInfo" << endl;
						fixPhase1Info_unfixedLeftReadHead.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedLeftReadHead, 
							annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
						//cout << "start to do fixPhase1_gapInfo" << endl;
						fixPhase1Info_unfixedLeftReadHead.fixPhase1_gapInfo(readInfo_unfixedLeftReadHead, indexInfo, Do_cirRNA,
							Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
						//cout << "start to initiatePeAlignInfo" << endl;
						peAlignInfo_unfixedLeftReadHead.initiatePeAlignInfo(fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor1,
							fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm1, fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor2,
							fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm2, indexInfo, true);
						//cout << "start to chooseBestAlignment_final_SE" << endl;
						peAlignInfo_unfixedLeftReadHead.chooseBestAlignment_final_SE();
						leftReadHeadUnfixed_fixed_bool 
							= peAlignInfo_unfixedLeftReadHead.mappedForFusionDetection_unfixedHead_uniqueCompleteMapped_SE_bool(
								min_fusion_distance, incompleteUniquePairedAlignmentInfo.returnChrNameInt(), 
								incompleteUniquePairedAlignmentInfo.returnStartPos_1(), indexInfo);
						if(leftReadHeadUnfixed_fixed_bool)
						{	
							leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool
								= peAlignInfo_unfixedLeftReadHead.returnOnlySeAlign_Nor_or_Rcm(); 
							leftReadHeadUnfixed_fixed_SAMstr 
								= peAlignInfo_unfixedLeftReadHead.returnSAMstr_mappedForFusion_headSoftClipUniquePaired_SE(
									readInfo_unfixedLeftReadHead, fasta_or_fastq_bool);
							//cout << "leftReadHeadUnfixed_fixed_SAMstr: " << leftReadHeadUnfixed_fixed_SAMstr << endl;
						}
					}

					if(rightReadTailUnfixed_bool)
					{
						//cout << "start to do rightReadTailUnfixed_bool" << endl;
						string tmpReadName_2
							= incompleteUniquePairedAlignmentInfo.returnReadName_2();
						//cout << "tmpReadName_2: " << tmpReadName_2 << endl; 
						string tmpUnfixedRightReadTail_readSeq
							= incompleteUniquePairedAlignmentInfo.returnUnfixedRightReadTail_readSeq();
						//cout << "tmpUnfixedRightReadTail_readSeq: " << tmpUnfixedRightReadTail_readSeq << endl;
						string tmpUnfixedRightReadTail_qualSeq
							= incompleteUniquePairedAlignmentInfo.returnUnfixedRightReadTail_qualSeq(fasta_or_fastq_bool);
						//cout << "tmpUnfixedRightReadTail_qualSeq: " << tmpUnfixedRightReadTail_qualSeq << endl;
						string tmpStr;
						//cout << "start to initiateReadInfo " << endl;
						readInfo_unfixedRightReadTail.initiateReadInfo(tmpReadName_2, tmpStr, tmpUnfixedRightReadTail_readSeq, tmpStr, 
							tmpUnfixedRightReadTail_qualSeq, tmpStr, fasta_or_fastq_bool, true);
						//cout << "start to do fixPhase1_segInfo" << endl;
						fixPhase1Info_unfixedRightReadTail.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
							preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray, 
							readInfo_unfixedRightReadTail, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
						//cout << "start to do fixPhase1_pathInfo" << endl;
						fixPhase1Info_unfixedRightReadTail.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedRightReadTail, 
							annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
						//cout << "start to do fixPhase1_gapInfo" << endl;
						fixPhase1Info_unfixedRightReadTail.fixPhase1_gapInfo(readInfo_unfixedRightReadTail, indexInfo, Do_cirRNA,
							Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
						//cout << "start to initiatePeAlignInfo" << endl;
						peAlignInfo_unfixedRightReadTail.initiatePeAlignInfo(fixPhase1Info_unfixedRightReadTail.pathInfo_Nor1,
							fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm1, fixPhase1Info_unfixedRightReadTail.pathInfo_Nor2,
							fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm2, indexInfo, true);
						//cout << "start to chooseBestAlignment_final_SE" << endl;
						peAlignInfo_unfixedRightReadTail.chooseBestAlignment_final_SE();
						rightReadTailUnfixed_fixed_bool 
							= peAlignInfo_unfixedRightReadTail.mappedForFusionDetection_unfixedTail_uniqueCompleteMapped_SE_bool(
								min_fusion_distance, incompleteUniquePairedAlignmentInfo.returnChrNameInt(), 
								incompleteUniquePairedAlignmentInfo.returnEndPos_2(), indexInfo);
						if(rightReadTailUnfixed_fixed_bool)
						{
							rightReadTailUnfixed_fixed_Nor_or_Rcm_bool
								= peAlignInfo_unfixedRightReadTail.returnOnlySeAlign_Nor_or_Rcm();	
							rightReadTailUnfixed_fixed_SAMstr 
								= peAlignInfo_unfixedRightReadTail.returnSAMstr_mappedForFusion_tailSoftClipUniquePaired_SE(
									readInfo_unfixedRightReadTail, fasta_or_fastq_bool);
							//cout << "rightReadTailUnfixed_fixed_SAMstr: " << rightReadTailUnfixed_fixed_SAMstr << endl;
						}
					}			
				}

				outputOriSamStrVec[tmpOpenMP] = "";			
				outputNonfusionSamStrVec[tmpOpenMP] = "";
				outputFusionSamStrVec[tmpOpenMP] = "";
				outputFusionSamStrVec_withBreakPoint[tmpOpenMP] = "";
				outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP] = "";
				//cout << "start to generte all results str ...." << endl;
				if(incompleteUniquePaired_bool)
				{
					if(leftReadHeadUnfixed_fixed_bool || rightReadTailUnfixed_fixed_bool)
					{
						string tmpSAMoutputStr = inputSamStr1Vec[tmpOpenMP];
						if(leftReadHeadUnfixed_fixed_bool)
						{
							if(leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool)
								tmpSAMoutputStr = tmpSAMoutputStr + "\tZF:Z:FUS_leftHead_for\n" 
									+ leftReadHeadUnfixed_fixed_SAMstr + "\tZF:Z:FUS_leftHead_for";
							else
								tmpSAMoutputStr = tmpSAMoutputStr + "\tZF:Z:FUS_leftHead_rev\n" 
									+ leftReadHeadUnfixed_fixed_SAMstr + "\tZF:Z:FUS_leftHead_rev";							
						}
						tmpSAMoutputStr += "\n";
						tmpSAMoutputStr += inputSamStr2Vec[tmpOpenMP];
						if(rightReadTailUnfixed_fixed_bool)
						{
							if(rightReadTailUnfixed_fixed_Nor_or_Rcm_bool)
								tmpSAMoutputStr = tmpSAMoutputStr + "\tZF:Z:FUS_rightTail_rev\n" 
									+ rightReadTailUnfixed_fixed_SAMstr + "\tZF:Z:FUS_rightTail_rev";
							else
								tmpSAMoutputStr = tmpSAMoutputStr + "\tZF:Z:FUS_rightTail_for\n" 
									+ rightReadTailUnfixed_fixed_SAMstr + "\tZF:Z:FUS_rightTail_for";							
						}
						outputFusionSamStrVec[tmpOpenMP] = tmpSAMoutputStr;

						string tmpLeftReadHeadBreakPoint;
						string tmpRightReadTailBreakPoint;
						string tmpLeftReadHeadBreakPoint_adjusted;
						string tmpRightReadTailBreakPoint_adjusted;	

						//cout << "start to generate breakPointStr" << endl;				
						if(leftReadHeadUnfixed_fixed_bool)
						{
							tmpLeftReadHeadBreakPoint 
								+= incompleteUniquePairedAlignmentInfo.returnReadName_1();
							tmpLeftReadHeadBreakPoint += "\t";
							string tmpLeftReadHead_chrName 
								= peAlignInfo_unfixedLeftReadHead.returnOnlySeAlign_chrNameStr();
							int tmpLeftReadHead_startMapPos
								= peAlignInfo_unfixedLeftReadHead.returnOnlySeAlign_startPos();
							int tmpLeftReadHead_endMapPos
								= peAlignInfo_unfixedLeftReadHead.returnOnlySeAlign_endPos();
							string tmpOriSAM_chrName
								= incompleteUniquePairedAlignmentInfo.returnChrName(indexInfo);
							int tmpOriSAM_startPos_1 
								= incompleteUniquePairedAlignmentInfo.returnStartPos_1();
							int tmpFusionAnchorLength_left 
							// this should the length of the last jumpcode in fixedLeftReadHead, 
							// for now = unfixedHeadLength_1() as only jumpcode M is accepted for this unfixedReadHead 
								= incompleteUniquePairedAlignmentInfo.returnUnfixedHeadLength_1();
							int tmpFusionAnchorLength_right 
								= incompleteUniquePairedAlignmentInfo.returnJumpCodeLength_1(1);
							if(leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool)
							{
								int canonicalAndStrandIdenfiedBreakPoint_left;
								string canonicalAndStrandIdenfiedChrName_left;
								int canonicalAndStrandIdenfiedBreakPoint_right;
								string canonicalAndStrandIdenfiedChrName_right;
								string canonicalAndStrandIdenfiedStrand_left;
								string canonicalAndStrandIdenfiedStrand_right;
								int canonicalAndStrandIdenfiedAnchorLength_left;
								int canonicalAndStrandIdenfiedAnchorLength_right;
								bool canonicalAndStrandIdenfiedOrNot 
									= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset(
										true, true, true, tmpLeftReadHead_chrName, tmpOriSAM_chrName,
										tmpLeftReadHead_endMapPos, tmpOriSAM_startPos_1,
										tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
										canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
										canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
										indexInfo, offsetForAdjustingFusionBreakPoint_max);
								if(canonicalAndStrandIdenfiedOrNot)
								{
									tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
										+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
										+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right);
								}
								else
								{
									tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
										+ tmpLeftReadHead_chrName + "\t" + tmpOriSAM_chrName + "\t"
										+ int_to_str(tmpLeftReadHead_endMapPos) + "\t"
										+ int_to_str(tmpOriSAM_startPos_1) + "\tN\tN\t"
										+ indexInfo->returnRawFusionJuncFlankString(
											tmpLeftReadHead_chrName, tmpOriSAM_chrName,
											tmpLeftReadHead_endMapPos, tmpOriSAM_startPos_1,
											true, true, indexInfo, true) 
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  
								}
								tmpLeftReadHeadBreakPoint = tmpLeftReadHeadBreakPoint
									+ tmpLeftReadHead_chrName + "\t" + tmpOriSAM_chrName + "\t"
									+ int_to_str(tmpLeftReadHead_endMapPos) + "\t"
									+ int_to_str(tmpOriSAM_startPos_1) + "\tFOR_B+FOR_A~REV_A\t"
									+ indexInfo->returnRawFusionJuncFlankString(
										tmpLeftReadHead_chrName, tmpOriSAM_chrName,
										tmpLeftReadHead_endMapPos, tmpOriSAM_startPos_1,
										true, true, indexInfo, true)
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  
							}
							else
							{
								int canonicalAndStrandIdenfiedBreakPoint_left;
								string canonicalAndStrandIdenfiedChrName_left;
								int canonicalAndStrandIdenfiedBreakPoint_right;
								string canonicalAndStrandIdenfiedChrName_right;
								string canonicalAndStrandIdenfiedStrand_left;
								string canonicalAndStrandIdenfiedStrand_right;
								int canonicalAndStrandIdenfiedAnchorLength_left;
								int canonicalAndStrandIdenfiedAnchorLength_right;
								bool canonicalAndStrandIdenfiedOrNot 
									= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset(
										true, false, true, tmpLeftReadHead_chrName, tmpOriSAM_chrName,
										tmpLeftReadHead_startMapPos, tmpOriSAM_startPos_1,
										tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
										canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
										canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,									
										indexInfo, offsetForAdjustingFusionBreakPoint_max);
								//cout << "canonicalAndStrandIdenfiedOrNot: " << canonicalAndStrandIdenfiedOrNot << endl;
								if(canonicalAndStrandIdenfiedOrNot)
								{
									tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
										+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
										+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right);
								}	
								else
								{
									tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
										+ tmpLeftReadHead_chrName + "\t" + tmpOriSAM_chrName + "\t"
										+ int_to_str(tmpLeftReadHead_startMapPos) + "\t"
										+ int_to_str(tmpOriSAM_startPos_1) + "\tN\tN\t"
										+ indexInfo->returnRawFusionJuncFlankString(
											tmpLeftReadHead_chrName, tmpOriSAM_chrName,
											tmpLeftReadHead_startMapPos, tmpOriSAM_startPos_1,
											false, true, indexInfo, true) 
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  
								}
								tmpLeftReadHeadBreakPoint = tmpLeftReadHeadBreakPoint
									+ tmpLeftReadHead_chrName + "\t" + tmpOriSAM_chrName + "\t"
									+ int_to_str(tmpLeftReadHead_startMapPos) + "\t"
									+ int_to_str(tmpOriSAM_startPos_1) + "\tREV_B+FOR_A~REV_A\t"
									+ indexInfo->returnRawFusionJuncFlankString(
										tmpLeftReadHead_chrName, tmpOriSAM_chrName,
										tmpLeftReadHead_startMapPos, tmpOriSAM_startPos_1,
										false, true, indexInfo, true) 
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  
							}
						}
						if(rightReadTailUnfixed_fixed_bool)
						{
							tmpRightReadTailBreakPoint 
								+= incompleteUniquePairedAlignmentInfo.returnReadName_2();
							tmpRightReadTailBreakPoint += "\t";
							string tmpOriSAM_chrName
								= incompleteUniquePairedAlignmentInfo.returnChrName(indexInfo);
							int tmpOriSAM_endPos_2 
								= incompleteUniquePairedAlignmentInfo.returnEndPos_2();
							string tmpRightReadTail_chrName
								= peAlignInfo_unfixedRightReadTail.returnOnlySeAlign_chrNameStr();
							int tmpRightReadTail_startPos
								= peAlignInfo_unfixedRightReadTail.returnOnlySeAlign_startPos();
							int tmpRightReadTail_endPos
								= peAlignInfo_unfixedRightReadTail.returnOnlySeAlign_endPos();
							int tmpFusionAnchorLength_left
								= incompleteUniquePairedAlignmentInfo.returnJumpCodeLength_rev_2(1);
							int tmpFusionAnchorLength_right
								= incompleteUniquePairedAlignmentInfo.returnUnfixedTailLength_2();
							if(rightReadTailUnfixed_fixed_Nor_or_Rcm_bool)
							{
								int canonicalAndStrandIdenfiedBreakPoint_left;
								string canonicalAndStrandIdenfiedChrName_left;
								int canonicalAndStrandIdenfiedBreakPoint_right;
								string canonicalAndStrandIdenfiedChrName_right;
								string canonicalAndStrandIdenfiedStrand_left;
								string canonicalAndStrandIdenfiedStrand_right;
								int canonicalAndStrandIdenfiedAnchorLength_left;
								int canonicalAndStrandIdenfiedAnchorLength_right;
								bool canonicalAndStrandIdenfiedOrNot 
									= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset(
										false, false, false, tmpOriSAM_chrName, tmpRightReadTail_chrName,
										tmpOriSAM_endPos_2, tmpRightReadTail_startPos,
										tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
										canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
										canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
										indexInfo, offsetForAdjustingFusionBreakPoint_max);						
								if(canonicalAndStrandIdenfiedOrNot)
								{
									tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
										+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
										+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right);
								}
								else
								{
									tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
										+ tmpOriSAM_chrName + "\t" + tmpRightReadTail_chrName + "\t" 
										+ int_to_str(tmpOriSAM_endPos_2) + "\t"
										+ int_to_str(tmpRightReadTail_startPos) + "\tN\tN\t" 
										+ indexInfo->returnRawFusionJuncFlankString(
											tmpRightReadTail_chrName, tmpOriSAM_chrName,
											tmpRightReadTail_startPos, tmpOriSAM_endPos_2,
											false, false, indexInfo, false)
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  								
								}
								tmpRightReadTailBreakPoint = tmpRightReadTailBreakPoint
									+ tmpOriSAM_chrName + "\t" + tmpRightReadTail_chrName + "\t" 
									+ int_to_str(tmpOriSAM_endPos_2) + "\t"
									+ int_to_str(tmpRightReadTail_startPos) + "\tFOR_A~REV_B+REV_A\t" 
									+ indexInfo->returnRawFusionJuncFlankString(
										tmpOriSAM_chrName, tmpRightReadTail_chrName, 
										tmpOriSAM_endPos_2, tmpRightReadTail_startPos,
										false, false, indexInfo, false)
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  
							}
							else
							{
								int canonicalAndStrandIdenfiedBreakPoint_left;
								string canonicalAndStrandIdenfiedChrName_left;
								int canonicalAndStrandIdenfiedBreakPoint_right;
								string canonicalAndStrandIdenfiedChrName_right;
								string canonicalAndStrandIdenfiedStrand_left;
								string canonicalAndStrandIdenfiedStrand_right;
								int canonicalAndStrandIdenfiedAnchorLength_left;
								int canonicalAndStrandIdenfiedAnchorLength_right;
								bool canonicalAndStrandIdenfiedOrNot 
									= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset(
										false, false, true, tmpOriSAM_chrName, tmpRightReadTail_chrName,
										tmpOriSAM_endPos_2, tmpRightReadTail_endPos,
										tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
										canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
										canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
										indexInfo, offsetForAdjustingFusionBreakPoint_max);
								if(canonicalAndStrandIdenfiedOrNot)
								{
									tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
										+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
										+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
										+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right);
								}
								else
								{
									tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
										+ tmpOriSAM_chrName + "\t" + tmpRightReadTail_chrName + "\t" 
										+ int_to_str(tmpOriSAM_endPos_2) + "\t"
										+ int_to_str(tmpRightReadTail_endPos) + "\tN\tN\t" 
										+ indexInfo->returnRawFusionJuncFlankString(
											tmpOriSAM_chrName, tmpRightReadTail_chrName,
											tmpOriSAM_endPos_2, tmpRightReadTail_endPos,
											false, true, indexInfo, false) 
										+ "\t" + int_to_str(tmpFusionAnchorLength_left)
										+ "\t" + int_to_str(tmpFusionAnchorLength_right);  	
								}
								tmpRightReadTailBreakPoint = tmpRightReadTailBreakPoint
										+ tmpOriSAM_chrName + "\t" + tmpRightReadTail_chrName + "\t" 
										+ int_to_str(tmpOriSAM_endPos_2) + "\t"
										+ int_to_str(tmpRightReadTail_endPos) + "\tFOR_A~FOR_B+REV_A\t"
									+ indexInfo->returnRawFusionJuncFlankString(
										tmpOriSAM_chrName, tmpRightReadTail_chrName, 
										tmpOriSAM_endPos_2, tmpRightReadTail_endPos, 
										false, true, indexInfo, false) 
									+ "\t" + int_to_str(tmpFusionAnchorLength_left)
									+ "\t" + int_to_str(tmpFusionAnchorLength_right);							
							}
						}

						//cout << "start to insert breakPointStr " << endl;
				 		if((leftReadHeadUnfixed_fixed_bool)&&(rightReadTailUnfixed_fixed_bool))
						{
							outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
								= tmpLeftReadHeadBreakPoint_adjusted + "\n" + tmpRightReadTailBreakPoint_adjusted;
							outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
								= tmpLeftReadHeadBreakPoint + "\n" + tmpRightReadTailBreakPoint;						
						}
						else if(leftReadHeadUnfixed_fixed_bool)
						{
							outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
								= tmpLeftReadHeadBreakPoint_adjusted;
							outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
								= tmpLeftReadHeadBreakPoint;
						}
						else if(rightReadTailUnfixed_fixed_bool)
						{
							outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
								= tmpRightReadTailBreakPoint_adjusted;
							outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
								= tmpRightReadTailBreakPoint;
						}
						else
						{}
						//cout << "end of inserting str" << endl;									
					}
					else
					{
						outputNonfusionSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
						//cout << "end of inserting str " << endl;
					}
				}
				else
				{
					outputOriSamStrVec[tmpOpenMP] = inputSamStr1Vec[tmpOpenMP] + "\n" + inputSamStr2Vec[tmpOpenMP];
					//cout << "end of inserting str" << endl;
				}
				fixPhase1Info_unfixedLeftReadHead.memoryFree();
				fixPhase1Info_unfixedRightReadTail.memoryFree();
				peAlignInfo_unfixedLeftReadHead.memoryFree();
				peAlignInfo_unfixedRightReadTail.memoryFree();					
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
				if(outputFusionSamStrVec_withBreakPoint[tmp] != "")
					breakPoint_ofs << outputFusionSamStrVec_withBreakPoint[tmp] << endl;
				if(outputFusionSamStrVec_withBreakPoint_adjusted[tmp] != "")
					breakPoint_adjusted_ofs << outputFusionSamStrVec_withBreakPoint_adjusted[tmp] << endl;			
			}
			nowtime = time(NULL);
			local = localtime(&nowtime);
			cout << endl << "[" << asctime(local) 
				<< "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
			incompleteUniquePairedAlignment_ifs.close();					
	   		oriSAM_ofs.close();
	   		nonFusionSAM_ofs.close();
	   		fusionSAM_ofs.close();
	   		breakPoint_ofs.close();
	   		breakPoint_adjusted_ofs.close();
		}
		incompleteUniquePairedAlignment_ifs.close();
	   	oriSAM_ofs.close();
	   	nonFusionSAM_ofs.close();
	   	fusionSAM_ofs.close();
	   	breakPoint_ofs.close();
	   	breakPoint_adjusted_ofs.close();		
	}

	void initiateInputOutputFilePath_generateFusionjuncFromAdjustedFusionBreakPointFile()
	{
		// input
		input_breakPoint_adjusted_generateFusionJuncFromAdjustedBreakPointFile
			= output_breakPoint_adjusted_globalMapOuterSoftClipUniquePairedSAM;
		// output
		output_rawFusionJunc_generateFusionJuncFromAdjustedBreakPointFile
			= output_folder_globalMapOuterSoftClipUniquePairedSAM + "/fusionJunc_raw.txt";
		//output_filteredFusionJunc_generateFusionJuncFromAdjustedBreakPointFile
		//	= output_folder_globalMapOuterSoftClipUniquePairedSAM + "/fusionJunc_filtered.txt";
	}

	void generateFusionJuncFromAdjustedFusionBreakPointFile(int threads_num, 
		bool fasta_or_fastq_bool, Index_Info* indexInfo, ofstream& log_ofs)
	{
		time_t nowtime;
		struct tm *local;

		int toCheckAnchorLengthMax = 30;
		bool checkAlterSpliceBool = true;
		checkAlterSpliceBool = false;

		ofstream rawFusionJunc_ofs(
			output_rawFusionJunc_generateFusionJuncFromAdjustedBreakPointFile.c_str());
		// ofstream filteredFusionJunc_ofs(
		// 	output_filteredFusionJunc_generateFusionJuncFromAdjustedBreakPointFile.c_str());
		ifstream fusionBreakPointResults_ifs(
			input_breakPoint_adjusted_generateFusionJuncFromAdjustedBreakPointFile.c_str());
	
		int chromNum = indexInfo->returnChromNum();
		vector< vector< FusionJuncPosPair2fusionInfoVecIndexMap > > fusionJuncMapVecVec;
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++ )
		{
			vector<FusionJuncPosPair2fusionInfoVecIndexMap> tmpFusionJuncMapVec;
			for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2 ++)
			{	
				FusionJuncPosPair2fusionInfoVecIndexMap tmpFusionJuncPosPair2fusionInfoVecIndexMap;
				tmpFusionJuncMapVec.push_back(tmpFusionJuncPosPair2fusionInfoVecIndexMap);
			}
			fusionJuncMapVecVec.push_back(tmpFusionJuncMapVec);
		}

		vector<int> fusionJuncSupNumVec;
		vector<string> fusionJuncStrandVec_1;
		vector<string> fusionJuncStrandVec_2;
		vector<string> fusionJuncFlankStringVec;
		vector<int> fusionJuncAnchorLengthVec_1;
		vector<int> fusionJuncAnchorLengthVec_2;

		while(!fusionBreakPointResults_ifs.eof())
		{
			string tmpFusionStr;
			getline(fusionBreakPointResults_ifs, tmpFusionStr);
			//cout << "tmpFusionStr: " << tmpFusionStr << endl;
			if((fusionBreakPointResults_ifs.eof())||(tmpFusionStr == ""))
				break;		
			int detectedChrNameInt_1, detectedChrNameInt_2;
			int detectedBreakPoint_1, detectedBreakPoint_2;
			string detectedFusionJuncStrand_1, detectedFusionJuncStrand_2;
			string detectedFusionJuncFlankString;
			int detectedFusionAnchorLength_1, detectedFusionAnchorLength_2;	
			extractFusionInfoFromAdjustedFusionBreakPointStr(
				detectedChrNameInt_1, detectedChrNameInt_2,
				detectedBreakPoint_1, detectedBreakPoint_2, 
				detectedFusionJuncStrand_1, detectedFusionJuncStrand_2,
				detectedFusionJuncFlankString,
				detectedFusionAnchorLength_1, detectedFusionAnchorLength_2,
				tmpFusionStr, indexInfo);
			//cout << "detectedChrNameInt_1: " << detectedChrNameInt_1 << endl;
			//cout << "detectedChrNameInt_2: " << detectedChrNameInt_2 << endl;
			//cout << "detectedBreakPoint_1: " << detectedBreakPoint_1 << endl;
			//cout << "detectedBreakPoint_2: " << detectedBreakPoint_2 << endl;
			insertUpdataFusionJuncMapVecVecAndFusionInfoVecWithDetectedFusionJunc(
				fusionJuncMapVecVec, 
				detectedChrNameInt_1, detectedChrNameInt_2,
				detectedBreakPoint_1, detectedBreakPoint_2,
				detectedFusionJuncStrand_1, detectedFusionJuncStrand_2,
				detectedFusionJuncFlankString,
				detectedFusionAnchorLength_1, detectedFusionAnchorLength_2,
				fusionJuncSupNumVec,
				fusionJuncStrandVec_1, fusionJuncStrandVec_2,
				fusionJuncFlankStringVec, fusionJuncAnchorLengthVec_1,
				fusionJuncAnchorLengthVec_2, indexInfo);
		}				

		cout << "start to output Fusion junction info ..." << endl;
		//output Fusion junction info
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
		{
			string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChr);
			for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2++)
			{
				string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChr2);
				for(FusionJuncPosPair2fusionInfoVecIndexMap::iterator tmp1stMapIter 
					= (fusionJuncMapVecVec[tmpChr])[tmpChr2].begin();
					tmp1stMapIter != (fusionJuncMapVecVec[tmpChr])[tmpChr2].end();
					tmp1stMapIter ++)
				{
					int tmpJunc_startPos = tmp1stMapIter->first;
					for(FusionJuncEndPos2fusionInfoVecIndexMap::iterator tmp2ndMapIter
						= (tmp1stMapIter->second).begin();
						tmp2ndMapIter != (tmp1stMapIter->second).end();
						tmp2ndMapIter ++)
					{
						int tmpJunc_endPos = tmp2ndMapIter->first;
						int tmpJunc_index = tmp2ndMapIter->second;
						int tmpJunc_supNum = fusionJuncSupNumVec[tmpJunc_index];
						string tmpJunc_strand_1 = fusionJuncStrandVec_1[tmpJunc_index];
						string tmpJunc_strand_2 = fusionJuncStrandVec_2[tmpJunc_index];
						string tmpJunc_flankString = fusionJuncFlankStringVec[tmpJunc_index];
						int tmpJunc_anchorLength_1 = fusionJuncAnchorLengthVec_1[tmpJunc_index];
						int tmpJunc_anchorLength_2 = fusionJuncAnchorLengthVec_2[tmpJunc_index];
						/*int validAnchorLength_1, validAnchorLength_2, penalty_1, penalty_2;
						bool tmpCheckMatchThroughSeqSimilarityBool 
							= checkMatchTroughSeqSimilarity(tmpJunc_strand_1, tmpJunc_strand_2, 
								tmpChr, tmpChr2, tmpJunc_startPos, tmpJunc_endPos,
								tmpJunc_anchorLength_1, tmpJunc_anchorLength_2, toCheckAnchorLengthMax, 
								indexInfo, validAnchorLength_1, validAnchorLength_2, penalty_1, penalty_2);

						int validAnchorLength_alterSplice_1, validAnchorLength_alterSplice_2,
							penalty_alterSplice_1, penalty_alterSplice_2;
						bool tmpCheckAlterSpliceAnchorSeqSimilarityBool;*/
						/*if(checkAlterSpliceBool)
						{	
							tmpCheckAlterSpliceAnchorSeqSimilarityBool = checkAlterSpliceAnchorSeqSimilarity(
								tmpJunc_strand_1, tmpJunc_strand_2, tmpChr, tmpChr2, tmpJunc_startPos, tmpJunc_endPos,
								tmpJunc_anchorLength_1, tmpJunc_anchorLength_2, toCheckAlterSpliceAnchorLengthMax,
								indexInfo, validAnchorLength_alterSplice_1, validAnchorLength_alterSplice_2,
								penalty_alterSplice_1, penalty_alterSplice_2);
						}*/

						rawFusionJunc_ofs << tmpChrNameStr_1 << "\t" << tmpChrNameStr_2 << "\t"
							<< tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" 
							<< tmpJunc_strand_1 << "\t" << tmpJunc_strand_2 << "\t" 
							<< tmpJunc_flankString << "\t" 
							<< tmpJunc_anchorLength_1 << "\t" << tmpJunc_anchorLength_2 << "\t"
							<< tmpJunc_supNum << "\t" 
							<< int_to_str(tmpJunc_supNum_encompassing) << endl;
						// if(tmpCheckMatchThroughSeqSimilarityBool)
						// {	
						// 	rawFusionJunc_ofs << validAnchorLength_2 << ":" << penalty_1 << "\t" 
						// 		<< validAnchorLength_1 << ":" << penalty_2 << "\t";
						// 	if((penalty_1 <= (validAnchorLength_2 / 4))||(penalty_2 <= (validAnchorLength_1 / 4)))
						// 		rawFusionJunc_ofs << "filterOut_matchThroughPenalty" << endl;
						// 	else 
						// 	{
						// 		filteredFusionJunc_ofs << tmpChrNameStr_1 << "\t" << tmpChrNameStr_2 << "\t"
						// 			<< tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" 
						// 			<< tmpJunc_strand_1 << "\t" << tmpJunc_strand_2 << "\t" 
						// 			<< tmpJunc_flankString << "\t" 
						// 			<< tmpJunc_anchorLength_1 << "\t" << tmpJunc_anchorLength_2 << "\t"
						// 			<< tmpJunc_supNum << "\t0" << endl; 								
						// 		rawFusionJunc_ofs << "pass" << endl;
						// 	}
						// }
						// else
						// 	rawFusionJunc_ofs << "N\tN" << "\tfilterOut_notCanonical" << endl;
					}
				}
			}
		}
		fusionBreakPointResults_ifs.close();
		rawFusionJunc_ofs.close();
		//filteredFusionJunc_ofs.close();
	}

	void extractFusionInfoFromAdjustedFusionBreakPointStr(
		int& detectedChrNameInt_1, int& detectedChrNameInt_2,
		int& detectedBreakPoint_1, int& detectedBreakPoint_2,
		string& detectedFusionJuncStrand_1, string& detectedFusionJuncStrand_2,
		string& detectedFusionJuncFlankString,
		int& detectedAnchorLength_1, int& detectedAnchorLength_2,
		string& tmpFusionBreakPointInfoStr, Index_Info* indexInfo)
	{
		vector<string> fusionJuncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 9; tmp++)
		{
			int tabLoc = tmpFusionBreakPointInfoStr.find("\t", startLoc);
			string tmpDetectedFusionJuncField = tmpFusionBreakPointInfoStr.substr(startLoc, tabLoc-startLoc);
			fusionJuncFieldVec.push_back(tmpDetectedFusionJuncField);
			startLoc = tabLoc + 1;
		}
		fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc));
		string tmpChrNameStr_1 = fusionJuncFieldVec[1];
		string tmpChrNameStr_2 = fusionJuncFieldVec[2];
		string tmpChrPosStr_1 = fusionJuncFieldVec[3];
		string tmpChrPosStr_2 = fusionJuncFieldVec[4];
		detectedChrNameInt_1 = indexInfo->convertStringToInt(fusionJuncFieldVec[1]);
		detectedChrNameInt_2 = indexInfo->convertStringToInt(fusionJuncFieldVec[2]);
		detectedBreakPoint_1 = atoi(tmpChrPosStr_1.c_str());
		detectedBreakPoint_2 = atoi(tmpChrPosStr_2.c_str());
		detectedFusionJuncStrand_1 = fusionJuncFieldVec[5];
		detectedFusionJuncStrand_2 = fusionJuncFieldVec[6];
		detectedFusionJuncFlankString = fusionJuncFieldVec[7];
		string tmpAnchorLengthStr_1 = fusionJuncFieldVec[8];
		string tmpAnchorLengthStr_2 = fusionJuncFieldVec[9];
		detectedAnchorLength_1 = atoi(tmpAnchorLengthStr_1.c_str());
		detectedAnchorLength_2 = atoi(tmpAnchorLengthStr_2.c_str());
	}

	void insertUpdataFusionJuncMapVecVecAndFusionInfoVecWithDetectedFusionJunc(
		vector < vector< FusionJuncPosPair2fusionInfoVecIndexMap > >& fusionJuncMapVecVec, 
		int detectedChrNameInt_1, int detectedChrNameInt_2,
		int detectedBreakPoint_1, int detectedBreakPoint_2,
		string& detectedFusionJuncStrand_1, string& detectedFusionJuncStrand_2,
		string& detectedFusionJuncFlankString,
		int detectedAnchorLength_1, int detectedAnchorLength_2,
		vector<int>& fusionJuncSupNumVec,
		vector<string>& fusionJuncStrandVec_1,
		vector<string>& fusionJuncStrandVec_2,
		vector<string>& fusionJuncFlankStringVec, 
		vector<int>& fusionJuncAnchorLengthVec_1,
		vector<int>& fusionJuncAnchorLengthVec_2,
		Index_Info* indexInfo)
	{
		FusionJuncPosPair2fusionInfoVecIndexMap::iterator tmp1stMapIter 
			= ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).find(detectedBreakPoint_1);
		if(tmp1stMapIter == ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).end())
		{
			int tmpFusionInfoJuncSize = fusionJuncSupNumVec.size();
			FusionJuncEndPos2fusionInfoVecIndexMap tmp2ndMap;
			tmp2ndMap.insert(pair<int,int>(detectedBreakPoint_2, tmpFusionInfoJuncSize));
			((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).insert(
				pair<int,FusionJuncEndPos2fusionInfoVecIndexMap>(detectedBreakPoint_1, tmp2ndMap));
			fusionJuncSupNumVec.push_back(1);
			fusionJuncStrandVec_1.push_back(detectedFusionJuncStrand_1);
			fusionJuncStrandVec_2.push_back(detectedFusionJuncStrand_2);
			fusionJuncFlankStringVec.push_back(detectedFusionJuncFlankString);
			fusionJuncAnchorLengthVec_1.push_back(detectedAnchorLength_1);
			fusionJuncAnchorLengthVec_2.push_back(detectedAnchorLength_2);
		}
		else
		{
			FusionJuncEndPos2fusionInfoVecIndexMap::iterator tmp2ndMapIter
				= (tmp1stMapIter->second).find(detectedBreakPoint_2);
			if(tmp2ndMapIter == (tmp1stMapIter->second).end())
			{	
				int tmpFusionInfoJuncSize = fusionJuncSupNumVec.size();
				(tmp1stMapIter->second).insert(pair<int,int>(detectedBreakPoint_2, tmpFusionInfoJuncSize));
				fusionJuncSupNumVec.push_back(1);
				fusionJuncStrandVec_1.push_back(detectedFusionJuncStrand_1);
				fusionJuncStrandVec_2.push_back(detectedFusionJuncStrand_2);
				fusionJuncFlankStringVec.push_back(detectedFusionJuncFlankString);			
				fusionJuncAnchorLengthVec_1.push_back(detectedAnchorLength_1);
				fusionJuncAnchorLengthVec_2.push_back(detectedAnchorLength_2);
			}
			else
			{
				int tmpIndexInFusionInfoVec = (tmp2ndMapIter->second);
				fusionJuncSupNumVec[tmpIndexInFusionInfoVec]++;
				int existingAnchorLength_1 = fusionJuncAnchorLengthVec_1[tmpIndexInFusionInfoVec];
				int existingAnchorLength_2 = fusionJuncAnchorLengthVec_2[tmpIndexInFusionInfoVec];
				if(detectedAnchorLength_1 > existingAnchorLength_1)
				{
					fusionJuncAnchorLengthVec_1[tmpIndexInFusionInfoVec] = detectedAnchorLength_1;
				}
				if(detectedAnchorLength_2 > existingAnchorLength_2)
				{
					fusionJuncAnchorLengthVec_2[tmpIndexInFusionInfoVec] = detectedAnchorLength_2;
				}
			}
		}	
	}
	/*
	bool checkMatchTroughSeqSimilarity(string& strand_1, string& strand_2, 
		int chrNameInt_1, int chrNameInt_2, int breakPointPos_1, int breakPointPos_2,
		int anchorLength_1, int anchorLength_2, int toCheckAnchorLengthMax, 
		Index_Info* indexInfo, int& validAnchorLength_1, int& validAnchorLength_2,
		int& penalty_1, int& penalty_2)
	{
		if((strand_1 == "N") || (strand_2 == "N"))
			return false;
		validAnchorLength_1 = anchorLength_1;
		validAnchorLength_2 = anchorLength_2;
		if(validAnchorLength_1 > toCheckAnchorLengthMax)
			validAnchorLength_1 = toCheckAnchorLengthMax;
		if(validAnchorLength_2 > toCheckAnchorLengthMax)
			validAnchorLength_2 = toCheckAnchorLengthMax;

		string matchThroughSeq_gene1;
		string anchorSeq_gene2;

		string matchThroughSeq_gene2;
		string anchorSeq_gene1;
		// cout << "validAnchorLength_1: " << validAnchorLength_1 << endl;
		// cout << "validAnchorLength_2: " << validAnchorLength_2 << endl;
		// cout << "strand_1: " << strand_1 << endl;
		// cout << "strand_2: " << strand_2 << endl;
		if((strand_1 == "+")&&(strand_2 == "+"))
		{
			matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1+1, validAnchorLength_2);
			anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2, validAnchorLength_2);
			matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_1 + 1, validAnchorLength_1);
			anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - validAnchorLength_1 + 1, validAnchorLength_1);
		}
		else if((strand_1 == "-")&&(strand_2 == "-"))
		{
			matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_2 + 1, validAnchorLength_2));
			anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - validAnchorLength_2 + 1, validAnchorLength_2));
			matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_1));
			anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1, validAnchorLength_1));
		}
		else if((strand_1 == "+")&&(strand_2 == "-"))
		{
			matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1+1, validAnchorLength_2);
			anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - validAnchorLength_2 + 1, validAnchorLength_2));
			matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_1));
			anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - validAnchorLength_1 + 1, validAnchorLength_1);
		}
		else
		{
			matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_2 + 1, validAnchorLength_2));
			anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2, validAnchorLength_2);
			matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_1 + 1, validAnchorLength_1);
			anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1, validAnchorLength_1));
		}
		//cout << "matchThroughSeq_gene1: " << matchThroughSeq_gene1 << endl;
		//cout << "matchThroughSeq_gene2: " << matchThroughSeq_gene2 << endl;
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1;
		tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(
			matchThroughSeq_gene1, anchorSeq_gene2);
		penalty_1 = tmpFixNWDPinfo_1.getPenalty();
		//delete tmpFixNWDPinfo_1;
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_2;
		tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(
			matchThroughSeq_gene2, anchorSeq_gene1);
		penalty_2 = tmpFixNWDPinfo_2.getPenalty();
		//delete tmpFixNWDPinfo_2;
		return true;
	}*/

	void fusionJuncFiltering_anchorSeqSimilarity(
		AlignInferJunctionHash_Info* alignInferJuncHashInfo, SJhash_Info* SJ,

		Index_Info* indexInfo)
	{

	}

	void checkAlterSpliceAnchorSeqSimilarity_fusionJunc(
		AlignInferJunctionHash_Info* alignInferJuncHashInfo, SJhash_Info* SJ,
		string& strand_1, string& strand_2, int chrNameInt_1, int chrNameInt_2, 
		int breakPointPos_1, int breakPointPos_2, int anchorLength_1, int anchorLength_2, 
		int toCheckAnchorLengthMax, Index_Info* indexInfo, 
		vector<int>& validAnchorLengthVec_alterSiteAtGene1toReplaceGene2, 
		vector<int>& validAnchorLengthVec_alterSiteAtGene2toReplaceGene1,
		vector<int>& penaltyVec_alterSiteAtGene1toReplaceGene2, 
		vector<int>& penaltyVec_alterSiteAtGene2toReplaceGene1,
		vector<int>& alterSiteVec_atGene1toReplaceGene2, 
		vector<int>& alterSiteVec_atGene2toReplaceGene1)
	{
		if(((strand_1 == "+")&&(strand_2 == "+"))
			||((strand_1 == "-")&&(strand_2 == "-")))
		{
			int tmpChrNameInt_doner;
			int tmpChrNameInt_acceptor;
			int tmpChrPos_doner; 
			int tmpChrPos_acceptor;
			int tmpAnchorSize_doner;
			int tmpAnchorSize_acceptor;

			if((strand_1 == "+")&&(strand_2 == "+"))
			{
				tmpChrNameInt_doner = chrNameInt_1;
				tmpChrNameInt_acceptor = chrNameInt_2;
				tmpChrPos_doner = breakPointPos_1;
				tmpChrPos_acceptor = breakPointPos_2;
				tmpAnchorSize_doner = anchorLength_1;
				tmpAnchorSize_acceptor = anchorLength_2;
			}
			else
			{
				tmpChrNameInt_doner = chrNameInt_2;
				tmpChrNameInt_acceptor = chrNameInt_1;
				tmpChrPos_doner = breakPointPos_2;
				tmpChrPos_acceptor = breakPointPos_1;
				tmpAnchorSize_doner = anchorLength_2;
				tmpAnchorSize_acceptor = anchorLength_1;
			}
			// cout << "tmpChrNameInt_doner: " << tmpChrNameInt_doner << endl;
			// cout << "tmpChrNameInt_acceptor: " << tmpChrNameInt_acceptor << endl;
			// cout << "tmpChrPos_doner: " << tmpChrPos_doner << endl;
			// cout << "tmpChrPos_acceptor: " << tmpChrPos_acceptor << endl;
			// cout << "tmpAnchorSize_doner: " << tmpAnchorSize_doner << endl;
			// cout << "tmpAnchorSize_acceptor: " << tmpAnchorSize_acceptor << endl; 
			// check alternative acceptor start pos vec
			vector<int> tmpAlterAcceptorSpliceSitePosVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChrNameInt_doner,
				tmpChrPos_doner, tmpAlterAcceptorSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
				if((tmpAlterAcceptorSpliceSitePos != tmpChrPos_acceptor)
					||(tmpChrNameInt_doner != tmpChrNameInt_acceptor))
				{
					int tmpIndexInAlignInferJuncVec
						= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
							tmpChrNameInt_doner, tmpChrPos_doner, tmpAlterAcceptorSpliceSitePos);
					if(tmpIndexInAlignInferJuncVec < 0)
					{
						cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
						exit(1);
					}
					else
					{
						int tmpAlterAcceptorSite_anchorLength 
							= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
								tmpIndexInAlignInferJuncVec); 
						int tmpValidAnchorLength_acceptor = minAmongThreeNum(
							tmpAnchorSize_acceptor, tmpAlterAcceptorSite_anchorLength, toCheckAnchorLengthMax);
						string tmpAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_acceptor,
							tmpChrPos_acceptor, tmpValidAnchorLength_acceptor);
						string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_doner,
							tmpAlterAcceptorSpliceSitePos, tmpValidAnchorLength_acceptor);
						FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
						tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(tmpAcceptorSeq, tmpAlterAcceptorSeq);
						int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
						if((strand_1 == "+")&&(strand_2 == "+"))
						{
							validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_acceptor);
							penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterAcceptor);
							alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterAcceptorSpliceSitePos);
						}
						else
						{
							validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_acceptor);
							penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterAcceptor);
							alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterAcceptorSpliceSitePos);						
						}
					}
				}
			}
			// check alternative doner end pos vec
			vector<int> tmpAlterDonerSpliceSitePosVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChrNameInt_acceptor,
				tmpChrPos_acceptor, tmpAlterDonerSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
				if((tmpAlterDonerSpliceSitePos != tmpChrPos_doner)
					||(tmpChrNameInt_doner != tmpChrNameInt_acceptor))
				{
					int tmpIndexInAlignInferJuncVec
						= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
							tmpChrNameInt_acceptor, tmpAlterDonerSpliceSitePos, 
							tmpChrPos_acceptor);
					if(tmpIndexInAlignInferJuncVec < 0)
					{
						cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
						exit(1);
					}
					else
					{
						int tmpAlterSplice_anchorLength 
							= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
								tmpIndexInAlignInferJuncVec); 
						int tmpValidAnchorLength_doner = minAmongThreeNum(
							tmpAnchorSize_doner, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);
						string tmpDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_doner,
							tmpChrPos_doner - tmpValidAnchorLength_doner + 1, tmpValidAnchorLength_doner);
						string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_acceptor,
							tmpAlterDonerSpliceSitePos - tmpValidAnchorLength_doner + 1, 
							tmpValidAnchorLength_doner);
						FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
						tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
							tmpDonerSeq, tmpAlterDonerSeq);
						int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
						if((strand_1 == "+")&&(strand_2 == "+"))
						{
							validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_doner);
							penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterDoner);
							alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterDonerSpliceSitePos);
						}
						else
						{
							validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_doner);
							penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterDoner);
							alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterDonerSpliceSitePos);					
						}
					}
				}
			}
		}
		else if((strand_1 == "+")&&(strand_2 == "-"))
		{
			// donerEndPos = breakPointPos_1
			vector<int> tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt_1,
				breakPointPos_1, tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2);		
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2 
					= tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2[tmp];
				int tmpIndexInAlignInferJuncVec 
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_1,
						breakPointPos_1, tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterSplice_anchorLength
						= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
							tmpIndexInAlignInferJuncVec);
					int tmpValidAnchorLength_acceptor = minAmongThreeNum(
						anchorLength_2, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

					string tmpAcceptorSeq = convertStringToReverseComplement(
						indexInfo->returnChromStrSubstr(chrNameInt_2,
							breakPointPos_2 - tmpValidAnchorLength_acceptor + 1, 
							tmpValidAnchorLength_acceptor));
					string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(
						chrNameInt_1, tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2, 
						tmpValidAnchorLength_acceptor);
					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
					tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(
						tmpAcceptorSeq, tmpAlterAcceptorSeq);
					int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
					validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_acceptor);
					penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterAcceptor);
					alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2);
				}	
			}		
			// donerEndPos = breakPointPos_2		
			vector<int> tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt_2,
				breakPointPos_2, tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1);		
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1.size(); tmp++)
			{	
				int tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1
					= tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1[tmp];
				int tmpIndexInAlignInferJuncVec 
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_2, 
						breakPointPos_2, tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterSplice_anchorLength
						= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
							tmpIndexInAlignInferJuncVec);
					int tmpValidAnchorLength_acceptor = minAmongThreeNum(
						anchorLength_1, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);
					string tmpAcceptorSeq = convertStringToReverseComplement(
						indexInfo->returnChromStrSubstr(chrNameInt_1,
							breakPointPos_1 - tmpValidAnchorLength_acceptor + 1, 
							tmpValidAnchorLength_acceptor));
					string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(
						chrNameInt_2, tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1, 
						tmpValidAnchorLength_acceptor);
					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
					tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(
						tmpAcceptorSeq, tmpAlterAcceptorSeq);
					int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
					validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_acceptor);
					penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterAcceptor);
					alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1);
				}
			}	
		}
		else // strand_1 == "-" && strand_2 == "+"
		{
			// acceptorStartPos = breakPointPos_1
			vector<int> tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt_1,
				breakPointPos_1, tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2
					= tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2[tmp];
				int tmpIndexInAlignInferJuncVec
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_1,
						tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2, breakPointPos_1);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterSplice_anchorLength
						= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
							tmpIndexInAlignInferJuncVec);
					int tmpValidAnchorLength_doner = minAmongThreeNum(
						anchorLength_2, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

					string tmpDonerSeq = convertStringToReverseComplement(
						indexInfo->returnChromStrSubstr(chrNameInt_2,
							breakPointPos_2, tmpValidAnchorLength_doner));
					string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(chrNameInt_1, 
						tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2 - tmpValidAnchorLength_doner + 1, 
						tmpValidAnchorLength_doner);
				
					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
					tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
						tmpDonerSeq, tmpAlterDonerSeq);
					int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
					validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_doner);
					penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterDoner);
					alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2);			
				}
			}
			// acceptorStartPos = breakPointPos_2
			vector<int> tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt_2,
				breakPointPos_2, tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1
					= tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1[tmp];
				int tmpIndexInAlignInferJuncVec
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_2,
						tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1, breakPointPos_2);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterSplice_anchorLength
						= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
							tmpIndexInAlignInferJuncVec);
					int tmpValidAnchorLength_doner = minAmongThreeNum(
						anchorLength_1, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

					string tmpDonerSeq = convertStringToReverseComplement(
						indexInfo->returnChromStrSubstr(chrNameInt_1,
							breakPointPos_1, tmpValidAnchorLength_doner));
					string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(chrNameInt_2,
						tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1 - tmpValidAnchorLength_doner + 1,
						tmpValidAnchorLength_doner);

					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
					tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
						tmpDonerSeq, tmpAlterDonerSeq);
					int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
					validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_doner);
					penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterDoner);
					alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1);
				}
			}		
		}
	}

	void checkMatchTroughSeqSimilarity_fusionJunc(
		string& strand_1, string& strand_2, int chrNameInt_1, int chrNameInt_2, 
		int breakPointPos_1, int breakPointPos_2, int anchorLength_1, int anchorLength_2, 
		int toCheckAnchorLengthMax, Index_Info* indexInfo, 
		int& validAnchorLength_matchThroughAtGene1, int& validAnchorLength_matchThroughAtGene2,
		int& penalty_matchThroughAtGene1, int& penalty_matchThroughAtGene2)
	{
		// if((strand_1 == "N") || (strand_2 == "N"))
		// 	return false;
		validAnchorLength_matchThroughAtGene1 = anchorLength_2;
		validAnchorLength_matchThroughAtGene2 = anchorLength_1;
		if(validAnchorLength_matchThroughAtGene1 > toCheckAnchorLengthMax)
			validAnchorLength_matchThroughAtGene1 = toCheckAnchorLengthMax;
		if(validAnchorLength_matchThroughAtGene2 > toCheckAnchorLengthMax)
			validAnchorLength_matchThroughAtGene2 = toCheckAnchorLengthMax;

		string matchThroughSeq_gene1;
		string anchorSeq_gene2;

		string matchThroughSeq_gene2;
		string anchorSeq_gene1;
		// cout << "validAnchorLength_1: " << validAnchorLength_1 << endl;
		// cout << "validAnchorLength_2: " << validAnchorLength_2 << endl;
		// cout << "strand_1: " << strand_1 << endl;
		// cout << "strand_2: " << strand_2 << endl;
		if((strand_1 == "+")&&(strand_2 == "+"))
		{
			matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1+1, validAnchorLength_matchThroughAtGene1);
			anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2, validAnchorLength_matchThroughAtGene1);
			matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_matchThroughAtGene2 + 1, 
				validAnchorLength_matchThroughAtGene2);
			anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - validAnchorLength_matchThroughAtGene2 + 1, 
				validAnchorLength_matchThroughAtGene2);
		}
		else if((strand_1 == "-")&&(strand_2 == "-"))
		{
			matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_matchThroughAtGene1 + 1, 
				validAnchorLength_matchThroughAtGene1));
			anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - validAnchorLength_matchThroughAtGene1 + 1, 
				validAnchorLength_matchThroughAtGene1));
			matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_matchThroughAtGene2));
			anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1, validAnchorLength_matchThroughAtGene2));
		}
		else if((strand_1 == "+")&&(strand_2 == "-"))
		{
			matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1+1, validAnchorLength_matchThroughAtGene1);
			anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - validAnchorLength_matchThroughAtGene1 + 1, 
				validAnchorLength_matchThroughAtGene1));
			matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_matchThroughAtGene2));
			anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - validAnchorLength_matchThroughAtGene2 + 1, 
				validAnchorLength_matchThroughAtGene2);
		}
		else
		{
			matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_matchThroughAtGene1 + 1, 
				validAnchorLength_matchThroughAtGene1));
			anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2, validAnchorLength_matchThroughAtGene1);
			matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
				chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_matchThroughAtGene2 + 1, 
				validAnchorLength_matchThroughAtGene2);
			anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				chrNameInt_1, breakPointPos_1, validAnchorLength_matchThroughAtGene2));
		}
		//cout << "matchThroughSeq_gene1: " << matchThroughSeq_gene1 << endl;
		//cout << "matchThroughSeq_gene2: " << matchThroughSeq_gene2 << endl;
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1;
		tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(
			matchThroughSeq_gene1, anchorSeq_gene2);
		penalty_matchThroughAtGene1 = tmpFixNWDPinfo_1.getPenalty();
		//delete tmpFixNWDPinfo_1;
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_2;
		tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(
			matchThroughSeq_gene2, anchorSeq_gene1);
		penalty_matchThroughAtGene2 = tmpFixNWDPinfo_2.getPenalty();
		//delete tmpFixNWDPinfo_2;
		//return true;
	}

	void extractFusionJuncInfoFromStr(
		int& tmpChrNameInt_1, int& tmpChrNameInt_2,
		int& tmpBreakPointPos_1, int& tmpBreakPointPos_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int& tmpAnchorSize_1, int& tmpAnchorSize_2,
		string& tmpFusionJuncStr, Index_Info* indexInfo)
	{
		vector<string> tmpFusionJuncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 9; tmp++)
		{
			int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
			string tmpFusionJuncField = tmpFusionJuncStr.substr(startLoc, tabLoc-startLoc);
			tmpFusionJuncFieldVec.push_back(tmpFusionJuncField);
			startLoc = tabLoc + 1;
		}
		tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[0]);
		tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[1]);	
		string tmpChrPosStr_1 = tmpFusionJuncFieldVec[2];
		string tmpChrPosStr_2 = tmpFusionJuncFieldVec[3];
		tmpBreakPointPos_1 = atoi(tmpChrPosStr_1.c_str());
		tmpBreakPointPos_2 = atoi(tmpChrPosStr_2.c_str());
		tmpStrand_1 = tmpFusionJuncFieldVec[4];
		tmpStrand_2 = tmpFusionJuncFieldVec[5];
		tmpAnchorSize_1 = atoi(tmpFusionJuncFieldVec[7].c_str());
		tmpAnchorSize_2 = atoi(tmpFusionJuncFieldVec[8].c_str());
	}

	void initiateInputOutputFilePath_remapPairedAlignment()
	{
		// input
		input_incompletePairedAlignmentPath_remapOuterSoftClipUniquePairedSAM
			= output_folder_globalMapOuterSoftClipUniquePairedSAM + "/nonFusion.sam";
		input_fusionJunc_remapOuterSoftClipUniquePairedSAM
			= output_filteredFusionJunc_generateFusionJuncFromAdjustedBreakPointFile;
		// output
		string output_folder_path_remap = output_folder_path + "/remap";
		output_folder_remapOuterSoftClipUniquePairedSAM	= output_folder_path_remap;
		cout << "creating results folder for remap...." << endl;
		string mkdir_remap = "mkdir -p " + output_folder_path_remap;
		system(mkdir_remap.c_str());
		output_oriSAM_remapOuterSoftClipUniquePairedSAM
			= output_folder_path_remap + "/kept.sam";
		output_nonFusionSAM_remapOuterSoftClipUniquePairedSAM
			= output_folder_path_remap + "/nonFusion.sam";
		output_fusionSAM_remapOuterSoftClipUniquePairedSAM
			= output_folder_path_remap + "/fusion.sam";
		output_breakPoint_remapOuterSoftClipUniquePairedSAM 
			= output_folder_path_remap + "/breakPoint.txt";
		output_updatedFusionJunc_remapOuterSoftClipUniquePairedSAM 
			= output_folder_path_remap + "/fusionJunc_updated.txt";
	}

	void remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint(
		int threads_num, bool fasta_or_fastq_bool, Index_Info* indexInfo,
		ofstream& log_ofs)
	{
		time_t nowtime;
		struct tm *local;

		int chromNum = indexInfo->returnChromNum();
	   	ofstream oriSAM_ofs(output_oriSAM_remapOuterSoftClipUniquePairedSAM.c_str());
   		ofstream nonFusionSAM_ofs(output_nonFusionSAM_remapOuterSoftClipUniquePairedSAM.c_str());
   		ofstream fusionSAM_ofs(output_fusionSAM_remapOuterSoftClipUniquePairedSAM.c_str());
   		ofstream fusionReadsWithBreakPoint_ofs(output_breakPoint_remapOuterSoftClipUniquePairedSAM.c_str());

		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... MPS fusion starts to do remapping ......" << endl << endl;  
		log_ofs << endl << "[" << asctime(local) << "... MPS fusion starts to do remapping ......" << endl << endl; 		
		cout << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash ......" << endl << endl;

		FusionBreakPointHash_Info* fusionBreakPointHashInfo_merged = new FusionBreakPointHash_Info();
		fusionBreakPointHashInfo_merged->initiateWithChromNum(chromNum);
		fusionBreakPointHashInfo_merged->generateFusionBreakPointHashInfo_fromFuionJuncFile(
			input_fusionJunc_remapOuterSoftClipUniquePairedSAM, indexInfo);

		cout << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHashInfoVec ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHashInfoVec ......" << endl << endl;	

		vector<FusionBreakPointHash_Info*> fusionBreakPointHashInfoVec;
		for(int tmp = 0; tmp < threads_num; tmp++)
		{
			FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo = new FusionBreakPointHash_Info();
			tmpFusionBreakPointHashInfo->initiateWithChromNum(chromNum);
			tmpFusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromFuionJuncFile(
				input_fusionJunc_remapOuterSoftClipUniquePairedSAM, indexInfo);
			tmpFusionBreakPointHashInfo->clearSupportNum();
			tmpFusionBreakPointHashInfo->clearSupportNum_encompassing();
			fusionBreakPointHashInfoVec.push_back(tmpFusionBreakPointHashInfo);
		}

		string inputIncompleteUniquePairedAlignmentPath 
			= input_incompletePairedAlignmentPath_remapOuterSoftClipUniquePairedSAM;
		ifstream incompleteUniquePairedAlignment_ifs(inputIncompleteUniquePairedAlignmentPath.c_str());
		int normalRecordNum = normalRecordNum_remap;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;

		vector<string> inputSamStr1Vec(normalRecordNum);
		vector<string> inputSamStr2Vec(normalRecordNum);
		vector<string> outputOriSamStrVec(normalRecordNum);
		vector<string> outputNonfusionSamStrVec(normalRecordNum);
		vector<string> outputFusionSamStrVec(normalRecordNum);
		vector<string> outputFusionBreakPointInfoStrVec(normalRecordNum);

		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl;

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
						outputFusionBreakPointInfoStrVec[tmpOpenMP] = tmpFusionBreakPointInfoStr;
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

		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl;
		cout << endl << "[" << asctime(local) << "... start to output updated fusionJunc ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... start to output updated fusionJunc ......" << endl << endl;
		
		// sum up support numbers and output
		for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
			fusionBreakPointHashInfo_merged->addSupportNumFromAnotherFusionBreakPointHashInfo(
				fusionBreakPointHashInfoVec[tmp]);
		fusionBreakPointHashInfo_merged->outputFusionBreakPointHashInfoStr(
			output_updatedFusionJunc_remapOuterSoftClipUniquePairedSAM, indexInfo);
		fusionBreakPointHashInfo_merged->memoryFree();
		delete fusionBreakPointHashInfo_merged;
		for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		{
			fusionBreakPointHashInfoVec[tmp]->memoryFree();
			delete fusionBreakPointHashInfoVec[tmp];
		}

	   	oriSAM_ofs.close();
   		nonFusionSAM_ofs.close();
   		fusionSAM_ofs.close();
   		fusionReadsWithBreakPoint_ofs.close();
   		incompleteUniquePairedAlignment_ifs.close();
	}

	void initiateInputOutputFilePath_countUniqueUnpairedReadsEncompassingFusionBreakPoint()
	{
		// output
		//string output_folder_countEncompassingReads;
		string input_incomplete_unpaired_sam = "/phase2_output/fixHeadTail_incomplete_unpair.sam";
		string input_complete_unpaired_sam = "/phase2_output/fixHeadTail_complete_unpair.sam";

		string output_folder_path_countEncompassingReads = output_folder_path + "/countEncompassingReads";
		output_folder_countEncompassingReads = output_folder_path_countEncompassingReads;
		
		cout << "creating results folder for countEncompassingReads...." << endl;
		string mkdir_countEncompassingReads = "mkdir -p " + output_folder_path_countEncompassingReads;
		system(mkdir_countEncompassingReads.c_str());
		string input_unpairedSAM_merged = output_folder_countEncompassingReads + "/input_merged.sam";
		string merge_unpairedSAM_cmd = "cat " + input_incomplete_unpaired_sam 
				+ " " + input_complete_unpaired_sam + " > " + input_unpairedSAM_merged;
		system(merge_unpairedSAM_cmd.c_str());

		output_oriSAM_countEncompassingReads = output_folder_countEncompassingReads + "/kept.sam";
		output_nonFusionSAM_countEncompassingReads = output_folder_countEncompassingReads + "/nonFusion.sam";
		output_fusionSAM_countEncompassingReads = output_folder_countEncompassingReads + "/fusion.sam";
		output_fusionReadsEncompassingBreakPoint_countEncompassingReads
			= output_folder_countEncompassingReads + "/breakPoint.txt";
		output_updatedFusionJunc_countEncompassingReads
			= output_folder_countEncompassingReads + "/fusionJunc_updated.txt";

		// input
		input_unpairedSAM_countEncompassingReads = input_unpairedSAM_merged;
		input_fusionJunc_countEncompassingReads 
			= output_updatedFusionJunc_remapOuterSoftClipUniquePairedSAM;
	}

	void countEncompassingReads(
		int threads_num, bool fasta_or_fastq_bool, Index_Info* indexInfo,
		ofstream& log_ofs)
	{
		time_t nowtime;
		struct tm *local;		

		ofstream oriSAM_ofs(output_oriSAM_countEncompassingReads.c_str());	
		ofstream nonFusionSAM_ofs(output_nonFusionSAM_countEncompassingReads.c_str());;
		ofstream fusionEncompassingSAM_ofs(output_fusionSAM_countEncompassingReads.c_str());
		ofstream fusionEncompassingReadsWithreakPoint_ofs(
			output_fusionReadsEncompassingBreakPoint_countEncompassingReads.c_str());		
	
		int chromNum = indexInfo->returnChromNum();

		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... start to generate FusionBreakPointHash_merged ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... start to generate fusionBreakPointHash_merged ......" << endl << endl; 

		FusionBreakPointHash_Info* fusionBreakPointHashInfo_merged
			= new FusionBreakPointHash_Info();
		fusionBreakPointHashInfo_merged->initiateWithChromNum(chromNum);
		//string inputFusionBreakPointPath = argv[2];
		fusionBreakPointHashInfo_merged->generateFusionBreakPointHashInfo_fromFuionJuncFile(
			input_fusionJunc_countEncompassingReads, indexInfo);

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
				input_fusionJunc_countEncompassingReads, indexInfo);
			tmpFusionBreakPointHashInfo->clearSupportNum();
			tmpFusionBreakPointHashInfo->clearSupportNum_encompassing();
			fusionBreakPointHashInfoVec.push_back(tmpFusionBreakPointHashInfo);
		}

		nowtime = time(NULL);
		//struct tm *local;
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... start to count encompassing unique unpaired reads ......" << endl << endl;
		log_ofs << endl << "[" << asctime(local) << "... start to count encompassing unique unpaired reads ......" << endl << endl;

		string inputUniqueUnpairedAlignmentPath = input_unpairedSAM_countEncompassingReads;
		ifstream uniqueUnpairedAlignments_ifs(inputUniqueUnpairedAlignmentPath.c_str());
		int normalRecordNum = normalRecordNum_countEncompassingReads;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;

		vector<string> inputSamStr1Vec(normalRecordNum);
		vector<string> inputSamStr2Vec(normalRecordNum);
		vector<string> outputOriSamStrVec(normalRecordNum);
		vector<string> outputNonfusionSamStrVec(normalRecordNum);
		vector<string> outputFusionEncompassingSamStrVec(normalRecordNum);
		vector<string> outputFusionEncompassingReadsWithBreakPointStrVec(normalRecordNum);

		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
		log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 

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
						fusionBreakPointHashInfoVec[threadNO]->addSupportNum_withIndexInBreakPointInfoVec_encompassing(1, fusionBreakPointEncompassed_found_index);		
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
			output_updatedFusionJunc_countEncompassingReads, indexInfo);

		uniqueUnpairedAlignments_ifs.close();
		fusionBreakPointHashInfo_merged->memoryFree();
		delete fusionBreakPointHashInfo_merged;
		for(int tmp = 0; tmp < fusionBreakPointHashInfoVec.size(); tmp++)
		{
			fusionBreakPointHashInfoVec[tmp]->memoryFree();
			delete fusionBreakPointHashInfoVec[tmp];
		}

		oriSAM_ofs.close();
		nonFusionSAM_ofs.close();
		fusionEncompassingSAM_ofs.close();
		fusionEncompassingReadsWithreakPoint_ofs.close();
		uniqueUnpairedAlignments_ifs.close();	
	}
};
#endif