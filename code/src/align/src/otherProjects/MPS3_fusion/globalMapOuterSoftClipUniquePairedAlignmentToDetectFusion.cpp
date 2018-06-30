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

using namespace std;

bool adjustFusionBreakPoint_searchForCanonicalOnesWithOffset_allowVaryOffsetAtDiffEnd(
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
	Index_Info* indexInfo, int offset_max, int& fusionCase, 
	bool allowInsertedBaseOrNotBool)
{
	int tmp_offset_max = offset_max;
	int shorterAnchorLength = leftAnchorLength;
	if(rightAnchorLength < leftAnchorLength)
		shorterAnchorLength = rightAnchorLength;
	if(tmp_offset_max > shorterAnchorLength - 1)
		tmp_offset_max = shorterAnchorLength - 1;
	//cout << "tmp_offset_max: " << tmp_offset_max << endl;
	int leftChrNameInt = indexInfo->convertStringToInt(leftChrNameStr);
	int rightChrNameInt = indexInfo->convertStringToInt(rightChrNameStr);
	//vector<int> candiOffsetVec;
	// candiOffsetVec.push_back(0);
	// for(int tmp = 0; tmp < tmp_offset_max; tmp++)
	// {	
	// 	candiOffsetVec.push_back(tmp+1);
	// 	candiOffsetVec.push_back(-1-tmp);
	// }
	vector< pair<int,int> > candiOffsetPairVec;
	candiOffsetPairVec.push_back(pair<int,int>(0,0));
	for(int tmp = 0; tmp < tmp_offset_max; tmp++)
	{		
		candiOffsetPairVec.push_back(pair<int,int>(tmp+1,tmp+1));
		candiOffsetPairVec.push_back(pair<int,int>(-1-tmp,-1-tmp));
	}
	if(allowInsertedBaseOrNotBool)
	{	
		for(int tmp = 1 - tmp_offset_max; tmp < tmp_offset_max; tmp++)
		{
			for(int tmp2 = 1 - tmp_offset_max; tmp2 < tmp_offset_max; tmp2++)
		 		candiOffsetPairVec.push_back(pair<int,int>(tmp,tmp2));
		}	
	}
	if(leftGene_NorOrRcm_bool && rightGene_NorOrRcm_bool)
	{
		for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetPairVec.size(); tmpOffsetIndex ++)
		{
			int tmpOffset_1 = candiOffsetPairVec[tmpOffsetIndex].first;
			int tmpOffset_2 = candiOffsetPairVec[tmpOffsetIndex].second;
			string tmpRawFusionJuncFlankString;
			bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
				leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset_1, rightBreakPointPos + tmpOffset_2,
				true, true, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
			if(!try_getRawFusionJuncFlankString_bool)
				continue;
			if(tmpRawFusionJuncFlankString == "GTAG")
			{
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset_1;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset_2;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
				canonicalAndStrandIdenfiedStrand_1 = "+";
				canonicalAndStrandIdenfiedStrand_2 = "+";
				canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset_1;
				canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset_2;
				fusionCase = 2;
				return true;
			}
			else if(tmpRawFusionJuncFlankString == "CTAC")
			{
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset_1;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset_2;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
				canonicalAndStrandIdenfiedStrand_1 = "-";
				canonicalAndStrandIdenfiedStrand_2 = "-";
				canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset_1;
				canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset_2;
				fusionCase = 5;
				return true;
			}
			else
			{}
		}
		return false;	
	}
	else if((!leftGene_NorOrRcm_bool) && (!rightGene_NorOrRcm_bool))
	{
		for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetPairVec.size(); tmpOffsetIndex ++)
		{
			int tmpOffset_1 = candiOffsetPairVec[tmpOffsetIndex].first;
			int tmpOffset_2 = candiOffsetPairVec[tmpOffsetIndex].second;
			string tmpRawFusionJuncFlankString;
			bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
				leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset_1, rightBreakPointPos + tmpOffset_2,
				false, false, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
			if(!try_getRawFusionJuncFlankString_bool)
				continue;
			if(tmpRawFusionJuncFlankString == "GTAG")
			{
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset_1;
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset_2;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
				canonicalAndStrandIdenfiedStrand_1 = "+";
				canonicalAndStrandIdenfiedStrand_2 = "+";
				canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset_1;
				canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset_2;
				fusionCase = 1;
				return true;
			}
			else if(tmpRawFusionJuncFlankString == "CTAC")
			{
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset_1;
				canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset_2;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
				canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
				canonicalAndStrandIdenfiedStrand_2 = "-";
				canonicalAndStrandIdenfiedStrand_1 = "-";
				canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset_1;
				canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset_2;	
				fusionCase = 4;
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
			for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetPairVec.size(); tmpOffsetIndex ++)
			{
				int tmpOffset_1 = candiOffsetPairVec[tmpOffsetIndex].first;
				int tmpOffset_2 = candiOffsetPairVec[tmpOffsetIndex].second;
				string tmpRawFusionJuncFlankString;
				bool try_getRawFusionJuncFlankString_bool  = indexInfo->try_returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset_1, rightBreakPointPos - tmpOffset_2,
					false, true, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
				if(!try_getRawFusionJuncFlankString_bool)
					continue;
				if(tmpRawFusionJuncFlankString == "GTCT")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos + tmpOffset_1;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos - tmpOffset_2;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset_1;
					canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset_2;
					canonicalAndStrandIdenfiedStrand_2 = "-";
					canonicalAndStrandIdenfiedStrand_1 = "+";
					fusionCase = 8;
					return true;
				}
				else if(tmpRawFusionJuncFlankString == "CTGT")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos + tmpOffset_1;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos - tmpOffset_2;					
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset_1;
					canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset_2;
					canonicalAndStrandIdenfiedStrand_2 = "-";
					canonicalAndStrandIdenfiedStrand_1 = "+";
					fusionCase = 7;
					return true;				
				}
				else
				{}
			}		
		}
		else
		{
			for(int tmpOffsetIndex = 0; tmpOffsetIndex < candiOffsetPairVec.size(); tmpOffsetIndex ++)
			{
				int tmpOffset_1 = candiOffsetPairVec[tmpOffsetIndex].first;
				int tmpOffset_2 = candiOffsetPairVec[tmpOffsetIndex].second;
				string tmpRawFusionJuncFlankString;
				bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt, leftBreakPointPos - tmpOffset_1, rightBreakPointPos + tmpOffset_2,
					false, true, (breakPointInEnd1OrEnd2_bool), tmpRawFusionJuncFlankString);
				if(!try_getRawFusionJuncFlankString_bool)
					continue;
				if(tmpRawFusionJuncFlankString == "AGAC")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = leftBreakPointPos - tmpOffset_1;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = rightBreakPointPos + tmpOffset_2;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_2 = leftAnchorLength + tmpOffset_1;
					canonicalAndStrandIdenfiedAnchorLength_1 = rightAnchorLength - tmpOffset_2;
					canonicalAndStrandIdenfiedStrand_2 = "+";
					canonicalAndStrandIdenfiedStrand_1 = "-";
					fusionCase = 11;
					return true;
				}
				else if(tmpRawFusionJuncFlankString == "ACAG")
				{
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_1 = leftBreakPointPos - tmpOffset_1;
					canonicalAndStrandIdenfiedCanonicalFusionBreakPoint_2 = rightBreakPointPos + tmpOffset_2;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_1 = leftChrNameStr;
					canonicalAndStrandIdenfiedCanonicalFusionGeneChrName_2 = rightChrNameStr;
					canonicalAndStrandIdenfiedAnchorLength_1 = leftAnchorLength + tmpOffset_1;
					canonicalAndStrandIdenfiedAnchorLength_2 = rightAnchorLength - tmpOffset_2;
					canonicalAndStrandIdenfiedStrand_2 = "+";
					canonicalAndStrandIdenfiedStrand_1 = "-";
					fusionCase = 10;
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
		cout << "error in adjustFusionBreakPoint_searchForCanonicalOnesWithOffset from ";
		cout << "globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion.cpp" << endl;
		return false;	
	}
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
	Index_Info* indexInfo, int offset_max, int& fusionCase)
{
	// cout << "start to do adjustFusionBreakPoint_searchForCanonicalOnesWithOffset" << endl;
	// cout << "breakPointInEnd1OrEnd2_bool: " << breakPointInEnd1OrEnd2_bool << endl;
	// cout << "leftGene_NorOrRcm_bool: " << leftGene_NorOrRcm_bool << endl;
	// cout << "rightGene_NorOrRcm_bool: " << rightGene_NorOrRcm_bool << endl;
	//return false;
	int tmp_offset_max = offset_max;
	int shorterAnchorLength = leftAnchorLength;
	if(rightAnchorLength < leftAnchorLength)
		shorterAnchorLength = rightAnchorLength;
	if(tmp_offset_max > shorterAnchorLength - 1)
		tmp_offset_max = shorterAnchorLength - 1;
	//cout << "tmp_offset_max: " << tmp_offset_max << endl;
	int leftChrNameInt = indexInfo->convertStringToInt(leftChrNameStr);
	int rightChrNameInt = indexInfo->convertStringToInt(rightChrNameStr);
	// cout << "leftChrNameInt: " << leftChrNameInt << endl;
	// cout << "rightChrNameInt: " << rightChrNameInt << endl;
	// cout << "leftBreakPointPos: " << leftBreakPointPos << endl;
	// cout << "rightBreakPointPos: " << rightBreakPointPos << endl;	
	// cout << "leftAnchorLength: " << leftAnchorLength << endl;
	// cout << "rightAnchorLength: " << rightAnchorLength << endl;
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
			string tmpRawFusionJuncFlankString;
			bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
				leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
				true, true, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
			if(!try_getRawFusionJuncFlankString_bool)
				continue;
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
				fusionCase = 2;
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
				fusionCase = 5;
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
			//cout << "leftBreakPointPos+tmpOffset: " << leftBreakPointPos + tmpOffset << endl;
			//cout << "rightBreakPointPos+tmpOffset: " << rightBreakPointPos + tmpOffset << endl;
			string tmpRawFusionJuncFlankString;
			bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
				leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
				false, false, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
			//cout << "try_getRawFusionJuncFlankString_bool: " << try_getRawFusionJuncFlankString_bool << endl;
			if(!try_getRawFusionJuncFlankString_bool)
				continue;
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
				fusionCase = 1;
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
				fusionCase = 4;
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
				string tmpRawFusionJuncFlankString;
				bool try_getRawFusionJuncFlankString_bool  = indexInfo->try_returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt, leftBreakPointPos + tmpOffset, rightBreakPointPos - tmpOffset,
					false, true, breakPointInEnd1OrEnd2_bool, tmpRawFusionJuncFlankString);
				if(!try_getRawFusionJuncFlankString_bool)
					continue;
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
					fusionCase = 8;
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
					fusionCase = 7;
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
				string tmpRawFusionJuncFlankString;
				bool try_getRawFusionJuncFlankString_bool = indexInfo->try_returnRawFusionJuncFlankString(
					leftChrNameInt, rightChrNameInt,
					leftBreakPointPos - tmpOffset, rightBreakPointPos + tmpOffset,
					false, true, (breakPointInEnd1OrEnd2_bool), tmpRawFusionJuncFlankString);
				if(!try_getRawFusionJuncFlankString_bool)
					continue;
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
					fusionCase = 11;
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
					fusionCase = 10;
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

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputIndexInfoPath ";
		cout << "inputIncompleteUniquePairedAlignmentPath ";
		cout << "outputFolder threads_num fasta_or_fastq min_fusion_distance allowInsertedBaseOrNot" << endl;
		exit(1);
	}
	string allowInsertedBaseOrNotStr = argv[7];
	bool allowInsertedBaseOrNotBool = false;
	if((allowInsertedBaseOrNotStr == "Y")||(allowInsertedBaseOrNotStr == "y")
		||(allowInsertedBaseOrNotStr == "Yes")||(allowInsertedBaseOrNotStr == "yes")
		||(allowInsertedBaseOrNotStr == "YES"))
		allowInsertedBaseOrNotBool = true;

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
	vector<string> outputFusionSamStrVec_withBreakPoint(normalRecordNum);
	vector<string> outputFusionSamStrVec_withBreakPoint_adjusted(normalRecordNum);
	vector<string> outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering(normalRecordNum);
	vector<string> outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq(normalRecordNum);
	vector<string> outputFusionSamStrVec_withBreakPoint_adjusted_mapRange(normalRecordNum);
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
			//cout << "start to initiate incompleteUniquePairedAlignmentInfo" << endl;
			bool incompleteUniquePaired_bool
				= incompleteUniquePairedAlignmentInfo.initiateWith2samStr_globalMapOuterUnfixedEnd2detectFusionBreakPoint(
					inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
			//cout << "incompleteUniquePaired_bool: " << incompleteUniquePaired_bool << endl;
			// if(incompleteUniquePaired_bool)
			// {	
			// 	cout << endl << "****************************************************" << endl;
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
					// cout << "readSeq_Nor1: " << endl << readInfo_unfixedLeftReadHead.returnReadSeq_1() << endl;
					// cout << "readSeq_Rcm1: " << endl << readInfo_unfixedLeftReadHead.returnRcmReadSeq_1() << endl;
					// cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info_unfixedLeftReadHead.coutSegInfo_Nor1(indexInfo);
					// cout << "segInfo_Rcm1: " << endl; fixPhase1Info_unfixedLeftReadHead.coutSegInfo_Rcm1(indexInfo);
					// cout << "segInfo_Nor2: " << endl; fixPhase1Info_unfixedLeftReadHead.coutSegInfo_Nor2(indexInfo);
					// cout << "segInfo_Rcm2: " << endl; fixPhase1Info_unfixedLeftReadHead.coutSegInfo_Rcm2(indexInfo);	
					//cout << "start to do fixPhase1_pathInfo" << endl;
					fixPhase1Info_unfixedLeftReadHead.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedLeftReadHead, 
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
					//cout << "getPossiPathInfo_Nor1: " << endl; fixPhase1Info_unfixedLeftReadHead.coutPossiPathStr_Nor1();
					//cout << "getPossiPathInfo_Rcm1: " << endl; fixPhase1Info_unfixedLeftReadHead.coutPossiPathStr_Rcm1();
					//cout << "start to do fixPhase1_gapInfo" << endl;
					fixPhase1Info_unfixedLeftReadHead.fixPhase1_gapInfo(readInfo_unfixedLeftReadHead, indexInfo, Do_cirRNA,
						Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
					//cout << "output map info: " << endl; fixPhase1Info_unfixedLeftReadHead.coutDebugInfo(readInfo_unfixedLeftReadHead, indexInfo, SE_or_PE_bool);
					//cout << "finish output map info ..." << endl;
					//cout << "start to initiatePeAlignInfo" << endl;
					peAlignInfo_unfixedLeftReadHead.initiatePeAlignInfo(fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor1,
						fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm1, fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor2,
						fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm2, indexInfo, true);
					//cout << endl << endl << "peAlignInfo: \n" << peAlignInfo_unfixedLeftReadHead.returnPeAlignInfoStr() << endl << endl;
					//cout << "start to do alignments filtering and selection" << endl;
					//cout << "start to chooseBestAlignment_final_SE" << endl;
					peAlignInfo_unfixedLeftReadHead.chooseBestAlignment_final_SE();
					leftReadHeadUnfixed_fixed_bool 
						= peAlignInfo_unfixedLeftReadHead.mappedForFusionDetection_unfixedHead_uniqueCompleteMapped_SE_bool(
							min_fusion_distance, incompleteUniquePairedAlignmentInfo.returnChrNameInt(), 
							incompleteUniquePairedAlignmentInfo.returnStartPos_1(), indexInfo);
					//cout << "leftReadHeadUnfixed_fixed_bool: " << leftReadHeadUnfixed_fixed_bool << endl;
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
					//cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info_unfixedRightReadTail.coutSegInfo_Nor1(indexInfo);
					//cout << "segInfo_Rcm1: " << endl; fixPhase1Info_unfixedRightReadTail.coutSegInfo_Rcm1(indexInfo);
					//cout << "segInfo_Nor2: " << endl; fixPhase1Info_unfixedRightReadTail.coutSegInfo_Nor2(indexInfo);
					//cout << "segInfo_Rcm2: " << endl; fixPhase1Info_unfixedRightReadTail.coutSegInfo_Rcm2(indexInfo);	
					//cout << "start to do fixPhase1_pathInfo" << endl;
					fixPhase1Info_unfixedRightReadTail.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedRightReadTail, 
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
					//cout << "start to do fixPhase1_gapInfo" << endl;
					fixPhase1Info_unfixedRightReadTail.fixPhase1_gapInfo(readInfo_unfixedRightReadTail, indexInfo, Do_cirRNA,
						Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
					//cout << "start to initiatePeAlignInfo" << endl;
					//cout << "output map info: " << endl; fixPhase1Info_unfixedRightReadTail.coutDebugInfo(readInfo_unfixedLeftReadHead, indexInfo, SE_or_PE_bool);
					peAlignInfo_unfixedRightReadTail.initiatePeAlignInfo(fixPhase1Info_unfixedRightReadTail.pathInfo_Nor1,
						fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm1, fixPhase1Info_unfixedRightReadTail.pathInfo_Nor2,
						fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm2, indexInfo, true);
					//cout << "start to chooseBestAlignment_final_SE" << endl;
					//cout << endl << endl << "peAlignInfo: \n" << peAlignInfo_unfixedRightReadTail.returnPeAlignInfoStr() << endl << endl;
					peAlignInfo_unfixedRightReadTail.chooseBestAlignment_final_SE();
					rightReadTailUnfixed_fixed_bool 
						= peAlignInfo_unfixedRightReadTail.mappedForFusionDetection_unfixedTail_uniqueCompleteMapped_SE_bool(
							min_fusion_distance, incompleteUniquePairedAlignmentInfo.returnChrNameInt(), 
							incompleteUniquePairedAlignmentInfo.returnEndPos_2(), indexInfo);
					//cout << "rightReadTailUnfixed_fixed_bool: " << rightReadTailUnfixed_fixed_bool << endl;
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
			outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP] = "";
			outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmpOpenMP] = "";
			outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmpOpenMP] = "";
			//cout << "start to generte all results str ...." << endl;
			if(incompleteUniquePaired_bool)
			{
				if(leftReadHeadUnfixed_fixed_bool || rightReadTailUnfixed_fixed_bool)
				{
					string tmpSAMoutputStr = inputSamStr1Vec[tmpOpenMP];
					if(leftReadHeadUnfixed_fixed_bool)
					{
						if(leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool)
							tmpSAMoutputStr = leftReadHeadUnfixed_fixed_SAMstr + "\tZF:Z:FUS_leftHead_for\n" 
								+ tmpSAMoutputStr + "\tZF:Z:FUS_leftHead_for";
						else
							tmpSAMoutputStr =leftReadHeadUnfixed_fixed_SAMstr + "\tZF:Z:FUS_leftHead_rev\n"  
								+ tmpSAMoutputStr + "\tZF:Z:FUS_leftHead_rev";							
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
					//cout << "tmpSAMoutputStr: " << tmpSAMoutputStr << endl;
					string tmpLeftReadHeadBreakPoint;
					string tmpRightReadTailBreakPoint;
					string tmpLeftReadHeadBreakPoint_adjusted;
					string tmpRightReadTailBreakPoint_adjusted;	
					string tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq;
					string tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq;
					string tmpLeftReadHeadBreakPoint_adjusted_mapRange;
					string tmpRightReadTailBreakPoint_adjusted_mapRange;
					 // cout << "start to generate breakPointStr" << endl;
					//cout << "leftReadHeadUnfixed_fixed_bool: " << leftReadHeadUnfixed_fixed_bool << endl;
					//cout << "rightReadTailUnfixed_fixed_bool: " << rightReadTailUnfixed_fixed_bool << endl;			
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
						int tmpOriSAM_endPos_2
							= incompleteUniquePairedAlignmentInfo.returnEndPos_2();
						int tmpFusionAnchorLength_left 
						// this should the length of the last jumpcode in fixedLeftReadHead, 
						// for now = unfixedHeadLength_1() as only jumpcode M is accepted for this unfixedReadHead 
							= incompleteUniquePairedAlignmentInfo.returnUnfixedHeadLength_1();
						int tmpFusionAnchorLength_right 
							= incompleteUniquePairedAlignmentInfo.returnJumpCodeLength_1(1);
						int tmpXM = incompleteUniquePairedAlignmentInfo.returnXM();
						//cout << "leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool: " << leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool << endl;
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
							int fusionCase;
							bool canonicalAndStrandIdenfiedOrNot 
								= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset_allowVaryOffsetAtDiffEnd(
									true, true, true, tmpLeftReadHead_chrName, tmpOriSAM_chrName,
									tmpLeftReadHead_endMapPos, tmpOriSAM_startPos_1,
									tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
									canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
									canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
									canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
									canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
									indexInfo, offsetForAdjustingFusionBreakPoint_max, fusionCase, allowInsertedBaseOrNotBool);
							//cout << "canonicalAndStrandIdenfiedOrNot: " << canonicalAndStrandIdenfiedOrNot << endl;	
							if(canonicalAndStrandIdenfiedOrNot)
							{
								tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
									+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
									+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right) + "\t" 
									+ int_to_str(tmpXM) + "\t" + int_to_str(fusionCase) + ",";
								tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedBreakPoint_left,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedChrName_right, 
										canonicalAndStrandIdenfiedBreakPoint_right, canonicalAndStrandIdenfiedStrand_right, 
										anchorSeqSize2output);
								int tmpMapRange_most_1, tmpMapRange_most_2;
								incompleteUniquePairedAlignmentInfo.returnFusionSpanningReadMapRange_canonical(
									fusionCase, tmpLeftReadHead_startMapPos, tmpLeftReadHead_endMapPos,
									tmpMapRange_most_1, tmpMapRange_most_2);
								tmpLeftReadHeadBreakPoint_adjusted_mapRange 
									= int_to_str(tmpMapRange_most_1) + "\t" + int_to_str(tmpMapRange_most_2);
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM) + "\t2,5,";
								tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(tmpLeftReadHead_chrName, tmpLeftReadHead_endMapPos, positiveStrand,
										tmpOriSAM_chrName, tmpOriSAM_startPos_1, positiveStrand, anchorSeqSize2output) + "\t"
										+ indexInfo->returnFusionAnchorSeq(tmpOriSAM_chrName, tmpOriSAM_startPos_1, negativeStrand,
											tmpLeftReadHead_chrName, tmpLeftReadHead_endMapPos, negativeStrand, anchorSeqSize2output);							
								tmpLeftReadHeadBreakPoint_adjusted_mapRange 
									= int_to_str(tmpLeftReadHead_startMapPos) + "\t" + int_to_str(tmpOriSAM_endPos_2);
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM);// + "\t2,5,";
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
							int fusionCase;
							bool canonicalAndStrandIdenfiedOrNot
								= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset_allowVaryOffsetAtDiffEnd(
									true, false, true, tmpLeftReadHead_chrName, tmpOriSAM_chrName,
									tmpLeftReadHead_startMapPos, tmpOriSAM_startPos_1,
									tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
									canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
									canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
									canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
									canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,									
									indexInfo, offsetForAdjustingFusionBreakPoint_max, fusionCase, allowInsertedBaseOrNotBool);
							//cout << "canonicalAndStrandIdenfiedOrNot: " << canonicalAndStrandIdenfiedOrNot << endl;
							if(canonicalAndStrandIdenfiedOrNot)
							{
								tmpLeftReadHeadBreakPoint_adjusted = tmpLeftReadHeadBreakPoint
									+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
									+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right) + "\t"
									+ int_to_str(tmpXM) + "\t" + int_to_str(fusionCase) + ",";
								tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedBreakPoint_left,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedChrName_right, 
										canonicalAndStrandIdenfiedBreakPoint_right, canonicalAndStrandIdenfiedStrand_right, 
										anchorSeqSize2output);
								int tmpMapRange_most_1, tmpMapRange_most_2;
								incompleteUniquePairedAlignmentInfo.returnFusionSpanningReadMapRange_canonical(
									fusionCase, tmpLeftReadHead_startMapPos, tmpLeftReadHead_endMapPos,
									tmpMapRange_most_1, tmpMapRange_most_2);
								tmpLeftReadHeadBreakPoint_adjusted_mapRange 
									= int_to_str(tmpMapRange_most_1) + "\t" + int_to_str(tmpMapRange_most_2);
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM) + "\t10,11,";  
								tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(tmpLeftReadHead_chrName, tmpLeftReadHead_startMapPos, negativeStrand,
										tmpOriSAM_chrName, tmpOriSAM_startPos_1, positiveStrand, anchorSeqSize2output) + "\t"
										+ indexInfo->returnFusionAnchorSeq(tmpOriSAM_chrName, tmpOriSAM_startPos_1, negativeStrand,
											tmpLeftReadHead_chrName, tmpLeftReadHead_startMapPos, positiveStrand, anchorSeqSize2output);
								tmpLeftReadHeadBreakPoint_adjusted_mapRange 
									= int_to_str(tmpLeftReadHead_endMapPos) + "\t" + int_to_str(tmpOriSAM_endPos_2);
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM);// + "\t10,11,";
						}
					}
					if(rightReadTailUnfixed_fixed_bool)
					{
						tmpRightReadTailBreakPoint 
							+= incompleteUniquePairedAlignmentInfo.returnReadName_2();
						tmpRightReadTailBreakPoint += "\t";
						string tmpOriSAM_chrName
							= incompleteUniquePairedAlignmentInfo.returnChrName(indexInfo);
						int tmpOriSAM_startPos_1 
							= incompleteUniquePairedAlignmentInfo.returnStartPos_1();
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
						int tmpXM = incompleteUniquePairedAlignmentInfo.returnXM();
						//cout << "rightReadTailUnfixed_fixed_Nor_or_Rcm_bool: " << rightReadTailUnfixed_fixed_Nor_or_Rcm_bool << endl;
						if(rightReadTailUnfixed_fixed_Nor_or_Rcm_bool) // right read rev -- rev case 1, 4
						{
							int canonicalAndStrandIdenfiedBreakPoint_left;
							string canonicalAndStrandIdenfiedChrName_left;
							int canonicalAndStrandIdenfiedBreakPoint_right;
							string canonicalAndStrandIdenfiedChrName_right;
							string canonicalAndStrandIdenfiedStrand_left;
							string canonicalAndStrandIdenfiedStrand_right;
							int canonicalAndStrandIdenfiedAnchorLength_left;
							int canonicalAndStrandIdenfiedAnchorLength_right;
							int fusionCase;
							bool canonicalAndStrandIdenfiedOrNot 
								= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset_allowVaryOffsetAtDiffEnd(
									false, false, false, tmpOriSAM_chrName, tmpRightReadTail_chrName,
									tmpOriSAM_endPos_2, tmpRightReadTail_startPos,
									tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
									canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
									canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
									canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
									canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
									indexInfo, offsetForAdjustingFusionBreakPoint_max, fusionCase, allowInsertedBaseOrNotBool);
							//cout << "canonicalAndStrandIdenfiedOrNot: " << canonicalAndStrandIdenfiedOrNot << endl;				
							if(canonicalAndStrandIdenfiedOrNot)
							{
								tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
									+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
									+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right) + "\t"
									+ int_to_str(tmpXM) + "\t" + int_to_str(fusionCase) + ",";
								tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedBreakPoint_left,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedChrName_right, 
										canonicalAndStrandIdenfiedBreakPoint_right, canonicalAndStrandIdenfiedStrand_right, 
										anchorSeqSize2output);
								int tmpMapRange_most_1, tmpMapRange_most_2;
								incompleteUniquePairedAlignmentInfo.returnFusionSpanningReadMapRange_canonical(
									fusionCase, tmpRightReadTail_startPos, tmpRightReadTail_endPos,
									tmpMapRange_most_1, tmpMapRange_most_2);
								tmpRightReadTailBreakPoint_adjusted_mapRange 
									= int_to_str(tmpMapRange_most_1) + "\t" + int_to_str(tmpMapRange_most_2);
							}
							else
							{
								tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
									+ tmpOriSAM_chrName + "\t" + tmpRightReadTail_chrName + "\t" 
									+ int_to_str(tmpOriSAM_endPos_2) + "\t"
									+ int_to_str(tmpRightReadTail_startPos) + "\tN\tN\t" 
									+ indexInfo->returnRawFusionJuncFlankString(
										tmpOriSAM_chrName, tmpRightReadTail_chrName, 
										tmpOriSAM_endPos_2, tmpRightReadTail_startPos,
										false, false, indexInfo, false)
									+ "\t" + int_to_str(tmpFusionAnchorLength_left)
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM) + "\t1,4,";
								tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(tmpOriSAM_chrName, tmpOriSAM_endPos_2, positiveStrand,
										tmpRightReadTail_chrName, tmpRightReadTail_startPos, positiveStrand, anchorSeqSize2output) + "\t"
										+ indexInfo->returnFusionAnchorSeq(tmpRightReadTail_chrName, tmpRightReadTail_startPos, negativeStrand,
											tmpOriSAM_chrName, tmpOriSAM_endPos_2, negativeStrand, anchorSeqSize2output);	
								tmpRightReadTailBreakPoint_adjusted_mapRange 
									= int_to_str(tmpOriSAM_startPos_1) + "\t" + int_to_str(tmpRightReadTail_endPos);
								//cout << "end of getting tmpRightReadTailBreakPoint_adjusted !" << endl;					
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM);
							//cout << "end of getting tmpRightReadTailBreakPoint !" << endl;
						}
						else  // right read rev -- for, case 7, 8
						{
							int canonicalAndStrandIdenfiedBreakPoint_left;
							string canonicalAndStrandIdenfiedChrName_left;
							int canonicalAndStrandIdenfiedBreakPoint_right;
							string canonicalAndStrandIdenfiedChrName_right;
							string canonicalAndStrandIdenfiedStrand_left;
							string canonicalAndStrandIdenfiedStrand_right;
							int canonicalAndStrandIdenfiedAnchorLength_left;
							int canonicalAndStrandIdenfiedAnchorLength_right;
							int fusionCase;
							bool canonicalAndStrandIdenfiedOrNot 
								= adjustFusionBreakPoint_searchForCanonicalOnesWithOffset_allowVaryOffsetAtDiffEnd(
									false, false, true, tmpOriSAM_chrName, tmpRightReadTail_chrName,
									tmpOriSAM_endPos_2, tmpRightReadTail_endPos,
									tmpFusionAnchorLength_left, tmpFusionAnchorLength_right,
									canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedChrName_right,
									canonicalAndStrandIdenfiedBreakPoint_left, canonicalAndStrandIdenfiedBreakPoint_right,
									canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedStrand_right,
									canonicalAndStrandIdenfiedAnchorLength_left, canonicalAndStrandIdenfiedAnchorLength_right,
									indexInfo, offsetForAdjustingFusionBreakPoint_max, fusionCase, allowInsertedBaseOrNotBool);
							//cout << "canonicalAndStrandIdenfiedOrNot: " << canonicalAndStrandIdenfiedOrNot << endl;	
							if(canonicalAndStrandIdenfiedOrNot)
							{
								tmpRightReadTailBreakPoint_adjusted = tmpRightReadTailBreakPoint
									+ canonicalAndStrandIdenfiedChrName_left + "\t" + canonicalAndStrandIdenfiedChrName_right + "\t" 
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedBreakPoint_right) + "\t"
									+ canonicalAndStrandIdenfiedStrand_left + "\t" + canonicalAndStrandIdenfiedStrand_right + "\tGTAG\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_left) + "\t"
									+ int_to_str(canonicalAndStrandIdenfiedAnchorLength_right) + "\t"
									+ int_to_str(tmpXM) + "\t" + int_to_str(fusionCase) + ",";
								tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(
										canonicalAndStrandIdenfiedChrName_left, canonicalAndStrandIdenfiedBreakPoint_left,
										canonicalAndStrandIdenfiedStrand_left, canonicalAndStrandIdenfiedChrName_right, 
										canonicalAndStrandIdenfiedBreakPoint_right, canonicalAndStrandIdenfiedStrand_right, 
										anchorSeqSize2output);
								int tmpMapRange_most_1, tmpMapRange_most_2;
								incompleteUniquePairedAlignmentInfo.returnFusionSpanningReadMapRange_canonical(
									fusionCase, tmpRightReadTail_startPos, tmpRightReadTail_endPos,
									tmpMapRange_most_1, tmpMapRange_most_2);
								tmpRightReadTailBreakPoint_adjusted_mapRange 
									= int_to_str(tmpMapRange_most_1) + "\t" + int_to_str(tmpMapRange_most_2);
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
									+ "\t" + int_to_str(tmpFusionAnchorLength_right)
									+ "\t" + int_to_str(tmpXM) + "\t7,8,"; 	
								tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq
									= indexInfo->returnFusionAnchorSeq(tmpOriSAM_chrName, tmpOriSAM_endPos_2, positiveStrand,
										tmpRightReadTail_chrName, tmpRightReadTail_endPos, negativeStrand, anchorSeqSize2output) + "\t"
										+ indexInfo->returnFusionAnchorSeq(tmpRightReadTail_chrName, tmpRightReadTail_endPos, positiveStrand,
											tmpOriSAM_chrName, tmpOriSAM_endPos_2, negativeStrand, anchorSeqSize2output);	
								tmpRightReadTailBreakPoint_adjusted_mapRange 
									= int_to_str(tmpOriSAM_startPos_1) + "\t" + int_to_str(tmpRightReadTail_startPos);
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
								+ "\t" + int_to_str(tmpFusionAnchorLength_right)
								+ "\t" + int_to_str(tmpXM);							
						}
					}

					//cout << "start to insert breakPointStr " << endl;
			 		if((leftReadHeadUnfixed_fixed_bool)&&(rightReadTailUnfixed_fixed_bool))
					{
						outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted + "\n" + tmpRightReadTailBreakPoint_adjusted;
						outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint + "\n" + tmpRightReadTailBreakPoint;
						
						outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]	
							= tmpLeftReadHeadBreakPoint_adjusted + "\n" + inputSamStr1Vec[tmpOpenMP] + "\n"
								+ inputSamStr2Vec[tmpOpenMP] + "\n" + leftReadHeadUnfixed_fixed_SAMstr;
	
						if(leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool)
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_leftHead_for";
						else
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_leftHead_rev";

						outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]		
							= outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+ "\n" + tmpRightReadTailBreakPoint_adjusted + "\n" + inputSamStr1Vec[tmpOpenMP] + "\n"
								+ inputSamStr2Vec[tmpOpenMP] + "\n" + rightReadTailUnfixed_fixed_SAMstr;
						
						if(rightReadTailUnfixed_fixed_Nor_or_Rcm_bool)
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_rightTail_rev";
						else
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_rightTail_for";

						outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted + "\n" + tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq + "\n"
								+ tmpRightReadTailBreakPoint_adjusted + "\n" + tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq;

						outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted + "\t" + tmpLeftReadHeadBreakPoint_adjusted_mapRange + "\n"
								+ tmpRightReadTailBreakPoint_adjusted + "\t" + tmpRightReadTailBreakPoint_adjusted_mapRange;
					}
					else if(leftReadHeadUnfixed_fixed_bool)
					{
						outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted;
						outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint;

						outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]	
							= tmpLeftReadHeadBreakPoint_adjusted + "\n" + inputSamStr1Vec[tmpOpenMP] + "\n"
								+ inputSamStr2Vec[tmpOpenMP] + "\n" + leftReadHeadUnfixed_fixed_SAMstr;
						if(leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool)
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_leftHead_for";
						else
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_leftHead_rev";

						outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted + "\n" + tmpLeftReadHeadBreakPoint_adjusted_fusionAnchorSeq;
						
						outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmpOpenMP]
							= tmpLeftReadHeadBreakPoint_adjusted + "\t" + tmpLeftReadHeadBreakPoint_adjusted_mapRange;
					}
					else if(rightReadTailUnfixed_fixed_bool)
					{
						outputFusionSamStrVec_withBreakPoint_adjusted[tmpOpenMP]
							= tmpRightReadTailBreakPoint_adjusted;
						outputFusionSamStrVec_withBreakPoint[tmpOpenMP]
							= tmpRightReadTailBreakPoint;

						outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]		
							= tmpRightReadTailBreakPoint_adjusted + "\n" + inputSamStr1Vec[tmpOpenMP] + "\n"
								+ inputSamStr2Vec[tmpOpenMP] + "\n" + rightReadTailUnfixed_fixed_SAMstr;
						if(rightReadTailUnfixed_fixed_Nor_or_Rcm_bool)
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_rightTail_rev";
						else
							outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmpOpenMP]
								+= "\tZF:Z:FUS_rightTail_for";

						outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmpOpenMP]
							= tmpRightReadTailBreakPoint_adjusted + "\n" + tmpRightReadTailBreakPoint_adjusted_fusionAnchorSeq;

						outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmpOpenMP]
							= tmpRightReadTailBreakPoint_adjusted + "\t" + tmpRightReadTailBreakPoint_adjusted_mapRange;
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
	
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{	
			if(outputOriSamStrVec[tmp] != "")
				oriSAM_ofs << outputOriSamStrVec[tmp] << endl;		
			if(outputNonfusionSamStrVec[tmp] != "")
				nonFusionSAM_ofs << outputNonfusionSamStrVec[tmp] << endl;
			if(outputFusionSamStrVec[tmp] != "")
				fusionSAM_ofs << outputFusionSamStrVec[tmp] << endl;
			if(outputFusionSamStrVec_withBreakPoint[tmp] != "")
				fusionSAM_withBreakPoint_ofs << outputFusionSamStrVec_withBreakPoint[tmp] << endl;
			//if(outputFusionSamStrVec_withBreakPoint_adjusted[tmp] != "")
			//	fusionSAM_withBreakPoint_adjusted_ofs << outputFusionSamStrVec_withBreakPoint_adjusted[tmp] << endl;
			if(outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmp] != "")
				fusionSAM_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering_ofs
					<< outputFusionSamStrVec_withBreakPoint_adjusted_withFusionSAMforLocalIndexMapFiltering[tmp] << endl;
			if(outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmp] != "")
				fusionSAM_withBreakPoint_adjusted_withFusionAnchorSeq_ofs 
					<< outputFusionSamStrVec_withBreakPoint_adjusted_fusionAnchorSeq[tmp] << endl;
			if(outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmp] != "")
				fusionSAM_withBreakPoint_adjusted_withMapRegion_ofs
					<< outputFusionSamStrVec_withBreakPoint_adjusted_mapRange[tmp] << endl;
		}
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