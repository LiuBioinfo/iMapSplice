// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef UNIQUEUNPAIREDALIGNMENT_INFO_H
#define UNIQUEUNPAIREDALIGNMENT_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "fusionBreakPointHash_info.h"

using namespace std;

class UniqueUnpairedAlignment_Info
{
private:
	bool forOrRevBool_end1;
	string readName_1;
	int chrNameInt_end1;
	int startPos_1;
	vector<Jump_Code> cigarStringJumpCodeVec_1;
	int endPos_1;
	//int unfixedHeadLen_1;
	//int unfixedTailLen_1;

	bool forOrRevBool_end2;
	string readName_2;
	int chrNameInt_end2;
	int startPos_2;
	int breakPointPos_end2;
	vector<Jump_Code> cigarStringJumpCodeVec_2;
	int endPos_2;
	//int unfixedHeadLen_2;
	//int unfixedTailLen_2;
public:
	UniqueUnpairedAlignment_Info()
	{}

	bool pairedMappedOrNot(int flag)
	{
		if((flag & 0x1) && (flag & 0x2))
			return true;
		else
			return false;
	}

	string returnReadName_withoutSignForEnd()
	{
		string tmpOriReadName = readName_1;
		int tmpOriReadNameLength = tmpOriReadName.length();
		if((tmpOriReadName.substr(tmpOriReadNameLength-2,2) == "/1")
			||(tmpOriReadName.substr(tmpOriReadNameLength-2,2) == "/2"))
			return tmpOriReadName.substr(0, tmpOriReadNameLength-2);
		else
			return tmpOriReadName;
	}

	int searchFusionBreakPointEncompassedAndAddSupport(
		FusionBreakPointHash_Info* fusionBreakPointHashInfo)
	{
		vector<int> breakPointPosVec_gene1;
		vector<int> breakPointPosVec_gene2;
		vector<int> mapPosDistanceSumVec;
		vector<bool> end1_gene1_or_gene2_boolVec;

		// cout << "forOrRevBool_end1: " << forOrRevBool_end1 << endl;
		// cout << "forOrRevBool_end2: " << forOrRevBool_end2 << endl;
		// cout << "startPos_1: " << startPos_1 << endl;
		// cout << "endPos_1: " << endPos_1 << endl;
		// cout << "startPos_2: " << startPos_2 << endl;
		// cout << "endPos_2: " << endPos_2 << endl;

		if(forOrRevBool_end1 && forOrRevBool_end2) // case 9, +-, for...for 
		{
			bool tmp_end1_gene1_or_gene2_bool;
			int tmpChrNameInt_gene1;
			int tmpChrNameInt_gene2;
			int tmpLeftMostPos_breakPointSearch_gene1;
			int tmpRightMostPos_breakPointSearch_gene1;
			int tmpLeftMostPos_breakPointSearch_gene2;
			int tmpRightMostPos_breakPointSearch_gene2;
			string tmpSpecificFusionStrand = "+-";
			// end1 -- gene1, end2 -- gene2
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end1;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end2;
			tmp_end1_gene1_or_gene2_bool = true;
			tmpChrNameInt_gene1 = chrNameInt_end1;
			tmpChrNameInt_gene2 = chrNameInt_end2;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_1;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_1 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_2;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_2 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;

			// cout << "tmpLeftMostPos_breakPointSearch_gene1: " << tmpLeftMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene1: " << tmpRightMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpLeftMostPos_breakPointSearch_gene2: " << tmpLeftMostPos_breakPointSearch_gene2 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene2: " << tmpRightMostPos_breakPointSearch_gene2 << endl;
			// cout << "start to do searchFusionBreakPointWithinRangeForBothTwoGenes" << endl;
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand, 
				tmpFoundFusionBreakPointPosVec_gene1_end1, tmpFoundFusionBreakPointPosVec_gene2_end2);
			// cout << "end of doing searchFusionBreakPointWithinRangeForBothTwoGenes ...." << endl;
			// cout << "found tmpFoundFusionBreakPointPosVec_gene1_end1.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene1_end1.size() << endl; 
			// cout << "found tmpFoundFusionBreakPointPosVec_gene2_end2.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene2_end2.size() << endl;
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end1.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end1 = tmpFoundFusionBreakPointPosVec_gene1_end1[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end2 = tmpFoundFusionBreakPointPosVec_gene2_end2[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end1);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end2);
				mapPosDistanceSumVec.push_back(tmpFoundFusionBreakPointPos_gene1_end1 - endPos_1
					+ tmpFoundFusionBreakPointPos_gene2_end2 - endPos_2);
				end1_gene1_or_gene2_boolVec.push_back(true);
			}

			// end2 -- gene1, end1 -- gene2
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end2;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end1;			
			tmp_end1_gene1_or_gene2_bool = false;			
			tmpChrNameInt_gene1 = chrNameInt_end2;
			tmpChrNameInt_gene2 = chrNameInt_end1;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_2;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_2 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_1;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_1 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			// cout << "tmpLeftMostPos_breakPointSearch_gene1: " << tmpLeftMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene1: " << tmpRightMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpLeftMostPos_breakPointSearch_gene2: " << tmpLeftMostPos_breakPointSearch_gene2 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene2: " << tmpRightMostPos_breakPointSearch_gene2 << endl;
			// cout << "start to do searchFusionBreakPointWithinRangeForBothTwoGenes" << endl;
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand, 
				tmpFoundFusionBreakPointPosVec_gene1_end2, tmpFoundFusionBreakPointPosVec_gene2_end1);
			// cout << "end of doing searchFusionBreakPointWithinRangeForBothTwoGenes ...." << endl;
			// cout << "found tmpFoundFusionBreakPointPosVec_gene1_end1.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene1_end1.size() << endl; 
			// cout << "found tmpFoundFusionBreakPointPosVec_gene2_end2.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene2_end2.size() << endl;
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end2.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end2 = tmpFoundFusionBreakPointPosVec_gene1_end2[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end1 = tmpFoundFusionBreakPointPosVec_gene2_end1[tmp];				
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end2);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end1);
				mapPosDistanceSumVec.push_back(tmpFoundFusionBreakPointPos_gene1_end2 - endPos_2
					+ tmpFoundFusionBreakPointPos_gene2_end1 - endPos_1);
				end1_gene1_or_gene2_boolVec.push_back(false);
			}
		}
		else if((!forOrRevBool_end1) && (!forOrRevBool_end2)) // case 12, strand = "-+", REV...REV
		{
			bool tmp_end1_gene1_or_gene2_bool;
			int tmpChrNameInt_gene1;
			int tmpChrNameInt_gene2;
			int tmpLeftMostPos_breakPointSearch_gene1;
			int tmpRightMostPos_breakPointSearch_gene1;
			int tmpLeftMostPos_breakPointSearch_gene2;
			int tmpRightMostPos_breakPointSearch_gene2;
			string tmpSpecificFusionStrand = "-+";
			// end1 -- gene1, end2 -- gene2
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end1;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end2;
			tmp_end1_gene1_or_gene2_bool = true;
			tmpChrNameInt_gene1 = chrNameInt_end1;
			tmpChrNameInt_gene2 = chrNameInt_end2;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_1 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_1;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_2 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_2;
			// cout << "tmpLeftMostPos_breakPointSearch_gene1: " << tmpLeftMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene1: " << tmpRightMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpLeftMostPos_breakPointSearch_gene2: " << tmpLeftMostPos_breakPointSearch_gene2 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene2: " << tmpRightMostPos_breakPointSearch_gene2 << endl;
			// cout << "start to do searchFusionBreakPointWithinRangeForBothTwoGenes" << endl;
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end1, tmpFoundFusionBreakPointPosVec_gene2_end2);			
			// cout << "end of doing searchFusionBreakPointWithinRangeForBothTwoGenes ...." << endl;
			// cout << "found tmpFoundFusionBreakPointPosVec_gene1_end1.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene1_end1.size() << endl; 
			// cout << "found tmpFoundFusionBreakPointPosVec_gene2_end2.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene2_end2.size() << endl;
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end1.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end1 = tmpFoundFusionBreakPointPosVec_gene1_end1[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end2 = tmpFoundFusionBreakPointPosVec_gene2_end2[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end1);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end2);
				mapPosDistanceSumVec.push_back(startPos_1 - tmpFoundFusionBreakPointPos_gene1_end1
					+ startPos_2 - tmpFoundFusionBreakPointPos_gene2_end2);
				end1_gene1_or_gene2_boolVec.push_back(true);
			}			

			// end2 -- gene1, end1 -- gene2
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end2;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end1;
			tmp_end1_gene1_or_gene2_bool = true;
			tmpChrNameInt_gene1 = chrNameInt_end2;
			tmpChrNameInt_gene2 = chrNameInt_end1;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_2 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_2;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_1 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_1;
			// cout << "tmpLeftMostPos_breakPointSearch_gene1: " << tmpLeftMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene1: " << tmpRightMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpLeftMostPos_breakPointSearch_gene2: " << tmpLeftMostPos_breakPointSearch_gene2 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene2: " << tmpRightMostPos_breakPointSearch_gene2 << endl;
			// cout << "start to do searchFusionBreakPointWithinRangeForBothTwoGenes" << endl;
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end2, tmpFoundFusionBreakPointPosVec_gene2_end1);
			// cout << "end of doing searchFusionBreakPointWithinRangeForBothTwoGenes ...." << endl;
			// cout << "found tmpFoundFusionBreakPointPosVec_gene1_end1.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene1_end1.size() << endl; 
			// cout << "found tmpFoundFusionBreakPointPosVec_gene2_end2.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene2_end2.size() << endl;
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end2.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end2 = tmpFoundFusionBreakPointPosVec_gene1_end2[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end1 = tmpFoundFusionBreakPointPosVec_gene2_end1[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end2);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end1);
				mapPosDistanceSumVec.push_back(startPos_2 - tmpFoundFusionBreakPointPos_gene1_end2
					+ startPos_1 - tmpFoundFusionBreakPointPos_gene2_end1);
				end1_gene1_or_gene2_boolVec.push_back(false);
			}			
		}
		else if(forOrRevBool_end1 && (!forOrRevBool_end2)) // case 3, case 6, end1 -- gene1 -- for, end2 -- gene2 -- rcm
		{
			bool tmp_end1_gene1_or_gene2_bool;
			int tmpChrNameInt_gene1;
			int tmpChrNameInt_gene2;
			int tmpLeftMostPos_breakPointSearch_gene1;
			int tmpRightMostPos_breakPointSearch_gene1;
			int tmpLeftMostPos_breakPointSearch_gene2;
			int tmpRightMostPos_breakPointSearch_gene2;
			string tmpSpecificFusionStrand;
			// end1 -- gene1 -- for, end2 -- gene2 -- rev, case 3
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end1;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end2;
			tmp_end1_gene1_or_gene2_bool = true;
			tmpChrNameInt_gene1 = chrNameInt_end1;
			tmpChrNameInt_gene2 = chrNameInt_end2;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_1;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_1 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;			
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_2 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_2;
			// cout << "tmpLeftMostPos_breakPointSearch_gene1: " << tmpLeftMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene1: " << tmpRightMostPos_breakPointSearch_gene1 << endl;
			// cout << "tmpLeftMostPos_breakPointSearch_gene2: " << tmpLeftMostPos_breakPointSearch_gene2 << endl;
			// cout << "tmpRightMostPos_breakPointSearch_gene2: " << tmpRightMostPos_breakPointSearch_gene2 << endl;
			// cout << "start to do searchFusionBreakPointWithinRangeForBothTwoGenes" << endl;
			tmpSpecificFusionStrand = "++";
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end1, tmpFoundFusionBreakPointPosVec_gene2_end2);
			// cout << "end of doing searchFusionBreakPointWithinRangeForBothTwoGenes ...." << endl;
			// cout << "found tmpFoundFusionBreakPointPosVec_gene1_end1.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene1_end1.size() << endl; 
			// cout << "found tmpFoundFusionBreakPointPosVec_gene2_end2.size(): " 
			// 	<< tmpFoundFusionBreakPointPosVec_gene2_end2.size() << endl;
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end1.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end1 = tmpFoundFusionBreakPointPosVec_gene1_end1[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end2 = tmpFoundFusionBreakPointPosVec_gene2_end2[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end1);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end2);
				mapPosDistanceSumVec.push_back(tmpFoundFusionBreakPointPos_gene1_end1 - endPos_1
					+ startPos_2 - tmpFoundFusionBreakPointPos_gene2_end2);
				end1_gene1_or_gene2_boolVec.push_back(true);
			}			

			// end1 -- gene2 -- for, end2 -- gene1 -- rev, case 6
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end2;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end1;
			tmp_end1_gene1_or_gene2_bool = false;
			tmpChrNameInt_gene1 = chrNameInt_end2;
			tmpChrNameInt_gene2 = chrNameInt_end1;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_2 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_2;			
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_1;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_1 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpSpecificFusionStrand = "--";
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end2, tmpFoundFusionBreakPointPosVec_gene2_end1);
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end2.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end2 = tmpFoundFusionBreakPointPosVec_gene1_end2[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end1 = tmpFoundFusionBreakPointPosVec_gene2_end1[tmp];				
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end2);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end1);
				mapPosDistanceSumVec.push_back(startPos_2 - tmpFoundFusionBreakPointPos_gene1_end2
					+ tmpFoundFusionBreakPointPos_gene2_end1 - endPos_1);
				end1_gene1_or_gene2_boolVec.push_back(false);
			}					
		}
		else // ((!forOrRevBool_end1) && forOrRevBool_end2), case 3, case 6
		{
			// end1 -- rev, end2 -- for
			bool tmp_end1_gene1_or_gene2_bool;
			int tmpChrNameInt_gene1;
			int tmpChrNameInt_gene2;
			int tmpLeftMostPos_breakPointSearch_gene1;
			int tmpRightMostPos_breakPointSearch_gene1;
			int tmpLeftMostPos_breakPointSearch_gene2;
			int tmpRightMostPos_breakPointSearch_gene2;
			string tmpSpecificFusionStrand;
			// case 3, end1 -- gene2 -- rev, end2 -- gene1 -- for
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end2;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end1;
			tmp_end1_gene1_or_gene2_bool = false;
			tmpChrNameInt_gene1 = chrNameInt_end2;
			tmpChrNameInt_gene2 = chrNameInt_end1;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_2;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_2 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_1 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_1;
			tmpSpecificFusionStrand = "++";
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end2, tmpFoundFusionBreakPointPosVec_gene2_end1);
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end2.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end2 = tmpFoundFusionBreakPointPosVec_gene1_end2[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end1 = tmpFoundFusionBreakPointPosVec_gene2_end1[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end2);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end1);
				mapPosDistanceSumVec.push_back(tmpFoundFusionBreakPointPos_gene1_end2 - endPos_2
					+ startPos_1 - tmpFoundFusionBreakPointPos_gene2_end1);
				end1_gene1_or_gene2_boolVec.push_back(false);				
			}			

			// case 6, end1 -- gene1 -- rev, end2 -- gene2 -- for
			vector<int> tmpFoundFusionBreakPointPosVec_gene1_end1;
			vector<int> tmpFoundFusionBreakPointPosVec_gene2_end2;
			tmp_end1_gene1_or_gene2_bool = true;
			tmpChrNameInt_gene1 = chrNameInt_end1;
			tmpChrNameInt_gene2 = chrNameInt_end2;
			tmpLeftMostPos_breakPointSearch_gene1 = startPos_1 - FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpRightMostPos_breakPointSearch_gene1 = endPos_1;
			tmpLeftMostPos_breakPointSearch_gene2 = startPos_2;
			tmpRightMostPos_breakPointSearch_gene2 = endPos_2 + FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX;
			tmpSpecificFusionStrand = "--";
			fusionBreakPointHashInfo->searchFusionBreakPointWithinRangeForBothTwoGenes(
				tmpChrNameInt_gene1, tmpLeftMostPos_breakPointSearch_gene1, 
				tmpRightMostPos_breakPointSearch_gene1,
				tmpChrNameInt_gene2, tmpLeftMostPos_breakPointSearch_gene2, 
				tmpRightMostPos_breakPointSearch_gene2, tmpSpecificFusionStrand,
				tmpFoundFusionBreakPointPosVec_gene1_end1, tmpFoundFusionBreakPointPosVec_gene2_end2);
			for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1_end1.size(); tmp++)
			{
				int tmpFoundFusionBreakPointPos_gene1_end1 = tmpFoundFusionBreakPointPosVec_gene1_end1[tmp];
				int tmpFoundFusionBreakPointPos_gene2_end2 = tmpFoundFusionBreakPointPosVec_gene2_end2[tmp];
				breakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPos_gene1_end1);
				breakPointPosVec_gene2.push_back(tmpFoundFusionBreakPointPos_gene2_end2);
				mapPosDistanceSumVec.push_back(tmpFoundFusionBreakPointPos_gene2_end2 - endPos_2
					+ startPos_1 - tmpFoundFusionBreakPointPos_gene1_end1);
				end1_gene1_or_gene2_boolVec.push_back(true);
			}
		}

		// select the most confident fusionBreakPoint with the smallest value of mapPosDistanceSumVec;
		if(breakPointPosVec_gene1.size() == 0)
			return -1;

		int tmpIndex_mostConfidentFusionBreakPoint = -1;
		int tmpLeastMapPosDistance_mostConfidentFusionBreakPoint = 10000000;
		for(int tmp = 0; tmp < breakPointPosVec_gene1.size(); tmp++)
		{
			int tmpMapPosDistance = mapPosDistanceSumVec[tmp];
			if(tmpMapPosDistance < tmpLeastMapPosDistance_mostConfidentFusionBreakPoint)
			{
				tmpIndex_mostConfidentFusionBreakPoint = tmp;
				tmpLeastMapPosDistance_mostConfidentFusionBreakPoint = tmpMapPosDistance;
			}
		}

		int breakPointPos_gene1_final = breakPointPosVec_gene1[tmpIndex_mostConfidentFusionBreakPoint];
		int breakPointPos_gene2_final = breakPointPosVec_gene2[tmpIndex_mostConfidentFusionBreakPoint];
		bool end1_gene1_or_gene2_bool = end1_gene1_or_gene2_boolVec[tmpIndex_mostConfidentFusionBreakPoint];
		int chrNameInt_gene1_final, chrNameInt_gene2_final;
		if(end1_gene1_or_gene2_bool)
		{	
			chrNameInt_gene1_final = chrNameInt_end1;
			chrNameInt_gene2_final = chrNameInt_end2;
		}
		else
		{
			chrNameInt_gene1_final = chrNameInt_end2;
			chrNameInt_gene2_final = chrNameInt_end1;
		}

		int finalIndexInFusionBreakPointHashInfo 
			= fusionBreakPointHashInfo->searchAndReturnIndexInBreakPointInfoVec_for(
				chrNameInt_gene1_final, chrNameInt_gene2_final, 
				breakPointPos_gene1_final, breakPointPos_gene2_final);
		return finalIndexInFusionBreakPointHashInfo;
	}

	bool initiateWith2samStr(
		string& tmpSamStr_1, string& tmpSamStr_2, Index_Info* indexInfo)
	{
		vector<string> samFieldVec_1;
		vector<string> samFieldVec_2;
		int startLoc = 0;
		for(int tmp = 0; tmp < 13; tmp++)
		{
			int tabLoc = tmpSamStr_1.find("\t", startLoc);
			string tmpSamField = tmpSamStr_1.substr(startLoc, tabLoc-startLoc);
			samFieldVec_1.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_1.push_back(tmpSamStr_1.substr(startLoc));	
		startLoc = 0;
		for(int tmp = 0; tmp < 13; tmp++)
		{
			int tabLoc = tmpSamStr_2.find("\t", startLoc);
			string tmpSamField = tmpSamStr_2.substr(startLoc, tabLoc-startLoc);
			samFieldVec_2.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_2.push_back(tmpSamStr_2.substr(startLoc));	

		string chrNameStr_1 = samFieldVec_1[2];
		string chrNameStr_2 = samFieldVec_2[2];
		if((chrNameStr_1 == "*")||(chrNameStr_2 == "*"))
			return false;

		string IHfield_1 = samFieldVec_1[12];
		string IHfield_2 = samFieldVec_2[12];
		string IHintStr_1 = IHfield_1.substr(5);
		string IHintStr_2 = IHfield_2.substr(5);
		int IHint_1 = atoi(IHintStr_1.c_str());
		int IHint_2 = atoi(IHintStr_2.c_str());
		if(!((IHint_1 == 1)&&(IHint_2 == 1)))
			return false;
		string flagStr_1 = samFieldVec_1[1];
		string flagStr_2 = samFieldVec_2[1];
		int flagInt_1 = atoi(flagStr_1.c_str());
		int flagInt_2 = atoi(flagStr_2.c_str());
		bool pairedMappedOrNot_1 = this->pairedMappedOrNot(flagInt_1);
		bool pairedMappedOrNot_2 = this->pairedMappedOrNot(flagInt_2);
		if(pairedMappedOrNot_1 || pairedMappedOrNot_2)
			return false;

		if(flagInt_1 & 0x10)
			forOrRevBool_end1 = false;
		else
			forOrRevBool_end1 = true;

		if(flagInt_2 & 0x10)
			forOrRevBool_end2 = false;
		else
			forOrRevBool_end2 = true;		

		readName_1 = samFieldVec_1[0];
		readName_2 = samFieldVec_2[0];
		chrNameInt_end1 = indexInfo->convertStringToInt(chrNameStr_1);
		chrNameInt_end2 = indexInfo->convertStringToInt(chrNameStr_2);
		string startPosStr_1 = samFieldVec_1[3];
		startPos_1 = atoi(startPosStr_1.c_str());
		string startPosStr_2 = samFieldVec_2[3];
		startPos_2 = atoi(startPosStr_2.c_str());
		string cigarString_1 = samFieldVec_1[5];
		string cigarString_2 = samFieldVec_2[5];
		this->cigarString2jumpCodeVec(cigarString_1, cigarStringJumpCodeVec_1);
		this->cigarString2jumpCodeVec(cigarString_2, cigarStringJumpCodeVec_2);
		endPos_1 = this->getEndPosOfSpecificJumpCode(
			startPos_1, cigarStringJumpCodeVec_1, cigarStringJumpCodeVec_1.size()-1);
		endPos_2 = this->getEndPosOfSpecificJumpCode(
			startPos_2, cigarStringJumpCodeVec_2, cigarStringJumpCodeVec_2.size()-1);

		// if(cigarStringJumpCodeVec_1[0].type == "S")
		// 	unfixedHeadLen_1 = cigarStringJumpCodeVec_1[0].len;
		// else
		// 	unfixedHeadLen_1 = 0;

		// if(cigarStringJumpCodeVec_2[0].type == "S")
		// 	unfixedHeadLen_2 = cigarStringJumpCodeVec_2[0].len;
		// else
		// 	unfixedHeadLen_2 = 0; 

		return true;
	}


	int getEndPosOfSpecificJumpCode(int startPos, 
		vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, 
		vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNIDX";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	string jumpCodeVec2cigarString(vector<Jump_Code>& jumpCodeVec)
	{
		string tmpCigarString;
		//cout << "********" << "jumpCodeVecSize: " << jumpCodeVec.size() << endl;
		for(int tmp = 0; tmp < jumpCodeVec.size(); tmp++)
		{
			tmpCigarString += jumpCodeVec[tmp].toString();
		}
		return tmpCigarString;
	}	
};
#endif