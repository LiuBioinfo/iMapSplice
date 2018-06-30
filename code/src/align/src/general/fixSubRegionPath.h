// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXSUBREGIONPATH_INFO_H
#define FIXSUBREGIONPATH_INFO_H

#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <map>
#include "index_info.h"
//#include "fixDoubleAnchor_annotation_info.h"

using namespace std;

class SubSegGroup_Info
{
private:
	int mapRegionIndex;
	vector< pair< int, vector<int> > > subSegGroup; 
	//int startSegToBuildPath_index_inSubSegGroup; // index of the first seg that MPS3 use to build Path from
public:	
	SubSegGroup_Info()
	{}

	int returnMapRegionIndex()
	{
		return mapRegionIndex;
	}

	int returnSubSegGroupSize()
	{
		return subSegGroup.size();
	}

	unsigned int returnUnsignedMinusVal(unsigned int a, unsigned int b)
	{
		if(a > b)
			return (a-b);
		else
			return (b-a);
	}

	int chooseCandiSegEndPosClosest2CertainPos(Seg_Info* segInfo,
		int candiSegGroup_index_inSubSegGroup, unsigned int genomicPos, int& tmpGapCase, Index_Info* indexInfo)
	{
		int tmpGroupIndexInOriSegInfo = subSegGroup[candiSegGroup_index_inSubSegGroup].first;
		int tmpGroupSegLength = segInfo->returnSegmentLength(tmpGroupIndexInOriSegInfo);
		int tmpBestCandiSegPos_index = -1;
		unsigned int tmpClosestDistance = 1000000;
		int checkSegVecSize = (subSegGroup[candiSegGroup_index_inSubSegGroup].second).size();
		int genomicPos_chr = indexInfo->getChr(genomicPos);
		unsigned int tmpBestCandiSegPos_endPos = -1;

		for(int tmp = 0; tmp < checkSegVecSize; tmp++)
		{
			int tmpCandiMapSeg_index = (subSegGroup[candiSegGroup_index_inSubSegGroup].second)[tmp];
			unsigned int tmpCandiMapSeg_startPos = segInfo->returnSegmentMapPos(tmpGroupIndexInOriSegInfo, tmpCandiMapSeg_index);
			unsigned int tmpCandiMapSeg_endPos = tmpCandiMapSeg_startPos + tmpGroupSegLength - 1;

			unsigned int tmpDistance = returnUnsignedMinusVal(tmpCandiMapSeg_endPos, genomicPos);
			if(tmpDistance < tmpClosestDistance)
			{
				tmpClosestDistance = tmpDistance;
				tmpBestCandiSegPos_index = tmpCandiMapSeg_index;
				tmpBestCandiSegPos_endPos = tmpCandiMapSeg_endPos;
			}
		}
		if(tmpClosestDistance > MAX_SPLICE_DISTANCE_PHASE1)
		{
			return -1;
		}
		else
		{
			int tmpBestCandiSegPos_endPos_chr = indexInfo->getChr(tmpBestCandiSegPos_endPos);
			if(genomicPos_chr != tmpBestCandiSegPos_endPos_chr)
				return -1;
			if(tmpClosestDistance == 0)
			{
				tmpGapCase = FIX_MATCH;
			}
			else if(tmpClosestDistance <= MAX_INSERTION_LENGTH) // MAX_INSERTION_LENGTH == MAX_DELETION_LENGTH
			{
				if(genomicPos > tmpBestCandiSegPos_endPos)
					tmpGapCase = FIX_DELETION;
				else
					tmpGapCase = FIX_INSERTION;
			}
			else if(genomicPos > tmpBestCandiSegPos_endPos)
			{
				tmpGapCase = FIX_SPLICE;
			}
			else
			{
				return -1;
			}
			return tmpBestCandiSegPos_index;
		}
	}

	int chooseCandiSegStartPosClosest2CertainPos(Seg_Info* segInfo,
		int candiSegGroup_index_inSubSegGroup, unsigned int genomicPos, 
		int& tmpGapCase, Index_Info* indexInfo) 
	{
		int tmpGroupIndexInOriSegInfo = subSegGroup[candiSegGroup_index_inSubSegGroup].first;
		//int tmpGroupSegLength = segInfo->returnSegmentLength(tmpGroupIndexInOriSegInfo);
		int tmpBestCandiSegPos_index_oriSegInfo = -1;
		unsigned int tmpClosestDistance = 1000000;
		int checkSegVecSize = (subSegGroup[candiSegGroup_index_inSubSegGroup].second).size();
		
		int genomicPos_chr = indexInfo->getChr(genomicPos);
		unsigned int tmpBestCandiSegPos_startPos = -1;

		for(int tmp = 0; tmp < checkSegVecSize; tmp++)
		{
			int tmpCandiMapSeg_index = (subSegGroup[candiSegGroup_index_inSubSegGroup].second)[tmp];
			unsigned int tmpCandiMapSeg_startPos = segInfo->returnSegmentMapPos(tmpGroupIndexInOriSegInfo, tmpCandiMapSeg_index);
			//unsigned int tmpCandiMapSeg_endPos = tmpCandiMapSeg_startPos + tmpGroupSegLength - 1;
			unsigned int tmpDistance = returnUnsignedMinusVal(tmpCandiMapSeg_startPos, genomicPos);
			if(tmpDistance < tmpClosestDistance)
			{
				tmpClosestDistance = tmpDistance;
				tmpBestCandiSegPos_index_oriSegInfo = tmpCandiMapSeg_index;
				tmpBestCandiSegPos_startPos = tmpCandiMapSeg_startPos;
			}
		}
		if(tmpClosestDistance > MAX_SPLICE_DISTANCE_PHASE1)
			return -1;
		else
		{	
			int tmpBestCandiSegPos_startPos_char = indexInfo->getChr(tmpBestCandiSegPos_startPos);
			if(genomicPos_chr != tmpBestCandiSegPos_endPos_chr)
				return -1;
			if(tmpClosestDistance == 0)
			{
				tmpGapCase = FIX_MATCH;
			}
			else if(tmpClosestDistance <= MAX_INSERTION_LENGTH)
			{
				if(genomicPos < tmpBestCandiSegPos_startPos)
					tmpGapCase = FIX_DELETION;
				else
					tmpGapCase = FIX_INSERTION;				
			}
			else if(genomicPos < tmpBestCandiSegPos_startPos)
			{
				tmpCapCase = FIX_SPLICE
			}
			else
			{
				return-1;
			}
			return tmpBestCandiSegPos_index_oriSegInfo;		
		}
	}

	int longestSeg_index_inSubSegGroup(Seg_Info* segInfo)
	{
		int tmp_seg_length_max = 10000000;
		int tmp_start_seg_index = 10000000;
		for(int tmp = 0; tmp < subSegGroup.size(); tmp ++)
		{
			int tmpGroupIndex = subSegGroup[tmp].first;
			int tmpSegGroupLength = segInfo.returnSegmentLength(tmpGroupIndex);
			if(tmpSegGroupLength < tmp_seg_length_max)
			{
				tmp_seg_length_max = tmpSegGroupLength;
				tmp_start_seg_index = tmp;
			}
		}		
		return tmp_start_seg_index;
	}	

	int returnSegGroupNO_inOriSegInfo(int tmpIndex)
	{
		return subSegGroup[tmpIndex].first;
	}

	int returnCandiMapSegNO_inOriSegInfo(int tmpIndex_1, int tmpIndex_2)
	{
		return (subSegGroup[tmpIndex_1].second)[tmpIndex_2];
	}

	int returnLastSegGroupNO()
	{
		int currentSegGroupSize = subSegGroup.size();
		return this->returnSegGroupNO_inOriSegInfo(currentSegGroupSize-1);
	}

	void pushBackNewSeg(int segGroupNO, int segCandiNO)
	{
		int lastSegGroupNO = this->returnLastSegGroupNO();
		if(segGroupNO == lastSegGroupNO)
		{
			int tmpSubSegGroupSize = subSegGroup.size();
			((subSegGroup[tmpSubSegGroupSize-1]).second).push_back(segCandiNO);
		}
		else if(segGroupNO > lastSegGroupNO)
		{
			vector<int> newSegCandiGroup;
			newSegCandiGroup.push_back(segCandiNO);
			subSegGroup.push_back(pair<int, vector<int> > (segGroupNO, newSegCandiGroup)); 
		}
		else
		{
			cout << " error in segGroupNO (> lastSegGroupNO)" << endl;
			exit(1);
		}
	}

	void initiateWith1stSeg(int segGroupNO, int segCandiNO)
	{
		vector<int> newSegCandiGroup;
		newSegCandiGroup.push_back(segCandiNO);
		subSegGroup.push_back(segCroupNO, newSegCandiGroup);
	}
};

class SubSegGroupPair_Info
{
private:
	int mapRegionIndex_min;
	int mapRegionIndex_max;
	bool strand;
	
	SubSegGroup_Info subSegGroupInfo_Nor;
	SubSegGroup_Info subSegGroupInfo_Rcm;

	//CandiPath_Info candiPathInfo_Nor;
	//CandiPath_Info candiPathInfo_Rcm;

public:
	SubSegGroupPair_Info(bool for_rev_bool, int tmpMapRegionIndex)
	{
		strand = for_rev_bool;
		mapRegionIndex_min = tmpMapRegionIndex;
		mapRegionIndex_max = tmpMapRegionIndex;
	}

	void pushBackNewSeg_NorOrRcm(int segGroupNO, int segCandiNO, bool Nor_or_Rcm_bool)
	{
		if(Nor_or_Rcm_bool)
			subSegGroupInfo_Nor.pushBackNewSeg(segGroupNO, segCandiNO)
		else
			subSegGroupInfo_Rcm.pushBackNewSeg(segGroupNO, segCandiNO)

	}

	void initiateWith1stSeg_NorOrRcm(int segGroupNO, int segCandiNO, bool Nor_or_Rcm_bool)
	{
		if(Nor_or_Rcm_bool)
			subSegGroupInfo_Nor.initiateWith1stSeg(segGroupNO, segCandiNO);
		else
			subSegGroupInfo_Rcm.initiateWith1stSeg(segGroupNO, segCandiNO);
	}

	int returnSubSegGroupSize_NorOrRcm(bool Nor_or_Rcm_bool)
	{
		if(Nor_or_Rcm_bool)
			return subSegGroupInfo_Nor.returnSubSegGroupSize();
		else
			return subSegGroupInfo_Rcm.returnSubSegGroupSize();
	}

	int returnSegGroupNO_inOriSegInfo_NorOrRcm(int tmpIndex, bool Nor_or_Rcm_bool)
	{
		if(Nor_or_Rcm_bool)
			return subSegGroupInfo_Nor.returnSegGroupNO_inOriSegInfo(tmpIndex);
		else
			return subSegGroupInfo_Rcm.returnSegGroupNO_inOriSegInfo(tmpIndex);
	}

	int chooseCandiSegStartPosClosest2CertainPos_NorOrRcm(Seg_Info* segInfo,
		int candiSegGroup_index_inSubSegGroup, unsigned int genomicPos, 
		int& tmpGapCase, Index_Info* indexInfo, bool Nor_or_Rcm_bool) 
	{
		if(Nor_or_Rcm_bool)
			subSegGroupInfo_Nor.chooseCandiSegStartPosClosest2CertainPos();
		else
			subSegGroupInfo_Rcm.chooseCandiSegStartPosClosest2CertainPos();
	}	
};	



/*
class CandiPath_Info
{
private:
	int chrInt;
	int mapStartPos;
	//int unfixedHeadLen;
	//int unfixedTailLen;
	vector< pair<int, int> > finalCandiPath; // < pair <segGroupNO, segCandiNO > >
	vector< int > finalCandiPath_gapCase;

	vector< pair<bool, vector<Jump_Code> > > finalCandiPathJumpCodeVec;
	vector< int > mismatchPosVec;
	vector< char > mismatchCharVec; 

	int LengthOfSeqPerMismatchAllowed;
public:

	CandiPath_Info()
	{
		LengthOfSeqPerMismatchAllowed = 8;
	}

	void generateCandiPath_uniqueOnly(
		SubSegGroup_Info* subSegGroupInfo, Seg_Info* segInfo, 
		Index_Info* indexInfo)
	{
		int firstSeg_groupNO = subSegGroupInfo->returnSegGroupNO_inOriSegInfo(0);
		int firstSeg_candiMapSegNO = subSegGroupIndo->returnCandiMapSegNO_inOriSegInfo(0, 0);
		finalCandiPath.push_back(pair<int,int> (firstSeg_groupNO, firstSeg_candiMapSegNO));

		unsigned int currentPathEndPos = segInfo->returnSegmentMapPos_end(firstSeg_groupNO, firstSeg_candiMapSegNO);
		int currentPathEndPos_inRead = segInfo->returnSegmentLocInRead(firstSeg_groupNO) + segInfo->returnSegmentLength(firstSeg_groupNO) - 1;

		int subSegGroupSize = subSegGroupInfo->returnSubSegGroupSize()
		for(int tmp = 1; tmp < subSegGroupSize; tmp++)
		{
			int tmpSegGroupNO_inOriSegInfo = subSegGroupInfo->returnSegGroupNO_inOriSegInfo(tmp);
			int tmpSegStartPos_inRead = segInfo->returnSegmentLocInRead(tmpSegGroupNO_inOriSegInfo);
			
			int tmpGapCase;
			int tmpBestCandiMapSeg_index_inOriSegInfo 
				= subSegGroupInfo->chooseCandiSegStartPosClosest2CertainPos(segInfo, tmp, 
					currentPathEndPos + tmpSegStartPos_inRead - currentPathEndPos_inRead, tmpGapCase);
			if(tmpBestCandiMapSeg_index_inOriSegInfo == -1)
			{}
			else
			{
				currentPathEndPos = segInfo->returnSegmentMapPos_end(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo);
				finalCandiPath.push_back(pair<int,int>(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo));
				finalCandiPath_gapCase.push_back(tmpGapCase);
			}
		}
		
		unsigned int mapStartPosInWholeGenome = segInfo->returnSegmentMapPos(firstSeg_groupNO, firstSeg_candiMapSegNO);
		chrInt = indexInfo->getChr(mapStartPosInWholeGenome);
		mapStartPos = returnChromPosWithWholeGenomePosAndChrInt(mapStartPosInWholeGenome, chrInt);
		//unfixedHeadLen = segInfo->returnSegmentLocInRead(firstSeg_groupNO)-1;
		//unfixedTailLen = segInfo->
	}

	void generateCandiPath_uniqueOnly_2TargetPathVec_seg(
		SubSegGroup_Info* subSegGroupInfo, Seg_Info* segInfo, 
		Index_Info* indexInfo,
		vector< pair<int, int> >& targetPathVec_seg,
		vector<int>& targetPathVec_seg_gapCase)
	{
		int firstSeg_groupNO = subSegGroupInfo->returnSegGroupNO_inOriSegInfo(0);
		int firstSeg_candiMapSegNO = subSegGroupIndo->returnCandiMapSegNO_inOriSegInfo(0, 0);
		finalCandiPath.push_back(pair<int,int> (firstSeg_groupNO, firstSeg_candiMapSegNO));

		unsigned int currentPathEndPos = segInfo->returnSegmentMapPos_end(firstSeg_groupNO, firstSeg_candiMapSegNO);
		int currentPathEndPos_inRead = segInfo->returnSegmentLocInRead(firstSeg_groupNO) + segInfo->returnSegmentLength(firstSeg_groupNO) - 1;

		int subSegGroupSize = subSegGroupInfo->returnSubSegGroupSize()
		for(int tmp = 1; tmp < subSegGroupSize; tmp++)
		{
			int tmpSegGroupNO_inOriSegInfo = subSegGroupInfo->returnSegGroupNO_inOriSegInfo(tmp);
			int tmpSegStartPos_inRead = segInfo->returnSegmentLocInRead(tmpSegGroupNO_inOriSegInfo);
			
			int tmpGapCase;
			int tmpBestCandiMapSeg_index_inOriSegInfo 
				= subSegGroupInfo->chooseCandiSegStartPosClosest2CertainPos(segInfo, tmp, 
					currentPathEndPos + tmpSegStartPos_inRead - currentPathEndPos_inRead, tmpGapCase);
			if(tmpBestCandiMapSeg_index_inOriSegInfo == -1)
			{}
			else
			{
				currentPathEndPos = segInfo->returnSegmentMapPos_end(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo);
				targetPathVec_seg.push_back(pair<int,int>(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo));
				targetPathVec_seg_gapCase.push_back(tmpGapCase);
			}
		}
	}

	void fixGapsInFinalCandiPath(Seg_Info* segInfo, Index_Info* indexInfo, 
		const string& readSeqInProcessing)
	{
		for(int tmp = 0; tmp < finalCandiPath.size()-1; tmp++)
		{
			vector<int> tmpGapMismatchPosVec;
			vector<char> tmpGapMismatchCharVec;
			vector<Jump_Code> tmpGapJumpCodeVec;
			int tmpGapCase = finalCandiPath_gapCase[tmp];

			int tmpSegGroupNO_1 = finalCandiPath[tmp].first;
			int tmpSegGroupNO_2 = finalCandiPath[tmp+1].first;
			int tmpCandiMapSegNO_1 = finalCandiPath[tmp].second;
			int tmpCandiMapSegNO_2 = finalCandiPath[tmp+1].second;
			int segLocInRead_1 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_1);
			int segLocInRead_2 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_2);
			int segmentLength_1 = segInfo->returnSegmentLength(tmpSegGroupNO_1);
			int segmentLength_2 = segInfo->returnSegmentLength(tmpSegGroupNO_2);
			unsigned int segPosInWholeGenome_1 = segInfo->returnSegmentMapPos(tmpSegGroupNO_1, tmpCandiMapSegNO_1);
			unsigned int segPosInWholeGenome_2 = segInfo->returnSegmentMapPos(tmpSegGroupNO_2, tmpCandiMapSegNO_2);
			int segPosInChr_1 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_1, chrInt);
			int segPosInChr_2 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_2, chrInt);

			bool fixDoubleAnchorGap_bool = this->fixDoubleAnchorGap(tmpGapMismatchPosVec, tmpGapMismatchCharVec,
				tmpGapJumpCodeVec, tmpGapCase, chrInt, segLocInRead_1, segLocInRead_2, 
				segmentLength_1, segmentLength_2, segPosInChr_1, segPosInChr_2,
				indexInfo, readSeqInProcessing);
			if(fixDoubleAnchorGap_bool)
			{
				finalCandiPathJumpCodeVec.push_back(pair<bool, vector<Jump_Code> >(true, tmpGapJumpCodeVec));
				this->pushBackMismatchPosVec();
				this->pushBackMismatchCharVec();
			}
			else
			{
				finalCandiPathJumpCodeVec.push_back(pair<bool, vector<Jump_Code> >(false, tmpGapJumpCodeVec));
			}
		}		
	}

	void fixGapsInFinalCandiPath2PathVec_seg(
		Seg_Info* segInfo, Index_Info* indexInfo, 
		const string& readSeqInProcessing)
	{
		for(int tmp = 0; tmp < finalCandiPath.size()-1; tmp++)
		{
			vector<int> tmpGapMismatchPosVec;
			vector<char> tmpGapMismatchCharVec;
			vector<Jump_Code> tmpGapJumpCodeVec;
			int tmpGapCase = finalCandiPath_gapCase[tmp];

			int tmpSegGroupNO_1 = finalCandiPath[tmp].first;
			int tmpSegGroupNO_2 = finalCandiPath[tmp+1].first;
			int tmpCandiMapSegNO_1 = finalCandiPath[tmp].second;
			int tmpCandiMapSegNO_2 = finalCandiPath[tmp+1].second;
			int segLocInRead_1 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_1);
			int segLocInRead_2 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_2);
			int segmentLength_1 = segInfo->returnSegmentLength(tmpSegGroupNO_1);
			int segmentLength_2 = segInfo->returnSegmentLength(tmpSegGroupNO_2);
			unsigned int segPosInWholeGenome_1 = segInfo->returnSegmentMapPos(tmpSegGroupNO_1, tmpCandiMapSegNO_1);
			unsigned int segPosInWholeGenome_2 = segInfo->returnSegmentMapPos(tmpSegGroupNO_2, tmpCandiMapSegNO_2);
			int segPosInChr_1 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_1, chrInt);
			int segPosInChr_2 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_2, chrInt);

			bool fixDoubleAnchorGap_bool = this->fixDoubleAnchorGap(tmpGapMismatchPosVec, tmpGapMismatchCharVec,
				tmpGapJumpCodeVec, tmpGapCase, chrInt, segLocInRead_1, segLocInRead_2, 
				segmentLength_1, segmentLength_2, segPosInChr_1, segPosInChr_2,
				indexInfo, readSeqInProcessing);
			if(fixDoubleAnchorGap_bool)
			{
				finalCandiPathJumpCodeVec.push_back(pair<bool, vector<Jump_Code> >(true, tmpGapJumpCodeVec));
				this->pushBackMismatchPosVec();
				this->pushBackMismatchCharVec();
			}
			else
			{
				finalCandiPathJumpCodeVec.push_back(pair<bool, vector<Jump_Code> >(false, tmpGapJumpCodeVec));
			}
		}		
	}	

	void fixDoubleAnchorGap(vector<int>& tmpMismatchPosVec, vector<char>& tmpMismatchCharVec, 
		vector<Jump_Code>& tmpJumpCodeVec,
		int tmpGapCase, int tmpChrInt, 
		int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2,
		int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		bool fixDoubleAnchor_bool = false;
		int extendBackNumMax = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		if(extendBackNumMax > segmentMapPos_2 - 1)
		{
			extendBackNumMax = segmentMapPos_2 - 1; 
		}
		int extendBackNum = this->extendBack(segmentLocInRead_2, readSeq_inProcess, 
			segmentMapPos_2, tmpChrInt, extendBackNumMax);
		
		segmentLocInRead_2 = segmentLocInRead_2 - extendBackNum;
		segmentLength_2 = segmentLength_2 + extendBackNum;
		segmentMapPos_2 = segmentMapPos_2 - extendBackNum;

		if(tmpGapCase = FIX_MATCH)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorMatch(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_INSERTION)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorInsertion(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_DELETION)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorDeletion(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_SPLICE)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorSplice(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else
		{
			cout << "error in tmpGapCase ..." << endl;
			exit(1);
		}
		return fixDoubleAnchor_bool;
	}

	bool fixDoubleAnchorInsertion(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{

	}

	bool fixDoubleAnchorMatch(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		int subSeq_toProcess_len = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		string readSubSeq_toProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeq_toProcess_len); 
		string chromSubSeq_toProcess = indexInfo->returnChromStrSubstr(chrInt, segmentMapPos_1 + segmentLength_1,
			subSeq_toProcess_len);

		int max_mismatch = subSeq_toProcess_len / LengthOfSeqPerMismatchAllowed + 1;

		FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
		bool scoreStringBool = fixMatchInfo->fixMatch(readSubSeqInProcess, chromSubSeqInProcess,
			max_mismatch, segmentLocInRead_1 + segmentLength_1);
		if(scoreStringBool)
		{
			Jump_Code matchJumpCode(//segmentLength_1 + 
				subSeq_toProcess_len + segmentLength_2, "M");

			fixMatchInfo->insertMismatchPosVec_fixMatchInfo(fixtmpGapMismatchPosVec);
			fixMatchInfo->insertMismatchCharVec_fixMatchInfo(fixMatchInfo, tmpGapMismatchCharVec);		
		}	
		delete fixMatchInfo;
		return fixDoubleAnchor_bool;
	}

	bool fixMatch_MismatchNum_MismatchPos_MismatchChar(const string& s1, 
		const string& s2, 
		int max_mismatch, int startPosInRead, vector<int>& tmpMismatchPosVec, 
		vector<char>& tmpMismatchCharVec)
	{
		bool match_bool = true;

		for(int tmp = 0; tmp < s1.length(); tmp++)
		{
			if(s1[tmp] != s2[tmp])
			{	
				int tmpMismatchPos = startPosInRead + tmp;
				mismatchPosInReadVec.push_back(tmpMismatchPos);
				mismatchCharVec.push_back(s2[tmp]);
				mismatchNum++;
				if(mismatchNum > max_mismatch)
				{
					match_bool = false; // not match
					break;
				}
			}
		}
		//num_mismatch = mismatchNum;
		return match_bool;
	}	

	void generateCandiPath_exhaust(SubSegGroup_Info* subSegGroupInfo, Seg_Info* segInfo)
	{
		//int segGroup_start2buildPath = subSegGroupInfo->longestSeg_index_inSubSegGroup(segInfo);
		//int candiMapSeg_start2buildPath = 0;  // fix me: there should be multiple choices when subSegGroup[segGroup_start2buildPath].size() > 1
											 // fix me: just for testing here with setting = 0 to check how much time can be saved.
		for(int tmpIndex_subSegGroup = 0; tmpIndex_subSegGroup < subSegGroup.size(); tmpIndex_subSegGroup)
		{
			int tmpSegGroupNO_inOriSegInfo = subSegGroupInfo->returnSegGroupNO_inSubSegGroup(tmpIndex_subSegGroup);


		}
	}
};
*/



class CandiPathPair_Info
{
private:
	int chrInt;

	//////////////   Nor   ////////////////////
	vector< pair<int, int> > finalCandiPath_Nor; // < pair <segGroupNO, segCandiNO > >
	vector< int > finalCandiPath_gapCase_Nor;

	vector< pair<bool, vector<Jump_Code> > > finalCandiPathJumpCodeVec_Nor;
	vector< vector<int> > mismatchPosVec_Nor;
	vector< vector<char> > mismatchCharVec_Nor; 

	//////////////   Rcm    ///////////////////
	vector< pair<int, int> > finalCandiPath_Rcm; // < pair <segGroupNO, segCandiNO > >
	vector< int > finalCandiPath_gapCase_Rcm;

	vector< pair<bool, vector<Jump_Code> > > finalCandiPathJumpCodeVec_Rcm;
	vector< vector<int> > mismatchPosVec_Rcm;
	vector< vector<char> > mismatchCharVec_Rcm; 


	int LengthOfSeqPerMismatchAllowed;	

public:
    CandiPathPair_Info()
    {
    	LengthOfSeqPerMismatchAllowed = 8;
    }


	bool fixGapsInPathPair(
		Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm, 
		Index_Info* indexInfo, bool Nor_or_Rcm_bool)
	{
		bool fixGapsInPathPair_bool = true;
		for(int tmp = 0; tmp < finalCandiPath_Nor.size()-1; tmp++)
		{
			vector<int> tmpGapMismatchPosVec;
			vector<char> tmpGapMismatchCharVec;
			vector<Jump_Code> tmpGapJumpCodeVec;
			int tmpGapCase = finalCandiPath_gapCase_Nor[tmp];

			int tmpSegGroupNO_1 = finalCandiPath_Nor[tmp].first;
			int tmpSegGroupNO_2 = finalCandiPath_Nor[tmp+1].first;
			int tmpCandiMapSegNO_1 = finalCandiPath_Nor[tmp].second;
			int tmpCandiMapSegNO_2 = finalCandiPath_Nor[tmp+1].second;
			int segLocInRead_1 = segInfo_Nor->returnSegmentLocInRead(tmpSegGroupNO_1);
			int segLocInRead_2 = segInfo_Nor->returnSegmentLocInRead(tmpSegGroupNO_2);
			int segmentLength_1 = segInfo_Nor->returnSegmentLength(tmpSegGroupNO_1);
			int segmentLength_2 = segInfo_Nor->returnSegmentLength(tmpSegGroupNO_2);
			unsigned int segPosInWholeGenome_1 = segInfo_Nor->returnSegmentMapPos(tmpSegGroupNO_1, tmpCandiMapSegNO_1);
			unsigned int segPosInWholeGenome_2 = segInfo_Nor->returnSegmentMapPos(tmpSegGroupNO_2, tmpCandiMapSegNO_2);
			int segPosInChr_1 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_1, chrInt);
			int segPosInChr_2 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_2, chrInt);

			bool fixDoubleAnchorGap_bool = this->fixDoubleAnchorGap(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec,
				tmpGapJumpCodeVec, tmpGapCase, chrInt, segLocInRead_1, segLocInRead_2, 
				segmentLength_1, segmentLength_2, segPosInChr_1, segPosInChr_2,
				indexInfo, readSeqInProcessing);
			if(fixDoubleAnchorGap_bool)
			{
				finalCandiPathJumpCodeVec_Nor.push_back(pair<bool, vector<Jump_Code> >(true, tmpGapJumpCodeVec));
			}
			else
			{
				finalCandiPathJumpCodeVec_Nor.push_back(pair<bool, vector<Jump_Code> >(false, tmpGapJumpCodeVec));				
				fixGapsInPathPair_bool = false;
			}
			mismatchPosVec_Nor.push_back(tmpGapMismatchPosVec);
			mismatchCharVec_Nor.push_back(tmpGapMismatchCharVec);
		}	

		for(int tmp = 0; tmp < finalCandiPath_Rcm.size()-1; tmp++) // gap # == segGroup # -1 
		{
			vector<int> tmpGapMismatchPosVec;
			vector<char> tmpGapMismatchCharVec;
			vector<Jump_Code> tmpGapJumpCodeVec;
			int tmpGapCase = finalCandiPath_gapCase_Rcm[tmp];

			int tmpSegGroupNO_1 = finalCandiPath_Rcm[tmp].first;
			int tmpSegGroupNO_2 = finalCandiPath_Rcmr[tmp+1].first;
			int tmpCandiMapSegNO_1 = finalCandiPath_Rcm[tmp].second;
			int tmpCandiMapSegNO_2 = finalCandiPath_Rcm[tmp+1].second;
			int segLocInRead_1 = segInfo_Rcm->returnSegmentLocInRead(tmpSegGroupNO_1);
			int segLocInRead_2 = segInfo_Rcm->returnSegmentLocInRead(tmpSegGroupNO_2);
			int segmentLength_1 = segInfo_Rcm->returnSegmentLength(tmpSegGroupNO_1);
			int segmentLength_2 = segInfo_Rcm->returnSegmentLength(tmpSegGroupNO_2);
			unsigned int segPosInWholeGenome_1 = segInfo_Rcm->returnSegmentMapPos(tmpSegGroupNO_1, tmpCandiMapSegNO_1);
			unsigned int segPosInWholeGenome_2 = segInfo_Rcm->returnSegmentMapPos(tmpSegGroupNO_2, tmpCandiMapSegNO_2);
			int segPosInChr_1 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_1, chrInt);
			int segPosInChr_2 = indexInfo->returnChromPosWithWholeGenomePosAndChrInt(segPosInWholeGenome_2, chrInt);

			bool fixDoubleAnchorGap_bool = this->fixDoubleAnchorGap(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec,
				tmpGapJumpCodeVec, tmpGapCase, chrInt, segLocInRead_1, segLocInRead_2, 
				segmentLength_1, segmentLength_2, segPosInChr_1, segPosInChr_2,
				indexInfo, readSeqInProcessing);
			if(fixDoubleAnchorGap_bool)
			{
				finalCandiPathJumpCodeVec_Rcm.push_back(pair<bool, vector<Jump_Code> >(true, tmpGapJumpCodeVec));
			}
			else
			{
				finalCandiPathJumpCodeVec_Rcm.push_back(pair<bool, vector<Jump_Code> >(false, tmpGapJumpCodeVec));				
				fixGapsInPathPair_bool = false;
			}
			mismatchPosVec_Rcm.push_back(tmpGapMismatchPosVec);
			mismatchCharVec_Rcm.push_back(tmpGapMismatchCharVec);
		}		
		return fixGapsInPathPair_bool;
	}

	void generateCandiPath_uniqueOnly_pair(
		SubSegGroupPair_Info* subSegGroupPairInfo, 
		Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm, 
		Index_Info* indexInfo)
	{
		// get Nor alignments 

		int firstSeg_groupNO = subSegGroupPairInfo->returnSegGroupNO_inOriSegInfo(0);
		int firstSeg_candiMapSegNO = subSegGroupPairInfo->returnCandiMapSegNO_inOriSegInfo(0, 0);
		finalCandiPath_Nor.push_back(pair<int,int> (firstSeg_groupNO_Nor, firstSeg_candiMapSegNO_Nor));

		unsigned int currentPathEndPos = segInfo_Nor->returnSegmentMapPos_end(firstSeg_groupNO, firstSeg_candiMapSegNO);
		int currentPathEndPos_inRead = segInfo_Nor->returnSegmentLocInRead(firstSeg_groupNO) 
			+ segInfo_Nor->returnSegmentLength(firstSeg_groupNO) - 1;

		int subSegGroupSize = subSegGroupPairInfo->returnSubSegGroupSize_NorOrRcm(true);
		for(int tmp = 1; tmp < subSegGroupSize; tmp++)
		{
			int tmpSegGroupNO_inOriSegInfo 
				= subSegGroupPairInfo->returnSegGroupNO_inOriSegInfo_NorOrRcm(tmp, true); // (index_subSegGroup, nor_or_rcm_bool)
			int tmpSegStartPos_inRead = segInfo_Nor->returnSegmentLocInRead(tmpSegGroupNO_inOriSegInfo);
			
			int tmpGapCase;
			int tmpBestCandiMapSeg_index_inOriSegInfo 
				= subSegGroupPairInfo->chooseCandiSegStartPosClosest2CertainPos_NorOrRcm(
					segInfo_Nor, tmp, 
					currentPathEndPos + tmpSegStartPos_inRead - currentPathEndPos_inRead, 
					tmpGapCase, true);
			if(tmpBestCandiMapSeg_index_inOriSegInfo == -1)
			{}
			else
			{
				currentPathEndPos = segInfo_Nor->returnSegmentMapPos_end(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo);
				finalCandiPath_Nor.push_back(pair<int,int>(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo));
				finalCandiPath_gapCase_Nor.push_back(tmpGapCase);
			}
		}

		// get Rcm alignments 
		subSegGroupSize = = subSegGroupPairInfo->returnSubSegGroupSize_NorOrRcm(false);
		for(int tmp = 0; tmp < subSegGroupSize; tmp++)
		{
			int tmpSegGroupNO_inOriSegInfo 
				= subSegGroupPairInfo->returnSegGroupNO_inOriSegInfo_NorOrRcm(tmp, false);
			int tmpSegStartPos_inRead = segInfo_Rcm->returnSegmentLocInRead(tmpSegGroupNO_inOriSegInfo);
			
			int tmpGapCase;
			int tmpBestCandiMapSeg_index_inOriSegInfo 
				= subSegGroupPairInfo->chooseCandiSegStartPosClosest2CertainPos_NorOrRcm(segInfo_Nor, tmp, 
					currentPathEndPos + tmpSegStartPos_inRead - currentPathEndPos_inRead, 
					tmpGapCase, false);
			if(tmpBestCandiMapSeg_index_inOriSegInfo == -1)
			{}
			else
			{
				currentPathEndPos = segInfo_Rcm->returnSegmentMapPos_end(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo);
				finalCandiPath_Rcm.push_back(pair<int,int>(tmpSegGroupNO_inOriSegInfo, tmpBestCandiMapSeg_index_inOriSegInfo));
				finalCandiPath_gapCase_Rcm.push_back(tmpGapCase);
			}
		}

		chrInt = this->checkTwoEndsChrConsistentOrNot_bool(indexInfo);
		if(chrInt < 0)
			return false;
	}

	void pushBack2TargetPathInfo(Path_Info* pathInfo_Nor, Path_Info* pathInfo_Rcm, int tmpPathNumIndex)
	{
		pathInfo_Nor->pushBack_CandiPathPairInfo(finalCandiPath_Nor, finalCandiPathJumpCodeVec_Nor.second,
			mismatchPosVec_Nor, mismatchCharVec_Nor, tmpPathNumIndex);
		pathInfo_Nor->pushBack_CandiPathPairInfo(finalCandiPath_Rcm, finalCandiPathJumpCodeVec_Rcm.second,
			mismatchPosVec_Rcm, mismatchCharVec_Rcm, tmpPathNumIndex);
	}

	void fixDoubleAnchorGap(
		vector<int>& tmpMismatchPosVec, vector<char>& tmpMismatchCharVec, 
		vector<Jump_Code>& tmpJumpCodeVec,
		int tmpGapCase, int tmpChrInt, 
		int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2,
		int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		bool fixDoubleAnchor_bool = false;
		int extendBackNumMax = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		if(extendBackNumMax > segmentMapPos_2 - 1)
		{
			extendBackNumMax = segmentMapPos_2 - 1; 
		}
		int extendBackNum = this->extendBack(segmentLocInRead_2, readSeq_inProcess, 
			segmentMapPos_2, tmpChrInt, extendBackNumMax);
		
		segmentLocInRead_2 = segmentLocInRead_2 - extendBackNum;
		segmentLength_2 = segmentLength_2 + extendBackNum;
		segmentMapPos_2 = segmentMapPos_2 - extendBackNum;

		if(tmpGapCase = FIX_MATCH)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorMatch(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_INSERTION)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorInsertion(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_DELETION)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorDeletion(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else if(tmpGapCase = FIX_SPLICE)
		{
			fixDoubleAnchor_bool = this->fixDoubleAnchorSplice(
				tmpGapMismatchPosVec, tmpGapMismatchCharVec, tmpGapJumpCodeVec, 
				chrInt, segmentLocInRead_1, segmentLocInRead_2, 
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2,
				indexInfo, readSeq_inProcess)
		}
		else
		{
			cout << "error in tmpGapCase ..." << endl;
			exit(1);
		}
		return fixDoubleAnchor_bool;
	}

	bool fixDoubleAnchorSplice(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, 
		vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		int tmpBuffer_left = FixSpliceBuffer;//+1;
		if(tmpBuffer_left > segmentLength_1 - 2) //anchor >= 2
		{
			tmpBuffer_left = segmentLength_1 - 2;
		}
		int tmpBuffer_right = FixSpliceBuffer;//+1;
		if(tmpBuffer_right > segmentLength_2 - 2)
		{
			tmpBuffer_right = segmentLength_2 - 2;
		}

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1 + tmpBuffer_left + tmpBuffer_right;
		int spliceJunctionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);

		//size_t prefix_length = 0;
		size_t max_double_splice_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;

		if((annotation_provided_bool && Do_annotation_only_bool) || !(annotation_provided_bool))
		{
			FixDoubleAnchor_Splice_Info* fixSpliceInfo = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo->detectBestSpliceSite(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
			
				if(splice_fixed && fixSpliceInfo->fixSpliceResultConfident())//((fixSpliceInfo->returnBestSplice_canonicalOrNot()) || (fixSpliceInfo->returnBestSplice_mismatchNum() <= 2) ) )
				{
					int prefix_length = fixSpliceInfo->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

					//fixGapVec[index_fixGapVec].first = true;
					//(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo->returnBestSplice_mismatchNum();
					//((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
					//((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
					//((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

					fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
					//this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
					//this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
					fixSpliceInfo->copyMismatchPos2TargetVec(tmpMismatchPosVec);
					fixSpliceInfo->copyMismatchPos2TargetVec(tmpMismatchCharVec);
					tmpGapJumpCodeVec.push_back(firstMatchJumpCode);
					tmpGapJumpCodeVec.push_back(spliceJumpCode);	
					tmpGapJumpCodeVec.push_back(secondMatchJumpCode);
					delete fixSpliceInfo;
					return true;
				}
				else
				{
					//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
					FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo = new FixDoubleAnchor_Splice_Complicate_Info();
					bool complicate_splice_fixed = fixComplicateSpliceInfo->detectComplicateSplice(
						segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
						segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
						readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
					//cout << "finish fixing complicated SJ ....." << endl;	
					bool complicated_splice_fixed_better_bool = false;
					if(complicate_splice_fixed)
					{
						complicated_splice_fixed_better_bool = fixComplicateSpliceInfo->compared2SpliceInfo(splice_fixed,
							fixSpliceInfo);
					}
					else
					{}

					if(complicated_splice_fixed_better_bool)
					{
						int prefix_match_length = fixComplicateSpliceInfo->return_prefix_match_length_best();
						int first_jumpCode_length = fixComplicateSpliceInfo->return_first_jumpCode_length_best();
						string first_jumpCode_type = fixComplicateSpliceInfo->return_first_jumpCode_type_best();
						int mid_match_length = fixComplicateSpliceInfo->return_mid_match_length_best();
						int second_jumpCode_length = fixComplicateSpliceInfo->return_second_jumpCode_length_best();
						string second_jumpCode_type = fixComplicateSpliceInfo->return_second_jumpCode_type_best();
						int suffix_match_length = fixComplicateSpliceInfo->return_suffix_match_length_best();
						
						Jump_Code prefixMatchJumpCode(prefix_match_length - tmpBuffer_left, "M");
						Jump_Code firstJumpCode(first_jumpCode_length, first_jumpCode_type);
						Jump_Code midMatchJumpCode(mid_match_length, "M");
						Jump_Code secondJumpCode(second_jumpCode_length, second_jumpCode_type); 				
						Jump_Code suffixMatchJumpCode(suffix_match_length - tmpBuffer_right + segmentLength_2, "M");

						//fixGapVec[index_fixGapVec].first = true;
						//(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo->return_mismatch_bestComplicatedSplice();
						//((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
						//((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
						//((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
						//((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
						//((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

						fixComplicateSpliceInfo->generateBestComplicateSJMismatchVec(
							segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
						//this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo, index_fixGapVec);
						//this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo, index_fixGapVec);
						fixComplicateSpliceInfo->copyMismatchPos2TargetVec(tmpMismatchPosVec);
						fixComplicateSpliceInfo->copyMismatchChar2TargetVec(tmpMismatchCharVec);
						tmpGapJumpCodeVec.push_back(prefixMatchJumpCode);	
						tmpGapJumpCodeVec.push_back(firstJumpCode);
						tmpGapJumpCodeVec.push_back(midMatchJumpCode);	
						tmpGapJumpCodeVec.push_back(secondJumpCode);	
						tmpGapJumpCodeVec.push_back(suffixMatchJumpCode);	

						delete fixSpliceInfo;
						delete fixComplicateSpliceInfo;
						return true;
					}
					else if(splice_fixed)
					{
						int prefix_length = fixSpliceInfo->returnBestSplice_prefixMatchLength();
						int firstMatchLength = //segmentLength_1 
							- tmpBuffer_left + prefix_length;
						int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
						//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code spliceJumpCode(spliceJunctionLength, "N");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

						//fixGapVec[index_fixGapVec].first = true;
						//(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo->returnBestSplice_mismatchNum();
						//((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						//((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						//((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						//this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
						//this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
						fixSpliceInfo->copyMismatchPos2TargetVec(tmpMismatchPosVec);
						fixSpliceInfo->copyMismatchPos2TargetVec(tmpMismatchCharVec);
						tmpGapJumpCodeVec.push_back(firstMatchJumpCode);
						tmpGapJumpCodeVec.push_back(spliceJumpCode);	
						tmpGapJumpCodeVec.push_back(secondMatchJumpCode);						
						
						delete fixSpliceInfo;
						delete fixComplicateSpliceInfo;
						return true;
					}
					else
					{
						delete fixSpliceInfo;
						delete fixComplicateSpliceInfo;
						return false;
						//fixGapVec[index_fixGapVec].first = false;
					}
					//delete fixComplicateSpliceInfo;
				}
		}/*
		else // annotation_provided_bool == true && Do_annotation_only_bool == false
		{ 
			bool annotated_splice_fixed_bool = false;
			FixDoubleAnchor_Splice_Info* fixSpliceInfo_annotated = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo_annotated->detectBestSpliceSite//_prefer_canonical_lessMismatch
				(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, //Do_annotation_only_bool, 
				true, annotationInfo);
			//if(FIX_INDEL_AROUND_SJ_BOOL)
			//{
				if(splice_fixed && 
					((fixSpliceInfo_annotated->returnBestSplice_canonicalOrNot()) 
						|| (fixSpliceInfo_annotated->returnBestSplice_mismatchNum() <= 2) ) )
				{
					int prefix_length = fixSpliceInfo_annotated->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

					fixGapVec[index_fixGapVec].first = true;					
					// fixed with annotated SJ
					annotated_splice_fixed_bool = true;

					(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
					((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
					((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
					((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

					//if(STORE_MISMATCH_POS)
					//{
					fixSpliceInfo_annotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
					this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
					//	if(STORE_MISMATCH_CHA)
					//	{
					this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
					//	}
					//}
				}
				else
				{
					/////  test fixing complicated SJ  ///////
					//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
					FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_annotated
						= new FixDoubleAnchor_Splice_Complicate_Info();
					bool complicate_splice_fixed = fixComplicateSpliceInfo_annotated->detectComplicateSplice(
						segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
						segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
						readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch,
						annotation_provided_bool, true, annotationInfo);
					//cout << "finish fixing complicated SJ ....." << endl;	
					bool complicated_splice_fixed_better_bool = false;
					if(complicate_splice_fixed)
					{
						complicated_splice_fixed_better_bool = fixComplicateSpliceInfo_annotated->compared2SpliceInfo(splice_fixed,
							fixSpliceInfo_annotated);							
					}
					else
					{}

					if(complicated_splice_fixed_better_bool)
					{
						int prefix_match_length = fixComplicateSpliceInfo_annotated->return_prefix_match_length_best();
						int first_jumpCode_length = fixComplicateSpliceInfo_annotated->return_first_jumpCode_length_best();
						string first_jumpCode_type = fixComplicateSpliceInfo_annotated->return_first_jumpCode_type_best();
						int mid_match_length = fixComplicateSpliceInfo_annotated->return_mid_match_length_best();
						int second_jumpCode_length = fixComplicateSpliceInfo_annotated->return_second_jumpCode_length_best();
						string second_jumpCode_type = fixComplicateSpliceInfo_annotated->return_second_jumpCode_type_best();
						int suffix_match_length = fixComplicateSpliceInfo_annotated->return_suffix_match_length_best();
						
						Jump_Code prefixMatchJumpCode(prefix_match_length - tmpBuffer_left, "M");
						Jump_Code firstJumpCode(first_jumpCode_length, first_jumpCode_type);
						Jump_Code midMatchJumpCode(mid_match_length, "M");
						Jump_Code secondJumpCode(second_jumpCode_length, second_jumpCode_type); 				
						Jump_Code suffixMatchJumpCode(suffix_match_length - tmpBuffer_right + segmentLength_2, "M");

						fixGapVec[index_fixGapVec].first = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;

						(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo_annotated->return_mismatch_bestComplicatedSplice();
						((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

						if(STORE_MISMATCH_POS)
						{
							fixComplicateSpliceInfo_annotated->generateBestComplicateSJMismatchVec(
								segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
							this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_annotated, index_fixGapVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_annotated, index_fixGapVec);
							}
						}		
					}
					else if(splice_fixed)
					{
						int prefix_length = fixSpliceInfo_annotated->returnBestSplice_prefixMatchLength();
						int firstMatchLength = //segmentLength_1 
							- tmpBuffer_left + prefix_length;
						int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
						//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code spliceJumpCode(spliceJunctionLength, "N");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

						fixGapVec[index_fixGapVec].first = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;			

						(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
						((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						if(STORE_MISMATCH_POS)
						{
							fixSpliceInfo_annotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
							}
						}
					}
					else
					{
						fixGapVec[index_fixGapVec].first = false;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = false;	
					}
					delete fixComplicateSpliceInfo_annotated;
				}
			delete fixSpliceInfo_annotated;

			if(annotated_splice_fixed_bool) // fixed with annotated SJ
			{}
			else
			{
				FixDoubleAnchor_Splice_Info* fixSpliceInfo_unannotated = new FixDoubleAnchor_Splice_Info();
				bool splice_fixed = fixSpliceInfo_unannotated->detectBestSpliceSite//_prefer_canonical_lessMismatch
					(
					segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
					segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
					readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
					false, false, annotationInfo);

					if(splice_fixed && 
						((fixSpliceInfo_unannotated->returnBestSplice_canonicalOrNot()) 
							|| (fixSpliceInfo_unannotated->returnBestSplice_mismatchNum() <= 2) ) )
					{
						int prefix_length = fixSpliceInfo_unannotated->returnBestSplice_prefixMatchLength();
						int firstMatchLength = //segmentLength_1 
							- tmpBuffer_left + prefix_length;
						int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
						//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code spliceJumpCode(spliceJunctionLength, "N");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

						fixGapVec[index_fixGapVec].first = true;					
						// fixed with unannotated SJ

						(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
						((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						if(STORE_MISMATCH_POS)
						{
							fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
							}
						}
					}
					else
					{
						/////  test fixing complicated SJ  ///////
						//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
						FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_unannotated
							= new FixDoubleAnchor_Splice_Complicate_Info();
						bool complicate_splice_fixed = fixComplicateSpliceInfo_unannotated->detectComplicateSplice(
							segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
							segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
							readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch,
							false, false, annotationInfo);
						//cout << "finish fixing complicated SJ ....." << endl;	
						bool complicated_splice_fixed_better_bool = false;
						if(complicate_splice_fixed)
						{
							complicated_splice_fixed_better_bool 
								= fixComplicateSpliceInfo_unannotated->compared2SpliceInfo(splice_fixed, fixSpliceInfo_unannotated);
						}
						else
						{}

						if(complicated_splice_fixed_better_bool)
						{
							int prefix_match_length = fixComplicateSpliceInfo_unannotated->return_prefix_match_length_best();
							int first_jumpCode_length = fixComplicateSpliceInfo_unannotated->return_first_jumpCode_length_best();
							string first_jumpCode_type = fixComplicateSpliceInfo_unannotated->return_first_jumpCode_type_best();
							int mid_match_length = fixComplicateSpliceInfo_unannotated->return_mid_match_length_best();
							int second_jumpCode_length = fixComplicateSpliceInfo_unannotated->return_second_jumpCode_length_best();
							string second_jumpCode_type = fixComplicateSpliceInfo_unannotated->return_second_jumpCode_type_best();
							int suffix_match_length = fixComplicateSpliceInfo_unannotated->return_suffix_match_length_best();
							
							Jump_Code prefixMatchJumpCode(prefix_match_length - tmpBuffer_left, "M");
							Jump_Code firstJumpCode(first_jumpCode_length, first_jumpCode_type);
							Jump_Code midMatchJumpCode(mid_match_length, "M");
							Jump_Code secondJumpCode(second_jumpCode_length, second_jumpCode_type); 				
							Jump_Code suffixMatchJumpCode(suffix_match_length - tmpBuffer_right + segmentLength_2, "M");

							fixGapVec[index_fixGapVec].first = true;
							// fixed with unannotated SJ

							(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo_unannotated->return_mismatch_bestComplicatedSplice();
							((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
							((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

							//if(STORE_MISMATCH_POS)
							//{
								fixComplicateSpliceInfo_unannotated->generateBestComplicateSJMismatchVec(
									segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
								this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_unannotated, index_fixGapVec);
								//if(STORE_MISMATCH_CHA)
								//{
									this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_unannotated, index_fixGapVec);
								//}
							//}		
						}
						else if(splice_fixed)
						{
							int prefix_length = fixSpliceInfo_unannotated->returnBestSplice_prefixMatchLength();
							int firstMatchLength = //segmentLength_1 
								- tmpBuffer_left + prefix_length;
							int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
							//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

							Jump_Code firstMatchJumpCode(firstMatchLength, "M");
							Jump_Code spliceJumpCode(spliceJunctionLength, "N");
							Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

							fixGapVec[index_fixGapVec].first = true;
							// fixed with unannotated SJ	

							(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
							((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
							((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

							fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
						}
						else
						{
							fixGapVec[index_fixGapVec].first = false;
							// fixed with annotated SJ
							//annotated_splice_fixed_bool = false;	
						}
						delete fixComplicateSpliceInfo_unannotated;
					}

				delete fixSpliceInfo_unannotated;
			}
		}*/
		return;
	}

	bool fixDoubleAnchorDeletion(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, 
		vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
			//+ buffer_left + buffer_right;
		int deletionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);

		if(subSeqLengthInProcess < 2)
		{
			Jump_Code deletionJumpCode(deletionLength, "D");
			Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess, "M");

			fixGapVec[index_fixGapVec].first = true;
			(fixGapVec[index_fixGapVec].second).first = subSeqLengthInProcess;
			((fixGapVec[index_fixGapVec].second).second).push_back(deletionJumpCode);
			((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);	

			if(STORE_MISMATCH_POS)
			{
				//delInfo->
				this->insertMismatchPosVec_fixDeletionInfo_SmallGap(index_fixGapVec, segmentLocInRead_1 + segmentLength_1,
					subSeqLengthInProcess);
				if(STORE_MISMATCH_CHA)
				{
					this->insertMismatchCharVec_fixDeletionInfo_SmallGap(index_fixGapVec, subSeqLengthInProcess, 
						chrNameInt, segmentMapPos_2 - subSeqLengthInProcess, indexInfo);
				}
			}
		}
		else
		{
			FixDoubleAnchor_Deletion_Info* delInfo = new FixDoubleAnchor_Deletion_Info();
			int tmp_toFix_deletion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_deletion_read_end = segmentLocInRead_2 - 1;
			//int subSeqLengthInProcess = tmp_toFix_deletion_read_end - tmp_toFix_deletion_read_start + 1;
			int tmp_toFix_deletion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_deletion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 1;
			bool deletion_fixed = delInfo->detectBestDeletion_lessMismatch(tmp_toFix_deletion_read_start, 
				tmp_toFix_deletion_read_end, tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, 
				readSeq_inProcess, indexInfo, chrNameInt, tmp_max_allowed_mismatchNum);
			if(deletion_fixed)
			{
				int firstMatchLength = //segmentLength_1 + 
						//prefix_length;
						delInfo->return_best_deletion_prefix_match_length();
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - firstMatchLength;
				int mismatch_bits = delInfo->return_best_deletion_mismatch();
				Jump_Code firstMatchJumpCode(firstMatchLength
					//- buffer_left
					, "M");
				Jump_Code deletionJumpCode(deletionLength, "D");
				Jump_Code secondMatchJumpCode(secondMatchLength
					//- buffer_right
					, "M");
				tmpGapJumpCodeVec.push_back(firstMatchJumpCode);
				tmpGapJumpCodeVec.push_back(deletionJumpCode);
				tmpGapJumpCodeVec.push_back(secondMatchJumpCode);

				delInfo->generateBestDeletionMismatchVec((segmentLocInRead_1 + segmentLength_1));
				delInfo->copyMismatchPos2TargetVec(tmpGapMismatchPosVec);
				delInfo->copyMismatchChar2TargetVec(tmpGapMismatchCharVec);
				delete insInfo;
				return true;
			}
			else
			{
				delete delInfo;
				return false;
			}
		}
	}

	bool fixDoubleAnchorInsertion(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, 
		vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int insertionLength = (segmentMapPos_1 - segmentLocInRead_1) - (segmentMapPos_2 - segmentLocInRead_2);

		if(subSeqLengthInProcess <= insertionLength)
		{
			int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - insertionLength;
			if(secondMatchLength > 0)
			{
				Jump_Code midInsertionJumpCode(insertionLength, "I");	
				Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess - insertionLength, "M");		
				
				tmpGapJumpCodeVec.push_back(midInsertionJumpCode);
				tmpGapJumpCodeVec.push_back(secondMatchJumpCode);
				return true;
				// no mismatch pos or char generated 		
			}
			else
			{
				return false;
			}
		}	 
		else
		{
			FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info();
			int tmp_toFix_insertion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_insertion_read_end = segmentLocInRead_2 - 1;
			int tmp_toFix_insertion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_insertion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;  
			// how to set the max-mismatch parameter

			bool insertion_fixed = insInfo->detectBestInsertion_lessMismatch(
				tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
				tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
				chromNameInt, tmp_max_allowed_mismatchNum);
			if(insertion_fixed)
			{
				int firstMatchLength = insInfo->return_best_insertion_prefix_match_length();
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - firstMatchLength - insertionLength;
				int mismatch_bits = insInfo->return_best_insertion_mismatch();
				if(secondMatchLength > 0)			
				{
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");

					tmpGapJumpCodeVec.push_back(firstMatchJumpCode);
					tmpGapJumpCodeVec.push_back(insertionJumpCode);		
					tmpGapJumpCodeVec.push_back(secondMatchJumpCode);
					insInfo->generateBestInsertionMismatchVec((segmentLocInRead_1 + segmentLength_1));
					insInfo->copyMismatchPos2TargetVec(tmpGapMismatchPosVec);
					insInfo->copyMismatchChar2TargetVec(tmpGapMismatchCharVec);
					delete insInfo;
					return true;
				}
				else
				{
					delete insInfo;
					return false;
				}
			}
			else
			{
				delete insInfo;
				return false;
			}
		}
	}

	bool fixDoubleAnchorMatch(
		vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec, 
		vector<Jump_Code>& tmpGapJumpCodeVec, 
		int chrInt, int segmentLocInRead_1, int segmentLocInRead_2, 
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		Index_Info* indexInfo, const string& readSeq_inProcess)
	{
		int subSeq_toProcess_len = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		if(subSeq_toProcess_len < 2)
		{
			Jump_Code matchJumpCode(subSeq_toProcess_len + segmentLength_2, "M");
			tmpGapJumpCodeVec.push_back(matchJumpCode);
			for(int tmp = 0; tmp < subSeq_toProcess_len; tmp++)
			{
				tmpGapMismatchPosVec.push_back(segmentLocInRead_1 + segmentLength_1 + tmp);
				tmpGapMismatchCharVec.push_back(indexInfo->getCharInChromosome(chrInt, segmentMapPos_1 + segmentLength_1 + tmp));
			}
			return true;
		}
		else
		{
			string readSubSeq_toProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeq_toProcess_len); 
			string chromSubSeq_toProcess = indexInfo->returnChromStrSubstr(chrInt, segmentMapPos_1 + segmentLength_1,
				subSeq_toProcess_len);

			int max_mismatch = subSeq_toProcess_len / LengthOfSeqPerMismatchAllowed + 1;

			FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
			bool scoreStringBool = fixMatchInfo->fixMatch(readSubSeqInProcess, chromSubSeqInProcess,
				max_mismatch, segmentLocInRead_1 + segmentLength_1);
			if(scoreStringBool)
			{
				Jump_Code matchJumpCode(//segmentLength_1 + 
					subSeq_toProcess_len + segmentLength_2, "M");
				tmpGapJumpCodeVec.push_back(matchJumpCode);
				fixMatchInfo->copyMismatchPos2TargetVec(tmpGapMismatchPosVec);
				fixMatchInfo->copyMismatchCahr2TargetVec(tmpGapMismatchCharVec);		
			}	
			delete fixMatchInfo;
		}
		return fixDoubleAnchor_bool;
	}	
};


class SubRegion_Info
{
private:
	vector< SubSegGroupPair_Info* > subSegGroupPairVec_Nor1Rcm2;
	vector< int > twoEndsSegExistsIndexVec_Nor1Rcm2; // < index in subSegGroupPairVec_Nor1Rcm2 >
	map<int, int> subSegGroupPairVec_Nor1Rcm2_indexMap;  // <mapRegionIndex, index_in_subSegGroupPairVec_NorRcm>
	vector< int > candiPathGeneratedPairIndexVec_Nor1Rcm2; // < index in subSegGroupPairVec_Nor1Rcm2 >
	vector< CandiPathPair_Info* > candiPathPairVec_Nor1Rcm2;
	vector< int > stitchedCandiPathPairIndexVec_Nor1Rcm2; // <index in candiPathpairVec_Nor1Rcm2>
 
	vector< SubSegGroupPair_Info* > subSegGroupPairVec_Nor2Rcm1;
	vector< int > twoEndsSegExistsIndexVec_Nor2Rcm1; // < index in subSegGroupPairVec_Nor2Rcm1 >
	map<int, int> subSegGroupPairVec_Nor2Rcm1_indexMap; // <mapRegionIndex, index_in_subSegGroupPairVec_NorRcm>
	vector< int > candiPathGeneratedPairIndexVec_Nor2Rcm1; // < index in subSegGroupPairVec_Nor2Rcm1 >
	vector< CandiPathPair_Info* > candiPathPairVec_Nor2Rcm1;
	vector< int > stitchedCandiPathPairIndexVec_Nor2Rcm1; // <index in candiPathpairVec_Nor2Rcm1>
 
	bool subSegGroupPairExists_bool;
public:
	SubRegion_Info()
	{
		SubRegionSize = 300000;
		subSegGroupPairExists_bool = false;
	}

	void generateCandiPath(Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2)
	{
		this->takeSegInfo(segInfo_Nor1, segInfo_Rcm1, segInfo_Nor2, segInfo_Rcm2);
		this->getCandiPathPair(segInfo_Nor1, segInfo_Rcm1, segInfo_Nor2, segInfo_Rcm2,
			indexInfo);
		this->stitchCandiPathPair(segInfo_Nor1, segInfo_Rcm1,
			segInfo_Nor2, segInfo_Rcm2, indexInfo);
		this->pushBackFixedPath2PathInfo(pathInfo_Nor1, pathInfo_Rcm1,
			pathInfo_Nor2, pathInfo_Rcm2, indexInfo);
	}

	void pushBackFixedPath2PathInfo(Path_Info* pathInfo_Nor1, Path_Info* pathInfo_Rcm1,
		Path_Info* pathInfo_Nor2, Path_Info* pathInfo_Rcm2)
	{
		for(int tmp = 0; tmp < stitchedCandiPathPairIndexVec_Nor1Rcm2.size(); tmp++)
		{
			int tmpIndex_inCandiPathPairVec_Nor1Rcm2 = stitchedCandiPathPairIndexVec_Nor1Rcm2[tmp];
			candiPathPairVec_Nor1Rcm2[tmpIndex_inCandiPathPairVec_Nor1Rcm2]->pushBack2TargetPathInfo(
				pathInfo_Nor1, pathInfo_Rcm2, tmp);
			tmp++;
		}
		for(int tmp = 0; tmp < stitchedCandiPathPairIndexVec_Nor2Rcm1.size(); tmp++)
		{
			int tmpIndex_inCandiPathPairVec_Nor2Rcm1 = stitchedCandiPathPairIndexVec_Nor2Rcm1[tmp];
			candiPathPairVec_Nor2Rcm1[tmpIndex_inCandiPathPairVec_Nor2Rcm1]->pushBack2TargetPathInfo(
				pathInfo_Nor2, pathInfo_Rcm1, tmp);
			tmp++;
		}	
	}

	void takeSegInfo(Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2, 
		bool normalMapMain, bool rcmMapMain, 
		bool normalMapMain_PE, bool rcmMapMain_PE) // fix me, for now, do not consider the segments in the borders
	{
		map<int,int>::iterator tmpMapIter;
		int tmpSize_subSegGroupPairVec_Nor1Rcm2 = 0;
		int tmpSize_subSegGroupPairVec_Nor2Rcm1 = 0;
		if(normalMapMain)
		{	
			for(int tmp = 0; tmp < segInfo_Nor1->returnSegmentNum(); tmp++)
			{
				for(int tmp2 = 0; tmp2 < segInfo_Nor1->returnSegmentAlignNUm(); tmp2++)
				{
					unsigned int tmpMapPos = segInfo_Nor1->returnSegmentMapPos(tmp, tmp2);
					int tmpRegionIndex = tmpMapPos / SubRegionSize;
					tmpMapIter = subSegGroupPairVec_Nor1Rcm2_indexMap.find(tmpRegionIndex);
					if(tmpMapIter != subSegGroupPairVec_Nor1Rcm2_indexMap.end()) // // fix me, for now, do not consider the segments in the borders
					{
						subSegGroupPairVec_Nor1Rcm2[tmpMapIter->second]->pushBackNewSeg(tmp, tmp2, true); // segGroupNO, segCandiNO, nor_or_rcm
					}
					else // can not find any existing subSegGroupPair, creat a new one
					{
						SubSegGroupPair_Info* newSubSegGroupPairInfo = new SubSegGroupPair_Info(true, tmpRegionIndex);// for_or_rev, regionIndex
						newSubSegGroupPairInfo->initiateWith1stSeg(tmp, tmp2, true); // segGroupNO, segCandiNO, nor_or_rcm
						subSegGroupPairVec_Nor1Rcm2.push_back(newSubSegGroupPairInfo);
						//twoEndsSegExistsVec_bool_Nor1Rcm2.push_back(false);
						subSegGroupPairVec_Nor1Rcm2_indexMap.insert(pair<int,int>(tmpRegionIndex, tmpSize_subSegGroupPairVec_Nor1Rcm2));
						tmpSize_subSegGroupPairVec_Nor1Rcm2 ++;
					}
				}
			}
		}
		if(normalMapMain_PE)
		{	
			for(int tmp = 0; tmp < segInfo_Nor2->returnSegmentNum(); tmp++)
			{
				for(int tmp2 = 0; tmp2 < segInfo_Nor2->returnSegmentAlignNUm(); tmp2++)
				{
					unsigned int tmpMapPos = segInfo_Nor2->returnSegmentMapPos(tmp, tmp2);
					int tmpRegionIndex = tmpMapPos / SubRegionSize;
					tmpMapIter = subSegGroupPairVec_Nor2Rcm1_indexMap.find(tmpRegionIndex);
					if(tmpMapIter != subSegGroupPairVec_Nor2Rcm1_indexMap.end()) // // fix me, for now, do not consider the segments in the borders
					{
						subSegGroupPairVec_Nor2Rcm1[tmpMapIter->second]->pushBackNewSeg(tmp, tmp2, true); // segGroupNO, segCandiNO, nor_or_rcm
					}
					else // can not find any existing subSegGroupPair, creat a new one
					{
						SubSegGroupPair_Info* newSubSegGroupPairInfo = new SubSegGroupPair_Info(false, tmpRegionIndex); // for_or_rev, regionIndex
						newSubSegGroupPairInfo->initiateWith1stSeg(tmp, tmp2, true); // segGroupNO, segCandiNO, nor_or_rcm
						subSegGroupPairVec_Nor2Rcm1.push_back(newSubSegGroupPairInfo);
						//twoEndsSegExistsVec_bool_Nor2Rcm1.push_back(false);
						subSegGroupPairVec_Nor2Rcm1_indexMap.insert(pair<int,int>(tmpRegionIndex, tmpSize_subSegGroupPairVec_Nor2Rcm1));
						tmpSize_subSegGroupPairVec_Nor2Rcm1 ++;
					}
				}
			}
		}
		for(int tmp = 0; tmp < segInfo_Rcm2->returnSegmentNum(); tmp++)
		{
			for(int tmp2 = 0; tmp2 < segInfo_Rcm2->returnSegmentAlignNUm(); tmp2++)
			{
				unsigned int tmpMapPos = segInfo_Rcm2->returnSegmentMapPos(tmp, tmp2);
				int tmpRegionIndex = tmpMapPos / SubRegionSize;
				tmpMapIter = subSegGroupPairVec_Nor1Rcm2_indexMap.find(tmpRegionIndex);
				if(tmpMapIter != subSegGroupPairVec_Nor1Rcm2_indexMap.end()) // // fix me, for now, do not consider the segments in the borders
				{
					subSegGroupPairVec_Nor1Rcm2[tmpMapIter->second]->pushBackNewSeg(tmp, tmp2, false); // segGroupNO, segCandiNO, nor_or_rcm
					//twoEndsSegExistsVec_bool_Nor1Rcm2[tmpMapIter->second] = true;
					twoEndsSegExistsIndexVec_Nor1Rcm2.push_back(tmpMapIter->second);
					subSegGroupPairExists_bool = true;
				}
				else // can not find any existing subSegGroupPair, creat a new one
				{
					SubSegGroupPair_Info* newSubSegGroupPairInfo = new SubSegGroupPair_Info(true, tmpRegionIndex); // for_or_rev, regionIndex
					newSubSegGroupPairInfo->initiateWith1stSeg(tmp, tmp2, false); // segGroupNO, segCandiNO, nor_or_rcm
					subSegGroupPairVec_Nor1Rcm2.push_back(newSubSegGroupPairInfo);
					//twoEndsSegExistsVec_bool_Nor1Rcm2.push_back(false);
					subSegGroupPairVec_Nor1Rcm2_indexMap.insert(pair<int,int>(tmpRegionIndex, tmpSize_subSegGroupPairVec_Nor1Rcm2));
					tmpSize_subSegGroupPairVec_Nor1Rcm2 ++;
				}
			}
		}	

		for(int tmp = 0; tmp < segInfo_Rcm1->returnSegmentNum(); tmp++)
		{
			for(int tmp2 = 0; tmp2 < segInfo_Rcm1->returnSegmentAlignNUm(); tmp2++)
			{
				unsigned int tmpMapPos = segInfo_Rcm1->returnSegmentMapPos(tmp, tmp2);
				int tmpRegionIndex = tmpMapPos / SubRegionSize;
				tmpMapIter = subSegGroupPairVec_Nor2Rcm1_indexMap.find(tmpRegionIndex);
				if(tmpMapIter != subSegGroupPairVec_Nor2Rcm1_indexMap.end()) // // fix me, for now, do not consider the segments in the borders
				{
					subSegGroupPairVec_Nor2Rcm1[tmpMapIter->second]->pushBackNewSeg(tmp, tmp2, false); // segGroupNO, segCandiNO, nor_or_rcm
					//twoEndsSegExistsVec_bool_Nor2Rcm1[tmpMapIter->second] = true;
					twoEndsSegExistsIndexVec_Nor2Rcm1.push_back(tmpMapIter->second);
					subSegGroupPairExists_bool = true;
				}
				else // can not find any existing subSegGroupPair, creat a new one
				{
					SubSegGroupPair_Info* newSubSegGroupPairInfo = new SubSegGroupPair_Info(false, tmpRegionIndex); // for_or_rev, regionIndex
					newSubSegGroupPairInfo->initiateWith1stSeg(tmp, tmp2, false); // segGroupNO, segCandiNO, nor_or_rcm
					subSegGroupPairVec_Nor2Rcm1.push_back(newSubSegGroupPairInfo);
					//twoEndsSegExistsVec_bool_Nor2Rcm1.push_back(false);
					subSegGroupPairVec_Nor2Rcm1_indexMap.insert(pair<int,int>(tmpRegionIndex, tmpSize_subSegGroupPairVec_Nor2Rcm1));
					tmpSize_subSegGroupPairVec_Nor2Rcm1 ++;
				}
			}
		}			
	}

	void getCandiPathPair(Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2, Index_Info* indexInfo)
	{

		for(int tmp = 0; tmp < twoEndsSegExistsIndexVec_Nor1Rcm2.size(); tmp++)
		{
			int tmpIndexInSubSegGroupPairVec_Nor1Rcm2 = twoEndsSegExistsIndexVec_Nor1Rcm2[tmp];
			CandiPathPair_Info* tmpCandiPathPairInfo = new candiPathPair_Info();

			bool pathGenerated_bool = tmpCandiPathPairInfo->generateCandiPath_uniqueOnly_pair(
					subSegGroupPairVec_Nor1Rcm2[tmpIndexInSubSegGroupPairVec_Nor1Rcm2], segInfo_Nor1,
					segInfo_Rcm2, indexInfo);
			if(pathGenerated_bool)
				candiPathGeneratedPairIndexVec_Nor1Rcm2.push_back(tmpIndexInSubSegGroupPairVec_Nor1Rcm2);
			else
			{}
		}
		for(int tmp = 0; tmp < twoEndsSegExistsIndexVec_Nor2Rcm1.size(); tmp++)
		{
			int tmpIndexInSubSegGroupPairVec_Nor2Rcm1 = twoEndsSegExistsIndexVec_Nor2Rcm1[tmp];
			CandiPathPair_Info* tmpCandiPathPairInfo = new candiPathPair_Info();

			bool pathGenerated_bool = tmpCandiPathPairInfo->generateCandiPath_uniqueOnly_pair(
				subSegGroupPairVec_Nor2Rcm1[tmpIndexInSubSegGroupPairVec_Nor2Rcm1], segInfo_Nor2,
				segInfo_Rcm1, indexInfo);
			if(pathGenerated_bool)
				candiPathGeneratedPairIndexVec_Nor2Rcm1.push_back(tmpIndexInSubSegGroupPairVec_Nor2Rcm1);
			else
			{}
		}
	}

	void stitchCandiPathPair(Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < candiPathGeneratedPairIndexVec_Nor1Rcm2.size(); tmp++)
		{	
			bool stitchedSuccess_bool = candiPathGeneratedPairIndexVec_Nor1Rcm2[tmp]->fixGapsInPathPair(
				segInfo_Nor1, segInfo_Rcm2, indexInfo);
			if(stitchedSuccess_bool)
			{
				stitchedCandiPathPairIndexVec_Nor1Rcm2.push_back(tmp);
			}
		}	
		for(int tmp = 0; tmp < candiPathGeneratedPairIndexVec_Nor2Rcm1.size(); tmp++)
		{
			bool stitchedSuccess_bool = candiPathGeneratedPairIndexVec_Nor2Rcm1[tmp]->fixGapsInPathPair(
				segInfo_Nor2, segInfo_Rcm1, indexInfo);
			if(stitchedSuccess_bool)
			{
				stitchedCandiPathPairIndexVec_Nor2Rcm1.push_back(tmp);
			}
		}
	}

};












#endif