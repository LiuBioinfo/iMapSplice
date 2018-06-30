// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENOMICREGION_INFO_H
#define GENOMICREGION_INFO_H

#include "../../../general/index_info.h"

using namespace std;

class GenomicRegion_Info
{
private:
	int index_inGenomicRegionVec;

	int chrNameInt;
	int chrStartPos;
	int chrEndPos;
	int chrPosRange;
public:
	GenomicRegion_Info()
	{}

	void initiate(int tmpIndexInGenomicRegionVec,
		int tmpChrNameInt, int tmpChrStartPos, int tmpChrEndPos)
	{
		index_inGenomicRegionVec = tmpIndexInGenomicRegionVec;
		chrNameInt = tmpChrNameInt;
		chrStartPos = tmpChrStartPos;
		chrEndPos = tmpChrEndPos;
		chrPosRange = tmpChrEndPos - tmpChrStartPos + 1;
	}

	int returnGenomicRegionIndex_genomicRegionInfo()
	{
		return index_inGenomicRegionVec;
	}

	int returnChrNameInt_genomicRegionInfo()
	{
		return chrNameInt;
	}

	int returnChrStartPos_genomicRegionInfo()
	{
		return chrStartPos;
	}

	int returnChrEndPos_genomicRegionInfo()
	{
		return chrEndPos;
	}
};


class GenomicRegionVec_Info
{
private:
	vector<GenomicRegion_Info*> genomicRegionVec;
	vector<int> genomicRegionIndexNumVec_eachChr;
	vector<int> genomicRegionIndexNumVec_eachChr_accumulative;
	int posRange_max_inGenomicRange;

public:
	GenomicRegionVec_Info()
	{
		posRange_max_inGenomicRange = 300000;
	}

	void initiate_withIndexInfo(Index_Info* indexInfo)
	{
		int accumulativeGenomicRegionIndexNum = 0;
		int accumulativeGenomicRegionIndexNum_chr = 0;
		int chrNum = indexInfo->returnChromNum();
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			int tmpChrLength = indexInfo->returnChromLength(tmp);
			int tmpChrRegionIndexNum = (tmpChrLength-1) / posRange_max_inGenomicRange + 1;
			accumulativeGenomicRegionIndexNum_chr += tmpChrRegionIndexNum;
			genomicRegionIndexNumVec_eachChr.push_back(tmpChrRegionIndexNum);
			genomicRegionIndexNumVec_eachChr_accumulative.push_back(
				accumulativeGenomicRegionIndexNum_chr);
			for(int tmp2 = 0; tmp2 < tmpChrRegionIndexNum-1; tmp2++)
			{
				int tmpChrNameInt = tmp;
				int tmpChrStartPos = tmp2 * posRange_max_inGenomicRange + 1;
				int tmpChrEndPos = (tmp2 + 1) * posRange_max_inGenomicRange; 
				GenomicRegion_Info* tmpGenomicRegionInfo
					= new GenomicRegion_Info();
				tmpGenomicRegionInfo->initiate(accumulativeGenomicRegionIndexNum,
					tmpChrNameInt, tmpChrStartPos, tmpChrEndPos);
				genomicRegionVec.push_back(tmpGenomicRegionInfo);
				accumulativeGenomicRegionIndexNum ++;
			}
			// last genomicRegionInfo in each Chr
			int tmpChrNameInt_lastInChr = tmp;
			int tmpChrStartPos_lastInChr 
				= (tmpChrRegionIndexNum-1) * posRange_max_inGenomicRange + 1;
			int tmpChrEndPos_lastInChr = tmpChrLength;
			GenomicRegion_Info* tmpGenomicRegionInfo_lastInChr
				= new GenomicRegion_Info();
			tmpGenomicRegionInfo_lastInChr->initiate(accumulativeGenomicRegionIndexNum,
				tmpChrNameInt_lastInChr, tmpChrStartPos_lastInChr, tmpChrEndPos_lastInChr);
			genomicRegionVec.push_back(tmpGenomicRegionInfo_lastInChr);
			accumulativeGenomicRegionIndexNum ++;
		}
	}

	int returnGenomicRegionIndex(int chrNameInt, int chrPos)
	{
		int tmpGenomicRegionInfoIndex_inChr 
			= (chrPos-1) / posRange_max_inGenomicRange;
		if(chrNameInt == 0)
			return tmpGenomicRegionInfoIndex_inChr;	
		else
		{
			return (tmpGenomicRegionInfoIndex_inChr 
				+ genomicRegionIndexNumVec_eachChr_accumulative[chrNameInt-1]);
		}
	}

	void generateFusionRegionIndexPairVec_unpairComplete(
		PE_Read_Alignment_Info* peAlignInfo,
		vector< pair<int,int> >& genomicRegionIndexPairVec,
		Index_Info* indexInfo)
	{
		vector< set<int> > alignInfoGenomicRegionIndexSetVec; 
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4;
			tmpAlignInfoType ++)
		{
			set<int> tmpAlignInfoGenomicRegionIndexSet;
			int tmpAlignInfoNum = peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType);
			for(int tmpAlignInfoIndex = 0; 
				tmpAlignInfoIndex < tmpAlignInfoNum; tmpAlignInfoIndex++)
			{
				string tmpAlignInfo_chrNameStr 
					= (peAlignInfo->returnAlignInfoInPeAlignInfo(tmpAlignInfoType,
						tmpAlignInfoIndex))->returnAlignChromName();
				int tmpAlignIfno_chrNameInt = indexInfo->convertStringToInt(
					tmpAlignInfo_chrNameStr);
				int tmpAlignInfo_chrMapPos 
					= (peAlignInfo->returnAlignInfoInPeAlignInfo(tmpAlignInfoType,
						tmpAlignInfoIndex))->returnAlignChromPos();
				int tmpAlignInfo_genomicRegionIndex
					= this->returnGenomicRegionIndex(tmpAlignIfno_chrNameInt,
						tmpAlignInfo_chrMapPos);
				tmpAlignInfoGenomicRegionIndexSet.insert(
					tmpAlignInfo_genomicRegionIndex);
			}
			alignInfoGenomicRegionIndexSetVec.push_back(tmpAlignInfoGenomicRegionIndexSet);
		}

		// generate genomicRegionIndexPairVec from Nor1Rcm2 pair
		for(set<int>::iterator tmpSetIter_1 = alignInfoGenomicRegionIndexSetVec[0].begin();
			tmpSetIter_1 != alignInfoGenomicRegionIndexSetVec[0].end(); tmpSetIter_1 ++)
		{
			int tmpAlignInfoGenomicRegionIndex_Nor1 = (*tmpSetIter_1);
			for(set<int>::iterator tmpSetIter_2 = alignInfoGenomicRegionIndexSetVec[3].begin();
				tmpSetIter_2 != alignInfoGenomicRegionIndexSetVec[3].end(); tmpSetIter_2 ++)			
			{
				int tmpAlignInfoGenomicRegionIndex_Rcm2 = (*tmpSetIter_2);
				genomicRegionIndexPairVec.push_back(pair<int,int>
					(tmpAlignInfoGenomicRegionIndex_Nor1, 
						tmpAlignInfoGenomicRegionIndex_Rcm2));
			}
		}

		// generate genomicRegionIndexPairVec from Nor2Rcm1 pair
		for(set<int>::iterator tmpSetIter_1 = alignInfoGenomicRegionIndexSetVec[2].begin();
			tmpSetIter_1 != alignInfoGenomicRegionIndexSetVec[2].end(); tmpSetIter_1 ++)
		{
			int tmpAlignInfoGenomicRegionIndex_Nor2 = (*tmpSetIter_1);
			for(set<int>::iterator tmpSetIter_2 = alignInfoGenomicRegionIndexSetVec[1].begin();
				tmpSetIter_2 != alignInfoGenomicRegionIndexSetVec[1].end(); tmpSetIter_2 ++)			
			{
				int tmpAlignInfoGenomicRegionIndex_Rcm1 = (*tmpSetIter_2);
				genomicRegionIndexPairVec.push_back(pair<int,int>
					(tmpAlignInfoGenomicRegionIndex_Nor2, 
						tmpAlignInfoGenomicRegionIndex_Rcm1));
			}
		}
	}

	int returnChrNameInt_genomicRegionInfoVec(int tmpIndex)
	{
		return genomicRegionVec[tmpIndex]->returnChrNameInt_genomicRegionInfo();
	}

	int returnChrStartPos_genomicRegionInfoVec(int tmpIndex)
	{
		return genomicRegionVec[tmpIndex]->returnChrStartPos_genomicRegionInfo();
	}

	int returnChrEndPos_genomicRegionInfoVec(int tmpIndex)
	{
		return genomicRegionVec[tmpIndex]->returnChrEndPos_genomicRegionInfo();
	}

	void memoryFree()
	{
		for(int tmp = 0; tmp < genomicRegionVec.size(); tmp++)
			delete genomicRegionVec[tmp];
	}

	void freeMemory()
	{
		for(int tmp = 0; tmp < genomicRegionVec.size(); tmp++)
			delete genomicRegionVec[tmp];
	}	
};

#endif