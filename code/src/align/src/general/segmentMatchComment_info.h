// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEGMENTMATCHCOMMENT_INFO_H
#define SEGMENTMATCHCOMMENT_INFO_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;

class SegmentMatchComment_Info
{
private:
	bool multiMapSeg_exists_bool;
	int multiMapSeg_maxLength_length;
	int multiMapSeg_maxLength_mapPosNum;
	bool multiMapSeg_maxLength_end1_or_end2_bool;
	bool multiMapSeg_maxLength_Nor_or_Rcm_bool;

public:
	SegmentMatchComment_Info()
	{
		multiMapSeg_exists_bool = false;
		multiMapSeg_maxLength_length = 0;
	}

	bool return_multiMapSeg_exists_bool()
	{
		return multiMapSeg_exists_bool;
	}

	int return_multiMapSeg_maxLength_length()
	{
		return multiMapSeg_maxLength_length;
	}

	void initiate_withSegInfo(
		bool mapBool_Nor1, Seg_Info* segInfo_Nor1,
		bool mapBool_Rcm1, Seg_Info* segInfo_Rcm1,
		bool mapBool_Nor2, Seg_Info* segInfo_Nor2,
		bool mapBool_Rcm2, Seg_Info* segInfo_Rcm2)
	{
		int tmpMultiMapSeg_maxLength_Nor1_length = 0;
		int tmpMultiMapSeg_maxLength_Rcm1_length = 0;
		int tmpMultiMapSeg_maxLength_Nor2_length = 0;
		int tmpMultiMapSeg_maxLength_Rcm2_length = 0;		

		int tmpMultiMapSeg_maxLength_Nor1_maxPosNum = 0;
		int tmpMultiMapSeg_maxLength_Rcm1_maxPosNum = 0;
		int tmpMultiMapSeg_maxLength_Nor2_maxPosNum = 0;
		int tmpMultiMapSeg_maxLength_Rcm2_maxPosNum = 0;

		checkMultiMapSeg(mapBool_Nor1, segInfo_Nor1,
			tmpMultiMapSeg_maxLength_Nor1_length,
			tmpMultiMapSeg_maxLength_Nor1_maxPosNum);
		checkMultiMapSeg(mapBool_Rcm1, segInfo_Rcm1,
			tmpMultiMapSeg_maxLength_Rcm1_length,
			tmpMultiMapSeg_maxLength_Rcm1_maxPosNum);
		checkMultiMapSeg(mapBool_Nor2, segInfo_Nor2,
			tmpMultiMapSeg_maxLength_Nor2_length,
			tmpMultiMapSeg_maxLength_Nor2_maxPosNum);
		checkMultiMapSeg(mapBool_Rcm2, segInfo_Rcm2,
			tmpMultiMapSeg_maxLength_Rcm2_length,
			tmpMultiMapSeg_maxLength_Rcm2_maxPosNum);

		if(tmpMultiMapSeg_maxLength_Nor1_length > multiMapSeg_maxLength_length)
		{
			multiMapSeg_exists_bool = true;
			multiMapSeg_maxLength_length = tmpMultiMapSeg_maxLength_Nor1_length;
			multiMapSeg_maxLength_mapPosNum = tmpMultiMapSeg_maxLength_Nor1_maxPosNum;
			multiMapSeg_maxLength_end1_or_end2_bool = true;
			multiMapSeg_maxLength_Nor_or_Rcm_bool = true;			
		}

		if(tmpMultiMapSeg_maxLength_Rcm1_length > multiMapSeg_maxLength_length)
		{
			multiMapSeg_exists_bool = true;
			multiMapSeg_maxLength_length = tmpMultiMapSeg_maxLength_Rcm1_length;
			multiMapSeg_maxLength_mapPosNum = tmpMultiMapSeg_maxLength_Rcm1_maxPosNum;
			multiMapSeg_maxLength_end1_or_end2_bool = true;
			multiMapSeg_maxLength_Nor_or_Rcm_bool = false;			
		}

		if(tmpMultiMapSeg_maxLength_Nor2_length > multiMapSeg_maxLength_length)
		{
			multiMapSeg_exists_bool = true;
			multiMapSeg_maxLength_length = tmpMultiMapSeg_maxLength_Nor2_length;
			multiMapSeg_maxLength_mapPosNum = tmpMultiMapSeg_maxLength_Nor2_maxPosNum;
			multiMapSeg_maxLength_end1_or_end2_bool = false;
			multiMapSeg_maxLength_Nor_or_Rcm_bool = true;			
		}

		if(tmpMultiMapSeg_maxLength_Rcm1_length > multiMapSeg_maxLength_length)
		{
			multiMapSeg_exists_bool = true;
			multiMapSeg_maxLength_length = tmpMultiMapSeg_maxLength_Rcm1_length;
			multiMapSeg_maxLength_mapPosNum = tmpMultiMapSeg_maxLength_Rcm1_maxPosNum;
			multiMapSeg_maxLength_end1_or_end2_bool = false;
			multiMapSeg_maxLength_Nor_or_Rcm_bool = false;			
		}		
	}

	void checkMultiMapSeg(bool mapMain_bool, Seg_Info* tmpSegInfo,
		int& tmpMultiMapSeg_maxLength_length, 
		int& tmpMultiMapSeg_maxLength_mapPosNum)
	{
		if(!mapMain_bool)
			return;
		int tmpCandiMultiMapSeg_maxLength_length = 0;
		int tmpCandiMultiMapSeg_maxLength_mapPosNum = 0;
		int tmpSegmentNum = tmpSegInfo->returnSegmentNum();
		for(int tmp = 0; tmp < tmpSegmentNum; tmp++)
		{
			int tmpSegLen = tmpSegInfo->returnSegmentLength(tmp);
			int tmpSegMapPosNum = tmpSegInfo->returnSegmentAlignNum(tmp);
			if((tmpSegLen > tmpCandiMultiMapSeg_maxLength_length)
				&&(tmpSegMapPosNum > 1))
			{
				tmpCandiMultiMapSeg_maxLength_length = tmpSegLen;
				tmpCandiMultiMapSeg_maxLength_mapPosNum = tmpSegMapPosNum;
			}
		}
		tmpMultiMapSeg_maxLength_length = tmpCandiMultiMapSeg_maxLength_length;
		tmpMultiMapSeg_maxLength_mapPosNum = tmpCandiMultiMapSeg_maxLength_mapPosNum;
	}
};

#endif