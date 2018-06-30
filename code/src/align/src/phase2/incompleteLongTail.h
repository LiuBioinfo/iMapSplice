// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INCOMPLETELONGTAIL_H
#define INCOMPLETELONGTAIL_H

#include <string>
#include <string.h>
//#include "splice_info.h"

using namespace std;

class Incomplete_Long_Tail
{
//public:
private:
	//unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedTailLength;
	int midPartLength;

	int midPartEndMapPosInChr; // map pos of 1st base in midPart
	

	string midPartMapChrName;
	
	int midPartMapChrInt;

	int secondLevelIndexNum;
	int chrPosStartIn2ndLevelIndex;

	int mapPosIntervalStart;
	int mapPosIntervalEnd;

	int incompleteTailAndMidPartLength;

	int midPartMapPosInChr;

	int readLength;

	int midPartLocInRead;

public:
	int returnIncompleteTailAndMidPartLength()
	{
		return incompleteTailAndMidPartLength;
	}
	int returnMidPartLength()
	{
		return midPartLength;
	}
	int returnSecondLevelIndexNum()
	{
		return secondLevelIndexNum;
	}
	string returnMidPartMapChrName()
	{
		return midPartMapChrName;
	}
	int returnMidPartMapPosInChr()
	{
		return midPartMapPosInChr;
	}
	int returnMidPartMapChrInt()
	{
		return midPartMapChrInt;
	}
	int returnMidPartLocInRead()
	{
		return midPartLocInRead;
	}
	int returnMapPosIntervalStart()
	{
		return mapPosIntervalStart;
	}
	int returnMapPosIntervalEnd()
	{
		return mapPosIntervalEnd;
	}
	int returnChrPosStartIn2ndLevelIndex()
	{
		return chrPosStartIn2ndLevelIndex;
	}
	string returnOtherCigarString()
	{
		return otherCigarString;
	}
	int returnUnfixedTailLength()
	{
		return unfixedTailLength;
	}
	Incomplete_Long_Tail()
	{

	}

	void getIncompleteLongTailInfoFromRecordWithAlignInfoType_new(
		PE_Read_Info& peReadInfo, int alignInfoType, 
		Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		if(alignInfoType <= 2)
			readLength = peReadInfo.returnReadSeqLength_1();
		else
			readLength = peReadInfo.returnReadSeqLength_2();

		midPartMapChrName = alignInfo->returnAlignChromName();
		midPartEndMapPosInChr = alignInfo->getEndMatchedPosInChr();
		

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();

		int cigarStringJumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[cigarStringJumpCodeVecSize-1].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[cigarStringJumpCodeVecSize-2].len;

		secondLevelIndexNum = indexInfo->getSecondLevelIndexFromChrAndPos(
			midPartMapChrInt, midPartEndMapPosInChr); // Xinan: need to debug

		chrPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
			midPartMapChrInt, secondLevelIndexNum, 1);

		mapPosIntervalStart = midPartEndMapPosInChr - MAX_INSERTION_LENGTH; 
		if(mapPosIntervalStart <= 0)
			mapPosIntervalStart = 1;

		mapPosIntervalEnd = midPartEndMapPosInChr + READ_ALIGN_AREA_LENGTH;

		if(mapPosIntervalEnd > ((
			//(indexInfo->chromLength)[midPartMapChrInt]
			indexInfo->returnChromLength(midPartMapChrInt)
			) - READ_LENGTH_MAX) ) //Xinan: to debug: 1000 should be the maximum read length
		{

			mapPosIntervalEnd = (
				//(indexInfo->chromLength)[midPartMapChrInt]
				indexInfo->returnChromLength(midPartMapChrInt)
				) - READ_LENGTH_MAX;
		}

		incompleteTailAndMidPartLength = unfixedTailLength + midPartLength;

		midPartMapPosInChr = midPartEndMapPosInChr - midPartLength + 1;

		midPartLocInRead = readLength - incompleteTailAndMidPartLength + 1;
	}

};

#endif