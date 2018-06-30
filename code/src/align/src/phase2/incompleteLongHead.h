// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INCOMPLETELONGHEAD_H
#define INCOMPLETELONGHEAD_H

#include <string>
#include <string.h>
//#include "splice_info.h"

using namespace std;

class Incomplete_Long_Head
{
//public:
private:
	//unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedHeadLength;
	int midPartLength;

	int midPartMapPosInChr; // map pos of 1st base in midPart
	string midPartMapChrName;
	int midPartMapChrInt;

	int secondLevelIndexNum;
	int chrPosStartIn2ndLevelIndex;

	int mapPosIntervalStart;
	int mapPosIntervalEnd;

	int incompleteHeadAndMidPartLength;

public:
	int returnUnfixedHeadLength()
	{
		return unfixedHeadLength;
	}
	int returnSecondLevelIndexNum()
	{
		return secondLevelIndexNum;
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
	string returnMidPartMapChrName()
	{
		return midPartMapChrName;
	}
	int returnMidPartLength()
	{
		return midPartLength;
	}
	int returnMidPartMapPosInChr()
	{
		return midPartMapPosInChr;
	}
	int returnMidPartMapChrInt()
	{
		return midPartMapChrInt;
	}
	int returnIncompleteHeadAndMidPartLength()
	{
		return incompleteHeadAndMidPartLength;
	}

	Incomplete_Long_Head()
	{

	}
	void getIncompleteLongHeadInfoFromRecordWithAlignInfoType_new(
		//PE_Read_Info* readInfo, int alignInfoType, 
		Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		midPartMapChrName = alignInfo->returnAlignChromName();
		midPartMapPosInChr = alignInfo->returnAlignChromPos();

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		secondLevelIndexNum = indexInfo->getSecondLevelIndexFromChrAndPos(
			midPartMapChrInt, midPartMapPosInChr); // Xinan: need to debug

		chrPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
			midPartMapChrInt, secondLevelIndexNum, 1);

		mapPosIntervalStart = midPartMapPosInChr - READ_ALIGN_AREA_LENGTH;
		if(mapPosIntervalStart <= 0)
			mapPosIntervalStart = 1;

		mapPosIntervalEnd = midPartMapPosInChr + MAX_DELETION_LENGTH;
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

		incompleteHeadAndMidPartLength = unfixedHeadLength + midPartLength;
	}

};

#endif