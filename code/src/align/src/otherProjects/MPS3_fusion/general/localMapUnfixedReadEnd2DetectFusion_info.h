// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef LOCALMAPUNFIXEDREADEND2DETECTFUSION_INFO_H
#define LOCALMAPUNFIXEDREADEND2DETECTFUSION_INFO_H

//#include "../../../general/localMapUnfixedReadEnd2DetectFusion_info.h"

using namespace std;

class LocalMapUnfixedReadEnd2DetectFusion_Info
{
private:
	int geneSizeMax;
public:
	LocalMapUnfixedReadEnd2DetectFusion_Info()
	{
		geneSizeMax = 500000;
	}

	bool mapSeq2genomeWithLeftOrRightMostMapPos_noUnfixedHeadOrTail(
		bool withLeftOrRightMostMapPosBool,
		bool unfixedHeadNotAllowedBool, bool unfixedTailNotAllowedBool,
		string& unfixedReadSubSeq, int leftOrRightMostMapPos, 
		int leftOrRightMostMapPos_chrNameInt,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo
		)
	{
		int unfixedReadSubSeqLength = unfixedReadSubSeq.length();
		char* unfixedReadSubSeqChar = const_cast<char*>(unfixedReadSubSeq.c_str());
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
		int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(
			leftOrRightMostMapPos_chrNameInt, leftOrRightMostMapPos) - 1;
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO + 1))
		{
			delete seg2ndOriInfo;
			return false;
		}
		bool mapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
			unfixedReadSubSeqChar,
			secondLevelSa[secondLevelIndexNO], 
			secondLevelLcpCompress[secondLevelIndexNO],
			secondLevelChildTab[secondLevelIndexNO],
			secondLevelChrom[secondLevelIndexNO], 
			secondLevelDetChild[secondLevelIndexNO],
			unfixedReadSubSeqLength, indexInfo);
		if(!mapBool)
		{
			delete seg2ndOriInfo;
			return false;
		}

		int mapInterval_mostLeft, mapInterval_mostRight;
		if(withLeftOrRightMostMapPosBool)
		{
			mapInterval_mostLeft = leftMostMapPos;
			mapInterval_mostRight = leftMostMapPos + geneSizeMax;
		}
		else
		{
			mapInterval_mostLeft = rightMostMapPos - geneSizeMax;
			mapInterval_mostRight = rightMostMapPos;
		}

		Seg_Info* newSegInfo = new Seg_Info(seg2ndOriInfo,
			mapInterval_mostLeft,
			mapInterval_mostRight,
			indexInfo->getChrPosFromSecondLevelIndexPos(
				leftOrRightMostMapPos_chrNameInt, secondLevelIndexNum + 1, 1),
			indexInfo, indexInfo->returnChrNameStr(leftOrRightMostMapPos_chrNameInt));

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);

		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg_localMapping2detectFusion(); //Xinan 10/08/15 to implement
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete newSegInfo;
			delete seg2ndOriInfo;
			return false;
		}
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			unfixedReadSubSeq, Do_extendHeadTail, annotation_provided_bool, 
			Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
		if(unfixedHeadNotAllowedBool)
			pathInfo->filterFinalPath_incompleteHead();//Xinan 10/08/15 to implement
		if(unfixedTailNotAllowedBool)
			pathInfo->filterFinalPath_incompleteTail();	
		
		bool fusionDetection_success_bool 
			= pathInfo->localMapping2detectFusion_successOrNot();//Xinan 10/08/15 to implement
		
		this->generateMappedReadEndInfoForFusionDetection();

		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete newSegInfo;
		delete seg2ndOriInfo;
		return fusionDetection_success_bool;
	}

	/*
	bool mapSeq2genomeWithLeftMostMapPos_noUnfixedTail(
		string& unfixedReadSubSeq, int leftMostMapPos, 
		int leftMostMapPos_chrNameInt,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo
		)
	{
		int unfixedReadSubSeqLength = unfixedReadSubSeq.length();
		char* unfixedReadSubSeqChar = const_cast<char*>(unfixedReadSubSeq.c_str());
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
		int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(
			leftMostMapPos_chrNameInt, leftMostMapPos) - 1;
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO + 1))
		{
			delete seg2ndOriInfo;
			return false;
		}
		bool mapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
			unfixedReadSubSeqChar,
			secondLevelSa[secondLevelIndexNO], 
			secondLevelLcpCompress[secondLevelIndexNO],
			secondLevelChildTab[secondLevelIndexNO],
			secondLevelChrom[secondLevelIndexNO], 
			secondLevelDetChild[secondLevelIndexNO],
			unfixedReadSubSeqLength, indexInfo);
		if(!mapBool)
		{
			delete seg2ndOriInfo;
			return false;
		}
		Seg_Info* newSegInfo = new Seg_Info(seg2ndOriInfo,
			leftMostMapPos, 
			leftMostMapPos + geneSizeMax,
			indexInfo->getChrPosFromSecondLevelIndexPos(
				midPartMapChrInt, secondLevelIndexNum + 1, 1),
			indexInfo, indexInfo->returnChrNameStr(leftMostMapPos_chrNameInt));

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);

		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg_localMapping2detectFusion(); //Xinan 10/08/15 to implement
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete newSegInfo;
			delete seg2ndOriInfo;
			return false;
		}
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			unfixedReadSubSeq, Do_extendHeadTail, annotation_provided_bool, 
			Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);		
		pathInfo->filterFinalPath_incompleteTail();//Xinan 10/08/15 to implement
		bool fusionDetection_success_bool 
			= pathInfo->localMapping2detectFusion_successOrNot();//Xinan 10/08/15 to implement
		
		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete newSegInfo;
		delete seg2ndOriInfo;
		return fusionDetection_success_bool;
	}*/

	bool mapSeq2genomeWithLeftMostMapPos_noUnfixedTail(
		string& unfixedReadSubSeq, int leftMostMapPos, 
		int leftMostMapPos_chrNameInt,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo
		)
	{
		bool unfixedTailNotAllowedBool = true;
		bool unfixedHeadNotAllowedBool = false;
		bool withLeftOrRightMostMapPosBool = true;
	}	

	bool mapSeq2genomeWithRightMostMapPos_noUnfixedHead(
		string& unfixedReadSubSeq, int leftMostMapPos, 
		int leftMostMapPos_chrNameInt, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo
		)
	{
		bool unfixedTailNotAllowedBool = true;
		bool unfixedHeadNotAllowedBool = false;
		bool withLeftOrRightMostMapPosBool = false;
	}
};

#endif