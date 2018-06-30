// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GROUPSEG_INFO_H
#define GROUPSEG_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include "enhanced_suffix_array_info.h"
#include "local_seg_info.h"
#include "missingLongSeg_info.h"

using namespace std;

class SegSubGroup_PE_Info
{
private:
	unsigned int wholeGenomePos_intervalStart;
	unsigned int wholeGenomePos_intervalEnd;

	Seg_Info* subSegInfo_Nor;
	Seg_Info* subSegInfo_Rcm;

	// sum of found seg length
	// int segLengthSum_Nor;
	// int segLengthSum_Rcm;
	// int segLengthSum;
	bool pairedSegFoundBool;

	Path_Info subPathInfo_Nor;
	Path_Info subPathInfo_Rcm;

	PE_Read_Alignment_Info peAlignInfo;
public:
	SegSubGroup_PE_Info()
	{
		subSegInfo_Nor = new Seg_Info();
		subSegInfo_Rcm = new Seg_Info();
		pairedSegFoundBool = false;
	}

	void insertFinalPairAlignInfo2targetPeAlignInfo_Nor1Rcm2Only(
		PE_Read_Alignment_Info& targetPeReadAlignInfo)
	{
		peAlignInfo.insertNor1Rcm2FinalPairAlignInfo2targetNor1Rcm2FinalPairAlignInfo(
			targetPeReadAlignInfo);
	}

	void insertFinalPairAlignInfo2targetPeAlignInfo_Nor2Rcm1Only(
		PE_Read_Alignment_Info& targetPeReadAlignInfo)
	{
		peAlignInfo.insertNor1Rcm2FinalPairAlignInfo2targetNor2Rcm1FinalPairAlignInfo(
			targetPeReadAlignInfo);
	}

	void initiatePeAlignInfo_OneDirOnly(Index_Info* indexInfo) // Nor1 -- Rcm2
	{
		peAlignInfo.initiatePeAlignInfo_Nor1Rcm2Only(subPathInfo_Nor, subPathInfo_Rcm, indexInfo);
	}

	void alignmentFilter_fixPhase1_SJpenalty_OneDirOnly(
		int readLength_1, int readLength_2)
	{
		peAlignInfo.alignmentFilter_fixPhase1_SJpenalty_Nor1Rcm2Only(
			readLength_1, readLength_2);
	}

	bool pairedOrNot_bool()
	{
		return pairedSegFoundBool;
	}

	string subSegInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to output subSegInfoStr in SegSubGroup_PE_Info" << endl;
		string tmpSegInfoStr = "SubSegInfo:\nNor:";
		//cout << "Nor: " << endl;
		tmpSegInfoStr = tmpSegInfoStr + subSegInfo_Nor->segInfoStr(indexInfo);
		tmpSegInfoStr = tmpSegInfoStr + "\nRcm:";
		//cout << "Rcm: " << endl;
		tmpSegInfoStr = tmpSegInfoStr + subSegInfo_Rcm->segInfoStr(indexInfo);
		//tmpSegInfoStr += "\n";
		return tmpSegInfoStr;
	}

	string subPathInfoStr()
	{
		string tmpPathInfoStr = "SubPathInfo:\nNor:";
		tmpPathInfoStr = tmpPathInfoStr + subPathInfo_Nor.possiPathStr();
		tmpPathInfoStr = tmpPathInfoStr + "\nRcm:";
		tmpPathInfoStr = tmpPathInfoStr + subPathInfo_Rcm.possiPathStr();
		return tmpPathInfoStr;
	}

	string subGapInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to generate subGapInfoStr ..." << endl;
		string tmpGapInfoStr = "SubGapInfo:\nNor:";
		//cout << "subGapInfo--Nor:" << endl;
		tmpGapInfoStr = tmpGapInfoStr 
			+ subPathInfo_Nor.fixedPathVecStr(indexInfo, subSegInfo_Nor);
		//cout << "fixedPathVecStr_Nor: " << endl << subPathInfo_Nor.fixedPathVecStr(indexInfo, subSegInfo_Nor) << endl;
		tmpGapInfoStr = tmpGapInfoStr
			+ subPathInfo_Nor.finalFixedPathStr(indexInfo);
		//cout << "finalFixedPathStr_Nor: " << endl << subPathInfo_Nor.finalFixedPathStr(indexInfo) << endl;
		tmpGapInfoStr = tmpGapInfoStr
			+ subPathInfo_Nor.getFixGapInfoStr();
		//cout << "getFixGapInfoStr_Nor: " << endl << subPathInfo_Nor.getFixGapInfoStr() << endl;
		tmpGapInfoStr = tmpGapInfoStr + "\nRcm:";
		//cout << "subGapInfo--Rcm:" << endl;
		tmpGapInfoStr = tmpGapInfoStr 
			+ subPathInfo_Rcm.fixedPathVecStr(indexInfo, subSegInfo_Rcm);
		//cout << "fixedPathVecStr_Rcm: " << endl << subPathInfo_Rcm.fixedPathVecStr(indexInfo, subSegInfo_Rcm) << endl;	
		tmpGapInfoStr = tmpGapInfoStr
			+ subPathInfo_Rcm.finalFixedPathStr(indexInfo);
		//cout << "finalFixedPathStr_Rcm: " << endl << subPathInfo_Rcm.finalFixedPathStr(indexInfo) << endl;
		tmpGapInfoStr = tmpGapInfoStr
			+ subPathInfo_Rcm.getFixGapInfoStr();
		//cout << "getFixGapInfoStr_Rcm: " << endl << subPathInfo_Rcm.getFixGapInfoStr() << endl;
		return tmpGapInfoStr;
	}

	void initiateSegSubGroupPEinfo_withNorSegRcmSeg(
		Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm, 
		unsigned int tmpSegMapPos, int tmpSegGroupIndex,
		unsigned int tmpSegMapPos_Rcm, int tmpSegGroupIndex_Rcm)
	{

	}

	void addRcmSeg2SegSubGroup()
	{}

	void initiateSegSubGroupPEinfo_withNorSeg(
		Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm,
		unsigned int firstSegMapPos_Nor, int firstSegGroupIndex_Nor)
	{
		// assign wholeGenomePos_interval
		if(firstSegMapPos_Nor > READ_ALIGN_AREA_LENGTH)
			wholeGenomePos_intervalStart = firstSegMapPos_Nor - READ_ALIGN_AREA_LENGTH;
		else
			wholeGenomePos_intervalStart = 1;
		wholeGenomePos_intervalEnd = firstSegMapPos_Nor + READ_ALIGN_AREA_LENGTH;

		// initiate segAlignNum_subGroup_Nor/Rcm
		int segInfo_segmentNum_Nor = segInfo_Nor->returnSegmentNum();
		int segInfo_segmentNum_Rcm = segInfo_Rcm->returnSegmentNum();
		subSegInfo_Nor->assignSegmentNum(segInfo_segmentNum_Nor);
		subSegInfo_Rcm->assignSegmentNum(segInfo_segmentNum_Rcm);
		subSegInfo_Nor->copySegmentLength(segInfo_Nor);
		subSegInfo_Rcm->copySegmentLength(segInfo_Rcm);
		subSegInfo_Nor->copySegmentLocInRead(segInfo_Nor);
		subSegInfo_Rcm->copySegmentLocInRead(segInfo_Rcm);
		subSegInfo_Nor->assignZeroAlignNumToAllSeg();
		subSegInfo_Rcm->assignZeroAlignNumToAllSeg();

		// initiate first seg in this subGroup -- Nor
		subSegInfo_Nor->assignSegmentAlignNum(firstSegGroupIndex_Nor, 1);
		subSegInfo_Nor->assignSegmentAlignLoc(firstSegGroupIndex_Nor, 0, firstSegMapPos_Nor);
		//segAlignNum_subGroup_Nor[firstSegGroupIndex_Nor] = 1;
		//segAlignLoc_subGroup_Nor[firstSegGroupIndex_Nor * CANDALILOC] 
		//	= firstSegMapPos_Nor;
	}

	bool checkNewNorSegMapPos_addIfFit(int tmpSegGroupIndex, unsigned int tmpSegMapPos)
	{
		if((tmpSegMapPos < wholeGenomePos_intervalStart)||
			(tmpSegMapPos > wholeGenomePos_intervalEnd))
			return false;
		else
		{
			// if(Nor_or_Rcm_bool) // Nor
			// {
				int tmpCurrentSegAlignNum 
					= subSegInfo_Nor->returnSegmentAlignNum(tmpSegGroupIndex);
				subSegInfo_Nor->assignSegmentAlignLoc(tmpSegGroupIndex, tmpCurrentSegAlignNum,
					tmpSegMapPos);
				subSegInfo_Nor->assignSegmentAlignNum(tmpSegGroupIndex, tmpCurrentSegAlignNum + 1);
			// }
			// else // Rcm
			// {
			// 	pairedSegFoundBool = true;
			// 	int tmpCurrentSegAlignNum 
			// 		= subSegInfo_Rcm->returnSegmentAlignNum(tmpSegGroupIndex);
			// 	subSegInfo_Rcm->assignSegmentAlignLoc(tmpSegGroupIndex, tmpCurrentSegAlignNum,
			// 		tmpSegMapPos);
			// 	subSegInfo_Rcm->assignSegmentAlignNum(tmpSegGroupIndex, tmpCurrentSegAlignNum + 1);				
			// }
			return true;
		}
	}	
	bool checkNewSegMapPos_addIfFit(bool Nor_or_Rcm_bool,
		int tmpSegGroupIndex, unsigned int tmpSegMapPos)
	{
		if((tmpSegMapPos < wholeGenomePos_intervalStart)||
			(tmpSegMapPos > wholeGenomePos_intervalEnd))
			return false;
		else
		{
			if(Nor_or_Rcm_bool) // Nor
			{
				int tmpCurrentSegAlignNum 
					= subSegInfo_Nor->returnSegmentAlignNum(tmpSegGroupIndex);
				subSegInfo_Nor->assignSegmentAlignLoc(tmpSegGroupIndex, tmpCurrentSegAlignNum,
					tmpSegMapPos);
				subSegInfo_Nor->assignSegmentAlignNum(tmpSegGroupIndex, tmpCurrentSegAlignNum + 1);
			}
			else // Rcm
			{
				pairedSegFoundBool = true;
				int tmpCurrentSegAlignNum 
					= subSegInfo_Rcm->returnSegmentAlignNum(tmpSegGroupIndex);
				subSegInfo_Rcm->assignSegmentAlignLoc(tmpSegGroupIndex, tmpCurrentSegAlignNum,
					tmpSegMapPos);
				subSegInfo_Rcm->assignSegmentAlignNum(tmpSegGroupIndex, tmpCurrentSegAlignNum + 1);				
			}
			return true;
		}
	}	

	// void getSegLengthSum_checkPairOrNot(Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm)
	// {
	// 	int segNum_Nor = segInfo_Nor->returnSegmentNum();
	// 	int segNum_Rcm = segInfo_Rcm->returnSegmentNum();
	// 	for(int tmp = 0; tmp < segNum_Nor; tmp++)
	// 	{
	// 		int tmpSegAlignNum = subSegInfo_Nor->returnSegmentAlignNum(tmp); //segAlignNum_subGroup_Nor[tmp];
	// 		if((tmpSegAlignNum > 0)&&(tmpSegAlignNum <= CANDALILOC))
	// 		{
	// 			int tmpSegLength = subSegInfo_Nor->returnSegmentLength(tmp);
	// 			segLengthSum_Nor += tmpSegLength;
	// 		}
	// 	}
	// 	for(int tmp = 0; tmp < segNum_Rcm; tmp++)
	// 	{
	// 		int tmpSegAlignNum = subSegInfo_Rcm->returnSegmentAlignNum(tmp);
	// 		if((tmpSegAlignNum > 0)&&(tmpSegAlignNum <= CANDALILOC))
	// 		{
	// 			int tmpSegLength = subSegInfo_Rcm->returnSegmentLength(tmp);
	// 			segLengthSum_Rcm += tmpSegLength;
	// 		}
	// 	}

	// 	segLengthSum = segLengthSum_Nor + segLengthSum_Rcm;
	// 	if((segLengthSum_Nor > 0)&&(segLengthSum_Rcm > 0))
	// 		pairedSegFoundBool = true;
	// 	else
	// 		pairedSegFoundBool = false;
	// }

	void getPossiPathFromSeg_fixGapAlso(Index_Info* indexInfo,
		const string& readSeq_Nor, const string& readSeq_Rcm,
		bool annotation_provided_bool,
		bool Do_annotation_only_bool,
		Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		subPathInfo_Nor.getPossiPathFromSeg_fixGapAlso(
			subSegInfo_Nor, indexInfo, readSeq_Nor, 
			annotation_provided_bool,
			Do_annotation_only_bool, annotationInfo, 
			spliceJunctionDistanceMax);
		subPathInfo_Rcm.getPossiPathFromSeg_fixGapAlso(
			subSegInfo_Rcm, indexInfo, readSeq_Rcm,
			annotation_provided_bool,
			Do_annotation_only_bool, annotationInfo,
			spliceJunctionDistanceMax); 
	}

	void fixAllPath_fixGapAlso(Index_Info* indexInfo,
		const string& readSeq_Nor, const string& readSeq_Rcm,
		bool Do_extendHeadTail)
	{
		//cout << "start to do fixAllPath_fixGapAlso " << endl;
		//cout << "start to do fixAllPath_fixGapAlso for Nor" << endl;
		subPathInfo_Nor.fixAllPath_fixGapAlso(
			subSegInfo_Nor, indexInfo, readSeq_Nor, 
			Do_extendHeadTail);
		//cout << "start to do fixAllPath_fixGapAlso for Rcm " << endl;
		subPathInfo_Rcm.fixAllPath_fixGapAlso(
			subSegInfo_Rcm, indexInfo, readSeq_Rcm,
			Do_extendHeadTail);
	}

	void freeMemory()
	{
		peAlignInfo.memoryFree();
		delete subSegInfo_Nor;
		delete subSegInfo_Rcm;
	}
};

class GroupSeg_Info_SingleDir
{
private:
	vector<SegSubGroup_PE_Info*> segSubGroupInfoVec;
public:
	GroupSeg_Info_SingleDir()
	{}

	void mergeSubGroupPeAlignInfo2targetPeAlignInfo_allSegSubGroup_Nor1Rcm2Only(
		PE_Read_Alignment_Info& targetPeReadAlignInfo)
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->insertFinalPairAlignInfo2targetPeAlignInfo_Nor1Rcm2Only(
					targetPeReadAlignInfo);
			}
		}
	}

	void mergeSubGroupPeAlignInfo2targetPeAlignInfo_allSegSubGroup_Nor2Rcm1Only(
		PE_Read_Alignment_Info& targetPeReadAlignInfo)
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->insertFinalPairAlignInfo2targetPeAlignInfo_Nor2Rcm1Only(
					targetPeReadAlignInfo);
			}
		}
	}

	string returnSubSegGroup_segInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to output subSegGroup_segInfoStr in GroupSeg_Info_SingleDir" << endl;
		string tmpStr = "SubSegInfo Vec Str:\n";
		//cout << "segSubGroupInfoVecSize: " << segSubGroupInfoVec.size() << endl;
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			tmpStr = tmpStr + "subGroup -- " + int_to_str(tmp+1) + "\n";
			//cout << "subGroup: " << tmp+1 << endl;
			tmpStr = tmpStr + segSubGroupInfoVec[tmp]->subSegInfoStr(indexInfo);
			tmpStr = tmpStr + "\n";
		}
		return tmpStr;
	}

	string returnSubSegGroup_pathInfoStr()
	{
		string tmpStr = "SubPathInfo Vec Str:\n";
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			tmpStr = tmpStr + "subGroup -- " + int_to_str(tmp+1) + "\n";
			tmpStr = tmpStr + segSubGroupInfoVec[tmp]->subPathInfoStr();
			tmpStr = tmpStr + "\n";
		}
		return tmpStr;
	}

	string returnSubSegGroup_gapInfoStr(Index_Info* indexInfo)
	{
		string tmpStr = "SubGapInfo Vec Str:\n";
		//cout << "subGapInfo Vec Str: " << endl;
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			//cout << "segSubGroupInfoVec index: tmp -- " << tmp << endl; 
			tmpStr = tmpStr + "subGroup -- " + int_to_str(tmp+1) + "\n";
			tmpStr = tmpStr + segSubGroupInfoVec[tmp]->subGapInfoStr(indexInfo);
			tmpStr = tmpStr + "\n";
		}
		return tmpStr;
	}	

	void initiateGroupSegInfo_withPeSegInfo(Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm)
	{
		int segNum_Nor = segInfo_Nor->returnSegmentNum();
		int segNum_Rcm = segInfo_Rcm->returnSegmentNum();
		for(int tmpSegGroupIndex = 0; tmpSegGroupIndex < segNum_Nor; tmpSegGroupIndex++)
		{
			int tmpSegGroup_alignNum = segInfo_Nor->returnSegmentAlignNum(tmpSegGroupIndex);
			int tmpSegGroup_length = segInfo_Nor->returnSegmentLength(tmpSegGroupIndex);
			if((tmpSegGroup_alignNum <= CANDALILOC)
				&&(tmpSegGroup_length >= LONG_SEG_LENGTH_THRESHOLD_PHASE1))
			{	
				for(int tmpSegCandiIndex = 0; tmpSegCandiIndex < tmpSegGroup_alignNum; tmpSegCandiIndex ++)
				{
					bool tmpSegMatchSomeSubGroup_bool = false;
					unsigned int tmpSegMapPos = segInfo_Nor->returnSegmentMapPos(tmpSegGroupIndex, tmpSegCandiIndex);
					int currentSegSubGroupNum = segSubGroupInfoVec.size();
					for(int tmpSubGroupIndex = 0; tmpSubGroupIndex < currentSegSubGroupNum; tmpSubGroupIndex++)
					{
						bool tmpSegMatchThisSubGroup_bool 
							= segSubGroupInfoVec[tmpSubGroupIndex]->checkNewSegMapPos_addIfFit(true,
								tmpSegGroupIndex, tmpSegMapPos);
						if(tmpSegMatchThisSubGroup_bool)
							tmpSegMatchSomeSubGroup_bool = true;
					}
					if(!tmpSegMatchSomeSubGroup_bool)
					{
						SegSubGroup_PE_Info* newSegSubGroupPEinfo = new SegSubGroup_PE_Info();
						newSegSubGroupPEinfo->initiateSegSubGroupPEinfo_withNorSeg(
							segInfo_Nor, segInfo_Rcm, tmpSegMapPos, tmpSegGroupIndex);
						segSubGroupInfoVec.push_back(newSegSubGroupPEinfo);
					}
				}
			}
		}
		for(int tmpSegGroupIndex = 0; tmpSegGroupIndex < segNum_Rcm; tmpSegGroupIndex++)
		{
			int tmpSegGroup_alignNum = segInfo_Rcm->returnSegmentAlignNum(tmpSegGroupIndex);
			int tmpSegGroup_length = segInfo_Rcm->returnSegmentLength(tmpSegGroupIndex);
			if((tmpSegGroup_alignNum <= CANDALILOC)
				&&(tmpSegGroup_length >= LONG_SEG_LENGTH_THRESHOLD_PHASE1))
			{	
				for(int tmpSegCandiIndex = 0; tmpSegCandiIndex < tmpSegGroup_alignNum; tmpSegCandiIndex ++)
				{
					bool tmpSegMatchSomeSubGroup_bool = false;
					unsigned int tmpSegMapPos = segInfo_Rcm->returnSegmentMapPos(tmpSegGroupIndex, tmpSegCandiIndex);
					int currentSegSubGroupNum = segSubGroupInfoVec.size();
					for(int tmpSubGroupIndex = 0; tmpSubGroupIndex < currentSegSubGroupNum; tmpSubGroupIndex++)
					{
						bool tmpSegMatchThisSubGroup_bool 
							= segSubGroupInfoVec[tmpSubGroupIndex]->checkNewSegMapPos_addIfFit(false,
								tmpSegGroupIndex, tmpSegMapPos);
						if(tmpSegMatchThisSubGroup_bool)
							tmpSegMatchSomeSubGroup_bool = true;
					}
					//if(!tmpSegMatchSomeSubGroup_bool)
					//{
						//SegSubGroup_PE_Info* newSegSubGroupPEinfo = new SegSubGroup_PE_Info();
						//newSegSubGroupPEinfo->initiateSegSubGroupPEinfo_withRcmSeg(
						//	segInfo_Nor, segInfo_Rcm, tmpSegMapPos, tmpSegGroupIndex);
						//segSubGroupInfoVec.push_back(newSegSubGroupPEinfo);
					//}
				}
			}
		}
	}

	void initiateGroupSegInfo_withPeSegInfo_onlyPaired(Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm)
	{
		int segNum_Nor = segInfo_Nor->returnSegmentNum();
		int segNum_Rcm = segInfo_Rcm->returnSegmentNum();
		for(int tmpSegGroupIndex = 0; tmpSegGroupIndex < segNum_Nor; tmpSegGroupIndex++)
		{
			int tmpSegGroup_alignNum = segInfo_Nor->returnSegmentAlignNum(tmpSegGroupIndex);
			int tmpSegGroup_length = segInfo_Nor->returnSegmentLength(tmpSegGroupIndex);
			if((tmpSegGroup_alignNum <= CANDALILOC)
				&&(tmpSegGroup_length >= LONG_SEG_LENGTH_THRESHOLD_PHASE1))
			{	
				for(int tmpSegCandiIndex = 0; tmpSegCandiIndex < tmpSegGroup_alignNum; tmpSegCandiIndex ++)
				{
					bool tmpSegMatchSomeSubGroup_bool = false;
					unsigned int tmpSegMapPos = segInfo_Nor->returnSegmentMapPos(tmpSegGroupIndex, tmpSegCandiIndex);
					int currentSegSubGroupNum = segSubGroupInfoVec.size();
					for(int tmpSubGroupIndex = 0; tmpSubGroupIndex < currentSegSubGroupNum; tmpSubGroupIndex++)
					{
						bool tmpSegMatchThisSubGroup_bool 
							= segSubGroupInfoVec[tmpSubGroupIndex]->checkNewNorSegMapPos_addIfFit(//true,
								tmpSegGroupIndex, tmpSegMapPos);
						if(tmpSegMatchThisSubGroup_bool)
							tmpSegMatchSomeSubGroup_bool = true;
					}
					if(!tmpSegMatchSomeSubGroup_bool)
					{
						// new Nor seg, candidate to initiate a new segSubGroup,
						// check Rcm seg, if found a paired one, then initiate the new segSubGroup,
						// and also continue search for all Rcm segs that can fit this new segSubGroup
						unsigned int tmpSegMapPos_intervalStart, tmpSegMapPos_intervalEnd;
						if(tmpSegMapPos > READ_ALIGN_AREA_LENGTH)
							tmpSegMapPos_intervalStart = tmpSegMapPos - READ_ALIGN_AREA_LENGTH;
						else
							tmpSegMapPos_intervalStart = 1;
						tmpSegMapPos_intervalEnd = tmpSegMapPos + READ_ALIGN_AREA_LENGTH;						
						bool newSegSubGroupAlreadyInitiatedBool = false;
						int newSegSubGroupIndex = -1;
						for(int tmpSegGroupIndex_Rcm = 0; tmpSegGroupIndex_Rcm < segNum_Rcm; tmpSegGroupIndex_Rcm++)
						{
							int tmpSegGroup_alignNum_Rcm = segInfo_Rcm->returnSegmentAlignNum(tmpSegGroupIndex_Rcm);
							int tmpSegGroup_length_Rcm = segInfo_Rcm->returnSegmentLength(tmpSegGroupIndex_Rcm);						
							if((tmpSegGroup_alignNum_Rcm <= CANDALILOC)
								&&(tmpSegGroup_length_Rcm >= LONG_SEG_LENGTH_THRESHOLD_PHASE1))							
							{
								for(int tmpSegCandiIndex_Rcm = 0; tmpSegCandiIndex_Rcm < tmpSegGroup_alignNum_Rcm; 
									tmpSegCandiIndex_Rcm ++)
								{	
									unsigned int tmpSegMapPos_Rcm 
										= segInfo_Rcm->returnSegmentMapPos(tmpSegGroupIndex_Rcm, tmpSegCandiIndex_Rcm);
									if((tmpSegMapPos_Rcm >= tmpSegMapPos_intervalStart)
										&&(tmpSegMapPos_Rcm <= tmpSegMapPos_intervalEnd))
									{
										if(newSegSubGroupAlreadyInitiatedBool) // newSegSubGroupAlreadyInititated
										{
											SegSubGroup_PE_Info* newSegSubGroupPEinfo = new SegSubGroup_PE_Info();
											newSegSubGroupPEinfo->initiateSegSubGroupPEinfo_withNorSegRcmSeg(
												segInfo_Nor, segInfo_Rcm, tmpSegMapPos, tmpSegGroupIndex,
												tmpSegMapPos_Rcm, tmpSegGroupIndex_Rcm);
											segSubGroupInfoVec.push_back(newSegSubGroupPEinfo);
											newSegSubGroupIndex = segSubGroupInfoVec.size()-1;
										}
										else // newSegSubGroup has not been inititated, newSegSubGroup to be initiated
										{
											segSubGroupInfoVec[newSegSubGroupIndex]->addRcmSeg2SegSubGroup();
										}
									}
								}							
							}
						}


					}
				}
			}
		}
	}

	void getPossiPathFromSeg_fixGapAlso_allSegSubGroup(Index_Info* indexInfo,
		const string& readSeq_Nor, const string& readSeq_Rcm,
		bool annotation_provided_bool, bool Do_annotation_only_bool,
		Annotation_Info* annotationInfo, int spliceJunctionDistanceMax)// Note: only fix paired subgroups
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->getPossiPathFromSeg_fixGapAlso(indexInfo,
					readSeq_Nor, readSeq_Rcm, annotation_provided_bool,
					Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
			}
		}
	}

	void fixAllPath_fixGapAlso_allSegSubGroup(Index_Info* indexInfo,
		const string& readSeq_Nor, const string& readSeq_Rcm,
		bool Do_extendHeadTail) // Note: only fix paired subgroups
	{
		//cout << "start to do fixAllPath_fixGapAlso_allSegSubGroup" << endl;
		//cout << "segSubGroupInfoVecSize: " << segSubGroupInfoVec.size() << endl;
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{		
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->fixAllPath_fixGapAlso(
					indexInfo, readSeq_Nor, readSeq_Rcm, Do_extendHeadTail);
			}
		}	
	}

	void initiatePeAlignInfo_allSegSubGroup_OneDirOnly(Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{		
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->initiatePeAlignInfo_OneDirOnly(
					indexInfo);
			}
		}			
	}

	void alignmentFilter_fixPhase1_SJpenalty_allSegSubGroup_OneDirOnly(
		int readLength_1, int readLength_2)
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{		
			bool tmpSegSubGroup_pairedOrNot_bool 
				= segSubGroupInfoVec[tmp]->pairedOrNot_bool();
			if(tmpSegSubGroup_pairedOrNot_bool)
			{	
				segSubGroupInfoVec[tmp]->alignmentFilter_fixPhase1_SJpenalty_OneDirOnly(
					readLength_1, readLength_2);
			}
		}				
	}

	void freeMemory()
	{
		for(int tmp = 0; tmp < segSubGroupInfoVec.size(); tmp++)
		{
			segSubGroupInfoVec[tmp]->freeMemory();
			delete segSubGroupInfoVec[tmp];
		}
	}
};

class GroupSeg_Info_BothDir
{
private:
	GroupSeg_Info_SingleDir groupSegInfo_Nor1Rcm2;
	GroupSeg_Info_SingleDir groupSegInfo_Nor2Rcm1;

	//PE_Read_Alignment_Info peReadAlignInfo;
public:
	GroupSeg_Info_BothDir()
	{}

	string returnBothDirSubSegGroup_segInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to output BothDirSubSegGroup_segInfo" << endl;
		string tmpStr = "Nor1Rcm2 groupSegInfo: \n";
		//cout << "start to output Nor1Rcm2 groupSegInfo" << endl;
		tmpStr = tmpStr 
			+ groupSegInfo_Nor1Rcm2.returnSubSegGroup_segInfoStr(indexInfo) + "\n";
		tmpStr = tmpStr + "Nor2Rcm1 groupSegInfo: \n";
		//cout << "start to output Nor2Rcm1 groupSegInfo" << endl;
		tmpStr = tmpStr
			+ groupSegInfo_Nor2Rcm1.returnSubSegGroup_segInfoStr(indexInfo);
		return tmpStr;
	}

	string returnBothDirSubSegGroup_pathInfoStr()
	{
		string tmpStr = "Nor1Rcm2 groupPathInfo: \n";
		tmpStr = tmpStr 
			+ groupSegInfo_Nor1Rcm2.returnSubSegGroup_pathInfoStr() + "\n";
		tmpStr = tmpStr + "Nor2Rcm1 groupPathInfo: \n";
		tmpStr = tmpStr
			+ groupSegInfo_Nor2Rcm1.returnSubSegGroup_pathInfoStr();
		return tmpStr;		
	}

	string returnBothDirSubSegGroup_gapInfoStr(Index_Info* indexInfo)
	{
		//cout << "Nor1Rcm2 group gap info " << endl;
		string tmpStr = "\nNor1Rcm2 groupGapInfo: \n";
		tmpStr = tmpStr 
			+ groupSegInfo_Nor1Rcm2.returnSubSegGroup_gapInfoStr(indexInfo) + "\n";
		//cout << "Nor2Rcm1 group gap info " << endl; 
		tmpStr = tmpStr + "\nNor2Rcm1 groupGapInfo: \n";
		tmpStr = tmpStr
			+ groupSegInfo_Nor2Rcm1.returnSubSegGroup_gapInfoStr(indexInfo);
		//cout << "end of generating gapInfoStr " << endl;
		return tmpStr;		
	}

	void initiateGroupSegInfo_withPeSegInfo_BothDir(
		Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2)
	{
		groupSegInfo_Nor1Rcm2.initiateGroupSegInfo_withPeSegInfo(
			segInfo_Nor1, segInfo_Rcm2);
		groupSegInfo_Nor2Rcm1.initiateGroupSegInfo_withPeSegInfo(
			segInfo_Nor2, segInfo_Rcm1);
	}

	void getPossiPathFromSeg_fixGapAlso_allSegSubGroup_BothDir(
		Index_Info* indexInfo, PE_Read_Info& peReadInfo,
		bool annotation_provided_bool, bool Do_annotation_only_bool,
		Annotation_Info* annotationInfo, int spliceJunctionDistanceMax)
	{
		groupSegInfo_Nor1Rcm2.getPossiPathFromSeg_fixGapAlso_allSegSubGroup(
			indexInfo, peReadInfo.returnReadSeq_1(), peReadInfo.returnRcmReadSeq_2(),
			annotation_provided_bool, Do_annotation_only_bool,
			annotationInfo, spliceJunctionDistanceMax);
		groupSegInfo_Nor2Rcm1.getPossiPathFromSeg_fixGapAlso_allSegSubGroup(
			indexInfo, peReadInfo.returnReadSeq_2(), peReadInfo.returnRcmReadSeq_1(),
			annotation_provided_bool, Do_annotation_only_bool,
			annotationInfo, spliceJunctionDistanceMax);
	}

	void fixAllPath_fixGapAlso_allSegSubGroup_BothDir(
		Index_Info* indexInfo, PE_Read_Info& peReadInfo, 
		bool Do_extendHeadTail)
	{
		groupSegInfo_Nor1Rcm2.fixAllPath_fixGapAlso_allSegSubGroup(
			indexInfo, peReadInfo.returnReadSeq_1(), peReadInfo.returnRcmReadSeq_2(),
			Do_extendHeadTail);
		groupSegInfo_Nor2Rcm1.fixAllPath_fixGapAlso_allSegSubGroup(
			indexInfo, peReadInfo.returnReadSeq_2(), peReadInfo.returnRcmReadSeq_1(),
			Do_extendHeadTail);		
	}

	void initiatePeAlignInfo_allSegSubGroup_BothDir(Index_Info* indexInfo)
	{
		groupSegInfo_Nor1Rcm2.initiatePeAlignInfo_allSegSubGroup_OneDirOnly(indexInfo);
		groupSegInfo_Nor2Rcm1.initiatePeAlignInfo_allSegSubGroup_OneDirOnly(indexInfo);
	}

	void alignmentFilter_fixPhase1_SJpenalty_allSegSubGroup_BothDir(
		int readLength_1, int readLength_2)
	{
		groupSegInfo_Nor1Rcm2.alignmentFilter_fixPhase1_SJpenalty_allSegSubGroup_OneDirOnly(
			readLength_1, readLength_2);
		groupSegInfo_Nor2Rcm1.alignmentFilter_fixPhase1_SJpenalty_allSegSubGroup_OneDirOnly(
			readLength_2, readLength_1);
	}

	void mergeSubGroupPeAlignInfo_allSegSubGroup_BothDir(PE_Read_Alignment_Info& peReadAlignInfo)
	{
		groupSegInfo_Nor1Rcm2.mergeSubGroupPeAlignInfo2targetPeAlignInfo_allSegSubGroup_Nor1Rcm2Only(
			peReadAlignInfo);

		groupSegInfo_Nor2Rcm1.mergeSubGroupPeAlignInfo2targetPeAlignInfo_allSegSubGroup_Nor2Rcm1Only(
			peReadAlignInfo);
	}

	void freeMemory()
	{
		groupSegInfo_Nor1Rcm2.freeMemory();
		groupSegInfo_Nor2Rcm1.freeMemory();
	}
};


#endif