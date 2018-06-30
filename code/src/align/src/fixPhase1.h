// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXPHASE1_H
#define FIXPHASE1_H

#include <string>
#include <string.h>
#include "general/otherFunc2.h"
#include "general/segmentMatchComment_info.h"
//#include "general/groupSeg_info.h"
//#include "general/fixSubRegionPath.h"
//#include "splice_info.h"

using namespace std;

class FixPhase1Info
{
//public:
private:
	Seg_Info* segInfo_Nor1;
	Seg_Info* segInfo_Rcm1;
	Seg_Info* segInfo_Nor2;
	Seg_Info* segInfo_Rcm2;

	bool normalMapMain;
	bool rcmMapMain;
	bool normalMapMain_PE;
	bool rcmMapMain_PE;

	unsigned int norValLength;// = 0;
	unsigned int norValLength_PE;// = 0;
	unsigned int rcmValLength;// = 0;
	unsigned int rcmValLength_PE;// = 0;
	unsigned int minValLengthToStitch;// = MIN_LENGTH_TO_STITCH;

	GroupSeg_Info_BothDir groupSegInfoBothDir;
public:
	Path_Info pathInfo_Nor1;
	Path_Info pathInfo_Rcm1;
	Path_Info pathInfo_Nor2;
	Path_Info pathInfo_Rcm2;		


	void initiateGroupSegInfo()
	{
		groupSegInfoBothDir.initiateGroupSegInfo_withPeSegInfo_BothDir(
			segInfo_Nor1, segInfo_Rcm1, segInfo_Nor2, segInfo_Rcm2);
	}

	void generate_segMapTranscriptIndexVec_segLengthVecVec(Seg_Info* segInfo_Nor, Seg_Info* segInfo_Rcm,
		vector<int>& segMapTranscriptIndexVec_final, vector< vector<int> >& segLengthVecVec, Index_Info* indexInfo)
	{
		//map<int,int> transcriptIndex2segMapTranscriptIndexInVecMap;
		//cout << "generate_segMapTranscriptIndexVec_segLengthVecVec starts ......" << endl;
		vector<int> segMapTranscriptIndexVec;
		vector< set<int> > segIndexSetVec_Nor;
		vector< set<int> > segIndexSetVec_Rcm;
		int segmentNum_Nor = segInfo_Nor->returnSegmentNum();
		//cout << "segmentNum_Nor: " << segmentNum_Nor << endl;
		for(int tmpSeg = 0; tmpSeg < segmentNum_Nor; tmpSeg++)
		{
			int tmpSegAlignNum = segInfo_Nor->returnSegmentAlignNum(tmpSeg);
			int tmpSegLength = segInfo_Nor->returnSegmentLength(tmpSeg);
			if((tmpSegAlignNum > SEGMENTNUM))//||(tmpSegLength < MIN_LENGTH_TO_STITCH))
				continue;
			for(int tmpAlignIndex = 0; tmpAlignIndex < tmpSegAlignNum; tmpAlignIndex++)
			{
				unsigned int tmpAlignPos = segInfo_Nor->returnSegmentMapPos(tmpSeg, tmpAlignIndex);
				int tmpTranscriptIndex = indexInfo->getChr(tmpAlignPos);
				bool foundInCurrentTranscriptIndexVec_bool = false;
				int currentSegMapTranscriptIndexVecSize = segMapTranscriptIndexVec.size();
				for(int tmpIndex = 0; tmpIndex < currentSegMapTranscriptIndexVecSize; tmpIndex ++)
				{
					int tmpTranscriptIndexInVec = segMapTranscriptIndexVec[tmpIndex];
					if(tmpTranscriptIndex == tmpTranscriptIndexInVec)
					{
						foundInCurrentTranscriptIndexVec_bool = true;
						segIndexSetVec_Nor[tmpIndex].insert(tmpSeg);
						break;
					}
				}
				if(!foundInCurrentTranscriptIndexVec_bool)
				{	
					segMapTranscriptIndexVec.push_back(tmpTranscriptIndex);
					set<int> tmpIntSet_Nor, tmpIntSet_Rcm;
					tmpIntSet_Nor.insert(tmpSeg);
					segIndexSetVec_Nor.push_back(tmpIntSet_Nor);
					segIndexSetVec_Rcm.push_back(tmpIntSet_Rcm);
				}
			}
		}
		int segmentNum_Rcm = segInfo_Rcm->returnSegmentNum();
		//cout << "segmentNum_Rcm: " << segmentNum_Rcm << endl;
		for(int tmpSeg = 0; tmpSeg < segmentNum_Rcm; tmpSeg++)
		{
			int tmpSegAlignNum = segInfo_Rcm->returnSegmentAlignNum(tmpSeg);
			int tmpSegLength = segInfo_Rcm->returnSegmentLength(tmpSeg);
			if((tmpSegAlignNum > SEGMENTNUM))//||(tmpSegLength < MIN_LENGTH_TO_STITCH))
				continue;
			for(int tmpAlignIndex = 0; tmpAlignIndex < tmpSegAlignNum; tmpAlignIndex++)
			{
				unsigned int tmpAlignPos = segInfo_Rcm->returnSegmentMapPos(tmpSeg, tmpAlignIndex);
				int tmpTranscriptIndex = indexInfo->getChr(tmpAlignPos);
				for(int tmpIndex = 0; tmpIndex < segMapTranscriptIndexVec.size(); tmpIndex ++)
				{
					int tmpTranscriptIndexInVec = segMapTranscriptIndexVec[tmpIndex];
					if(tmpTranscriptIndex == tmpTranscriptIndexInVec)
					{
						segIndexSetVec_Rcm[tmpIndex].insert(tmpSeg);
						break;
					}
				}
			}	
		}

		int segMapTranscriptIndexVecSize = segMapTranscriptIndexVec.size();
		//cout << "segMapTranscriptIndexVecSize: " << segMapTranscriptIndexVecSize << endl;
		for(int tmp = 0; tmp < segMapTranscriptIndexVecSize; tmp++)
		{
			int setSize_Nor = segIndexSetVec_Nor[tmp].size();
			int setSize_Rcm = segIndexSetVec_Rcm[tmp].size();
			if((setSize_Nor > 0)&&(setSize_Rcm > 0))
			{
				int tmpTranscriptIndex = segMapTranscriptIndexVec[tmp];
				segMapTranscriptIndexVec_final.push_back(tmpTranscriptIndex);
				vector<int> tmpSegLengthVec;
				for(set<int>::iterator tmpIntSetIter = segIndexSetVec_Nor[tmp].begin();
					tmpIntSetIter != segIndexSetVec_Nor[tmp].end(); tmpIntSetIter ++)
				{
					int tmpIndexInSeg_Nor = (*tmpIntSetIter);
					int tmpSegLength_Nor = segInfo_Nor->returnSegmentLength(tmpIndexInSeg_Nor);
					tmpSegLengthVec.push_back(tmpSegLength_Nor);
				}
				for(set<int>::iterator tmpIntSetIter = segIndexSetVec_Rcm[tmp].begin();
					tmpIntSetIter != segIndexSetVec_Rcm[tmp].end(); tmpIntSetIter ++)
				{
					int tmpIndexInSeg_Rcm = (*tmpIntSetIter);
					int tmpSegLength_Rcm = segInfo_Rcm->returnSegmentLength(tmpIndexInSeg_Rcm);
					tmpSegLengthVec.push_back(tmpSegLength_Rcm);
				}
				segLengthVecVec.push_back(tmpSegLengthVec);				
			}
		}		
	}

	void doGeneCount_segMapOnTranscript(Gene_Count_Vec* geneCountVecInfo, int index_countVec, 
		Transcript2geneMap_Info* transcript2geneMapInfo, int tmpReadLength_twoEnds, Index_Info* indexInfo)
	{
		int tmpValidMapScore = (tmpReadLength_twoEnds/25) * 625;
		//cout << "tmpValidMapScore: " << tmpValidMapScore << endl;
		vector<int> segMapTranscriptIndexVec_Nor1Rcm2;
		vector< vector<int> > segLengthVecVec_Nor1Rcm2;		
		vector<int> segMapTranscriptIndexVec_Nor2Rcm1;
		vector< vector<int> > segLengthVecVec_Nor2Rcm1;
		this->generate_segMapTranscriptIndexVec_segLengthVecVec(segInfo_Nor1, segInfo_Rcm2,
			segMapTranscriptIndexVec_Nor1Rcm2, segLengthVecVec_Nor1Rcm2, indexInfo);
		//cout << "segMapTranscriptIndexVecSize_Nor1Rcm2: " << segMapTranscriptIndexVec_Nor1Rcm2.size() << endl;
		this->generate_segMapTranscriptIndexVec_segLengthVecVec(segInfo_Nor2, segInfo_Rcm1,
			segMapTranscriptIndexVec_Nor2Rcm1, segLengthVecVec_Nor2Rcm1, indexInfo);		
		//cout << "segMapTranscriptIndexVecSize_Nor2Rcm1: " << segMapTranscriptIndexVec_Nor2Rcm1.size() << endl;
		int tmpScoreMax = -1;
		vector<int> segMapScoreVec_Nor1Rcm2;
		int transcriptIndexVecSize_Nor1Rcm2 = segMapTranscriptIndexVec_Nor1Rcm2.size();
		for(int tmp = 0; tmp < transcriptIndexVecSize_Nor1Rcm2; tmp++)
		{
			int tmpScore = 0;
			int tmpSegNum = segLengthVecVec_Nor1Rcm2[tmp].size();
			for(int tmpSeg = 0; tmpSeg < tmpSegNum; tmpSeg++)
			{
				int tmpSegScore = ((segLengthVecVec_Nor1Rcm2[tmp])[tmpSeg]) * ((segLengthVecVec_Nor1Rcm2[tmp])[tmpSeg]);
				tmpScore += tmpSegScore;
			}
			//cout << "tmpScore: " << tmpScore << endl;
			segMapScoreVec_Nor1Rcm2.push_back(tmpScore);
			if(tmpScore > tmpScoreMax)
				tmpScoreMax = tmpScore;
		}
		
		vector<int> segMapScoreVec_Nor2Rcm1;
		int transcriptIndexVecSize_Nor2Rcm1 = segMapTranscriptIndexVec_Nor2Rcm1.size();
		for(int tmp = 0; tmp < transcriptIndexVecSize_Nor2Rcm1; tmp++)
		{
			int tmpScore = 0;
			int tmpSegNum = segLengthVecVec_Nor2Rcm1[tmp].size();
			for(int tmpSeg = 0; tmpSeg < tmpSegNum; tmpSeg++)
			{
				int tmpSegScore = ((segLengthVecVec_Nor2Rcm1[tmp])[tmpSeg]) * ((segLengthVecVec_Nor2Rcm1[tmp])[tmpSeg]);
				tmpScore += tmpSegScore;
			}
			//cout << "tmpScore: " << tmpScore << endl;
			segMapScoreVec_Nor2Rcm1.push_back(tmpScore);
			if(tmpScore > tmpScoreMax)
				tmpScoreMax = tmpScore;
		}

		set<int> geneIndexSet;
		if(tmpScoreMax >= tmpValidMapScore)
		{	
			for(int tmp = 0; tmp < transcriptIndexVecSize_Nor1Rcm2; tmp++)
			{
				int tmpScore = segMapScoreVec_Nor1Rcm2[tmp];
				if(tmpScore >= tmpScoreMax)
				{
					int tmpTranscriptIndex = segMapTranscriptIndexVec_Nor1Rcm2[tmp];
					int tmpGeneIndex = transcript2geneMapInfo->getGeneIndexFromTranscriptIndex(tmpTranscriptIndex);
					geneIndexSet.insert(tmpGeneIndex);
				}
			}
			for(int tmp = 0; tmp < transcriptIndexVecSize_Nor2Rcm1; tmp++)
			{
				int tmpScore = segMapScoreVec_Nor2Rcm1[tmp];
				if(tmpScore >= tmpScoreMax)
				{
					int tmpTranscriptIndex = segMapTranscriptIndexVec_Nor2Rcm1[tmp];
					int tmpGeneIndex = transcript2geneMapInfo->getGeneIndexFromTranscriptIndex(tmpTranscriptIndex);
					geneIndexSet.insert(tmpGeneIndex);
				}
			}
			int geneIndexNum = geneIndexSet.size();
			if(geneIndexNum == 1)
			{
				geneCountVecInfo->uniqueMappedReadCountIncrement_inOneCountInfo(index_countVec);
				int geneIndex = *(geneIndexSet.begin());
				geneCountVecInfo->addNewReadCount_inOneCountInfo(geneIndex, 1.0, index_countVec);
			}
			else
				geneCountVecInfo->multiMappedReadCountIncrement_inOneCountInfo(index_countVec);
		}
	}

	void fixPhase1_generateSegmentMatchCommentInfo(SegmentMatchComment_Info& tmpSegMatchCommentInfo)
	{
		tmpSegMatchCommentInfo.initiate_withSegInfo(
			normalMapMain, segInfo_Nor1, rcmMapMain, segInfo_Rcm1,
			normalMapMain_PE, segInfo_Nor2, rcmMapMain_PE, segInfo_Rcm2);
	}

	void fixPhase1_pathInfo_withGroupSegInfo(Index_Info* indexInfo,
		PE_Read_Info& peReadInfo, bool annotation_provided_bool,
		bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		groupSegInfoBothDir.getPossiPathFromSeg_fixGapAlso_allSegSubGroup_BothDir(
			indexInfo, peReadInfo,
			annotation_provided_bool, Do_annotation_only_bool,
			annotationInfo, spliceJunctionDistanceMax);
	}

	void fixPhase1_gapInfo_withGroupSegInfo(Index_Info* indexInfo,
		PE_Read_Info& peReadInfo, bool Do_extendHeadTail)
	{
		groupSegInfoBothDir.fixAllPath_fixGapAlso_allSegSubGroup_BothDir(
			indexInfo, peReadInfo, Do_extendHeadTail);
	}

	void fixPhase1_initiatePeAlignInfo_withGroupSegInfo(Index_Info* indexInfo)
	{
		groupSegInfoBothDir.initiatePeAlignInfo_allSegSubGroup_BothDir(indexInfo);
	}

	void alignmentFilter_fixPhase1_SJpenalty_withGroupSegInfo(
		int readLength_1, int readLength_2)
	{
		groupSegInfoBothDir.alignmentFilter_fixPhase1_SJpenalty_allSegSubGroup_BothDir(
			readLength_1, readLength_2);
	}

	void mergeSubGroupPeAlignInfo2targetPeAlignInfo_withGroupSegInfo(
		PE_Read_Alignment_Info& peAlignInfo)
	{
		groupSegInfoBothDir.mergeSubGroupPeAlignInfo_allSegSubGroup_BothDir(
			peAlignInfo);
	}

	void freeMemory_groupSegInfo()
	{
		groupSegInfoBothDir.freeMemory();
	}

	Seg_Info* returnSegInfo_Nor1()
	{
		return segInfo_Nor1;
	}
	Seg_Info* returnSegInfo_Rcm1()
	{
		return segInfo_Rcm1;
	}
	Seg_Info* returnSegInfo_Nor2()
	{
		return segInfo_Nor2;
	}
	Seg_Info* returnSegInfo_Rcm2()
	{
		return segInfo_Rcm2;
	}

	void memoryFree()
	{
		pathInfo_Nor1.memoryFree();
		pathInfo_Rcm1.memoryFree();
		pathInfo_Nor2.memoryFree();
		pathInfo_Rcm2.memoryFree();
		delete(segInfo_Nor1); delete(segInfo_Nor2); 
		delete(segInfo_Rcm1); delete(segInfo_Rcm2);		
	}

	bool perfectMatch(PE_Read_Info& peReadInfo)
	{
		if(
			( 
			(segInfo_Nor1->returnSegmentNum() == 1)
			&&(segInfo_Rcm2->returnSegmentNum() == 1)
			&&(segInfo_Nor1->returnSegmentLength(0) == peReadInfo.returnReadLength_end1())
			&&(segInfo_Rcm2->returnSegmentLength(0) == peReadInfo.returnReadLength_end2()) )
			||
			( 
			(segInfo_Nor2->returnSegmentNum() == 1)
			&&(segInfo_Rcm1->returnSegmentNum() == 1)
			&&(segInfo_Nor2->returnSegmentLength(0) == peReadInfo.returnReadLength_end1())
			&&(segInfo_Rcm1->returnSegmentLength(0) == peReadInfo.returnReadLength_end2()) )
			)
		{
			return true;
		}
		else 
			return false;
	}

	FixPhase1Info()
	{
		norValLength = 0;
		norValLength_PE = 0;
		rcmValLength = 0;
		rcmValLength_PE = 0;
		minValLengthToStitch = MIN_LENGTH_TO_STITCH;

	   	segInfo_Nor1 = new Seg_Info();		
	    segInfo_Rcm1 = new Seg_Info();
		segInfo_Nor2 = new Seg_Info();
		segInfo_Rcm2 = new Seg_Info();		
	}

	void filterSegInfo(
		Seg_Info* oldSegInfo_Nor1, Seg_Info* oldSegInfo_Rcm1, 
		Seg_Info* oldSegInfo_Nor2, Seg_Info* oldSegInfo_Rcm2, 
		Seg_Info* newSegInfo_Nor1, Seg_Info* newSegInfo_Rcm1, 
		Seg_Info* newSegInfo_Nor2, Seg_Info* newSegInfo_Rcm2, 
		bool& validSegInfoBool_Nor1, bool& validSegInfoBool_Rcm1, 
		bool& validSegInfoBool_Nor2, bool& validSegInfoBool_Rcm2) // decrease # of segments, with paired Information,  
	{
		bool pairAlignmentsExist_Nor1Rcm2_bool = false;
		if(validSegInfoBool_Nor1 && validSegInfoBool_Rcm2)
		{
			//do pairing
			pairAlignmentsExist_Nor1Rcm2_bool = this->pairingSegInfoWithAnother(
				oldSegInfo_Nor1, oldSegInfo_Rcm2, newSegInfo_Nor1, newSegInfo_Rcm2); // if paired seg found, return true;
		}

		bool pairAlignmentsExist_Nor2Rcm1_bool = false;
		if(validSegInfoBool_Nor2 && validSegInfoBool_Rcm1)
		{
			//do pairing
			pairAlignmentsExist_Nor2Rcm1_bool = this->pairingSegInfoWithAnother(
				oldSegInfo_Nor2, oldSegInfo_Rcm1, newSegInfo_Nor2, newSegInfo_Rcm1); // if paired Seg found, return true;
		}
	
		// if(pairAlignmentsExist_Nor1Rcm2_bool && pairAlignmentsExist_Nor2Rcm1_bool)
		// {

		// }

	}
	bool pairingSegInfoWithAnother(Seg_Info* oldSegInfo_1, Seg_Info* oldSegInfo_2, 
		Seg_Info* newSegInfo_1, Seg_Info* newSegInfo_2) // if paired seg found, return true; if not, return false
	{
		bool pairedSegFound_bool = false;
		map<int, pair< vector< pair<int, int> >, vector< pair<int, int> > > > segInfoGenomeRegionIndexMap;
		this->insertOldSegInfo2GenomeRegionIndexMap(oldSegInfo_1, oldSegInfo_2, segInfoGenomeRegionIndexMap);		

		return pairedSegFound_bool;
	}

	bool insertOldSegInfo2GenomeRegionIndexMap(Seg_Info* oldSegInfo_1, Seg_Info* oldSegInfo_2,
		map<int, pair< vector< pair<int, int> >, vector< pair<int, int> > > >& genomeRegionIndexMap)
	{
		//insertOldSegInfo_1
		if(oldSegInfo_1->returnRepeatRegion_index() >= 0)
		{
			return false;
		}
		// for(int tmp = 0; tmp < oldSegInfo_1)
		// {

		// }
		return true;
	}
	void outputGenomeRegionIndexMap2NewSegInfo(
		map<int, pair< vector< pair<int, int> >, vector< pair<int, int> > > >& genomeRegionIndexMap,
		Seg_Info* newSegInfo_1, Seg_Info* newSegInfo_2)
	{

	}

	void determineReadAlignDirection(bool& pairAlignmentsExist_Nor1Rcm2_bool, 
		bool& pairAlignmentsExist_Nor2Rcm1_bool,
		Seg_Info* segInfo_Nor1, Seg_Info* segInfo_Rcm1,
		Seg_Info* segInfo_Nor2, Seg_Info* segInfo_Rcm2
		) // choose the best direction, if forward, pairAlignmentsExist_Nor1Rcm2_bool = true;
		// if reverse pairAlignemntsExist_Nor2Rcm1_bool = true; if both directions, both true.
	{
		if(pairAlignmentsExist_Nor1Rcm2_bool && pairAlignmentsExist_Nor2Rcm1_bool)
		{

		}
		else
		{
			return;
		}
	}



	void fixPhase1_segInfo(
		unsigned int* sa, BYTE* lcpCompress,
		unsigned int* childTab, char* chrom,
		BYTE* verifyChild, Index_Info* indexInfo,
		int* preIndexMapLengthArray, unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray, PE_Read_Info& readInfo, 
		RepeatRegion_Info* repeatRegionInfo, bool checkQualSeqForReadSegSeq, bool SE_or_PE_bool)
	{
		//normalMapMain = false;
		/****************** the positive strand of end1 starts .************************/
		//cout << endl << "## do segment mapping for Nor_1 ##" << endl;
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());
		Seg_Info oriSegInfo_Nor1;// = new Seg_Info();
		MissingLongSeg_Info missingLongSegInfo_Nor1;
		normalMapMain = oriSegInfo_Nor1.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read, sa, lcpCompress, childTab, 
			chrom, &norValLength, verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_1()), repeatRegionInfo,
			missingLongSegInfo_Nor1);
		//normalMapMain = false;
		//cout << "Nor1 segInfo: " << endl << oriSegInfo_Nor1.segInfoStr(indexInfo) << endl;
		// cout << "Nor1 missingSegInfo: " << endl;// << missingLongSegInfo_Nor1->segInfoStr(indexInfo) << endl;
		if(normalMapMain)
		{
			if(missingLongSegInfo_Nor1.returnMissingSegGroupNum() > 0)
			{
				//cout << "missingLongSeg exists" << endl;
				normalMapMain = segInfo_Nor1->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Nor1, missingLongSegInfo_Nor1);
			}
			else 
			{
				(*segInfo_Nor1) = (oriSegInfo_Nor1);
			}

			// fix me: test: checkQualSeqForReadSegSeq
			if(checkQualSeqForReadSegSeq)
			{
				segInfo_Nor1->filterLowQualitySeg(readInfo, 1);
			}
		}
		/****************** the positive strand of end1 ends .************************/
		//cout << "Nor1 segInfo: " << endl << segInfo_Nor1->segInfoStr(indexInfo) << endl;
		//cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
	
		//rcmMapMain = false;
		/****************** the other strand of end1 starts .************************/
		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
	    Seg_Info oriSegInfo_Rcm1;// = new Seg_Info();
	    MissingLongSeg_Info missingLongSegInfo_Rcm1;
	    //rcmMapMain = false;
	   	rcmMapMain = oriSegInfo_Rcm1.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_RC, sa, lcpCompress, childTab,
			chrom, &rcmValLength, verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_1()), repeatRegionInfo,
			missingLongSegInfo_Rcm1);		
	   	
	 	//cout << "Rcm1 segInfo: " << endl << oriSegInfo_Rcm1.segInfoStr(indexInfo) << endl;
		// cout << "Rcm1 missingSegInfo: " << endl;// << missingLongSegInfo_Rcm1->segInfoStr(indexInfo) << endl;
		if(rcmMapMain)	
		{	
			if(missingLongSegInfo_Rcm1.returnMissingSegGroupNum() > 0)
			{
				//cout << "missingLongSeg exists" << endl;
				rcmMapMain = segInfo_Rcm1->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Rcm1, missingLongSegInfo_Rcm1);
			}	
			else
			{
				(*segInfo_Rcm1) = oriSegInfo_Rcm1;
			}

			// fix me: test: checkQualSeqForReadSegSeq
			if(checkQualSeqForReadSegSeq)
			{
				segInfo_Rcm1->filterLowQualitySeg(readInfo, 2);
			}
		}
		//cout << "Rcm1 segInfo: " << endl << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
		//cout << endl << "## do segment mapping for Nor_2 ##" << endl;
		/**************************** the other strand of end1 ends .*******************/


		if(!SE_or_PE_bool)
		{
			//normalMapMain_PE = false;
			/****************** the positive strand of end2 starts .************************/
			char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
			Seg_Info oriSegInfo_Nor2;// = new Seg_Info();
			MissingLongSeg_Info missingLongSegInfo_Nor2;
			//normalMapMain_PE = false;
			normalMapMain_PE = oriSegInfo_Nor2.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_PE, sa, lcpCompress, childTab, 
				chrom, &norValLength_PE, verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_2()), repeatRegionInfo,
				missingLongSegInfo_Nor2);
			
			//cout << "Nor2 segInfo: " << endl << oriSegInfo_Nor2.segInfoStr(indexInfo) << endl;
			//cout << "Nor2 missingSegInfo: " << endl;// << missingLongSegInfo_Nor2->segInfoStr(indexInfo) << endl;
			if(normalMapMain_PE)
			{	
				if(missingLongSegInfo_Nor2.returnMissingSegGroupNum() > 0)
				{
					//cout << "missingLongSeg exists" << endl;
				 	normalMapMain_PE = segInfo_Nor2->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Nor2, missingLongSegInfo_Nor2);
				}
				else
				{	
					(*segInfo_Nor2) = oriSegInfo_Nor2;
				}

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo_Nor2->filterLowQualitySeg(readInfo, 3);
				}
			}
			/****************** the positive strand of end2 ends .************************/

			//rcmMapMain_PE = false;
			/****************** the other strand of end2 starts .************************/
			//cout << "Nor2 segInfo: " << endl << segInfo_Nor2->segInfoStr(indexInfo) << endl;
			//cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
			char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
			Seg_Info oriSegInfo_Rcm2;// = new Seg_Info();
			MissingLongSeg_Info missingLongSegInfo_Rcm2;
			rcmMapMain_PE = oriSegInfo_Rcm2.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_RC_PE, sa, lcpCompress, childTab, 
				chrom, &rcmValLength_PE, verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_2()), repeatRegionInfo,
				missingLongSegInfo_Rcm2);	
			//rcmMapMain_PE = false;
			//cout << "Rcm2 oriSegInfo: " << endl << oriSegInfo_Rcm2.segInfoStr(indexInfo) << endl;
			//cout << "Rcm2 missingSegInfo: " << endl;// << missingLongSegInfo_Rcm2->segInfoStr(indexInfo) << endl;
			if(rcmMapMain_PE)
			{
				if(missingLongSegInfo_Rcm2.returnMissingSegGroupNum() > 0)
				{
					//cout << "missingLongSeg exists" << endl;
				 	rcmMapMain_PE = segInfo_Rcm2->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Rcm2, missingLongSegInfo_Rcm2);
				}
				else
				{
					(*segInfo_Rcm2) = oriSegInfo_Rcm2;
				}

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo_Rcm2->filterLowQualitySeg(readInfo, 4);
				}			
			}
			/****************** the other strand of end2 ends .************************/
			//cout << "Rcm2 segInfo: " << endl << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
		}
	}

	void fixPhase1_segInfo(
		unsigned int* sa, BYTE* lcpCompress,
		unsigned int* childTab, char* chrom,
		BYTE* verifyChild, Index_Info* indexInfo,
		int* preIndexMapLengthArray, unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray, PE_Read_Info& readInfo, 
		RepeatRegion_Info* repeatRegionInfo, bool checkQualSeqForReadSegSeq)
	{
		fixPhase1_segInfo_PE(
			sa, lcpCompress, childTab, chrom,
			verifyChild, indexInfo,
			preIndexMapLengthArray, preIndexIntervalStartArray,
			preIndexIntervalEndArray, readInfo, 
			repeatRegionInfo, checkQualSeqForReadSegSeq);
	}	

	void fixPhase1_segInfo_PE(
		unsigned int* sa, BYTE* lcpCompress,
		unsigned int* childTab, char* chrom,
		BYTE* verifyChild, Index_Info* indexInfo,
		int* preIndexMapLengthArray, unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray, PE_Read_Info& readInfo, 
		RepeatRegion_Info* repeatRegionInfo, bool checkQualSeqForReadSegSeq)
	{
		//cout << endl << "## do segment mapping for Nor_1 ##" << endl;
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());
		Seg_Info oriSegInfo_Nor1;// = new Seg_Info();
		MissingLongSeg_Info missingLongSegInfo_Nor1;
		normalMapMain = oriSegInfo_Nor1.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read, sa, lcpCompress, childTab, 
			chrom, &norValLength, verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_1()), repeatRegionInfo,
			missingLongSegInfo_Nor1);
	
		// cout << "Nor1 segInfo: " << endl << oriSegInfo_Nor1.segInfoStr(indexInfo) << endl;
		// cout << "Nor1 missingSegInfo: " << endl << missingLongSegInfo_Nor1->segInfoStr(indexInfo) << endl;
		if(normalMapMain)
		{
			if(missingLongSegInfo_Nor1.returnMissingSegGroupNum() > 0)
			{
				//cout << "missingLongSeg exists" << endl;
				normalMapMain = segInfo_Nor1->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Nor1, missingLongSegInfo_Nor1);
			}
			else 
			{
				(*segInfo_Nor1) = (oriSegInfo_Nor1);
			}

			// fix me: test: checkQualSeqForReadSegSeq
			if(checkQualSeqForReadSegSeq)
			{
				segInfo_Nor1->filterLowQualitySeg(readInfo, 1);
			}
		}
		//cout << "Nor1 segInfo: " << endl << segInfo_Nor1->segInfoStr(indexInfo) << endl;
		//cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
	    Seg_Info oriSegInfo_Rcm1;// = new Seg_Info();
	    MissingLongSeg_Info missingLongSegInfo_Rcm1;
	   	rcmMapMain = oriSegInfo_Rcm1.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_RC, sa, lcpCompress, childTab,
			chrom, &rcmValLength, verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMapLengthArray,
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_1()), repeatRegionInfo,
			missingLongSegInfo_Rcm1);		

	 	// cout << "Rcm1 segInfo: " << endl << oriSegInfo_Rcm1.segInfoStr(indexInfo) << endl;
		// cout << "Rcm1 missingSegInfo: " << endl << missingLongSegInfo_Rcm1->segInfoStr(indexInfo) << endl;
		if(rcmMapMain)	
		{	
			if(missingLongSegInfo_Rcm1.returnMissingSegGroupNum() > 0)
			{
				//cout << "missingLongSeg exists" << endl;
				rcmMapMain = segInfo_Rcm1->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Rcm1, missingLongSegInfo_Rcm1);
			}	
			else
			{
				(*segInfo_Rcm1) = oriSegInfo_Rcm1;
			}

			// fix me: test: checkQualSeqForReadSegSeq
			if(checkQualSeqForReadSegSeq)
			{
				segInfo_Rcm1->filterLowQualitySeg(readInfo, 2);
			}
		}
		//cout << "Rcm1 segInfo: " << endl << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
		//cout << endl << "## do segment mapping for Nor_2 ##" << endl;
	
		//if(!SE_or_PE_bool)
		//{
			char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
			Seg_Info oriSegInfo_Nor2;// = new Seg_Info();
			MissingLongSeg_Info missingLongSegInfo_Nor2;
			normalMapMain_PE = oriSegInfo_Nor2.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_PE, sa, lcpCompress, childTab, 
				chrom, &norValLength_PE, verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_2()), repeatRegionInfo,
				missingLongSegInfo_Nor2);
			
			//cout << "Nor2 segInfo: " << endl << oriSegInfo_Nor2.segInfoStr(indexInfo) << endl;
			//cout << "Nor2 missingSegInfo: " << endl << missingLongSegInfo_Nor2->segInfoStr(indexInfo) << endl;
			if(normalMapMain_PE)
			{	
				if(missingLongSegInfo_Nor2.returnMissingSegGroupNum() > 0)
				{
					//cout << "missingLongSeg exists" << endl;
				 	normalMapMain_PE = segInfo_Nor2->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Nor2, missingLongSegInfo_Nor2);
				}
				else
				{	
					(*segInfo_Nor2) = oriSegInfo_Nor2;
				}

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo_Nor2->filterLowQualitySeg(readInfo, 3);
				}
			}
			//cout << "Nor2 segInfo: " << endl << segInfo_Nor2->segInfoStr(indexInfo) << endl;
			//cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
			char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
			Seg_Info oriSegInfo_Rcm2;// = new Seg_Info();
			MissingLongSeg_Info missingLongSegInfo_Rcm2;
			rcmMapMain_PE = oriSegInfo_Rcm2.mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(read_RC_PE, sa, lcpCompress, childTab, 
				chrom, &rcmValLength_PE, verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMapLengthArray,
				preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_2()), repeatRegionInfo,
				missingLongSegInfo_Rcm2);	
		
			//cout << "Rcm2 oriSegInfo: " << endl << oriSegInfo_Rcm2.segInfoStr(indexInfo) << endl;
			//cout << "Rcm2 missingSegInfo: " << endl << missingLongSegInfo_Rcm2->segInfoStr(indexInfo) << endl;
			if(rcmMapMain_PE)
			{
				if(missingLongSegInfo_Rcm2.returnMissingSegGroupNum() > 0)
				{
					//cout << "missingLongSeg exists" << endl;
				 	rcmMapMain_PE = segInfo_Rcm2->combineOriSegInfoWithMissingLongSegInfo(oriSegInfo_Rcm2, missingLongSegInfo_Rcm2);
				}
				else
				{
					(*segInfo_Rcm2) = oriSegInfo_Rcm2;
				}

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo_Rcm2->filterLowQualitySeg(readInfo, 4);
				}			
			}		
			//cout << "Rcm2 segInfo: " << endl << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
		//}
	}
	#ifdef PERSONALIZED_CHR_SEQ
	void fixPhase1_segInfo_map2syntheticSNPtransSeq(
		unsigned int* sa, BYTE* lcpCompress, unsigned int* childTab, char* chrom,
		BYTE* verifyChild, Index_Info* indexInfo, PE_Read_Info& readInfo)
	{
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());
		normalMapMain = segInfo_Nor1->greedyMapWithoutPreIndexArray(read, sa, lcpCompress, childTab, chrom, 
			verifyChild, readInfo.returnReadLength_end1(), indexInfo);

		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
		rcmMapMain = segInfo_Rcm1->greedyMapWithoutPreIndexArray(read_RC, sa, lcpCompress, childTab, chrom, 
			verifyChild, readInfo.returnReadLength_end1(), indexInfo);

		char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
		normalMapMain_PE = segInfo_Nor2->greedyMapWithoutPreIndexArray(read_PE, sa, lcpCompress, childTab, chrom, 
			verifyChild, readInfo.returnReadLength_end2(), indexInfo);

		char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
		rcmMapMain_PE = segInfo_Rcm2->greedyMapWithoutPreIndexArray(read_RC_PE, sa, lcpCompress, childTab, chrom, 
			verifyChild, readInfo.returnReadLength_end2(), indexInfo);
	}
	
	void fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, PE_Read_Info& readInfo, Index_Info* indexInfo, int SNPlocInSyntheticSNPseq)
	{
		//cout << "start to do fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo( " << endl;
		Seg_Info segInfo2snpSeq_Nor1, segInfo2snpSeq_Rcm1, segInfo2snpSeq_Nor2, segInfo2snpSeq_Rcm2;
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());		
		bool tmpSnpSeqMap_Nor1 = segInfo2snpSeq_Nor1.greedyMapWithoutPreIndexArray(read, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end1(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Nor1: " << endl << segInfo2snpSeq_Nor1.segInfoStr(indexInfo_snp) << endl;
		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
		bool tmpSnpSeqMap_Rcm1 = segInfo2snpSeq_Rcm1.greedyMapWithoutPreIndexArray(read_RC, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end1(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Rcm1: " << endl << segInfo2snpSeq_Rcm1.segInfoStr(indexInfo_snp) << endl;
		char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
		bool tmpSnpSeqMap_Nor2 = segInfo2snpSeq_Nor2.greedyMapWithoutPreIndexArray(read_PE, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end2(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Nor2: " << endl << segInfo2snpSeq_Nor2.segInfoStr(indexInfo_snp) << endl;
		char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
		bool tmpSnpSeqMap_Rcm2 = segInfo2snpSeq_Rcm2.greedyMapWithoutPreIndexArray(read_RC_PE, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end2(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Rcm2: " << endl << segInfo2snpSeq_Rcm2.segInfoStr(indexInfo_snp) << endl;
		//cout << "start to update oriSegInfo ..." << endl;
		//cout << "start to update Nor_1" << endl;
		if(normalMapMain && tmpSnpSeqMap_Nor1)
			segInfo_Nor1->update_includeSNPseqMapSegInfo(segInfo2snpSeq_Nor1, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		//cout << "start to update Rcm_1" << endl;
		if(rcmMapMain && tmpSnpSeqMap_Rcm1)
			segInfo_Rcm1->update_includeSNPseqMapSegInfo(segInfo2snpSeq_Rcm1, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		//cout << "start to update Nor_2" << endl;
		if(normalMapMain_PE && tmpSnpSeqMap_Nor2)
			segInfo_Nor2->update_includeSNPseqMapSegInfo(segInfo2snpSeq_Nor2, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		//cout << "start to update Rcm_2" << endl;
		if(rcmMapMain_PE && tmpSnpSeqMap_Rcm2)
			segInfo_Rcm2->update_includeSNPseqMapSegInfo(segInfo2snpSeq_Rcm2, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
	}

	void fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo_varySNPmer(
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, PE_Read_Info& readInfo, Index_Info* indexInfo)
	{
		//cout << "start to do fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo( " << endl;
		Seg_Info segInfo2snpSeq_Nor1, segInfo2snpSeq_Rcm1, segInfo2snpSeq_Nor2, segInfo2snpSeq_Rcm2;
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());		
		bool tmpSnpSeqMap_Nor1 = segInfo2snpSeq_Nor1.greedyMapWithoutPreIndexArray(read, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end1(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Nor1: " << endl << segInfo2snpSeq_Nor1.segInfoStr(indexInfo_snp) << endl;
		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
		bool tmpSnpSeqMap_Rcm1 = segInfo2snpSeq_Rcm1.greedyMapWithoutPreIndexArray(read_RC, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end1(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Rcm1: " << endl << segInfo2snpSeq_Rcm1.segInfoStr(indexInfo_snp) << endl;
		char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
		bool tmpSnpSeqMap_Nor2 = segInfo2snpSeq_Nor2.greedyMapWithoutPreIndexArray(read_PE, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end2(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Nor2: " << endl << segInfo2snpSeq_Nor2.segInfoStr(indexInfo_snp) << endl;
		char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
		bool tmpSnpSeqMap_Rcm2 = segInfo2snpSeq_Rcm2.greedyMapWithoutPreIndexArray(read_RC_PE, sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
			verifyChild_snp, readInfo.returnReadLength_end2(), indexInfo_snp);
		//cout << "segInfo2snpSeq_Rcm2: " << endl << segInfo2snpSeq_Rcm2.segInfoStr(indexInfo_snp) << endl;
		//cout << "start to update oriSegInfo ..." << endl;
		//cout << "start to update Nor_1" << endl;
		if(normalMapMain && tmpSnpSeqMap_Nor1)
			segInfo_Nor1->update_includeSNPseqMapSegInfo_varySNPmer(segInfo2snpSeq_Nor1, indexInfo_snp, indexInfo);
		//cout << "start to update Rcm_1" << endl;
		if(rcmMapMain && tmpSnpSeqMap_Rcm1)
			segInfo_Rcm1->update_includeSNPseqMapSegInfo_varySNPmer(segInfo2snpSeq_Rcm1, indexInfo_snp, indexInfo);
		//cout << "start to update Nor_2" << endl;
		if(normalMapMain_PE && tmpSnpSeqMap_Nor2)
			segInfo_Nor2->update_includeSNPseqMapSegInfo_varySNPmer(segInfo2snpSeq_Nor2, indexInfo_snp, indexInfo);
		//cout << "start to update Rcm_2" << endl;
		if(rcmMapMain_PE && tmpSnpSeqMap_Rcm2)
			segInfo_Rcm2->update_includeSNPseqMapSegInfo_varySNPmer(segInfo2snpSeq_Rcm2, indexInfo_snp, indexInfo);
	}

	void fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(unsigned int* sa_snp, BYTE* lcpCompress_snp, 
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		PE_Read_Info& readInfo, Index_Info* indexInfo, int SNPlocInSyntheticSNPseq)
	{
		// cout << "###############################" << endl;
		// cout << "start to do fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo( " << endl;
		// cout << "start to do update_targetMap2SNPseq( for Nor1" << endl;
		if(normalMapMain)
			segInfo_Nor1->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
				indexInfo_snp, readInfo.returnReadSeq_1(), indexInfo, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Rcm1" << endl;
		if(rcmMapMain)
			segInfo_Rcm1->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
				indexInfo_snp, readInfo.returnRcmReadSeq_1(), indexInfo, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Nor2" << endl;
		if(normalMapMain_PE)
			segInfo_Nor2->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
				indexInfo_snp, readInfo.returnReadSeq_2(), indexInfo, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Rcm2" << endl;
		if(rcmMapMain_PE)
			segInfo_Rcm2->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
				indexInfo_snp, readInfo.returnRcmReadSeq_2(), indexInfo, SNPlocInSyntheticSNPseq);
	}

	void fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo_varySNPmer(unsigned int* sa_snp, BYTE* lcpCompress_snp, 
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		PE_Read_Info& readInfo, Index_Info* indexInfo)//, int SNPlocInSyntheticSNPseq)
	{
		// cout << "###############################" << endl;
		// cout << "start to do fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo( " << endl;
		// cout << "start to do update_targetMap2SNPseq( for Nor1" << endl;
		if(normalMapMain)
			segInfo_Nor1->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
				verifyChild_snp, indexInfo_snp, readInfo.returnReadSeq_1(), indexInfo);//, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Rcm1" << endl;
		if(rcmMapMain)
			segInfo_Rcm1->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
				verifyChild_snp, indexInfo_snp, readInfo.returnRcmReadSeq_1(), indexInfo);//, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Nor2" << endl;
		if(normalMapMain_PE)
			segInfo_Nor2->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
				verifyChild_snp, indexInfo_snp, readInfo.returnReadSeq_2(), indexInfo);//, SNPlocInSyntheticSNPseq);
		//cout << "start to do update_targetMap2SNPseq( for Rcm2" << endl;
		if(rcmMapMain_PE)
			segInfo_Rcm2->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
				verifyChild_snp, indexInfo_snp, readInfo.returnRcmReadSeq_2(), indexInfo);//, SNPlocInSyntheticSNPseq);
	}	

	#endif
	void fixPhase1_pathInfo(bool Do_cirRNA, Index_Info* indexInfo, PE_Read_Info& readInfo,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool SE_or_PE_bool
		)
	{
		//cout << "start to generate Path ..." << endl;
		if(!Do_cirRNA)
		{
			bool tooManyValSegs_Nor1_skip = tooManyValSegs_bool_phase1(
				segInfo_Nor1, readInfo.returnReadLength_end1(), normalMapMain);
			bool tooManyValSegs_Rcm1_skip = tooManyValSegs_bool_phase1(
				segInfo_Rcm1, readInfo.returnReadLength_end1(), rcmMapMain);
			//cout << "tooManyValSegs_Nor1_skip: " << tooManyValSegs_Nor1_skip << endl;
			//cout << "tooManyValSegs_Rcm1_skip: " << tooManyValSegs_Rcm1_skip << endl;
			if(tooManyValSegs_Nor1_skip)
			{
				normalMapMain = false;
			}
			if(tooManyValSegs_Rcm1_skip)
			{
				rcmMapMain = false;
			}
			//cout << endl << "*** getPossiPathFromSeg for nor_1" << endl << endl;
			if(normalMapMain)
			{	
				pathInfo_Nor1.getPossiPathFromSeg_fixGapAlso(segInfo_Nor1, indexInfo, readInfo.returnReadSeq_1(),
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
			}
			//cout << endl << "*** getPossiPathFromSeg for rcm_1" << endl << endl;
			if(rcmMapMain)
			{
				pathInfo_Rcm1.getPossiPathFromSeg_fixGapAlso(segInfo_Rcm1, indexInfo, readInfo.returnRcmReadSeq_1(),
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
			}	


			if(!SE_or_PE_bool)
			{	
				bool tooManyValSegs_Nor2_skip = tooManyValSegs_bool_phase1(
					segInfo_Nor2, readInfo.returnReadLength_end2(), normalMapMain_PE);
				bool tooManyValSegs_Rcm2_skip = tooManyValSegs_bool_phase1(
					segInfo_Rcm2, readInfo.returnReadLength_end2(), rcmMapMain_PE);

				if(tooManyValSegs_Nor2_skip)
				{
					normalMapMain_PE = false;
				}
				if(tooManyValSegs_Rcm2_skip)
				{
					rcmMapMain_PE = false;
				}					
				//cout << endl << "*** getPossiPathFromSeg for nor_2" << endl << endl;
				if(normalMapMain_PE)
				{
					pathInfo_Nor2.getPossiPathFromSeg_fixGapAlso(segInfo_Nor2, indexInfo, readInfo.returnReadSeq_2(),
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
				}
				//cout << endl << "*** getPossiPathFromSeg for rcm_2" << endl << endl;
				if(rcmMapMain_PE)
				{
					pathInfo_Rcm2.getPossiPathFromSeg_fixGapAlso(segInfo_Rcm2, indexInfo, readInfo.returnRcmReadSeq_2(),
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
				}
			}
			//cout << endl << "*** finish getting all possible paths from seg !" << endl << endl;
		}
		else
		{
			if(normalMapMain)
				pathInfo_Nor1.getPossiPathFromSeg_cirRNA(segInfo_Nor1);
			if(rcmMapMain)
				pathInfo_Rcm1.getPossiPathFromSeg_cirRNA(segInfo_Rcm1);
			if(!SE_or_PE_bool)
			{	
				if(normalMapMain_PE)
					pathInfo_Nor2.getPossiPathFromSeg_cirRNA(segInfo_Nor2);
				if(rcmMapMain_PE)
					pathInfo_Rcm2.getPossiPathFromSeg_cirRNA(segInfo_Rcm2);			
			}
		}
	}

	void fixPhase1_gapInfo(PE_Read_Info& readInfo, Index_Info* indexInfo, 
		bool Do_cirRNA, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		bool SE_or_PE_bool)
	{
		//cout << "start to fix gapInfo_Nor1 ..." << endl;
		if(!Do_cirRNA)
		{
			//cout << "*** start to fix gapInfo_Nor1 ... " << endl;
			pathInfo_Nor1.fixAllPath_fixGapAlso(segInfo_Nor1, indexInfo, readInfo.returnReadSeq_1(), Do_extendHeadTail);
			//cout << "*** start to fix gapInfo_Rcm1 ... " << endl;			
			pathInfo_Rcm1.fixAllPath_fixGapAlso(segInfo_Rcm1, indexInfo, readInfo.returnRcmReadSeq_1(), Do_extendHeadTail);
			if(!SE_or_PE_bool)
			{	
				//cout << "*** start to fix gapInfo_Nor2 ... " << endl;			
				pathInfo_Nor2.fixAllPath_fixGapAlso(segInfo_Nor2, indexInfo, readInfo.returnReadSeq_2(), Do_extendHeadTail);
				//cout << "*** start to fix gapInfo_Rcm2 ... " << endl;			
				pathInfo_Rcm2.fixAllPath_fixGapAlso(segInfo_Rcm2, indexInfo, readInfo.returnRcmReadSeq_2(), Do_extendHeadTail);
			}
		}
		else
		{
			pathInfo_Nor1.fixAllPath_cirRNA(segInfo_Nor1, indexInfo, readInfo.returnReadSeq_1(), Do_extendHeadTail, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
			pathInfo_Rcm1.fixAllPath_cirRNA(segInfo_Rcm1, indexInfo, readInfo.returnRcmReadSeq_1(), Do_extendHeadTail, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
			if(!SE_or_PE_bool)
			{	
				pathInfo_Nor2.fixAllPath_cirRNA(segInfo_Nor2, indexInfo, readInfo.returnReadSeq_2(), Do_extendHeadTail, 
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
				pathInfo_Rcm2.fixAllPath_cirRNA(segInfo_Rcm2, indexInfo, readInfo.returnRcmReadSeq_2(), Do_extendHeadTail, 
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo);			
			}
		}
	}

	void coutSegInfo_Nor1(Index_Info* indexInfo)
	{
		cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
	}
	void coutSegInfo_Rcm1(Index_Info* indexInfo)
	{
		cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
	}
	void coutSegInfo_Nor2(Index_Info* indexInfo)
	{
		cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
	}
	void coutSegInfo_Rcm2(Index_Info* indexInfo)
	{
		cout << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
	}

	void coutPossiPathStr_Nor1()
	{
		cout << pathInfo_Nor1.possiPathStr() << endl;
	}

	void coutPossiPathStr_Rcm1()
	{
		cout << pathInfo_Rcm1.possiPathStr() << endl;
	}

	void coutDebugInfo(PE_Read_Info& readInfo, Index_Info* indexInfo, bool SE_or_PE_bool)
	{
			cout << endl << "##### readName_1: " << readInfo.returnReadName_1() << " #####" << endl;
			if(!SE_or_PE_bool)
			{	
				cout << "##### readName_2: " << readInfo.returnReadName_2() << " #####"<< endl;
			}	
			cout << endl << "## do segment mapping for Nor_1 ##" << endl;
			cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
			cout << "repeatRegion_index: " << segInfo_Nor1->returnRepeatRegion_index() << endl;
			cout << pathInfo_Nor1.possiPathStr() << endl;
			cout << pathInfo_Nor1.fixedPathVecStr(indexInfo, segInfo_Nor1) << endl;
			cout << pathInfo_Nor1.finalFixedPathStr(indexInfo) << endl;
			cout << pathInfo_Nor1.getFixGapInfoStr() << endl;

		    cout << endl << "## do segment mapping for Rcm_1 ##" << endl;
		   	cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
		   	cout << "repeatRegion_index: " << segInfo_Rcm1->returnRepeatRegion_index() << endl;
			cout << pathInfo_Rcm1.possiPathStr() << endl;
			cout << pathInfo_Rcm1.fixedPathVecStr(indexInfo, segInfo_Rcm1) << endl;
			cout << pathInfo_Rcm1.finalFixedPathStr(indexInfo) << endl;
			cout << pathInfo_Rcm1.getFixGapInfoStr() << endl;
			if(!SE_or_PE_bool)
			{	
				cout << endl << "## do segment mapping for Nor_2 ##" << endl;
				cout << segInfo_Nor2->segInfoStr(indexInfo) << endl;
				cout << "repeatRegion_index: " << segInfo_Nor2->returnRepeatRegion_index() << endl;
				cout << pathInfo_Nor2.possiPathStr() << endl;
				cout << pathInfo_Nor2.fixedPathVecStr(indexInfo, segInfo_Nor2) << endl;
				cout << pathInfo_Nor2.finalFixedPathStr(indexInfo) << endl;
				cout << pathInfo_Nor2.getFixGapInfoStr() << endl;

				cout << endl << "## do segment mapping for Rcm_2 ##" << endl;
				cout << segInfo_Rcm2->segInfoStr(indexInfo) << endl;
				cout << "repeatRegion_index: " << segInfo_Rcm2->returnRepeatRegion_index() << endl;
				cout << pathInfo_Rcm2.possiPathStr() << endl;
				cout << pathInfo_Rcm2.fixedPathVecStr(indexInfo, segInfo_Rcm2) << endl;
				cout << pathInfo_Rcm2.finalFixedPathStr(indexInfo) << endl;
				cout << pathInfo_Rcm2.getFixGapInfoStr() << endl;
			}
	}

	void coutDebugInfo_SE(PE_Read_Info& readInfo, Index_Info* indexInfo)
	{
		cout << endl << "##### readName: " << readInfo.returnReadName_1() << " #####" << endl;
		cout << endl << "## do segment mapping for Nor ##" << endl;
		cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
		cout << pathInfo_Nor1.finalFixedPathStr(indexInfo) << endl;
		cout << endl << "## do segment mapping for Rcm ##" << endl;
		cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
		cout << pathInfo_Rcm1.finalFixedPathStr(indexInfo) << endl;
	}

	void coutFixedPathInfo_SE(PE_Read_Info& readInfo, Index_Info* indexInfo)
	{
		cout << endl << "##### readName: " << readInfo.returnReadName_1() << " #####" << endl;
		cout << endl << "## do segment mapping for Nor ##" << endl;
		//cout << segInfo_Nor1->segInfoStr(indexInfo) << endl;
		cout << pathInfo_Nor1.finalFixedPathStr(indexInfo) << endl;
		cout << endl << "## do segment mapping for Rcm ##" << endl;
		//cout << segInfo_Rcm1->segInfoStr(indexInfo) << endl;
		cout << pathInfo_Rcm1.finalFixedPathStr(indexInfo) << endl;
	}	

	string returnSubGroup_segInfoStr(Index_Info* indexInfo)
	{
		return groupSegInfoBothDir.returnBothDirSubSegGroup_segInfoStr(indexInfo);
	}

	string returnSubGroup_pathInfoStr()
	{
		return groupSegInfoBothDir.returnBothDirSubSegGroup_pathInfoStr();
	}

	string returnSubGroup_gapInfoStr(Index_Info* indexInfo)
	{
		return groupSegInfoBothDir.returnBothDirSubSegGroup_gapInfoStr(indexInfo);
	}
};

#endif