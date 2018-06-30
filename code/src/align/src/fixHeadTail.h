// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXHEADTAIL_H
#define FIXHEADTAIL_H

#include <string>
#include <string.h>

using namespace std;

class FixHeadTailInfo
{
public:

	FixHeadTailInfo()
	{

	}

	int returnMaxSpliceDistanceTargetMapping(int singleAnchorLength)
	{
		if(singleAnchorLength < 10)
			return MAX_SPLICE_DISTANCE_TARGETMAPPING_SHORTSINGLEANCHOR;
		else
			return MAX_SPLICE_DISTANCE_TARGETMAPPING_LONGSINGLEANCHOR;
	}

	void fixHead_shortAnchorRemappingOnly_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignInfo)
	{
		//cout << "start fixHead function ..." << endl;
		
		Unfixed_Head unfixedHeadInfo;
		//unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType_new(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
		unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType_new(peReadInfo, tmpAlignInfoType, tmpAlignInfo, indexInfo);
		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);
		//cout << endl << "original tmpAlignInfo: " << endl << tmpAlignInfo->returnAlignInfoStr() << endl;				
		///////////////////////////////////////////////////////////////////////////////////////////	
		/////////////////////////  try remapping with splice junction hash ////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start SJsearchInSJhash !" << endl;
		bool spliceJunctionFoundInHash;

		if(spliceJunctionHashExists)
			spliceJunctionFoundInHash
				= unfixedHeadInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
		else
			spliceJunctionFoundInHash = false;

		//fix me
		//spliceJunctionFoundInHash = false;
		//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
		if(spliceJunctionFoundInHash)
		{
			int tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[0].first;
			int tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[0].second;
			int tmpFirstMatchLength = tmpDonerEndPosInRead;
			//cout << "tmpDonerEndPosInRead: " << tmpDonerEndPosInRead << endl;
			//cout << "tmpSpliceJunctionDistance: " << tmpSpliceJunctionDistance << endl;
			//cout << "tmpFirstMatchLength: " << endl;
			int newMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[0];

			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			#ifdef DETECT_CIRCULAR_RNA
			newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
				tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, 
				newMismatchToAddInHead,
				(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
				(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[0],
				tmpAlignInfo);
			#else
			newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
				tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, 
				newMismatchToAddInHead,
				(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
				(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[0],
				tmpAlignInfo);
			#endif
			//cout << endl << "newTmpAlignInfo: " << endl << newTmpAlignInfo->returnAlignInfoStr() << endl;
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmpSJposVec = 1; tmpSJposVec < (unfixedHeadInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
			{
				int tmpNewMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

				tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].first;
				tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].second;
				tmpFirstMatchLength = tmpDonerEndPosInRead;	

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
					tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInHead,
					(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
					(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec],
					tmpAlignInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
					tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInHead,
					(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
					(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec],
					tmpAlignInfo);
				#endif
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
				//cout << endl << "newTmpAlignInfo: " << endl << newTmpAlignInfo->returnAlignInfoStr() << endl;
			}
			//continue;
			return;
		}
		else // no candidate SJ found in hash
		{}
	}
	
	void fixHead_shortAnchorRemappingOnly_withAlignInfer(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJinfo, 
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignInfo,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{
		//cout << "start fixHead_shortAnchorRemappingOnly_withAlignInfer ..." << endl;
		
		Unfixed_Head unfixedHeadInfo;
		//unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType_new(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
		unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType_new(peReadInfo, tmpAlignInfoType, tmpAlignInfo, indexInfo);
		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);
						
		///////////////////////////////////////////////////////////////////////////////////////////	
		/////////////////////////  try remapping with splice junction hash ////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start SJsearchInSJhash !" << endl;
		bool spliceJunctionFoundInHash;
		vector< vector<Jump_Code> > inferedPathJumpCodeVecVec;
		vector< vector<int> > inferedPathMismatchPosVecVec;
		vector< vector<char> > inferedPathMismatchCharVecVec;		
		if(spliceJunctionHashExists)
		{
			spliceJunctionFoundInHash
				//= unfixedHeadInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
				= unfixedHeadInfo.SJsearchInSJhash_areaStringHash_withAlignInferJuncHash(
					SJinfo, alignInferJunctionHashInfo, 
					readSeqWithDirection, indexInfo, SJinfo->areaSize,
					inferedPathJumpCodeVecVec,
					inferedPathMismatchPosVecVec,
					inferedPathMismatchCharVecVec);
		}
		else
		{
			spliceJunctionFoundInHash = false;
		}

		//fix me
		//spliceJunctionFoundInHash = false;
		// cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
		// cout << "inferedPathJumpCodeVecVecSize: " << inferedPathJumpCodeVecVec.size() << endl;
		// cout << "inferedPathMismatchPosVecVecSize: " << inferedPathMismatchPosVecVec.size() << endl;
		// cout << "inferedPathMismatchCharVecVecSize: " << inferedPathMismatchCharVecVec.size() << endl;
		// cout << ">>>>>>>>>>>>>>>>>>>>" << endl 
		// 	<< " spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;

		if(spliceJunctionFoundInHash)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			// for(int tmp = 0; tmp < inferedPathJumpCodeVecVec[0].size(); tmp++)
			// {
			// 	cout << "tmpJumpCodeType: " << (inferedPathJumpCodeVecVec[0])[tmp].type;
			// 	cout << "  tmpJumpCodeLength: " << (inferedPathJumpCodeVecVec[0])[tmp].len << endl;
			// }
			newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_withAlignInfer(
				inferedPathJumpCodeVecVec[0], 
				inferedPathMismatchPosVecVec[0],
				inferedPathMismatchCharVecVec[0],
				indexInfo, tmpAlignInfo);

			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmp = 1; tmp < inferedPathJumpCodeVecVec.size(); tmp++)
			{
				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_withAlignInfer(
					inferedPathJumpCodeVecVec[tmp], 
					inferedPathMismatchPosVecVec[tmp],
					inferedPathMismatchCharVecVec[tmp],
					indexInfo, tmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
			//continue;
			return;
		}
		else // no candidate SJ found in hash
		{}
	}

	#ifdef PERSONALIZED_CHR_SEQ
	void fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo, 
		Alignment_Info* tmpAlignmentInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, int SNPlocInSyntheticSNPseq
		)
	{
		//cout << "start to fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap(" << endl;
		Incomplete_Long_Head* incompleteHeadInfo = new Incomplete_Long_Head();
		incompleteHeadInfo->getIncompleteLongHeadInfoFromRecordWithAlignInfoType_new(
			tmpAlignmentInfo, indexInfo);		
		string incompleteLongHeadSeq = peReadInfo.returnIncompleteLongHeadSeq(
			tmpAlignInfoType, incompleteHeadInfo->returnUnfixedHeadLength());
		char* incompleteLongHeadSeqChar = const_cast<char*>(incompleteLongHeadSeq.c_str());
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
		int secondLevelIndexNO = incompleteHeadInfo->returnSecondLevelIndexNum() - 1;
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO + 1))
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
			//continue;
		}	
		//cout << "start to do incompleteMap " << endl;
		bool incompleteMapBool;
		//cout << "secondLevelIndexNO: " << secondLevelIndexNO << endl;
		//cout << "incompleteHeadInfo_unfixedHeadLength: " << incompleteHeadInfo->returnUnfixedHeadLength() << endl;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongHeadSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteHeadInfo->returnUnfixedHeadLength(), indexInfo);
		//cout << "incompleteMapBool: " << incompleteMapBool << endl; 
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
			incompleteHeadInfo->returnMapPosIntervalStart(),
			incompleteHeadInfo->returnMapPosIntervalEnd(), 
			incompleteHeadInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteHeadInfo->returnMidPartMapChrName());
		//cout << "#########################" << endl;
		//cout << "beforeUpdating segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;

		#ifdef PERSONALIZED_CHR_SEQ
		Seg_Info segInfo_SNPseq;
		int tmpIncompleteSeqLength = incompleteLongHeadSeq.length();
		bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(incompleteLongHeadSeqChar,
			sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
		if(tmpMapWithSnpSeqIndexBool)
			segInfo->update_includeSNPseqMapSegInfo(segInfo_SNPseq, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		//cout << "afterUpdating with greedyMapping segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;
		segInfo->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
			indexInfo_snp, incompleteLongHeadSeq, indexInfo, SNPlocInSyntheticSNPseq);
		//cout << "afterUpdating with targetMapping segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;
		#endif

		if(segInfo->returnSegmentNum() >= SEGMENTNUM)
		{
			delete seg2ndOriInfo;
			delete segInfo;
			delete incompleteHeadInfo;
			return;
		}

		segInfo->addMidPartSeg_incompleteHead(incompleteHeadInfo->returnMidPartLength(), 
			(incompleteHeadInfo->returnUnfixedHeadLength() + 1), 
			incompleteHeadInfo->returnMidPartMapChrInt(), incompleteHeadInfo->returnMidPartMapPosInChr(), indexInfo);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}
		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg_incompleteHead(segInfo, spliceJunctionDistanceMax); // Fix me: why not replace it by getPossiPathFromSeg_incompleteHead  
	
		pathInfo->filterPath_incompleteHead(segInfo); // filter out those path without the midPartSeg
	
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongHeadSeq(
				tmpAlignInfoType, incompleteHeadInfo->returnIncompleteHeadAndMidPartLength()), 
			incompleteHeadInfo->returnIncompleteHeadAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		//cout << "start to get peAlignInfo" << endl;
		peAlignInfo->incompleteHead_replaceAndPushBackPathInfo2PeAlignInfo_new(pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteHeadInfo->returnMidPartLength(), indexInfo, tmpAlignmentInfo);

		//cout << "peAlignInfo: " << endl << peAlignInfo->getTmpAlignInfo(
		//	(peReadInfo.readInfo_pe1).readName, (peReadInfo.readInfo_pe2).readName, 
		//	(peReadInfo.readInfo_pe1).readSeq, (peReadInfo.readInfo_pe2).readSeq, "*", "*");

		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete seg2ndOriInfo;
		delete incompleteHeadInfo;
		return;
	}

	#endif
	
	void fixHead_shortAnchorGreedyMappingOnly(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo, 
		Alignment_Info* tmpAlignmentInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq)
	{
		//cout << "start to fixHead_incomplete" << endl;
		Incomplete_Long_Head* incompleteHeadInfo = new Incomplete_Long_Head();		
		incompleteHeadInfo->getIncompleteLongHeadInfoFromRecordWithAlignInfoType_new(
			tmpAlignmentInfo, indexInfo);		
		string incompleteLongHeadSeq = peReadInfo.returnIncompleteLongHeadSeq(
			tmpAlignInfoType, incompleteHeadInfo->returnUnfixedHeadLength());		
		char* incompleteLongHeadSeqChar = const_cast<char*>(incompleteLongHeadSeq.c_str());		
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();		
		int secondLevelIndexNO = incompleteHeadInfo->returnSecondLevelIndexNum() - 1;							
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO + 1))
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
			//continue;
		}	
		//cout << "start to do incompleteMapBool " << endl;
		//cout << "secondLevelIndexNO: " << secondLevelIndexNO << endl;
		bool incompleteMapBool;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongHeadSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteHeadInfo->returnUnfixedHeadLength(), indexInfo);

		//cout << "incompleteMapBool: " << incompleteMapBool << endl;
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
			incompleteHeadInfo->returnMapPosIntervalStart(),
			incompleteHeadInfo->returnMapPosIntervalEnd(), 
			incompleteHeadInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteHeadInfo->returnMidPartMapChrName());

		if(segInfo->returnSegmentNum() >= SEGMENTNUM)
		{
			delete seg2ndOriInfo;
			delete segInfo;
			delete incompleteHeadInfo;
			return;
		}

		segInfo->addMidPartSeg_incompleteHead(incompleteHeadInfo->returnMidPartLength(), 
			(incompleteHeadInfo->returnUnfixedHeadLength() + 1), 
			incompleteHeadInfo->returnMidPartMapChrInt(), incompleteHeadInfo->returnMidPartMapPosInChr(), indexInfo);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}
		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg_incompleteHead(segInfo, spliceJunctionDistanceMax); // Fix me: why not replace it by getPossiPathFromSeg_incompleteHead  
	
		pathInfo->filterPath_incompleteHead(segInfo); // filter out those path without the midPartSeg
	
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongHeadSeq(
				tmpAlignInfoType, incompleteHeadInfo->returnIncompleteHeadAndMidPartLength()), 
			incompleteHeadInfo->returnIncompleteHeadAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		//cout << "start to get peAlignInfo" << endl;
		peAlignInfo->incompleteHead_replaceAndPushBackPathInfo2PeAlignInfo_new(pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteHeadInfo->returnMidPartLength(), indexInfo, tmpAlignmentInfo);

		//cout << "peAlignInfo: " << endl << peAlignInfo->getTmpAlignInfo(
		//	(peReadInfo.readInfo_pe1).readName, (peReadInfo.readInfo_pe2).readName, 
		//	(peReadInfo.readInfo_pe1).readSeq, (peReadInfo.readInfo_pe2).readSeq, "*", "*");

		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete seg2ndOriInfo;
		delete incompleteHeadInfo;
		return;
	}

	void fixHead_shortAnchorSJ_remapping_targetMapping_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo,
		bool checkQualSeqForShortAnchorSeqToTargetMap)
	{
		//  cout << "****************************************************" << endl 
		//  	<< "start to do fixHead_shortAnchorSJ_remapping_targetMapping_new ..." << endl
		//  	<< "****************************************************" << endl;
		// cout << "start fixHead_shortAnchorSJ_remapping_targetMapping function ..." << endl;
		Unfixed_Head unfixedHeadInfo;
		unfixedHeadInfo.getUnfixedHeadInfoFromRecordWithAlignInfoType_new(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);

		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);
						
						///////////////////////////////////////////////////////////////////////////////////////////	
						/////////////////////////  try remapping with splice junction hash ////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////	
						//cout << "start SJsearchInSJhash !" << endl;
						bool spliceJunctionFoundInHash;

						if(spliceJunctionHashExists)
						{
							/////////////////////////////////////////////////////////////////////////////
							////////////////////////      string hash      //////////////////////////////
							/////////////////////////////////////////////////////////////////////////////
							spliceJunctionFoundInHash
								= unfixedHeadInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
						}
						else
						{
							spliceJunctionFoundInHash = false;
						}
						//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
						// fix me
						//spliceJunctionFoundInHash = false;

						if(spliceJunctionFoundInHash)
						{

							int tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[0].first;
							int tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[0].second;
							int tmpFirstMatchLength = tmpDonerEndPosInRead;

							int newMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[0];

							Alignment_Info* newTmpAlignInfo = new Alignment_Info();
							#ifdef DETECT_CIRCULAR_RNA
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
								tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, 
								newMismatchToAddInHead,
								(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
								(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[0], tmpAlignmentInfo);
							#else
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
								tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, newMismatchToAddInHead,
								(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
								(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[0], tmpAlignmentInfo);
							#endif
							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
							peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							for(int tmpSJposVec = 1; tmpSJposVec < (unfixedHeadInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
							{
								tmpDonerEndPosInRead = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].first;
								tmpSpliceJunctionDistance = (unfixedHeadInfo.SJposFromRemappingVec)[tmpSJposVec].second;
								tmpFirstMatchLength = tmpDonerEndPosInRead;	

								int tmpNewMismatchToAddInHead = (unfixedHeadInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

								Alignment_Info* newTmpAlignInfo = new Alignment_Info();
								#ifdef DETECT_CIRCULAR_RNA
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
									tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInHead,
									(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
									(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], tmpAlignmentInfo);								
								#else
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
									tmpFirstMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInHead,
									(unfixedHeadInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
									(unfixedHeadInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], tmpAlignmentInfo);
								#endif
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
							//continue;
							//return;
						}

						
						///////////////////////////////////////////////////////////////////////////////////////////	
						///////////////////////// remapping with splice junction hash failed //////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////	
						//cout << "start to getPossibleSJpos " << endl;
						unfixedHeadInfo.getPossibleSJpos(readSeqWithDirection, indexInfo);
						//cout << "finish getting Possible SJ pos" << endl;
						if((unfixedHeadInfo.possiGTAGpos).size() + (unfixedHeadInfo.possiCTACpos).size() == 0)
						{
						 	//cout << "no possible SJpos found !" << endl;
						 	return;
						}
						// cout << "output possiGTAGpos positions: " << endl;
						// for(int tmp = 0; tmp < unfixedHeadInfo.possiGTAGpos.size(); tmp++)
						// {
						// 	cout << unfixedHeadInfo.possiGTAGpos[tmp] << endl;
						// }
						// //cout << "output possiCTACpos positions: " << endl;
						// for(int tmp = 0; tmp < unfixedHeadInfo.possiCTACpos.size(); tmp++)
						// {
						// 	cout << unfixedHeadInfo.possiCTACpos[tmp] << endl;
						// }		
						// cout << "end of getPossibleSJpos ..." << endl;
			    		///////////////////////////////////////////////////////////////////////////////////////////
			    		///////////////////////////////   check possible SJs    ///////////////////////////////////
						int midPartMapPosSecondLevelIndexNO 
							= indexInfo->getSecondLevelIndexFromChrAndPos(unfixedHeadInfo.returnMidPartMapChrInt(), 
								unfixedHeadInfo.returnMidPartMapPosInChr());
						midPartMapPosSecondLevelIndexNO --;

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		

						if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(midPartMapPosSecondLevelIndexNO+1))
						{
							// cout << "invalid index ..." << endl;
							// cout << "unfixedHeadInfo.returnMidPartMapChrInt(): " << unfixedHeadInfo.returnMidPartMapChrInt() << endl;
							// cout << "unfixedHeadInfo.returnMidPartMapPosInChr(): " << unfixedHeadInfo.returnMidPartMapPosInChr() << endl;
							return;
						}

						//cout << "2ndLevel Index: " << midPartMapPosSecondLevelIndexNO << endl;		
			    		unsigned int midPartMapPosForLongHeadInSecondLevelIndex = unfixedHeadInfo.returnMidPartMapPosInChr() - 
							((unfixedHeadInfo.returnMidPartMapPosInChr())/(
								//indexInfo->secondLevelIndexNormalSize
								indexInfo->returnSecondLevelIndexNormalSize()
								))*(indexInfo->returnSecondLevelIndexNormalSize());

						//   check GTAG splice junctions
						//cout << "unfixedHeadInfo.possiGTAGpos.size(): " << (unfixedHeadInfo.possiGTAGpos).size() << endl;
						for(int tmp = 0; tmp < (unfixedHeadInfo.possiGTAGpos).size(); tmp++)
						{
							int tmpSJposInRead = unfixedHeadInfo.possiGTAGpos[tmp];

							if(tmpSJposInRead-1 < min_anchor_length)
								continue;
							if(checkQualSeqForShortAnchorSeqToTargetMap) // check short anchor seq confident or not
							{
								bool confidenceInShortAnchorSeq2TargetMapping_bool 
									= peReadInfo.checkConfidenceInShortAnchorHeadSeq(tmpAlignInfoType, tmpSJposInRead-1);
								if(!confidenceInShortAnchorSeq2TargetMapping_bool)
									continue;
							}	
							string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
							string targetMappingStr = tmpShortAnchorStr + "GT";
							//cout << "targetStrLen: " << tmpShortAnchorStr.length() << endl;
							char* headChar = const_cast<char*>(targetMappingStr.c_str());
							int targetMappingNum = 0;
							unsigned int targetMappingLoc[100];
						
							unsigned int finalMidPartMappingPos
								= unfixedHeadInfo.returnMidPartMapPosInWholeGenome() + tmpSJposInRead - unfixedHeadInfo.returnUnfixedHeadLength() - 1;
							bool headSegMapMain;
							
							headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.returnMidPartMapPosInWholeGenome(), 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							//cout << "headSegMapMain: " << headSegMapMain << endl;
							int tmp_MaxSpliceDistance_TargetMapping = this->returnMaxSpliceDistanceTargetMapping(targetMappingStr.length()-2);
							if(headSegMapMain)
							{
								int tmpMinSpliceDistance = 300001;
								int tmpMinSpliceDistance_abs = 300001;
								for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
								{
									int tmpSpliceDistance = finalMidPartMappingPos - (targetMappingLoc[tmp2] + unfixedHeadInfo.possiGTAGpos[tmp] - 2) - 1;
									
									int tmpSpliceDistance_abs = tmpSpliceDistance;
									if(tmpSpliceDistance < 0)
									{
										tmpSpliceDistance_abs = 0 - tmpSpliceDistance;
									}

									if((tmpSpliceDistance < tmp_MaxSpliceDistance_TargetMapping) 
										&& (tmpSpliceDistance > -4))
									{
										if(tmpSpliceDistance_abs < tmpMinSpliceDistance_abs)
										{
											tmpMinSpliceDistance_abs = tmpSpliceDistance_abs;
											tmpMinSpliceDistance = tmpSpliceDistance;
										}
									}						
								}	
								if(tmpMinSpliceDistance_abs < tmp_MaxSpliceDistance_TargetMapping)
								{
									(unfixedHeadInfo.GTAGsjPos).push_back(pair<int,int>(unfixedHeadInfo.possiGTAGpos[tmp], tmpMinSpliceDistance));
									(unfixedHeadInfo.GTAGsjPos_mismatch).push_back(unfixedHeadInfo.possiGTAGpos_mismatch[tmp]);
									//if(STORE_MISMATCH_POS)
									//{
										(unfixedHeadInfo.GTAGsjPos_mismatchPos).push_back(unfixedHeadInfo.possiGTAGpos_mismatchPos[tmp]);
										//if(STORE_MISMATCH_CHA)
										//{
											(unfixedHeadInfo.GTAGsjPos_mismatchChar).push_back(unfixedHeadInfo.possiGTAGpos_mismatchChar[tmp]);
										//}
									//}
								}				
							}

						}
						//   check CTAC splice junctions
						//cout << "unfixedHeadInfo.possiCTACpos.size(): " << (unfixedHeadInfo.possiCTACpos).size() << endl;
						for(int tmp = 0; tmp < (unfixedHeadInfo.possiCTACpos).size(); tmp++)
						{
							int tmpSJposInRead = unfixedHeadInfo.possiCTACpos[tmp];
						
							if(tmpSJposInRead-1 < min_anchor_length)
								continue;
							if(checkQualSeqForShortAnchorSeqToTargetMap) // check short anchor seq confident or not
							{	
								bool confidenceInShortAnchorSeq2TargetMapping_bool 
									= peReadInfo.checkConfidenceInShortAnchorHeadSeq(tmpAlignInfoType, tmpSJposInRead-1);
								if(!confidenceInShortAnchorSeq2TargetMapping_bool)
									continue;
							}
							string tmpShortAnchorStr = readSeqWithDirection.substr(0, tmpSJposInRead-1);
							string targetMappingStr = tmpShortAnchorStr + "CT";
							//cout << "targetStrLen: " << tmpShortAnchorStr.length() << endl;
							char* headChar = const_cast<char*>(targetMappingStr.c_str());
							int targetMappingNum = 0;
							unsigned int targetMappingLoc[100];
							
							unsigned int finalMidPartMappingPos
								= unfixedHeadInfo.returnMidPartMapPosInWholeGenome() + tmpSJposInRead - unfixedHeadInfo.returnUnfixedHeadLength() - 1;
							

							bool headSegMapMain;

							headSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(headChar,
									secondLevelSa[midPartMapPosSecondLevelIndexNO], 
									secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
									secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
									secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
									secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
									targetMappingStr.length(), 3000000, unfixedHeadInfo.returnMidPartMapPosInWholeGenome(), 
									midPartMapPosForLongHeadInSecondLevelIndex, &targetMappingNum, targetMappingLoc,
									indexInfo
									);
							//cout << "headSegMapMain: " << headSegMapMain << endl;
							int tmp_MaxSpliceDistance_TargetMapping = this->returnMaxSpliceDistanceTargetMapping(targetMappingStr.length()-2);
							if(headSegMapMain)
							{
								int tmpMinSpliceDistance = 300001;
								int tmpMinSpliceDistance_abs = 300001;
								for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
								{
									int tmpSpliceDistance = finalMidPartMappingPos - (targetMappingLoc[tmp2] + unfixedHeadInfo.possiCTACpos[tmp] - 2) - 1;
									
									int tmpSpliceDistance_abs = tmpSpliceDistance;
									if(tmpSpliceDistance < 0)
									{
										tmpSpliceDistance_abs = 0 - tmpSpliceDistance;
									}

									if((tmpSpliceDistance < tmp_MaxSpliceDistance_TargetMapping) 
										&& (tmpSpliceDistance > -4))
									{
										if(tmpSpliceDistance_abs < tmpMinSpliceDistance_abs)
										{
											tmpMinSpliceDistance_abs = tmpSpliceDistance_abs;
											tmpMinSpliceDistance = tmpSpliceDistance;
										}
									}						
								}		
								if(tmpMinSpliceDistance_abs < tmp_MaxSpliceDistance_TargetMapping)
								{
									(unfixedHeadInfo.CTACsjPos).push_back(pair<int,int>(unfixedHeadInfo.possiCTACpos[tmp], tmpMinSpliceDistance));
									(unfixedHeadInfo.CTACsjPos_mismatch).push_back(unfixedHeadInfo.possiCTACpos_mismatch[tmp]);
									//if(STORE_MISMATCH_POS)
									//{
										(unfixedHeadInfo.CTACsjPos_mismatchPos).push_back(unfixedHeadInfo.possiCTACpos_mismatchPos[tmp]);
										//if(STORE_MISMATCH_CHA)
										//{
											(unfixedHeadInfo.CTACsjPos_mismatchChar).push_back(unfixedHeadInfo.possiCTACpos_mismatchChar[tmp]);
										//}
									//}
								}			
							}

						}
						//cout << "unfixedHeadInfo.GTAGsjPos.size(): " << (unfixedHeadInfo.GTAGsjPos).size() << endl;
						//cout << "unfixedHeadInfo.CTACsjPos.size(): " << (unfixedHeadInfo.CTACsjPos).size() << endl;
						if((unfixedHeadInfo.GTAGsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();

							int tmpNewMismatchToAddInHead = unfixedHeadInfo.GTAGsjPos_mismatch[0];
							#ifdef DETECT_CIRCULAR_RNA
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
								(unfixedHeadInfo.GTAGsjPos[0]).first - 1, (unfixedHeadInfo.GTAGsjPos[0]).second, 
								indexInfo, tmpNewMismatchToAddInHead,
								unfixedHeadInfo.GTAGsjPos_mismatchPos[0], unfixedHeadInfo.GTAGsjPos_mismatchChar[0], tmpAlignmentInfo);
							#else
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
								(unfixedHeadInfo.GTAGsjPos[0]).first - 1, (unfixedHeadInfo.GTAGsjPos[0]).second, 
								indexInfo, tmpNewMismatchToAddInHead,
								unfixedHeadInfo.GTAGsjPos_mismatchPos[0], unfixedHeadInfo.GTAGsjPos_mismatchChar[0], tmpAlignmentInfo);
							#endif
							//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
							peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);

							for(int tmp = 1; tmp < (unfixedHeadInfo.GTAGsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.GTAGsjPos_mismatch[tmp]; //0;

								Alignment_Info* newTmpAlignInfo = new Alignment_Info();  
								#ifdef DETECT_CIRCULAR_RNA
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
									(unfixedHeadInfo.GTAGsjPos[tmp]).first - 1, (unfixedHeadInfo.GTAGsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.GTAGsjPos_mismatchPos[tmp], unfixedHeadInfo.GTAGsjPos_mismatchChar[tmp], tmpAlignmentInfo);								
								#else
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
									(unfixedHeadInfo.GTAGsjPos[tmp]).first - 1, (unfixedHeadInfo.GTAGsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.GTAGsjPos_mismatchPos[tmp], unfixedHeadInfo.GTAGsjPos_mismatchChar[tmp], tmpAlignmentInfo);
								#endif
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
							
							for(int tmp = 0; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[tmp]; //0;

								Alignment_Info* newTmpAlignInfo = new Alignment_Info(); 
								#ifdef DETECT_CIRCULAR_RNA
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.CTACsjPos_mismatchPos[tmp], unfixedHeadInfo.CTACsjPos_mismatchChar[tmp], tmpAlignmentInfo);
								#else
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.CTACsjPos_mismatchPos[tmp], unfixedHeadInfo.CTACsjPos_mismatchChar[tmp], tmpAlignmentInfo);
								#endif
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}
						}
						else if((unfixedHeadInfo.CTACsjPos).size() > 0)
						{
							Alignment_Info* newTmpAlignInfo = new Alignment_Info();

							int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[0]; //0;
							#ifdef DETECT_CIRCULAR_RNA
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
								(unfixedHeadInfo.CTACsjPos[0]).first - 1, (unfixedHeadInfo.CTACsjPos[0]).second, 
								indexInfo, tmpNewMismatchToAddInHead,
								unfixedHeadInfo.CTACsjPos_mismatchPos[0], unfixedHeadInfo.CTACsjPos_mismatchChar[0], tmpAlignmentInfo);
							#else
							newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
								(unfixedHeadInfo.CTACsjPos[0]).first - 1, (unfixedHeadInfo.CTACsjPos[0]).second, 
								indexInfo, tmpNewMismatchToAddInHead,
								unfixedHeadInfo.CTACsjPos_mismatchPos[0], unfixedHeadInfo.CTACsjPos_mismatchChar[0], tmpAlignmentInfo);
							#endif							
							//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
							// peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
							peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							for(int tmp = 1; tmp < (unfixedHeadInfo.CTACsjPos).size(); tmp++)
							{
								int tmpNewMismatchToAddInHead = unfixedHeadInfo.CTACsjPos_mismatch[tmp];

								Alignment_Info* newTmpAlignInfo = new Alignment_Info();
								#ifdef DETECT_CIRCULAR_RNA
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.CTACsjPos_mismatchPos[tmp], unfixedHeadInfo.CTACsjPos_mismatchChar[tmp], tmpAlignmentInfo);
								#else
								newTmpAlignInfo->newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
									(unfixedHeadInfo.CTACsjPos[tmp]).first - 1, (unfixedHeadInfo.CTACsjPos[tmp]).second, 
									indexInfo, tmpNewMismatchToAddInHead,
									unfixedHeadInfo.CTACsjPos_mismatchPos[tmp], unfixedHeadInfo.CTACsjPos_mismatchChar[tmp], tmpAlignmentInfo);
								#endif
								//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
								peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
							}				
						}
						else
						{}

						// cout << "start to get unfixedHeadInfo.GTAGsjPos " << endl;
						// cout << "output GTAG positions: " << endl;
						// for(int tmp = 0; tmp < unfixedHeadInfo.GTAGsjPos.size(); tmp++)
						// {
						// 	cout << (unfixedHeadInfo.GTAGsjPos[tmp]).first << "," << (unfixedHeadInfo.GTAGsjPos[tmp]).second << endl;
						// }
						// cout << "output CTAC positions: " << endl;
						// for(int tmp = 0; tmp < unfixedHeadInfo.CTACsjPos.size(); tmp++)
						// {
						// 	cout << (unfixedHeadInfo.CTACsjPos[tmp]).first << "," << (unfixedHeadInfo.CTACsjPos[tmp]).second << endl;
						// }		
						// cout << "end of get unfixedHeadInfo.CTACsjPos ..." << endl;


						// cout //<< "****************************************************" << endl 
						// 	<< "end of fixHead_shortAnchorSJ_remapping_targetMapping_new ..." << endl
						// 	<< "****************************************************" << endl ;					
	} 

	void fixHead_extend2end_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		int headSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;

		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);		
		int midPartMapPosInChr = tmpAlignmentInfo->returnAlignChromPos();
		
		if(midPartMapPosInChr - headSeq_len < 1) // start pos < 0
			return;
		
		vector<int> tmpMismatchPosVec;
		vector<char> tmpMismatchCharVec;		
		bool matchBool;
		int extension_length = extensionBackwards_errorTolerated(midPartMapChrInt, 
			midPartMapPosInChr - 1, indexInfo, 
			peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType),
			headSeq_len, 1, tmpMismatchPosVec, tmpMismatchCharVec, MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH, MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH);
		if(extension_length == 0)
			matchBool = false;
		else
			matchBool = true;
		if(matchBool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new(
					indexInfo, headSeq_len, extension_length, tmpMismatchPosVec, tmpMismatchCharVec, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, 
				tmpAlignInfoType, tmpIndex_peAlignInfo);			
		}
		else
		{}
		return;
	}

	void fixHead_extend2end_new_fixIndelBesideReadEnd(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo, 
		int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		//cout << "fixHead_extend2end_new_fixIndelBesideReadEnd starts ......" << endl;
		int headSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;

		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);		
		int midPartMapPosInChr = tmpAlignmentInfo->returnAlignChromPos();
		
		//cout << "headSeq_len: " << headSeq_len << endl;
		//cout << "midPartMapChrName: " << midPartMapChrName << endl;
		//cout << "midPartMapPosInChr: " << midPartMapPosInChr << endl;

		if(midPartMapPosInChr - headSeq_len - 3 < 1) // start pos < 0
			return;

		if(headSeq_len < 3)
			return;

		string first3baseSeqInRead = (peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)).substr(0,3);
		int correspondingMapPos_firstBaseInRead = midPartMapPosInChr - headSeq_len;
		int possibleIndelLen[6] = {1,-1,2,-2,3,-3}; // <0 means deletion; >0 means insertion, in this method, max_ins/del <= 3
		
		//cout << "first3baseSeqInRead: " << first3baseSeqInRead << endl;
		//cout << "correspondingMapPos_firstBaseInRead: " << correspondingMapPos_firstBaseInRead << endl;

		vector<int> validIndelLenVec;
		
		for(int tmp = 0; tmp < 6; tmp++)
		{
			if((3 + possibleIndelLen[tmp]) <= headSeq_len)
			{
				//cout << "tmpFirstBaseMapPos: " << correspondingMapPos_firstBaseInRead + possibleIndelLen[tmp] << endl;
				int tmpFirstBaseMapPos = correspondingMapPos_firstBaseInRead + possibleIndelLen[tmp];
				// if(tmpFirstBaseMapPos <= 0)
				// 	continue;
				string tmpCorrespondingFirst3baseSeqInChr 
					= indexInfo->returnChromStrSubstr(midPartMapChrInt, tmpFirstBaseMapPos, 3);
				//cout << "tmp: " << tmp << endl;
				//cout << "tmpIndel: " << possibleIndelLen[tmp] << endl;
				//cout << "targetSeq: " << tmpCorrespondingFirst3baseSeqInChr << endl;
				if(tmpCorrespondingFirst3baseSeqInChr == first3baseSeqInRead)
				{
					//cout << "match !" << endl;
					validIndelLenVec.push_back(possibleIndelLen[tmp]);
				}
				else
				{
					//cout << "unmatch !" << endl;
				}
			}
		}
		bool fixIndel_bool = false;
		vector<Jump_Code> fixedIndelJumpCodeVec;
		int fixedIndelChromMapPos;
		for(int tmp = 0; tmp < validIndelLenVec.size(); tmp++)
		{
			int tmpIndelLen = validIndelLenVec[tmp];	
			//cout << "tmpIndelLen: " << tmpIndelLen << endl;
			if((tmpIndelLen) < 0) // deletion
			{
				//cout << "start to try deletion ...." << endl;
				int deletionLength = 0 - tmpIndelLen;
				if(headSeq_len-3 < 2)
				{
					int firstMatchLength = 3;
					int secondMatchLength = headSeq_len - 3;					
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code deletionJumpCode(deletionLength, "D");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
					fixedIndelJumpCodeVec.push_back(deletionJumpCode);
					fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
					fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
					fixIndel_bool = true;
					break;	
				}
				else
				{
					FixDoubleAnchor_Deletion_Info* delInfo = new FixDoubleAnchor_Deletion_Info();
					int tmp_toFix_deletion_read_start = 4;
					int tmp_toFix_deletion_read_end = headSeq_len+1;
					int subSeqLengthInProcess = tmp_toFix_deletion_read_end - tmp_toFix_deletion_read_start + 1;
					int tmp_toFix_deletion_chrom_start = midPartMapPosInChr - headSeq_len + tmpIndelLen + 3;
					int tmp_toFix_deletion_chrom_end = midPartMapPosInChr;
					int tmp_max_allowed_mismatchNum = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed_INDEL_READ_END; // add mismatches into alignInfo
					bool deletion_fixed = delInfo->detectBestDeletion_lessMismatch(tmp_toFix_deletion_read_start, 
						tmp_toFix_deletion_read_end, tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, 
						(peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)), indexInfo,
						midPartMapChrInt, tmp_max_allowed_mismatchNum);
					if(deletion_fixed)
					{
						int firstMatchLength = delInfo->return_best_deletion_prefix_match_length() + 3;
						int secondMatchLength = headSeq_len - firstMatchLength;
						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code deletionJumpCode(deletionLength, "D");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M");
						fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
						fixedIndelJumpCodeVec.push_back(deletionJumpCode);
						fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
						fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
						fixIndel_bool = true;
						delete delInfo;
						break;
					}
					else
					{
						delete delInfo;
					}
				}
			}
			else // insertion
			{
				//cout << "start to try insertion ...." << endl;
				if(headSeq_len-3 < tmpIndelLen)
				{
					int firstMatchLength = 3;
					int insertionLength = tmpIndelLen;
					int secondMatchLength = headSeq_len - 3- insertionLength;
					
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
					fixedIndelJumpCodeVec.push_back(insertionJumpCode);
					fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
					fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
					fixIndel_bool = true;
					break;					
				}
				else
				{	
					//cout << "fixInsInfo starts ...." << endl;
					FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info();
					int tmp_toFix_insertion_read_start = 4;
					int tmp_toFix_insertion_read_end = headSeq_len + 1;
					int tmp_toFix_insertion_chrom_start = midPartMapPosInChr - headSeq_len + tmpIndelLen + 3;
					int tmp_toFix_insertion_chrom_end = midPartMapPosInChr;
					int tmp_max_allowed_mismatchNum = headSeq_len/LengthOfSeqPerMismatchAllowed_INDEL_READ_END; // add mismatches into alignInfo, 
					//cout << "tmp_toFix_insertion_read_start: " << tmp_toFix_insertion_read_start << endl;
					//cout << "tmp_toFix_insertion_read_end: " << tmp_toFix_insertion_read_end << endl;
					//cout << "tmp_toFix_insertion_chrom_start: " << tmp_toFix_insertion_chrom_start << endl;
					//cout << "tmp_toFix_insertion_chrom_end: " << tmp_toFix_insertion_chrom_end << endl;
					bool fixInsertionBool = insInfo->detectBestInsertion_lessMismatch(
						tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
						tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, 
						(peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)), indexInfo,
						midPartMapChrInt, tmp_max_allowed_mismatchNum);
					//cout << "fixInsertionBool: " << fixInsertionBool << endl;
					int firstMatchLength = insInfo->return_best_insertion_prefix_match_length() + 3;
					int insertionLength = tmpIndelLen;
					int secondMatchLength = headSeq_len - firstMatchLength - insertionLength;
					//cout << "firstMatchLength: " << firstMatchLength << endl;
					//cout << "insertionLength: " << insertionLength << endl;
					//cout << "secondMatchLength: " << secondMatchLength << endl;
					if(fixInsertionBool)
					{
						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code insertionJumpCode(insertionLength, "I");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M");
						fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
						fixedIndelJumpCodeVec.push_back(insertionJumpCode);
						fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
						fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
						fixIndel_bool = true;
						delete insInfo;
						break;
					}
					else
					{
						delete insInfo;
					}
				}
			}
		}
		
		//cout << "fixIndel_bool: " << fixIndel_bool << endl;

		if(fixIndel_bool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2Head_fixIndelBesideReadEnd(
					indexInfo, fixedIndelJumpCodeVec, 
					fixedIndelChromMapPos, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, 
				tmpAlignInfoType, tmpIndex_peAlignInfo);				
		}
		return;
	}

	void fixHead_extend2end_new_finalStepForAligner(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo, 
		int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		int headSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);		
		int midPartMapPosInChr = tmpAlignmentInfo->returnAlignChromPos();
		if(midPartMapPosInChr - headSeq_len < 1) // start pos < 0
			return;		
		vector<int> tmpMismatchPosVec;
		vector<char> tmpMismatchCharVec;		
		bool matchBool;
		int extension_length = extensionBackwards_errorTolerated_finalStepForAligner(midPartMapChrInt, 
			midPartMapPosInChr - 1, indexInfo, 
			peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType),
			headSeq_len, 1, tmpMismatchPosVec, tmpMismatchCharVec, MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH, MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH);
		//cout << "final extension_length: " << extension_length << endl;
		if(extension_length == 0)
			matchBool = false;
		else
			matchBool = true;
		//cout << "matchBool: " << matchBool << endl;
		if(matchBool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new(
					indexInfo, headSeq_len, extension_length, tmpMismatchPosVec, tmpMismatchCharVec, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, 
				tmpAlignInfoType, tmpIndex_peAlignInfo);			
		}
		else
		{}
		return;
	}

	void fixHead_extend2end_nwdpBesideReadEnd(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo, 
		int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		
	}

	void fixTail_shortAnchorRemappingOnly_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		int readLength = peReadInfo.returnReadLength(tmpAlignInfoType);	

		//cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
		Unfixed_Tail unfixedTailInfo;
		unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, 
			tmpAlignmentInfo, indexInfo);
						
		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);

		///////////////////////////////////////////////////////////////////////////////////////////	
		/////////////////////////  try remapping with splice junction hash ////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start SJsearchInSJhash " << endl;
		bool spliceJunctionFoundInHash;
		
		if(spliceJunctionHashExists)
		{
			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////     String hash     //////////////////////////////					
			spliceJunctionFoundInHash 
				= unfixedTailInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
		}
		else
		{							
			spliceJunctionFoundInHash = false;
		}
		// fix me
		//spliceJunctionFoundInHash = false;
		//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
		if(spliceJunctionFoundInHash)
		{
			int tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[0].first;
			int tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[0].second;
			int tmpLastMatchLength = readLength - tmpDonerEndPosInRead;
			//cout << "*******************" << endl;
			//cout << "tmpDonerEndPosInRead: " << tmpDonerEndPosInRead << endl;
			//cout << "tmpSpliceJunctionDistance: " << tmpSpliceJunctionDistance << endl;
			//cout << "tmpLastMatchLength: " << tmpLastMatchLength << endl;
			int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[0];

			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			#ifdef DETECT_CIRCULAR_RNA
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
					tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
					(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
					(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[0], readLength, tmpAlignmentInfo);
			#else
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
					tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
					(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
					(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[0], readLength, tmpAlignmentInfo);
			#endif
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo );
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmpSJposVec = 1; tmpSJposVec < (unfixedTailInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
			{
				int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

				tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].first;
				tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].second;
				tmpLastMatchLength = readLength - tmpDonerEndPosInRead;		

			//cout << "*******************" << endl;
			//cout << "tmpDonerEndPosInRead: " << tmpDonerEndPosInRead << endl;
			//cout << "tmpSpliceJunctionDistance: " << tmpSpliceJunctionDistance << endl;
			//cout << "tmpLastMatchLength: " << tmpLastMatchLength << endl;

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
						tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
						(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
						(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], readLength, tmpAlignmentInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
						tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
						(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
						(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], readLength, tmpAlignmentInfo);
				#endif
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
			return;
		}
		else // can not find candidate SJ in hash
		{}
	} 

	void fixTail_shortAnchorRemappingOnly_withAlignInfer(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJinfo, 
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignInfo,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{
		int readLength = peReadInfo.returnReadLength(tmpAlignInfoType);	

		//cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
		Unfixed_Tail unfixedTailInfo;
		unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, 
			tmpAlignInfo, indexInfo);
						
		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);

		///////////////////////////////////////////////////////////////////////////////////////////	
		/////////////////////////  try remapping with splice junction hash ////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start SJsearchInSJhash for unfixedTail " << endl;
		bool spliceJunctionFoundInHash;
		vector< vector<Jump_Code> > inferedPathJumpCodeVecVec;
		vector< vector<int> > inferedPathMismatchPosVecVec;
		vector< vector<char> > inferedPathMismatchCharVecVec;		
		if(spliceJunctionHashExists)
		{
			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////     String hash     //////////////////////////////					
			spliceJunctionFoundInHash 
				= unfixedTailInfo.SJsearchInSJhash_areaStringHash_withAlignInferJuncHash(
					SJinfo, alignInferJunctionHashInfo,
					readSeqWithDirection, indexInfo, SJinfo->areaSize,
					inferedPathJumpCodeVecVec,
					inferedPathMismatchPosVecVec,
					inferedPathMismatchCharVecVec);
		}
		else
		{							
			spliceJunctionFoundInHash = false;
		}
		// fix me
		//spliceJunctionFoundInHash = false;
		//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
		if(spliceJunctionFoundInHash)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_withAlignInfer(
				inferedPathJumpCodeVecVec[0], 
				inferedPathMismatchPosVecVec[0],
				inferedPathMismatchCharVecVec[0],
				indexInfo, tmpAlignInfo);
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo );
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);

			for(int tmp = 1; tmp < inferedPathJumpCodeVecVec.size(); tmp++)
			{
				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_withAlignInfer(
					inferedPathJumpCodeVecVec[tmp], 
					inferedPathMismatchPosVecVec[tmp],
					inferedPathMismatchCharVecVec[tmp],
					indexInfo, tmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
			//continue;
			return;
		}
		else // can not find candidate SJ in hash
		{}
	} 

	#ifdef PERSONALIZED_CHR_SEQ
	void fixTail_shortAnchorGreedyMappingOnly_includeSNPseqMap(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo, 
		bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, int SNPlocInSyntheticSNPseq	
		)
	{
		Incomplete_Long_Tail* incompleteTailInfo = new Incomplete_Long_Tail();

		incompleteTailInfo->getIncompleteLongTailInfoFromRecordWithAlignInfoType_new(
			peReadInfo, tmpAlignInfoType,
			tmpAlignmentInfo, indexInfo);

		string incompleteLongTailSeq = peReadInfo.returnIncompleteLongTailSeq(
			tmpAlignInfoType, incompleteTailInfo->returnUnfixedTailLength());
	
		char* incompleteLongTailSeqChar = const_cast<char*>(incompleteLongTailSeq.c_str());

		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

		int secondLevelIndexNO = incompleteTailInfo->returnSecondLevelIndexNum() - 1;

		if(
			indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1)
			)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
			//continue;
		}		
		bool incompleteMapBool;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongTailSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteTailInfo->returnUnfixedTailLength(), indexInfo);
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		Seg_Info* segInfo_old = new Seg_Info(seg2ndOriInfo, 
			incompleteTailInfo->returnMapPosIntervalStart(),
			incompleteTailInfo->returnMapPosIntervalEnd(), 
			incompleteTailInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteTailInfo->returnMidPartMapChrName());
		
		#ifdef PERSONALIZED_CHR_SEQ
		Seg_Info segInfo_SNPseq;
		int tmpIncompleteSeqLength = incompleteLongTailSeq.length();
		bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(incompleteLongTailSeqChar, 
			sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
		if(tmpMapWithSnpSeqIndexBool)
			segInfo_old->update_includeSNPseqMapSegInfo(segInfo_SNPseq, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		
		segInfo_old->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
			indexInfo_snp, incompleteLongTailSeq, indexInfo, SNPlocInSyntheticSNPseq);
		#endif

		if(segInfo_old->returnSegmentNum() >= SEGMENTNUM)
		{
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}	
	
		Seg_Info* segInfo = new Seg_Info();

		segInfo->addMidPartSeg_incompleteTail(incompleteTailInfo->returnMidPartLength(), 
			incompleteTailInfo->returnMidPartLocInRead(), 
			incompleteTailInfo->returnMidPartMapChrInt(), 
			incompleteTailInfo->returnMidPartMapPosInChr(), indexInfo, segInfo_old);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}

		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
	
		Path_Info* pathInfo = new Path_Info();
	
		pathInfo->getPossiPathFromSeg_incompleteTail(segInfo, spliceJunctionDistanceMax);  

		pathInfo->filterPath_incompleteTail();//segInfo); // filter out those path without the midPartSeg

		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		//cout << "start to fix gaps" << endl;
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongTailSeq(
				tmpAlignInfoType, incompleteTailInfo->returnIncompleteTailAndMidPartLength()), 
			incompleteTailInfo->returnIncompleteTailAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		peAlignInfo->incompleteTail_replaceAndPushBackPathInfo2PeAlignInfo_new(
			pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteTailInfo->returnMidPartLength(), 
			indexInfo, tmpAlignmentInfo, incompleteTailInfo->returnMidPartLocInRead());

		// stop4
		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete segInfo_old;
		delete seg2ndOriInfo;
		delete incompleteTailInfo;
		return;
	}
	#endif

	void fixTail_shortAnchorGreedyMappingOnly(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo, 
		bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq)
	{
		Incomplete_Long_Tail* incompleteTailInfo = new Incomplete_Long_Tail();

		incompleteTailInfo->getIncompleteLongTailInfoFromRecordWithAlignInfoType_new(
			peReadInfo, tmpAlignInfoType,
			tmpAlignmentInfo, indexInfo);

		string incompleteLongTailSeq = peReadInfo.returnIncompleteLongTailSeq(
			tmpAlignInfoType, incompleteTailInfo->returnUnfixedTailLength());
	
		char* incompleteLongTailSeqChar = const_cast<char*>(incompleteLongTailSeq.c_str());

		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

		int secondLevelIndexNO = incompleteTailInfo->returnSecondLevelIndexNum() - 1;

		if(
			indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1)
			)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
			//continue;
		}		
		bool incompleteMapBool;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongTailSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteTailInfo->returnUnfixedTailLength(), indexInfo);
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		Seg_Info* segInfo_old = new Seg_Info(seg2ndOriInfo, 
			incompleteTailInfo->returnMapPosIntervalStart(),
			incompleteTailInfo->returnMapPosIntervalEnd(), 
			incompleteTailInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteTailInfo->returnMidPartMapChrName());
		
		if(segInfo_old->returnSegmentNum() >= SEGMENTNUM)
		{
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}	
	
		Seg_Info* segInfo = new Seg_Info();

		segInfo->addMidPartSeg_incompleteTail(incompleteTailInfo->returnMidPartLength(), 
			incompleteTailInfo->returnMidPartLocInRead(), 
			incompleteTailInfo->returnMidPartMapChrInt(), 
			incompleteTailInfo->returnMidPartMapPosInChr(), indexInfo, segInfo_old);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}

		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
	
		Path_Info* pathInfo = new Path_Info();
	
		pathInfo->getPossiPathFromSeg_incompleteTail(segInfo, spliceJunctionDistanceMax);  

		pathInfo->filterPath_incompleteTail();//segInfo); // filter out those path without the midPartSeg

		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		//cout << "start to fix gaps" << endl;
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongTailSeq(
				tmpAlignInfoType, incompleteTailInfo->returnIncompleteTailAndMidPartLength()), 
			incompleteTailInfo->returnIncompleteTailAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		peAlignInfo->incompleteTail_replaceAndPushBackPathInfo2PeAlignInfo_new(
			pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteTailInfo->returnMidPartLength(), 
			indexInfo, tmpAlignmentInfo, incompleteTailInfo->returnMidPartLocInRead());

		// stop4
		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete segInfo_old;
		delete seg2ndOriInfo;
		delete incompleteTailInfo;
		return;
	}

	void fixTail_shortAnchorSJ_remapping_targetMapping_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo,
		bool checkQualSeqForShortAnchorSeqToTargetMap)
	{
		// cout << "****************************************************" << endl 
		// 	<< "start to do fixTail_shortAnchorSJ_remapping_targetMapping_new ..." << endl;
		// 	//<< "****************************************************" << endl;

		int readLength = peReadInfo.returnReadLength(tmpAlignInfoType);	

		//cout << "start to getUnfixedTailInfoFromRecordWithAlignInfoType !" << endl;
		Unfixed_Tail unfixedTailInfo;
		unfixedTailInfo.getUnfixedTailInfoFromRecordWithAlignInfoType(peReadInfo, tmpAlignInfoType, tmpAlignmentInfo, indexInfo);
						
		string readSeqWithDirection = peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType);

		///////////////////////////////////////////////////////////////////////////////////////////	
		/////////////////////////  try remapping with splice junction hash ////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start SJsearchInSJhash " << endl;
		bool spliceJunctionFoundInHash;

		if(spliceJunctionHashExists)
		{
			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////     String hash     //////////////////////////////
			/////////////////////////////////////////////////////////////////////////////						
			spliceJunctionFoundInHash 
				= unfixedTailInfo.SJsearchInSJhash_areaStringHash(SJ, readSeqWithDirection, indexInfo, SJ->areaSize);
		}
		else
		{							
			spliceJunctionFoundInHash = false;
		}
		// fix me
		//spliceJunctionFoundInHash = false;
		//cout << "spliceJunctionFoundInHash: " << spliceJunctionFoundInHash << endl;
		if(spliceJunctionFoundInHash)
		{
			//if((unfixedTailInfo.SJposFromRemapping).size() == 1)
			int tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[0].first;
			int tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[0].second;
			int tmpLastMatchLength = readLength - tmpDonerEndPosInRead;

			int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[0];

			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			#ifdef DETECT_CIRCULAR_RNA
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
				tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
				(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
				(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[0], readLength, tmpAlignmentInfo);
			#else
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
				tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
				(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[0], 
				(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[0], readLength, tmpAlignmentInfo);			
			#endif
			//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo );
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmpSJposVec = 1; tmpSJposVec < (unfixedTailInfo.SJposFromRemappingVec).size(); tmpSJposVec++)
			{
				int tmpNewMismatchToAddInTail = (unfixedTailInfo.SJposFromRemappingVec_mismatch)[tmpSJposVec];

				tmpDonerEndPosInRead = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].first;
				tmpSpliceJunctionDistance = (unfixedTailInfo.SJposFromRemappingVec)[tmpSJposVec].second;
				tmpLastMatchLength = readLength - tmpDonerEndPosInRead;		

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
					tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
					(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
					(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], readLength, tmpAlignmentInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
					tmpLastMatchLength, tmpSpliceJunctionDistance, indexInfo, tmpNewMismatchToAddInTail,
					(unfixedTailInfo.SJposFromRemappingVec_mismatchPosVec)[tmpSJposVec], 
					(unfixedTailInfo.SJposFromRemappingVec_mismatchCharVec)[tmpSJposVec], readLength, tmpAlignmentInfo);				
				#endif
				//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
			//return;
			//continue;
		}
						
		///////////////////////////////////////////////////////////////////////////////////////////	
		///////////////////////// remapping with splice junction hash failed //////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////	
		//cout << "start to get possible SJ pos" << endl; 
		unfixedTailInfo.getPossibleSJpos(readSeqWithDirection, indexInfo->returnChromString(), indexInfo);
		if((unfixedTailInfo.possiGTAGpos).size() + (unfixedTailInfo.possiCTACpos).size() == 0)
		{
			//cout << "no possible SJpos found !" << endl;
			return;
		}
		//cout << "output possiGTAGpos positions: " << endl;
		// for(int tmp = 0; tmp < unfixedTailInfo.possiGTAGpos.size(); tmp++)
		// {
		// 	cout << unfixedTailInfo.possiGTAGpos[tmp] << endl;
		// }
		// //cout << "output possiCTACpos positions: " << endl;
		// for(int tmp = 0; tmp < unfixedTailInfo.possiGTAGpos.size(); tmp++)
		// {
		// 	cout << unfixedTailInfo.possiCTACpos[tmp] << endl;
		// }		
		//cout << "end of getPossibleSJpos ..." << endl;

		//cout << "finish getting possible SJ pos " << endl;
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////   check possible SJs    ///////////////////////////////////			

		int midPartMapPosSecondLevelIndexNO
			= indexInfo->getSecondLevelIndexFromChrAndPos(unfixedTailInfo.returnMidPartMapChrInt(), 
					unfixedTailInfo.returnMidPartMapPosInChr());
		midPartMapPosSecondLevelIndexNO --;
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(midPartMapPosSecondLevelIndexNO+1))
		{
			// cout << "invalid index ...." << endl;
			// cout << "unfixedTailInfo.returnMidPartMapChrInt(): " << unfixedTailInfo.returnMidPartMapChrInt() << endl;
			// cout << "unfixedTailInfo.returnMidPartMapPosInChr(): " << unfixedTailInfo.returnMidPartMapPosInChr() << endl;
			return;
		}

		unsigned int midPartMapPosForLongTailInSecondLevelIndex = unfixedTailInfo.returnMidPartMapPosInChr() - 
			((unfixedTailInfo.returnMidPartMapPosInChr())/(indexInfo->returnSecondLevelIndexNormalSize()))*(indexInfo->returnSecondLevelIndexNormalSize());

		//cout << "start to check GTAG splice junctions " << endl;
		for(int tmp = 0; tmp < (unfixedTailInfo.possiGTAGpos).size(); tmp++)
		{
			int tmpSJposInRead = unfixedTailInfo.possiGTAGpos[tmp];

			//cout << "SJposInRead: " << tmpSJposInRead << endl;
			if(readLength - tmpSJposInRead + 1 < min_anchor_length)
				continue;
			if(checkQualSeqForShortAnchorSeqToTargetMap) // check short anchor seq confident or not
			{	
				bool confidenceInShortAnchorSeq2TargetMapping_bool 
					= peReadInfo.checkConfidenceInShortAnchorTailSeq(tmpAlignInfoType, readLength - tmpSJposInRead + 1);
				if(!confidenceInShortAnchorSeq2TargetMapping_bool)
					continue;
			}
			string tmpShortAnchorStr 
				= readSeqWithDirection.substr(tmpSJposInRead-1, readLength - tmpSJposInRead + 1);

			string targetMappingStr = "AG" + tmpShortAnchorStr;
			//cout << "targetStr: AG-" << tmpShortAnchorStr << endl;
			//cout << "targetStrLen: " << tmpShortAnchorStr.length() << endl;
			char* tailChar = const_cast<char*>(targetMappingStr.c_str());
			int targetMappingNum = 0;
			unsigned int targetMappingLoc[100];
						
			unsigned int finalMidPartMappingPos
				= unfixedTailInfo.returnMidPartMapPosInWholeGenome() + tmpSJposInRead - readLength + unfixedTailInfo.returnUnfixedTailLength() - 1;
							
			bool tailSegMapMain;
			//cout << "start to do mapMainSecondLevelForTargetMapping_compressedIndex ..." << endl;
			tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
				secondLevelSa[midPartMapPosSecondLevelIndexNO],
				secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
				secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
				secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
				secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
				targetMappingStr.length(), 3000000, unfixedTailInfo.returnMidPartMapPosInWholeGenome(), 
				midPartMapPosForLongTailInSecondLevelIndex, &targetMappingNum, targetMappingLoc, indexInfo);								
			//cout << "tailSegMapMain: " << tailSegMapMain << endl;
			int tmp_MaxSpliceDistance_TargetMapping = this->returnMaxSpliceDistanceTargetMapping(targetMappingStr.length()-2);
			//cout << "finish to do mapMainSecondLevelForTargetMapping_compressedIndex ..." << endl;
			if(tailSegMapMain)// select the nearest short anchor location
			{
				int tmpMinSpliceDistance = 300001;
				int tmpMinSpliceDistance_abs = 300001;
				for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
				{
					int tmpSpliceDistance = targetMappingLoc[tmp2] - finalMidPartMappingPos + 1; 
									
					int tmpSpliceDistance_abs = tmpSpliceDistance;
					if(tmpSpliceDistance < 0)
					{
						tmpSpliceDistance_abs = 0 - tmpSpliceDistance;
					}

					if((tmpSpliceDistance < tmp_MaxSpliceDistance_TargetMapping) 
							&& (tmpSpliceDistance > -4))
					{
						if(tmpSpliceDistance_abs < tmpMinSpliceDistance_abs)
						{
							tmpMinSpliceDistance_abs = tmpSpliceDistance_abs;
							tmpMinSpliceDistance = tmpSpliceDistance;
						}
					}
				}
				if(tmpMinSpliceDistance_abs < tmp_MaxSpliceDistance_TargetMapping)
				{
					(unfixedTailInfo.GTAGsjPos).push_back(pair<int,int>(unfixedTailInfo.possiGTAGpos[tmp], tmpMinSpliceDistance));
					(unfixedTailInfo.GTAGsjPos_mismatch).push_back(unfixedTailInfo.possiGTAGpos_mismatch[tmp]);
					//if(STORE_MISMATCH_POS)
					//{
						(unfixedTailInfo.GTAGsjPos_mismatchPos).push_back(unfixedTailInfo.possiGTAGpos_mismatchPos[tmp]);
						//if(STORE_MISMATCH_CHA)
						//{
							(unfixedTailInfo.GTAGsjPos_mismatchChar).push_back(unfixedTailInfo.possiGTAGpos_mismatchChar[tmp]);
						//}
					//}
				}
			}
		}		
		//   check CTAC splice junctions
		//cout << "start to check CTAC splice junctions " << endl;
		for(int tmp = 0; tmp < (unfixedTailInfo.possiCTACpos).size(); tmp++)
		{
			int tmpSJposInRead = unfixedTailInfo.possiCTACpos[tmp];

			//cout << "SJposInRead: " << tmpSJposInRead << endl;
			if(readLength - tmpSJposInRead + 1 < min_anchor_length)
				continue;
			if(checkQualSeqForShortAnchorSeqToTargetMap) // check short anchor seq confident or not
			{	
				bool confidenceInShortAnchorSeq2TargetMapping_bool 
					= peReadInfo.checkConfidenceInShortAnchorTailSeq(tmpAlignInfoType, readLength - tmpSJposInRead + 1);
				if(!confidenceInShortAnchorSeq2TargetMapping_bool)
					continue;
			}
			string tmpShortAnchorStr 
				= readSeqWithDirection.substr(tmpSJposInRead-1, readLength - tmpSJposInRead + 1);

			string targetMappingStr = "AC" + tmpShortAnchorStr;
			//cout << "targetStr: AC-" << tmpShortAnchorStr << endl;
			//cout << "targetStrLen: " << tmpShortAnchorStr.length() << endl;
			char* tailChar = const_cast<char*>(targetMappingStr.c_str());
			int targetMappingNum = 0;
			unsigned int targetMappingLoc[100];
						
			unsigned int finalMidPartMappingPos
				= unfixedTailInfo.returnMidPartMapPosInWholeGenome() + tmpSJposInRead - readLength + unfixedTailInfo.returnUnfixedTailLength() - 1;

			bool tailSegMapMain;
			//cout << "start to do mapMainSecondLevelForTargetMapping_compressedIndex ..." << endl;
			tailSegMapMain = mapMainSecondLevelForTargetMapping_compressedIndex(tailChar, 
					secondLevelSa[midPartMapPosSecondLevelIndexNO],
					secondLevelLcpCompress[midPartMapPosSecondLevelIndexNO], 
					secondLevelChildTab[midPartMapPosSecondLevelIndexNO],
					secondLevelDetChild[midPartMapPosSecondLevelIndexNO],
					secondLevelChrom[midPartMapPosSecondLevelIndexNO], 
					targetMappingStr.length(), 3000000, unfixedTailInfo.returnMidPartMapPosInWholeGenome(), 
					midPartMapPosForLongTailInSecondLevelIndex, 
					&targetMappingNum, targetMappingLoc, indexInfo);								
			//cout << "tailSegMapMain: " << tailSegMapMain << endl;
			int tmp_MaxSpliceDistance_TargetMapping = this->returnMaxSpliceDistanceTargetMapping(targetMappingStr.length()-2);
			//cout << "finish do mapping" << endl;
			//cout << "finish to do mapMainSecondLevelForTargetMapping_compressedIndex ..." << endl;
			if(tailSegMapMain)
			{
				int tmpMinSpliceDistance = 300001;
				int tmpMinSpliceDistance_abs = 300001;
				for(int tmp2 = 0; tmp2 < targetMappingNum; tmp2++)
				{
					int tmpSpliceDistance = targetMappingLoc[tmp2] - finalMidPartMappingPos + 1; 
									
					int tmpSpliceDistance_abs = tmpSpliceDistance;
					if(tmpSpliceDistance < 0)
					{
						tmpSpliceDistance_abs = 0 - tmpSpliceDistance;
					}

					if((tmpSpliceDistance < tmp_MaxSpliceDistance_TargetMapping) 
						&& (tmpSpliceDistance > -4))
					{
						if(tmpSpliceDistance_abs < tmpMinSpliceDistance_abs)
						{
							tmpMinSpliceDistance_abs = tmpSpliceDistance_abs;
							tmpMinSpliceDistance = tmpSpliceDistance;
						}
					}
				}
				if(tmpMinSpliceDistance_abs < tmp_MaxSpliceDistance_TargetMapping)
				{
					(unfixedTailInfo.CTACsjPos).push_back(pair<int,int>(unfixedTailInfo.possiCTACpos[tmp], tmpMinSpliceDistance));
					(unfixedTailInfo.CTACsjPos_mismatch).push_back(unfixedTailInfo.possiCTACpos_mismatch[tmp]);
					//if(STORE_MISMATCH_POS)
					//{
						(unfixedTailInfo.CTACsjPos_mismatchPos).push_back(unfixedTailInfo.possiCTACpos_mismatchPos[tmp]);
						//if(STORE_MISMATCH_CHA)
						//{
							(unfixedTailInfo.CTACsjPos_mismatchChar).push_back(unfixedTailInfo.possiCTACpos_mismatchChar[tmp]);
						//}
					//}
				}
			}
		}
		//cout << "generate new GTAG sj alignmentInfo and new CTAC sj alignmentInfo" << endl;
		if((unfixedTailInfo.GTAGsjPos).size() > 0)
		{
			int tmpNewMismatchToAddInTail = (unfixedTailInfo.GTAGsjPos_mismatch)[0];

			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			#ifdef DETECT_CIRCULAR_RNA
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
				readLength - (unfixedTailInfo.GTAGsjPos[0]).first + 1, 
				(unfixedTailInfo.GTAGsjPos[0]).second, indexInfo, 
				tmpNewMismatchToAddInTail,
				unfixedTailInfo.GTAGsjPos_mismatchPos[0], 
				unfixedTailInfo.GTAGsjPos_mismatchChar[0], 
				readLength, tmpAlignmentInfo);
			#else
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
				readLength - (unfixedTailInfo.GTAGsjPos[0]).first + 1, 
				(unfixedTailInfo.GTAGsjPos[0]).second, indexInfo, 
				tmpNewMismatchToAddInTail,
				unfixedTailInfo.GTAGsjPos_mismatchPos[0], 
				unfixedTailInfo.GTAGsjPos_mismatchChar[0], 
				readLength, tmpAlignmentInfo);
			#endif
			//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmp = 1; tmp < (unfixedTailInfo.GTAGsjPos).size(); tmp++)
			{
				int tmpNewMismatchToAddInTail = unfixedTailInfo.GTAGsjPos_mismatch[tmp];

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
					readLength - (unfixedTailInfo.GTAGsjPos[tmp]).first + 1, 
					(unfixedTailInfo.GTAGsjPos[tmp]).second, indexInfo, 
					tmpNewMismatchToAddInTail,
					unfixedTailInfo.GTAGsjPos_mismatchPos[tmp], 
					unfixedTailInfo.GTAGsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
					readLength - (unfixedTailInfo.GTAGsjPos[tmp]).first + 1, 
					(unfixedTailInfo.GTAGsjPos[tmp]).second, indexInfo, 
					tmpNewMismatchToAddInTail,
					unfixedTailInfo.GTAGsjPos_mismatchPos[tmp], 
					unfixedTailInfo.GTAGsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#endif
					//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
							
			for(int tmp = 0; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
			{
				int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[tmp];//0;

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
					readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, 
					(unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo,
					tmpNewMismatchToAddInTail,
					unfixedTailInfo.CTACsjPos_mismatchPos[tmp], 
					unfixedTailInfo.CTACsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
					readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, 
					(unfixedTailInfo.CTACsjPos[tmp]).second, indexInfo,
					tmpNewMismatchToAddInTail,
					unfixedTailInfo.CTACsjPos_mismatchPos[tmp], 
					unfixedTailInfo.CTACsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#endif
				//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}
		}
		else if((unfixedTailInfo.CTACsjPos).size() > 0)
		{
			int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[0];

			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			#ifdef DETECT_CIRCULAR_RNA
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
				readLength - (unfixedTailInfo.CTACsjPos[0]).first + 1, 
				(unfixedTailInfo.CTACsjPos[0]).second, 
				indexInfo, tmpNewMismatchToAddInTail,
				unfixedTailInfo.CTACsjPos_mismatchPos[0], 
				unfixedTailInfo.CTACsjPos_mismatchChar[0],
				readLength, tmpAlignmentInfo);
			#else
			newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
				readLength - (unfixedTailInfo.CTACsjPos[0]).first + 1, 
				(unfixedTailInfo.CTACsjPos[0]).second, 
				indexInfo, tmpNewMismatchToAddInTail,
				unfixedTailInfo.CTACsjPos_mismatchPos[0], 
				unfixedTailInfo.CTACsjPos_mismatchChar[0],
				readLength, tmpAlignmentInfo);
			#endif
			//(peAlignInfo->norAlignmentInfo_PE_1)[tmpAlignmentNO] = newTmpAlignInfo;
			//peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
			peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			for(int tmp = 1; tmp < (unfixedTailInfo.CTACsjPos).size(); tmp++)
			{
				int tmpNewMismatchToAddInTail = unfixedTailInfo.CTACsjPos_mismatch[tmp]; //0;

				Alignment_Info* newTmpAlignInfo = new Alignment_Info();
				#ifdef DETECT_CIRCULAR_RNA
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
					readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, 
					(unfixedTailInfo.CTACsjPos[tmp]).second, 
					indexInfo, tmpNewMismatchToAddInTail,
					unfixedTailInfo.CTACsjPos_mismatchPos[tmp], 
					unfixedTailInfo.CTACsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#else
				newTmpAlignInfo->newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
					readLength - (unfixedTailInfo.CTACsjPos[tmp]).first + 1, 
					(unfixedTailInfo.CTACsjPos[tmp]).second, 
					indexInfo, tmpNewMismatchToAddInTail,
					unfixedTailInfo.CTACsjPos_mismatchPos[tmp], 
					unfixedTailInfo.CTACsjPos_mismatchChar[tmp],
					readLength, tmpAlignmentInfo);
				#endif
				//(peAlignInfo->norAlignmentInfo_PE_1).push_back(newTmpAlignInfo);
				peAlignInfo->fixHeadTail_addNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType);
			}				
		}
		else
		{
		}

		// cout << "start to get unfixedTailInfo.GTAGsjPos " << endl;
		// cout << "output GTAG positions: " << endl;
		// for(int tmp = 0; tmp < unfixedTailInfo.GTAGsjPos.size(); tmp++)
		// {
		// 	cout << (unfixedTailInfo.GTAGsjPos[tmp]).first << "," << (unfixedTailInfo.GTAGsjPos[tmp]).second << endl;
		// }
		// cout << "output CTAC positions: " << endl;
		// for(int tmp = 0; tmp < unfixedTailInfo.CTACsjPos.size(); tmp++)
		// {
		// 	cout << (unfixedTailInfo.CTACsjPos[tmp]).first << "," << (unfixedTailInfo.CTACsjPos[tmp]).second << endl;
		// }		
		// cout << "end of get unfixedTailInfo.CTACsjPos ..." << endl;

		// cout //<< "****************************************************" << endl 
		// 	<< "end of fixTail_shortAnchorSJ_remapping_targetMapping_new ..." << endl
		// 	<< "****************************************************" << endl ;			
	}

	void fixTail_extend2end_new(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
		int tailSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].len;

		//string readSeq_tail = peReadInfo.returnIncompleteLongTailSeq(tmpAlignInfoType, tailSeq_len);
		int readSeq_len = peReadInfo.returnReadLength(tmpAlignInfoType);

		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);
		int midPartEndMapPosInChr = tmpAlignmentInfo->getEndMatchedPosInChr();

		if(midPartEndMapPosInChr + tailSeq_len > indexInfo->returnChromLength(midPartMapChrInt)) // endPos (if match) > chromEndPos
			return;

		vector<int> tmpMismatchPosVec;
		vector<char> tmpMismatchCharVec;
		bool matchBool;
		int extension_length = extensionForwards_errorTolerated(midPartMapChrInt,
			midPartEndMapPosInChr + 1, indexInfo, 
			peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType),
			readSeq_len - tailSeq_len + 1, readSeq_len, 
			tmpMismatchPosVec, tmpMismatchCharVec, 5, 8);
	
		if(extension_length == 0)
			matchBool = false;
		else
			matchBool = true;
		if(matchBool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2TailSoftClippingAlignment_errorTolerated_new(
				indexInfo, tailSeq_len, extension_length, tmpMismatchPosVec, tmpMismatchCharVec, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
		}
		else
		{}
	}

	void fixTail_extend2end_new_finalStepForAligner(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
		int tailSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].len;

		//string readSeq_tail = peReadInfo.returnIncompleteLongTailSeq(tmpAlignInfoType, tailSeq_len);
		int readSeq_len = peReadInfo.returnReadLength(tmpAlignInfoType);

		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);
		int midPartEndMapPosInChr = tmpAlignmentInfo->getEndMatchedPosInChr();

		if(midPartEndMapPosInChr + tailSeq_len > indexInfo->returnChromLength(midPartMapChrInt)) // endPos (if match) > chromEndPos
			return;

		vector<int> tmpMismatchPosVec;
		vector<char> tmpMismatchCharVec;
		bool matchBool;
		int extension_length = extensionForwards_errorTolerated_finalStepForAligner(midPartMapChrInt,
			midPartEndMapPosInChr + 1, indexInfo, 
			peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType),
			readSeq_len - tailSeq_len + 1, readSeq_len, 
			tmpMismatchPosVec, tmpMismatchCharVec, 5, 8);
	
		if(extension_length == 0)
			matchBool = false;
		else
			matchBool = true;
		if(matchBool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2TailSoftClippingAlignment_errorTolerated_new(
				indexInfo, tailSeq_len, extension_length, tmpMismatchPosVec, tmpMismatchCharVec, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, tmpAlignInfoType, tmpIndex_peAlignInfo);
		}
		else
		{}
	}

	void fixTail_extend2end_new_fixIndelBesideReadEnd(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo, 
		int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{
		//cout << "fixTail_extend2end_new_fixIndelBesideReadEnd starts ......" << endl;
		int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
		int tailSeq_len = (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].len;

		int readSeq_len = peReadInfo.returnReadLength(tmpAlignInfoType);

		string midPartMapChrName = tmpAlignmentInfo->returnAlignChromName();
		int midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);	
		int midPartEndMapPosInChr = tmpAlignmentInfo->getEndMatchedPosInChr();

		int tmpChromLength = indexInfo->returnChromLength(midPartMapChrInt);
		if(midPartEndMapPosInChr + tailSeq_len + 3> tmpChromLength) // endPos (if match) > chromEndPos
			return;

		if(tailSeq_len < 3)
			return;

		string last3baseSeqInRead = (peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)).substr(readSeq_len-3,3);
		int correspondingMapPos_lastBaseInRead = midPartEndMapPosInChr + tailSeq_len;
		int possibleIndelLen[6] = {1,-1,2,-2,3,-3}; // <0 means deletion; >0 means insertion, in this method, max_ins/del <= 3
		
		//cout << "last3baseSeqInRead: " << last3baseSeqInRead << endl;
		//cout << "correspondingMapPos_lastBaseInRead: " << correspondingMapPos_lastBaseInRead << endl;

		vector<int> validIndelLenVec;
		
		for(int tmp = 0; tmp < 6; tmp++)
		{
			if((3 + possibleIndelLen[tmp]) <= tailSeq_len)
			{
				//cout << "tmplastBaseMapPos: " << correspondingMapPos_lastBaseInRead + possibleIndelLen[tmp] << endl;
				// if(correspondingMapPos_lastBaseInRead - possibleIndelLen[tmp] > tmpChromLength)
				// 	continue;
				string tmpCorrespondingLast3baseSeqInChr 
					= indexInfo->returnChromStrSubstr(midPartMapChrInt, correspondingMapPos_lastBaseInRead - possibleIndelLen[tmp] - 2, 3);
				//cout << "tmp: " << tmp << endl;
				//cout << "tmpIndel: " << possibleIndelLen[tmp] << endl;
				//cout << "targetSeq: " << tmpCorrespondingLast3baseSeqInChr << endl;
				if(tmpCorrespondingLast3baseSeqInChr == last3baseSeqInRead)
				{
					//cout << "match !" << endl;
					validIndelLenVec.push_back(possibleIndelLen[tmp]);
				}
				else
				{
					//cout << "unmatch !" << endl;
				}
			}
		}
		bool fixIndel_bool = false;
		vector<Jump_Code> fixedIndelJumpCodeVec;
		//int fixedIndelChromMapPos;
		for(int tmp = 0; tmp < validIndelLenVec.size(); tmp++)
		{
			int tmpIndelLen = validIndelLenVec[tmp];	
			//cout << "tmpIndelLen: " << tmpIndelLen << endl;
			if((tmpIndelLen) < 0) // deletion
			{
				//cout << "start to try deletion ...." << endl;
				int deletionLength = 0 - tmpIndelLen;
				if(tailSeq_len-3 < 2)
				{
					int secondMatchLength = 3;
					int firstMatchLength = tailSeq_len - 3;					
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code deletionJumpCode(deletionLength, "D");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
					fixedIndelJumpCodeVec.push_back(deletionJumpCode);
					fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
					//fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
					fixIndel_bool = true;
					break;	
				}
				else
				{
					FixDoubleAnchor_Deletion_Info* delInfo = new FixDoubleAnchor_Deletion_Info();
					int tmp_toFix_deletion_read_start = readSeq_len - tailSeq_len;
					int tmp_toFix_deletion_read_end = readSeq_len -3;
					//int subSeqLengthInProcess = tmp_toFix_deletion_read_end - tmp_toFix_deletion_read_start + 1;
					int tmp_toFix_deletion_chrom_start = midPartEndMapPosInChr;
					int tmp_toFix_deletion_chrom_end = correspondingMapPos_lastBaseInRead - 3 + deletionLength;
					int tmp_max_allowed_mismatchNum = 0;//subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 1;
					bool deletion_fixed = delInfo->detectBestDeletion_lessMismatch(tmp_toFix_deletion_read_start, 
						tmp_toFix_deletion_read_end, tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, 
						(peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)), indexInfo,
						midPartMapChrInt, tmp_max_allowed_mismatchNum);
					if(deletion_fixed)
					{
						int firstMatchLength = delInfo->return_best_deletion_prefix_match_length() - 1;
						int secondMatchLength = tailSeq_len - firstMatchLength;
						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code deletionJumpCode(deletionLength, "D");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M");
						fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
						fixedIndelJumpCodeVec.push_back(deletionJumpCode);
						fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
						//fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
						fixIndel_bool = true;
						delete delInfo;
						break;
					}
					else
					{
						delete delInfo;
					}
				}
			}
			else // insertion
			{
				//cout << "start to try insertion ...." << endl;
				if(tailSeq_len-3 < tmpIndelLen)
				{
					//cout << "no need to call insInfo" << endl;
					int secondMatchLength = 3;
					int insertionLength = tmpIndelLen;
					int firstMatchLength = tailSeq_len - 3- insertionLength;
					
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
					fixedIndelJumpCodeVec.push_back(insertionJumpCode);
					fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
					//fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
					fixIndel_bool = true;
					break;					
				}
				else
				{	
					//cout << "call FixDoubleAnchor_Insertion_Info to fix insertion" << endl;
					FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info();
					int tmp_toFix_insertion_read_start = readSeq_len - tailSeq_len;
					int tmp_toFix_insertion_read_end = readSeq_len -3;
					int tmp_toFix_insertion_chrom_start = midPartEndMapPosInChr; //midPartMapPosInChr - headSeq_len + tmpIndelLen + 3;
					int tmp_toFix_insertion_chrom_end = correspondingMapPos_lastBaseInRead - tmpIndelLen - 3;
					int tmp_max_allowed_mismatchNum = 0; // no mismatch allowed, 
					bool fixInsertionBool = insInfo->detectBestInsertion_lessMismatch(
						tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
						tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, 
						(peReadInfo.returnReadSeqInDirection_alignInfoType(tmpAlignInfoType)), indexInfo,
						midPartMapChrInt, tmp_max_allowed_mismatchNum);
					int firstMatchLength = insInfo->return_best_insertion_prefix_match_length() - 1;
					int insertionLength = tmpIndelLen;
					int secondMatchLength = tailSeq_len - firstMatchLength - insertionLength;
					
					if(fixInsertionBool)
					{	
						if(secondMatchLength > 0)
						{
							Jump_Code firstMatchJumpCode(firstMatchLength, "M");
							Jump_Code insertionJumpCode(insertionLength, "I");
							Jump_Code secondMatchJumpCode(secondMatchLength, "M");
							fixedIndelJumpCodeVec.push_back(firstMatchJumpCode);
							fixedIndelJumpCodeVec.push_back(insertionJumpCode);
							fixedIndelJumpCodeVec.push_back(secondMatchJumpCode);
							//fixedIndelChromMapPos = midPartMapPosInChr - headSeq_len + tmpIndelLen;
							fixIndel_bool = true;
							delete insInfo;
							break;
						}
						else
						{
							delete insInfo;
						}
					}
					else
					{
						delete insInfo;
					}
				}
			}
		}
		
		//cout << "fixIndel_bool: " << fixIndel_bool << endl;

		if(fixIndel_bool)
		{
			Alignment_Info* newTmpAlignInfo = new Alignment_Info();
			newTmpAlignInfo->newAlignInfoAfterExtension2Tail_fixIndelBesideReadEnd(
					indexInfo, fixedIndelJumpCodeVec, tmpAlignmentInfo);
			peAlignInfo->fixHeadTail_replaceWithNewAlignInfo(newTmpAlignInfo, 
				tmpAlignInfoType, tmpIndex_peAlignInfo);				
		}
		return;
	}

	void fixTail_extend2end_nwdpBesideReadEnd(PE_Read_Info& peReadInfo,
		PE_Read_Alignment_Info* peAlignInfo, Index_Info* indexInfo,
		int tmpAlignInfoType, int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo)
	{

	}

	void fixHeadTail_areaAndStringHash_new_remappingOnly(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, 
		SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, bool Do_extendHeadTail, bool SE_or_PE_bool)
	{

		//cout << "fixHeadTail_areaAndStringHash_new_remappingOnly ..." << endl;
		int alignType_max = 4; 
		if(SE_or_PE_bool)
			alignType_max = 2;

		Alignment_Info* tmpAlignmentInfo;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}
				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				//cout << "unfixedHeadLength: " << unfixedHeadLength << endl;
				this->fixHead_shortAnchorRemappingOnly_new(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);
				tmpAlignmentInfo = NULL;
			}	
		}

		//cout << endl << "after fixing Head /////////////" << endl;
		//cout << "peAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
		//cout << "end of fixing head peAlignInfo ////////////////////" << endl << endl;
		/////////////////// fix tail //////////////////////////////
		//cout << "start to fix tail" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
				//cout << "... to fixTail in the specific alignInfo" << endl;
				this->fixTail_shortAnchorRemappingOnly_new(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
				tmpAlignmentInfo = NULL;
			}
		}
	}
	#ifdef PERSONALIZED_CHR_SEQ
	void fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq, bool SE_or_PE_bool,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, int SNPlocInSyntheticSNPseq)
	{
		//cout << "start to do fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap( " << endl;
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;

		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}

				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;

				this->fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
						spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, 
						sa_snp, lcpCompress_snp, childTab_snp, chrom_snp,
						verifyChild_snp, indexInfo_snp, SNPlocInSyntheticSNPseq
						);

				tmpAlignmentInfo = NULL;
			}	
		}
		
		/////////////////// fix tail //////////////////////////////
		//cout << "start to fix tail " << endl;		
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
						
				this->fixTail_shortAnchorGreedyMappingOnly_includeSNPseqMap(peReadInfo, peAlignInfo, SJ, 
							secondLevelChrom,
							secondLevelSa,
							secondLevelLcpCompress,
							secondLevelChildTab,
							secondLevelDetChild,
							spliceJunctionHashExists,
							indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
							annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
							spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, 
							sa_snp, lcpCompress_snp, childTab_snp, chrom_snp,
							verifyChild_snp, indexInfo_snp, SNPlocInSyntheticSNPseq);
				tmpAlignmentInfo = NULL;
			}
		}		
	}
	#endif

	void fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq, bool SE_or_PE_bool)
	{
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;

		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}

				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;

				this->fixHead_shortAnchorGreedyMappingOnly(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
						spliceJunctionDistanceMax, checkQualSeqForReadSegSeq);

				tmpAlignmentInfo = NULL;
			}	
		}
		
		//cout << "start to fix tail " << endl;
		/////////////////// fix tail //////////////////////////////
				
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
						
				this->fixTail_shortAnchorGreedyMappingOnly(peReadInfo, peAlignInfo, SJ, 
							secondLevelChrom,
							secondLevelSa,
							secondLevelLcpCompress,
							secondLevelChildTab,
							secondLevelDetChild,
							spliceJunctionHashExists,
							indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
							annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
							spliceJunctionDistanceMax, checkQualSeqForReadSegSeq);
				tmpAlignmentInfo = NULL;
			}
		}		
	}	
	
	void fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, 
		SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool checkQualSeqForShortAnchorSeqToTargetMap, bool SE_or_PE_bool)
	{
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;
		//cout << "********************************************************************" << endl;
		//cout << "start to fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping" << endl;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;
			//cout << "readLength: " << readLength << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}
				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				//cout << "unfixedHeadLength: " << unfixedHeadLength << endl;
				//if(unfixedHeadLength < SHORT_LONG_END_THRESHOLD)
				//{	
					this->fixHead_shortAnchorSJ_remapping_targetMapping_new(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo,
						checkQualSeqForShortAnchorSeqToTargetMap);	
				tmpAlignmentInfo = NULL;	
			}	
		}
		//cout << "start to fix tail " << endl;
		/////////////////// fix tail //////////////////////////////
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "alignInfoType: " << tmpAlignInfoType << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			//cout << "tmpVecSize: " << tmpVecSize << endl;
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: " << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
				//if((tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].len < SHORT_LONG_END_THRESHOLD )
				//{
					this->fixTail_shortAnchorSJ_remapping_targetMapping_new(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo,
						checkQualSeqForShortAnchorSeqToTargetMap);	
				tmpAlignmentInfo = NULL;
			}
		}
		// cout << "end of fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping" << endl;
		// cout << "********************************************************************" << endl;
	}

	void fixHeadTail_extend2end(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, bool SE_or_PE_bool)
	{
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;

				//cout << "start fixHead_extend2end ... " << endl;
				//////////////////// fix head //////////////////////////////
				//cout << "start to fix head" << endl;
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
				{
					int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
						tmpAlignmentNO++)
					{

						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
						{
							continue;
						}
						//cout << "Head need to extend2end -- tmpAlignInfoType: " << tmpAlignInfoType << " tmpAlignmentNO: " << tmpAlignmentNO << endl;
						//int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;

						this->fixHead_extend2end_new(peReadInfo, peAlignInfo,
								indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
						tmpAlignmentInfo = NULL;
					}	
				}
				/////////////////// fix tail //////////////////////////////
				//cout << "start fixTail_extend2end ... " << endl;
				for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
				{
					int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
					for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
						tmpAlignmentNO++)
					{
						tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
						int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
						if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
						{
							continue;
						}
						//cout << "Tail need to extend2end -- tmpAlignInfoType: " << tmpAlignInfoType << " tmpAlignmentNO: " << tmpAlignmentNO << endl;
						this->fixTail_extend2end_new(peReadInfo, peAlignInfo,
							indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
						tmpAlignmentInfo = NULL;
					}
				}		
	}

	void fixHeadTail_shortAnchorRemappingOnly_withAlignInfer(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, 
		SJhash_Info* SJ, 
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{

		//cout << endl << "fixHeadTail_shortAnchorRemappingOnly_withAlignInfer ..." << endl;
		int alignType_max = 4; 
		//if(SE_or_PE_bool)
		//	alignType_max = 2;

		Alignment_Info* tmpAlignmentInfo;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "tmpAlignInfoType: in head ..." << tmpAlignInfoType << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: in head ..." << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}
				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				//cout << "unfixedHeadLength: " << unfixedHeadLength << endl;
				this->fixHead_shortAnchorRemappingOnly_withAlignInfer(
						peReadInfo, peAlignInfo, SJ, 
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, 
						tmpAlignmentInfo, alignInferJunctionHashInfo);

				tmpAlignmentInfo = NULL;
			}	
		}
		/////////////////// fix tail //////////////////////////////
		//cout << "start to fix tail" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			//cout << "tmpAlignInfoType: in tail ..." << tmpAlignInfoType << endl;
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				//cout << "tmpAlignmentNO: in tail ..." << tmpAlignmentNO << endl;
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
				
				this->fixTail_shortAnchorRemappingOnly_withAlignInfer(
						peReadInfo, peAlignInfo, SJ, 
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, 
						tmpAlignmentInfo, alignInferJunctionHashInfo);	
				tmpAlignmentInfo = NULL;
			}
		}
	}

	void fixHeadTail_extend2end_finalStepForAligner(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, bool SE_or_PE_bool)
	{
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;
		//cout << "start fixHead_extend2end ... " << endl;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
				{
					continue;
				}
				//cout << "Head need to extend2end -- tmpAlignInfoType: " << tmpAlignInfoType << " tmpAlignmentNO: " << tmpAlignmentNO << endl;
				//int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				this->fixHead_extend2end_new_finalStepForAligner(peReadInfo, peAlignInfo,
					indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
				tmpAlignmentInfo = NULL;
			}	
		}
		/////////////////// fix tail //////////////////////////////
		//cout << "start fixTail_extend2end ... " << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
				{
					continue;
				}
				//cout << "Tail need to extend2end -- tmpAlignInfoType: " << tmpAlignInfoType << " tmpAlignmentNO: " << tmpAlignmentNO << endl;
				this->fixTail_extend2end_new_finalStepForAligner(peReadInfo, peAlignInfo,
					indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
				tmpAlignmentInfo = NULL;
			}
		}		
	}

	void fixHeadTail_extend2end_fixIndel(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		Index_Info* indexInfo, bool SE_or_PE_bool)
	{
		//cout << "start to do fixHeadTail_extend2end_fixIndel ......" << endl;
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;		

		Alignment_Info* tmpAlignmentInfo;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix indel beside head ......" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if(((tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S")
					|| ((tmpAlignmentInfo->cigarStringJumpCode)[0].len < 3) )
				{
					continue;
				}
				//cout << "Head need to check indel -- tmpAlignInfoType: " << tmpAlignInfoType 
				//	<< " tmpAlignmentNO: " << tmpAlignmentNO << endl;
				//int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				this->fixHead_extend2end_new_fixIndelBesideReadEnd(peReadInfo, peAlignInfo,
					indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
				tmpAlignmentInfo = NULL;
			}	
		}
		/////////////////// fix tail //////////////////////////////
		//cout << "start to fix indel beside tail ......" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( ((tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S")
					||((tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].len < 3) )
				{
					continue;
				}
				//cout << "Tail need to extend2end -- tmpAlignInfoType: " << tmpAlignInfoType 
				//	<< " tmpAlignmentNO: " << tmpAlignmentNO << endl;
				this->fixTail_extend2end_new_fixIndelBesideReadEnd(peReadInfo, peAlignInfo,
					indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo);	
				tmpAlignmentInfo = NULL;
			}
		}
	}

	int checkTwoStringMatchOrNot(const string& string_1, const string& string_2, int maxMismatch)
	{
		//bool matchBool = true;
		int mismatchNum = 0;
		for(int tmp = 0; tmp < string_1.length(); tmp++)
		{
			if(string_1[tmp] != string_2[tmp])
			{	mismatchNum++;
				if(mismatchNum > maxMismatch)
				{
					return -1; // not match
				}
			}
		}
		return mismatchNum;
	}

	#ifdef VARY_SNP_MER
	void fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap_varySNPmer(
		PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq, bool SE_or_PE_bool,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp)//, int SNPlocInSyntheticSNPseq)
	{
		//cout << "start to do fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap( " << endl;
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;
		Alignment_Info* tmpAlignmentInfo;
		//////////////////// fix head //////////////////////////////
		//cout << "start to fix head" << endl;
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				if( (tmpAlignmentInfo->cigarStringJumpCode)[0].type != "S" )
					continue;
				int unfixedHeadLength = (tmpAlignmentInfo->cigarStringJumpCode)[0].len;
				this->fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap_varySNPmer(peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
						spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, 
						sa_snp, lcpCompress_snp, childTab_snp, chrom_snp,
						verifyChild_snp, indexInfo_snp);//, SNPlocInSyntheticSNPseq);
				tmpAlignmentInfo = NULL;
			}	
		}
		
		/////////////////// fix tail //////////////////////////////
		//cout << "start to fix tail " << endl;		
		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= alignType_max; tmpAlignInfoType++)
		{
			int tmpVecSize = (peAlignInfo->fixHeadTail_getAlignInfoVecSize(tmpAlignInfoType));
			for(int tmpAlignmentNO = 0; tmpAlignmentNO < tmpVecSize; 
				tmpAlignmentNO++)
			{
				tmpAlignmentInfo = (peAlignInfo->fixHeadTail_getAlignInfo(tmpAlignInfoType, tmpAlignmentNO));
				int cigarStringJumpCodeSize = (tmpAlignmentInfo->cigarStringJumpCode).size();
				if( (tmpAlignmentInfo->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type != "S" )
					continue;		
				this->fixTail_shortAnchorGreedyMappingOnly_includeSNPseqMap_varySNPmer(peReadInfo, peAlignInfo, SJ, 
							secondLevelChrom,
							secondLevelSa,
							secondLevelLcpCompress,
							secondLevelChildTab,
							secondLevelDetChild,
							spliceJunctionHashExists,
							indexInfo, tmpAlignInfoType, tmpAlignmentNO, tmpAlignmentInfo, Do_extendHeadTail,
							annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
							spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, 
							sa_snp, lcpCompress_snp, childTab_snp, chrom_snp,
							verifyChild_snp, indexInfo_snp);//, SNPlocInSyntheticSNPseq);
				tmpAlignmentInfo = NULL;
			}
		}		
	}	
	#endif

	#ifdef VARY_SNP_MER
	void fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap_varySNPmer(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo, 
		Alignment_Info* tmpAlignmentInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp//, int SNPlocInSyntheticSNPseq
		)
	{
		//cout << "start to fixHead_shortAnchorGreedyMappingOnly_includeSNPseqMap(" << endl;
		Incomplete_Long_Head* incompleteHeadInfo = new Incomplete_Long_Head();
		incompleteHeadInfo->getIncompleteLongHeadInfoFromRecordWithAlignInfoType_new(
			tmpAlignmentInfo, indexInfo);		
		string incompleteLongHeadSeq = peReadInfo.returnIncompleteLongHeadSeq(
			tmpAlignInfoType, incompleteHeadInfo->returnUnfixedHeadLength());
		char* incompleteLongHeadSeqChar = const_cast<char*>(incompleteLongHeadSeq.c_str());
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
		int secondLevelIndexNO = incompleteHeadInfo->returnSecondLevelIndexNum() - 1;
		if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO + 1))
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
			//continue;
		}	
		//cout << "start to do incompleteMap " << endl;
		bool incompleteMapBool;
		//cout << "secondLevelIndexNO: " << secondLevelIndexNO << endl;
		//cout << "incompleteHeadInfo_unfixedHeadLength: " << incompleteHeadInfo->returnUnfixedHeadLength() << endl;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongHeadSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteHeadInfo->returnUnfixedHeadLength(), indexInfo);
		//cout << "incompleteMapBool: " << incompleteMapBool << endl; 
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
			incompleteHeadInfo->returnMapPosIntervalStart(),
			incompleteHeadInfo->returnMapPosIntervalEnd(), 
			incompleteHeadInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteHeadInfo->returnMidPartMapChrName());
		//cout << "#########################" << endl;
		//cout << "beforeUpdating segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;

		//#ifdef PERSONALIZED_CHR_SEQ
		Seg_Info segInfo_SNPseq;
		int tmpIncompleteSeqLength = incompleteLongHeadSeq.length();
		bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(incompleteLongHeadSeqChar,
			sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
		if(tmpMapWithSnpSeqIndexBool)
			segInfo->update_includeSNPseqMapSegInfo_varySNPmer(segInfo_SNPseq, indexInfo_snp, indexInfo);
		//cout << "afterUpdating with greedyMapping segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;
		segInfo->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
			indexInfo_snp, incompleteLongHeadSeq, indexInfo);
		//cout << "afterUpdating with targetMapping segInfo: " << endl << segInfo->segInfoStr(indexInfo) ;
		//#endif

		if(segInfo->returnSegmentNum() >= SEGMENTNUM)
		{
			delete seg2ndOriInfo;
			delete segInfo;
			delete incompleteHeadInfo;
			return;
		}

		segInfo->addMidPartSeg_incompleteHead(incompleteHeadInfo->returnMidPartLength(), 
			(incompleteHeadInfo->returnUnfixedHeadLength() + 1), 
			incompleteHeadInfo->returnMidPartMapChrInt(), incompleteHeadInfo->returnMidPartMapPosInChr(), indexInfo);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}
		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg_incompleteHead(segInfo, spliceJunctionDistanceMax); // Fix me: why not replace it by getPossiPathFromSeg_incompleteHead  
	
		pathInfo->filterPath_incompleteHead(segInfo); // filter out those path without the midPartSeg
	
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete seg2ndOriInfo;
			delete incompleteHeadInfo;
			return;
		}

		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongHeadSeq(
				tmpAlignInfoType, incompleteHeadInfo->returnIncompleteHeadAndMidPartLength()), 
			incompleteHeadInfo->returnIncompleteHeadAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		//cout << "start to get peAlignInfo" << endl;
		peAlignInfo->incompleteHead_replaceAndPushBackPathInfo2PeAlignInfo_new(pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteHeadInfo->returnMidPartLength(), indexInfo, tmpAlignmentInfo);

		//cout << "peAlignInfo: " << endl << peAlignInfo->getTmpAlignInfo(
		//	(peReadInfo.readInfo_pe1).readName, (peReadInfo.readInfo_pe2).readName, 
		//	(peReadInfo.readInfo_pe1).readSeq, (peReadInfo.readInfo_pe2).readSeq, "*", "*");

		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete seg2ndOriInfo;
		delete incompleteHeadInfo;
		return;
	}
	#endif

	#ifdef VARY_SNP_MER
	void fixTail_shortAnchorGreedyMappingOnly_includeSNPseqMap_varySNPmer(PE_Read_Info& peReadInfo, 
		PE_Read_Alignment_Info* peAlignInfo, SJhash_Info* SJ, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		bool spliceJunctionHashExists,
		Index_Info* indexInfo, int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* tmpAlignmentInfo, 
		bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
		int spliceJunctionDistanceMax, bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp//, int SNPlocInSyntheticSNPseq	
		)
	{
		Incomplete_Long_Tail* incompleteTailInfo = new Incomplete_Long_Tail();

		incompleteTailInfo->getIncompleteLongTailInfoFromRecordWithAlignInfoType_new(
			peReadInfo, tmpAlignInfoType,
			tmpAlignmentInfo, indexInfo);

		string incompleteLongTailSeq = peReadInfo.returnIncompleteLongTailSeq(
			tmpAlignInfoType, incompleteTailInfo->returnUnfixedTailLength());
	
		char* incompleteLongTailSeqChar = const_cast<char*>(incompleteLongTailSeq.c_str());

		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

		int secondLevelIndexNO = incompleteTailInfo->returnSecondLevelIndexNum() - 1;

		if(
			indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1)
			)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
			//continue;
		}		
		bool incompleteMapBool;
		incompleteMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									incompleteLongTailSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									incompleteTailInfo->returnUnfixedTailLength(), indexInfo);
		if(!incompleteMapBool)
		{
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		Seg_Info* segInfo_old = new Seg_Info(seg2ndOriInfo, 
			incompleteTailInfo->returnMapPosIntervalStart(),
			incompleteTailInfo->returnMapPosIntervalEnd(), 
			incompleteTailInfo->returnChrPosStartIn2ndLevelIndex(),
			indexInfo, incompleteTailInfo->returnMidPartMapChrName());
		
		//#ifdef PERSONALIZED_CHR_SEQ
		Seg_Info segInfo_SNPseq;
		int tmpIncompleteSeqLength = incompleteLongTailSeq.length();
		bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(incompleteLongTailSeqChar, 
			sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
		//if(tmpMapWithSnpSeqIndexBool)
		//	segInfo_old->update_includeSNPseqMapSegInfo(segInfo_SNPseq, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
		
		//segInfo_old->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
		//	indexInfo_snp, incompleteLongTailSeq, indexInfo, SNPlocInSyntheticSNPseq);

		if(tmpMapWithSnpSeqIndexBool)
			segInfo_old->update_includeSNPseqMapSegInfo_varySNPmer(segInfo_SNPseq, indexInfo_snp, indexInfo);
		
		segInfo_old->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, 
			chrom_snp, verifyChild_snp, indexInfo_snp, incompleteLongTailSeq, indexInfo);		
		//#endif

		if(segInfo_old->returnSegmentNum() >= SEGMENTNUM)
		{
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}	
	
		Seg_Info* segInfo = new Seg_Info();

		segInfo->addMidPartSeg_incompleteTail(incompleteTailInfo->returnMidPartLength(), 
			incompleteTailInfo->returnMidPartLocInRead(), 
			incompleteTailInfo->returnMidPartMapChrInt(), 
			incompleteTailInfo->returnMidPartMapPosInChr(), indexInfo, segInfo_old);

		// fix me: test: checkQualSeqForReadSegSeq
		if(checkQualSeqForReadSegSeq)
		{
			segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
		}

		//cout << "segInfo: " << endl << segInfo->segInfoStr(indexInfo);

		segInfo->assignLongSegMinLength(CONFIDENT_SEG_LENGTH_FIX_LONG_END);
	
		Path_Info* pathInfo = new Path_Info();
	
		pathInfo->getPossiPathFromSeg_incompleteTail(segInfo, spliceJunctionDistanceMax);  

		pathInfo->filterPath_incompleteTail();//segInfo); // filter out those path without the midPartSeg

		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete pathInfo;
			delete segInfo;
			delete segInfo_old;
			delete seg2ndOriInfo;
			delete incompleteTailInfo;
			return;
		}
		//cout << "start to fix gaps" << endl;
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath_phase2(pathInfo, segInfo, indexInfo,
			peReadInfo.returnIncompleteLongTailSeq(
				tmpAlignInfoType, incompleteTailInfo->returnIncompleteTailAndMidPartLength()), 
			incompleteTailInfo->returnIncompleteTailAndMidPartLength(), Do_extendHeadTail,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		peAlignInfo->incompleteTail_replaceAndPushBackPathInfo2PeAlignInfo_new(
			pathInfo, tmpAlignInfoType, 
			tmpIndex_peAlignInfo, incompleteTailInfo->returnMidPartLength(), 
			indexInfo, tmpAlignmentInfo, incompleteTailInfo->returnMidPartLocInRead());

		// stop4
		delete gapInfo;
		pathInfo->memoryFree();
		delete pathInfo;
		delete segInfo;
		delete segInfo_old;
		delete seg2ndOriInfo;
		delete incompleteTailInfo;
		return;
	}
	#endif	
};

#endif