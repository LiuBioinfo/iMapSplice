// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PEALIGN_INFO_H
#define PEALIGN_INFO_H

#include <set>
#include "transcript_count_vec.h"
#include "transcript2geneMap_info.h"

using namespace std;

class PE_Read_Alignment_Info
{
//public:
private:
	double highestPairAlignmentScore;// = 0;

	int repeatRegion_index_Nor1;
	int repeatRegion_index_Rcm1;
	int repeatRegion_index_Nor2;
	int repeatRegion_index_Rcm2;

	vector<Alignment_Info*> norAlignmentInfo_PE_1;
	vector<Alignment_Info*> rcmAlignmentInfo_PE_1;
	vector<Alignment_Info*> norAlignmentInfo_PE_2;
	vector<Alignment_Info*> rcmAlignmentInfo_PE_2;

	// ****************** needed by PE *********************** // 
	vector< pair< int, vector<int> > > oriAlignPair_Nor1Rcm2;
	vector< pair< int, vector<int> > > oriAlignPair_Nor2Rcm1;

	vector< pair<int, int> > oriAlignPair_Nor1Rcm2_new;
	vector< string > oriAlignPairNew_chrNameVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_startMapPosVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_endMapPosVec_Nor1Rcm2;
	vector< int > oriAlignPairNew_mappedLength_Nor1Rcm2;
	vector< int > oriAlignPairNew_mismatchNum_Nor1Rcm2;
	vector< int > oriAlignPairNew_pairDistance_Nor1Rcm2;
	vector< int > oriAlignPairNew_insertionLength_Nor1Rcm2;
	vector< int > oriAlignPairNew_deletionLength_Nor1Rcm2;
	vector< int > oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2;
	vector< double > oriAlignPairNew_score_Nor1Rcm2;

	vector< pair<int, int> > oriAlignPair_Nor2Rcm1_new; 
	vector< string > oriAlignPairNew_chrNameVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_startMapPosVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_endMapPosVec_Nor2Rcm1;
	vector< int > oriAlignPairNew_mappedLength_Nor2Rcm1;
	vector< int > oriAlignPairNew_mismatchNum_Nor2Rcm1;
	vector< int > oriAlignPairNew_pairDistance_Nor2Rcm1;
	vector< int > oriAlignPairNew_insertionLength_Nor2Rcm1;
	vector< int > oriAlignPairNew_deletionLength_Nor2Rcm1;
	vector< int > oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1;
	vector< double > oriAlignPairNew_score_Nor2Rcm1;

	vector< vector< int > > oriAlignPairGroupedByRegion_Nor1Rcm2;
	vector< vector< int > > oriAlignPairGroupedByRegion_Nor2Rcm1;

	vector< pair< int, int > > finalAlignPair_Nor1Rcm2;
	vector< pair< int, int > > finalAlignPair_Nor2Rcm1;

	/************************************************************/
	vector<int> unpairedSEalignVec_final_Nor_end1;
	vector<int> unpairedSEalignVec_final_Rcm_end1;
	vector<int> unpairedSEalignVec_final_Nor_end2;
	vector<int> unpairedSEalignVec_final_Rcm_end2;	
	/************************************************************/
	//vector<bool> otherEndUnmappedBoolVec;
	// ******************************************************* //
	// ****************** needed by SE *********************** //  
	vector< string > seAlign_chrNameVec_Nor;
	vector< int > seAlign_startMapPosVec_Nor;
	vector< int > seAlign_endMapPosVec_Nor;
	vector< int > seAlign_mappedLengthVec_Nor;
	vector< int > seAlign_mismatchNumVec_Nor;
	vector< int > seAlign_insertionLengthVec_Nor;
	vector< int > seAlign_deletionLengthVec_Nor;	
	vector< int > seAlign_SJconfidenceLevel_Nor;
	vector< double > seAlign_scoreVec_Nor;

	vector< string > seAlign_chrNameVec_Rcm;
	vector< int > seAlign_startMapPosVec_Rcm;
	vector< int > seAlign_endMapPosVec_Rcm;
	vector< int > seAlign_mappedLengthVec_Rcm;
	vector< int > seAlign_mismatchNumVec_Rcm;
	vector< int > seAlign_insertionLengthVec_Rcm;
	vector< int > seAlign_deletionLengthVec_Rcm;	
	vector< int > seAlign_SJconfidenceLevel_Rcm;
	vector< double > seAlign_scoreVec_Rcm;

	vector< int > seAlignVec_final_Nor;
	vector< int > seAlignVec_final_Rcm;
 	// ******************************************************* //	
public:
	void cout_finalAlignPair()
	{
		cout << "## finalAlignPair_Nor1Rcm2: ##" << endl;
		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp++)
		{
			int tmpNorIndex = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcmIndex = finalAlignPair_Nor1Rcm2[tmp].second;
			cout << tmpNorIndex << "," << tmpRcmIndex << ";\t";
		}
		cout << endl << "## finalAlignPair_Nor2Rcm1: ##" << endl;
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp++)
		{
			int tmpNorIndex = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcmIndex = finalAlignPair_Nor2Rcm1[tmp].second;
			cout << tmpNorIndex << "," << tmpRcmIndex << ";\t";
		}
		cout << endl;
	}

	void cout_seAlignVec_final()
	{
		//cout << endl << "" << endl;
		cout << "## seAlignVec_final_Nor: ##" << endl;
		for(int tmp = 0; tmp < seAlignVec_final_Nor.size(); tmp++)
			cout << seAlignVec_final_Nor[tmp] + 1 << ";\t";
		cout << endl;
		cout << "## seAlignVec_final_Rcm: ##" << endl;
		for(int tmp = 0; tmp < seAlignVec_final_Rcm.size(); tmp++)
			cout << seAlignVec_final_Rcm[tmp] + 1 << ";\t";
		cout << endl;
	}

	int getMaxOverlapLengthInAlignmentPair(int readLength_1, int readLength_2)
	{
		//cout << "start to do getMaxOverlapLengthInAlignmentPair...." << endl;
		vector< vector<bool> > coveredBaseVecVec_read1;
		vector< vector<bool> > coveredBaseVecVec_read2;
		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp ++)
		{
			//cout << "tmp in finalAlignPair_Nor1Rcm2: " << tmp << endl;
			vector<bool> tmpCoveredBaseVec_read1;
			vector<bool> tmpCoveredBaseVec_read2;
			int tmpNor1_index = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcm2_index = finalAlignPair_Nor1Rcm2[tmp].second;
			//cout << "tmpNor1_index: " << tmpNor1_index << endl;
			//cout << "tmpRcm2_index: " << tmpRcm2_index << endl;
			norAlignmentInfo_PE_1[tmpNor1_index]->generateCoveredBaseVec(readLength_1,
				tmpCoveredBaseVec_read1, true);
			rcmAlignmentInfo_PE_2[tmpRcm2_index]->generateCoveredBaseVec(readLength_2,
				tmpCoveredBaseVec_read2, false);
			coveredBaseVecVec_read1.push_back(tmpCoveredBaseVec_read1);
			coveredBaseVecVec_read2.push_back(tmpCoveredBaseVec_read2);
		}
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp ++)
		{
			//cout << "tmp in finalAlignPair_Nor1Rcm2: " << tmp << endl;
			vector<bool> tmpCoveredBaseVec_read1;
			vector<bool> tmpCoveredBaseVec_read2;
			int tmpNor2_index = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcm1_index = finalAlignPair_Nor2Rcm1[tmp].second;
			//cout << "tmpNor2_index: " << tmpNor2_index << endl;
			//cout << "tmpRcm1_index: " << tmpRcm1_index << endl;
			norAlignmentInfo_PE_2[tmpNor2_index]->generateCoveredBaseVec(readLength_2,
				tmpCoveredBaseVec_read2, true);
			rcmAlignmentInfo_PE_1[tmpRcm1_index]->generateCoveredBaseVec(readLength_1,
				tmpCoveredBaseVec_read1, false);
			coveredBaseVecVec_read1.push_back(tmpCoveredBaseVec_read1);
			coveredBaseVecVec_read2.push_back(tmpCoveredBaseVec_read2);
		}

		int tmpMaxNum_coveredBase = 0;
		int coveredBaseVecVec_read1_size = coveredBaseVecVec_read1.size();
		for(int tmp = 0; tmp < coveredBaseVecVec_read1_size; tmp++)
		{
			for(int tmp2 = tmp+1; tmp2 < coveredBaseVecVec_read1_size; tmp2++)
			{
				//cout << "tmp: " << tmp << endl;
				//cout << "tmp2: " << tmp2 << endl;
				int tmpCoveredBaseNum = this->getOverlappedBaseNum(
					coveredBaseVecVec_read1[tmp], coveredBaseVecVec_read2[tmp],
					coveredBaseVecVec_read1[tmp2], coveredBaseVecVec_read2[tmp2]);
				//cout << "tmpCoveredBaseNum " << tmpCoveredBaseNum << endl;
				if(tmpCoveredBaseNum > tmpMaxNum_coveredBase)
					tmpMaxNum_coveredBase = tmpCoveredBaseNum;
			}
		}
		return tmpMaxNum_coveredBase;
	}

	int getOverlappedBaseNum(vector<bool>& tmpCoveredBase_read1_1, vector<bool>& tmpCoveredBase_read2_1,
		vector<bool>& tmpCoveredBase_read1_2, vector<bool>& tmpCoveredBase_read2_2)
	{
		int tmpOverlapBaseNum = 0;
		for(int tmp = 0; tmp < tmpCoveredBase_read1_1.size(); tmp++)
		{
			bool coveredOrNot_read1_1_bool = tmpCoveredBase_read1_1[tmp];
			bool coveredOrNot_read1_2_bool = tmpCoveredBase_read1_2[tmp];
			if(coveredOrNot_read1_1_bool && coveredOrNot_read1_2_bool)
				tmpOverlapBaseNum ++;
		}
		for(int tmp = 0; tmp < tmpCoveredBase_read2_1.size(); tmp++)
		{
			bool coveredOrNot_read2_1_bool = tmpCoveredBase_read2_1[tmp];
			bool coveredOrNot_read2_2_bool = tmpCoveredBase_read2_2[tmp];
			if(coveredOrNot_read2_1_bool && coveredOrNot_read2_2_bool)
				tmpOverlapBaseNum ++;
		}
		return tmpOverlapBaseNum;
	}

	bool partialAlignmentExists()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			bool tmpAlignmentIncompleteBool = norAlignmentInfo_PE_1[tmp]->alignmentIncompleteBool();
			if(tmpAlignmentIncompleteBool)
				return true;
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			bool tmpAlignmentIncompleteBool = rcmAlignmentInfo_PE_1[tmp]->alignmentIncompleteBool();
			if(tmpAlignmentIncompleteBool)
				return true;			
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			bool tmpAlignmentIncompleteBool = norAlignmentInfo_PE_2[tmp]->alignmentIncompleteBool();
			if(tmpAlignmentIncompleteBool)
				return true;			
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			bool tmpAlignmentIncompleteBool = rcmAlignmentInfo_PE_2[tmp]->alignmentIncompleteBool();
			if(tmpAlignmentIncompleteBool)
				return true;			
		}
		return false;
	}

	bool returnOnlySeAlign_Nor_or_Rcm()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		return Nor_or_Rcm_bool;
	}

	string returnOnlySeAlign_chrNameStr()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromName();
		}
	}

	int returnOnlySeAlign_startPos()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromPos();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromPos();
		}
	}

	int returnOnlySeAlign_endPos()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnEndMatchedPosInChr();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnEndMatchedPosInChr();
		}		
	}

	bool mappedForFusionDetection_unfixedHead_uniqueCompleteMapped_SE_bool(
		int min_fusion_distance, int chrNameIntInOriSAM, int breakPointInOriSAM, Index_Info* indexInfo)
	{
		if(seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size() != 1)
			return false;
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool) // case 2,5
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			bool tmpAlignIncompleteBool 
				= norAlignmentInfo_PE_1[tmpNor_index]->alignmentIncompleteBool();
			if(tmpAlignIncompleteBool)
			{
				bool tmpUnfixedTailExistBool = norAlignmentInfo_PE_1[tmpNor_index]->unfixedTailExistsBool();
				if(tmpUnfixedTailExistBool)
				{
					int tmpUnfixedTailLength = norAlignmentInfo_PE_1[tmpNor_index]->unfixedTailLength();
					if(tmpUnfixedTailLength >= GlobalMapForFusionDetectionUnfixedHeadTailBuffer)
						return false;
					else
					{
						string tmpChrNameStrInUnfixedHead
							= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
						int tmpChrNameIntInUnfixedHead = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedHead);
						if(tmpChrNameIntInUnfixedHead != chrNameIntInOriSAM)
							return true;
						int tmpBreakPointInUnfixedHead 
							= norAlignmentInfo_PE_1[tmpNor_index]->returnEndMatchedPosInChr();
						//int tmpFusionDistance = tmpBreakPointInUnfixedHead - breakPointInOriSAM;
						int tmpFusionDistance = breakPointInOriSAM - tmpBreakPointInUnfixedHead - 1;
						//if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
						if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
							return true;
						else
							return false;						
					}
				}
				else
				{
					string tmpChrNameStrInUnfixedHead
						= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
					int tmpChrNameIntInUnfixedHead = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedHead);
					if(tmpChrNameIntInUnfixedHead != chrNameIntInOriSAM)
						return true;
					int tmpBreakPointInUnfixedHead 
						= norAlignmentInfo_PE_1[tmpNor_index]->returnEndMatchedPosInChr();
					//int tmpFusionDistance = tmpBreakPointInUnfixedHead - breakPointInOriSAM;
					int tmpFusionDistance = breakPointInOriSAM - tmpBreakPointInUnfixedHead - 1;
					//if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
						return true;
					else
						return false;					
				}
			}
			else
			{	
				string tmpChrNameStrInUnfixedHead
					= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
				int tmpChrNameIntInUnfixedHead = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedHead);
				if(tmpChrNameIntInUnfixedHead != chrNameIntInOriSAM)
					return true;
				int tmpBreakPointInUnfixedHead 
					= norAlignmentInfo_PE_1[tmpNor_index]->returnEndMatchedPosInChr();
				//int tmpFusionDistance = tmpBreakPointInUnfixedHead - breakPointInOriSAM;
				int tmpFusionDistance = breakPointInOriSAM - tmpBreakPointInUnfixedHead - 1;
				//if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
				if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
					return true;
				else
					return false;
			}
		}
		else // rev -- unfixed head;  for -- fixed tail --end_1, case 10,11
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			bool tmpAlignIncompleteBool 
				= rcmAlignmentInfo_PE_1[tmpRcm_index]->alignmentIncompleteBool();
			if(tmpAlignIncompleteBool)
			{
				bool tmpUnfixedHeadExistBool = rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedHeadExistsBool();
				if(tmpUnfixedHeadExistBool)
				{
					int tmpUnfixedHeadLength = rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedHeadLength();
					if(tmpUnfixedHeadLength >= GlobalMapForFusionDetectionUnfixedHeadTailBuffer)
						return false;
					else
						return true;
				}
				else
					return true;					
			}
			else
			{
				return true;
				//string tmpChrNameStrInUnfixedHead
				//	= rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromName();
				//int tmpChrNameIntInUnfixedHead = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedHead);
				//if(tmpChrNameIntInUnfixedHead != chrNameIntInOriSAM)
				//	return true;
				//int tmpBreakPointInUnfixedHead 
				//	= rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromPos();
				//int tmpFusionDistance = tmpBreakPointInUnfixedHead - breakPointInOriSAM;				
				// if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
				// 	return true;
				// else
				// 	return false;
			}
		}		
	}

	int return_onlySeAlignInfo_firstMatchLen()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnFirstMatchLen();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnFirstMatchLen();
		}		
	}

	int return_onlySeAlignInfo_lastMatchLen()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnLastMatchLen();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnLastMatchLen();
		}		
	}

	int return_onlySeAlignInfo_headSoftClipLen()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnHeadSoftClipLen();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnHeadSoftClipLen();
		}		
	}

	int return_onlySeAlignInfo_tailSoftClipLen()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnTailSoftClipLen();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnTailSoftClipLen();
		}		
	}

	int return_onlySeAlignInfo_chrNameInt(Index_Info* indexInfo)
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return indexInfo->convertStringToInt(norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName());
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return indexInfo->convertStringToInt(rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromName());
		}			
	}

	int return_onlySeAlignInfo_startMatchedPosInChr()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromPos();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromPos();
		}			
	}	

	int return_onlySeAlignInfo_endMatchedPosInChr()
	{
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			return norAlignmentInfo_PE_1[tmpNor_index]->returnEndMatchedPosInChr();
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			return rcmAlignmentInfo_PE_1[tmpRcm_index]->returnEndMatchedPosInChr();
		}			
	}

	bool mappedForFusionDetection_unfixedTail_uniqueCompleteMapped_SE_bool(
		int min_fusion_distance, int chrNameIntInOriSAM, int breakPointInOriSAM, Index_Info* indexInfo)
	{
		if(seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size() != 1)
			return false;
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool) // case -- 1,4
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			bool tmpAlignIncompleteBool 
				= norAlignmentInfo_PE_1[tmpNor_index]->alignmentIncompleteBool();
			if(tmpAlignIncompleteBool)
			{
				bool tmpUnfixedHeadExistBool = norAlignmentInfo_PE_1[tmpNor_index]->unfixedHeadExistsBool();
				if(tmpUnfixedHeadExistBool)
				{
					int tmpUnfixedHeadLength = norAlignmentInfo_PE_1[tmpNor_index]->unfixedHeadLength();
					if(tmpUnfixedHeadLength >= GlobalMapForFusionDetectionUnfixedHeadTailBuffer)
						return false;
					else
					{
						string tmpChrNameStrInUnfixedTail
							= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
						int tmpChrNameIntInUnfixedTail = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedTail);
						if(tmpChrNameIntInUnfixedTail != chrNameIntInOriSAM)
							return true;
						int tmpBreakPointInUnfixedTail 
							= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromPos();
						int tmpFusionDistance = tmpBreakPointInUnfixedTail - breakPointInOriSAM - 1;
						if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
							return true;
						else
							return false;							
					}
				}
				else
				{
					string tmpChrNameStrInUnfixedTail
						= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
					int tmpChrNameIntInUnfixedTail = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedTail);
					if(tmpChrNameIntInUnfixedTail != chrNameIntInOriSAM)
						return true;
					int tmpBreakPointInUnfixedTail 
						= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromPos();
					int tmpFusionDistance = tmpBreakPointInUnfixedTail - breakPointInOriSAM - 1;
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
						return true;
					else
						return false;					
				}
			}
			else
			{	
				string tmpChrNameStrInUnfixedTail
					= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromName();
				int tmpChrNameIntInUnfixedTail = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedTail);
				if(tmpChrNameIntInUnfixedTail != chrNameIntInOriSAM)
					return true;
				int tmpBreakPointInUnfixedTail 
					= norAlignmentInfo_PE_1[tmpNor_index]->returnAlignChromPos();
				int tmpFusionDistance = tmpBreakPointInUnfixedTail - breakPointInOriSAM - 1;
				if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance < 0))
					return true;
				else
					return false;
			}
		}
		else // case 7,8
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			bool tmpAlignIncompleteBool 
				= rcmAlignmentInfo_PE_1[tmpRcm_index]->alignmentIncompleteBool();
			if(tmpAlignIncompleteBool)
			{
				bool tmpUnfixedTailExistBool = rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedTailExistsBool();
				if(tmpUnfixedTailExistBool)
				{
					int tmpUnfixedTailLength = rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedTailLength();
					if(tmpUnfixedTailLength >= GlobalMapForFusionDetectionUnfixedHeadTailBuffer)
						return false;
					else
						return true;
				}
				else
					return true;		
			}
			else
			{
				return true;
				// string tmpChrNameStrInUnfixedTail
				// 	= rcmAlignmentInfo_PE_1[tmpRcm_index]->returnAlignChromName();
				// int tmpChrNameIntInUnfixedTail = indexInfo->convertStringToInt(tmpChrNameStrInUnfixedTail);
				// if(tmpChrNameIntInUnfixedTail != chrNameIntInOriSAM)
				// 	return true;
				// int tmpBreakPointInUnfixedTail 
				// 	= rcmAlignmentInfo_PE_1[tmpRcm_index]->returnEndMatchedPosInChr();
				// int tmpFusionDistance = tmpBreakPointInUnfixedTail - breakPointInOriSAM;
				// if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
				// 	return true;
				// else
				// 	return false;
			}
		}		
	}

	/*
	bool mappedForFusionDetection_outerSoftClipUniquePaired_SE_bool(
		bool headOrTail_inOriSAM_bool, int min_fusion_distance,
		int startPos_or_endPos_inOriSAM)
	{
		if(seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size() != 1)
			return false;
		bool Nor_or_Rcm_bool = (seAlignVec_final_Nor.size() == 1);
		if(Nor_or_Rcm_bool)
		{
			int tmpNor_index = seAlignVec_final_Nor[0];
			if(headOrTail_inOriSAM_bool)
			{
				bool unfixedTailExistsBool_tmpOuterSoftClipSeq 
					= norAlignmentInfo_PE_1[tmpNor_index]->unfixedTailExistsBool();
				if(unfixedTailExistsBool_tmpOuterSoftClipSeq)
					return false;
				else
				{
					int startPosInOriSAM = startPos_or_endPos_inOriSAM;
					int endPosInMappedSoftClipHead = this->returnOnlySeAlign_endPos();
					int tmpFusionDistance = startPosInOriSAM - endPosInMappedSoftClipHead - 1;
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
						return true;
					else
						return false;
				}
			}	
			else
			{
				bool unfixedHeadExistsBool_tmpOuterSoftClipSeq 
					= norAlignmentInfo_PE_1[tmpNor_index]->unfixedHeadExistsBool();
				if(unfixedHeadExistsBool_tmpOuterSoftClipSeq)
					return false;
				else
				{
					int endPosInOriSAM = startPos_or_endPos_inOriSAM;
					int startPosInMappedSoftClipHead = this->returnOnlySeAlign_startPos();
					int tmpFusionDistance = startPosInMappedSoftClipHead - endPosInOriSAM - 1;
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
						return true;
					else
						return false;	
				}
			}
		}
		else
		{
			int tmpRcm_index = seAlignVec_final_Rcm[0];
			if(headOrTail_inOriSAM_bool)
			{
				bool unfixedTailExistsBool_tmpOuterSoftClipSeq 
					= rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedTailExistsBool();
				if(unfixedTailExistsBool_tmpOuterSoftClipSeq)
					return false;
				else
				{	
					int startPosInOriSAM = startPos_or_endPos_inOriSAM;
					int endPosInMappedSoftClipHead = this->returnOnlySeAlign_endPos();
					int tmpFusionDistance = startPosInOriSAM - endPosInMappedSoftClipHead - 1;
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
						return true;
					else
						return false;
				}
			}	
			else
			{
				bool unfixedHeadExistsBool_tmpOuterSoftClipSeq 
					= rcmAlignmentInfo_PE_1[tmpRcm_index]->unfixedHeadExistsBool();
				if(unfixedHeadExistsBool_tmpOuterSoftClipSeq)
					return false;
				else
				{
					int endPosInOriSAM = startPos_or_endPos_inOriSAM;
					int startPosInMappedSoftClipHead = this->returnOnlySeAlign_startPos();
					int tmpFusionDistance = startPosInMappedSoftClipHead - endPosInOriSAM - 1;
					if((tmpFusionDistance >= min_fusion_distance)||(tmpFusionDistance <= 0 - min_fusion_distance))
						return true;
					else
						return false;	
				}	
			}
		}
	}*/

	string returnSAMstr_mappedForFusion_headSoftClipUniquePaired_SE(
		PE_Read_Info& seReadInfo, bool FastaOrFastq)
	{
		string tmpStr = getSAMformatForFinalAlignment_SE(seReadInfo, FastaOrFastq);
		return tmpStr;
	}

	string returnSAMstr_mappedForFusion_tailSoftClipUniquePaired_SE(
		PE_Read_Info& seReadInfo, bool FastaOrFastq)
	{
		string tmpStr = getSAMformatForFinalAlignment_SE(seReadInfo, FastaOrFastq);
		return tmpStr;
	}

	string returnSAMstr_mappedForFusion_outerSoftClipUniquePaired_SE(
		PE_Read_Info& seReadInfo, bool FastaOrFastq)
	{
		string tmpStr = getSAMformatForFinalAlignment_SE(seReadInfo, FastaOrFastq);
		return tmpStr;
	}

	int return_seAlignVec_final_Nor_size()
	{
		return seAlignVec_final_Nor.size();
	}

	int return_seAlignVec_final_Rcm_size()
	{
		return seAlignVec_final_Rcm.size();
	}

	int return_seAlignVec_final_size()
	{
		return seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size();
	}

	string returnSeAlignVec_final()
	{
		string tmpStr = "seAlignVec_final_Nor: ";
		for(int tmp = 0; tmp < seAlignVec_final_Nor.size(); tmp++)
		{
			tmpStr += "\t";
			tmpStr += int_to_str(seAlignVec_final_Nor[tmp]);
		}
		tmpStr += "\n";
		tmpStr += "seAlignVec_final_Rcm: ";
		for(int tmp = 0; tmp < seAlignVec_final_Rcm.size(); tmp++)
		{
			tmpStr += "\t";
			tmpStr += int_to_str(seAlignVec_final_Rcm[tmp]);
		}
		return tmpStr;
	}

	void pushBackNewAlignInfo_Nor1(Alignment_Info* newAlignInfo)
	{
		norAlignmentInfo_PE_1.push_back(newAlignInfo);
	}

	void pushBackNewAlignInfo_Rcm1(Alignment_Info* newAlignInfo)
	{
		rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
	}

	void pushBackNewAlignInfo_Nor2(Alignment_Info* newAlignInfo)
	{
		norAlignmentInfo_PE_2.push_back(newAlignInfo);
	}

	void pushBackNewAlignInfo_Rcm2(Alignment_Info* newAlignInfo)
	{
		rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
	}	

	void pushBackNewFinalPairAlignInfoIndex_Nor1Rcm2(int a, int b)
	{
		finalAlignPair_Nor1Rcm2.push_back(pair<int,int>(a,b));
	}

	void pushBackNewFinalPairAlignInfoIndex_Nor2Rcm1(int a, int b)
	{
		finalAlignPair_Nor2Rcm1.push_back(pair<int,int>(a,b));
	}

	double returnHighestPairAlignmentScore()
	{
		return highestPairAlignmentScore;
	}

	int returnFinalAlignPairSize_Nor1Rcm2()
	{
		return finalAlignPair_Nor1Rcm2.size();
	}

	int returnFinalAlignPairSize_Nor2Rcm1()
	{
		return finalAlignPair_Nor2Rcm1.size();
	}

	void insertNor1Rcm2FinalPairAlignInfo2targetNor1Rcm2FinalPairAlignInfo(
		PE_Read_Alignment_Info&	targetPeReadAlignInfo)
	{
		// update highest score
		double highestScore_targetPeReadAlignInfo = targetPeReadAlignInfo.returnHighestPairAlignmentScore();
		if(highestPairAlignmentScore > highestScore_targetPeReadAlignInfo)
			targetPeReadAlignInfo.assignHighestPairAlignmentScore(highestPairAlignmentScore);
		// insert finalPairAlignInfo
		int currentSize_targetNor1Rcm2FinalPairAlignInfoVec 
			= targetPeReadAlignInfo.returnFinalAlignPairSize_Nor1Rcm2();
		int size_finalPairAlignInfoVec_Nor1Rcm2
			= finalAlignPair_Nor1Rcm2.size();
		for(int tmp = 0; tmp < size_finalPairAlignInfoVec_Nor1Rcm2; tmp++)
		{
			int alignInfoIndex_Nor1 = finalAlignPair_Nor1Rcm2[tmp].first;
			int alignInfoIndex_Rcm2 = finalAlignPair_Nor1Rcm2[tmp].second;
			targetPeReadAlignInfo.pushBackNewAlignInfo_Nor1(norAlignmentInfo_PE_1[alignInfoIndex_Nor1]);
			targetPeReadAlignInfo.pushBackNewAlignInfo_Rcm2(rcmAlignmentInfo_PE_2[alignInfoIndex_Rcm2]);
			targetPeReadAlignInfo.pushBackNewFinalPairAlignInfoIndex_Nor1Rcm2(
				currentSize_targetNor1Rcm2FinalPairAlignInfoVec, currentSize_targetNor1Rcm2FinalPairAlignInfoVec);
			currentSize_targetNor1Rcm2FinalPairAlignInfoVec ++;
		}
	}

	void insertNor1Rcm2FinalPairAlignInfo2targetNor2Rcm1FinalPairAlignInfo(
		PE_Read_Alignment_Info&	targetPeReadAlignInfo)
	{
		// update highest score
		double highestScore_targetPeReadAlignInfo = targetPeReadAlignInfo.returnHighestPairAlignmentScore();
		if(highestPairAlignmentScore > highestScore_targetPeReadAlignInfo)
			targetPeReadAlignInfo.assignHighestPairAlignmentScore(highestPairAlignmentScore);
		// insert finalPairAlignInfo
		int currentSize_targetNor2Rcm1FinalPairAlignInfoVec 
			= targetPeReadAlignInfo.returnFinalAlignPairSize_Nor2Rcm1();
		int size_finalPairAlignInfoVec_Nor1Rcm2
			= finalAlignPair_Nor1Rcm2.size();
		for(int tmp = 0; tmp < size_finalPairAlignInfoVec_Nor1Rcm2; tmp++)
		{
			int alignInfoIndex_Nor1 = finalAlignPair_Nor1Rcm2[tmp].first;
			int alignInfoIndex_Rcm2 = finalAlignPair_Nor1Rcm2[tmp].second;
			targetPeReadAlignInfo.pushBackNewAlignInfo_Nor2(norAlignmentInfo_PE_1[alignInfoIndex_Nor1]);
			targetPeReadAlignInfo.pushBackNewAlignInfo_Rcm1(rcmAlignmentInfo_PE_2[alignInfoIndex_Rcm2]);
			targetPeReadAlignInfo.pushBackNewFinalPairAlignInfoIndex_Nor2Rcm1(
				currentSize_targetNor2Rcm1FinalPairAlignInfoVec, currentSize_targetNor2Rcm1FinalPairAlignInfoVec);
			currentSize_targetNor2Rcm1FinalPairAlignInfoVec ++;
		}			
	}

	void alignmentFilter_fixPhase1_SE(int readLength_SE)
	{
		double min_map_score_keptAlignmentPair = ((readLength_SE)/100) * Min_Map_Score_keptAlignmentPair_PerHundredBases;

		this->getMismatchNumForEveryAlignment_SE();
		this->getEndMatchPosForEveryAlignment_SE();		
		this->generateAlignConfidenceMetric_SE();
		if((norAlignmentInfo_PE_1.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
			||(rcmAlignmentInfo_PE_1.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES))
		{
			return;
		}
		this->getScoreForEachAlignment_SE();
		if(norAlignmentInfo_PE_1.size() + rcmAlignmentInfo_PE_1.size() == 0)
		{
			return;
		}

		vector<int> finalAlign_Nor_tmp; 
		vector<int> finalAlign_Rcm_tmp;		
		vector< double > finalAlign_Nor_tmp_score;  
		vector< double > finalAlign_Rcm_tmp_score;	
	
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp ++)
		{		
			double tmpScore = seAlign_scoreVec_Nor[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlign_Nor_tmp.push_back(tmp);
				finalAlign_Nor_tmp_score.push_back(tmpScore);
			}
		}	
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp ++)
		{		
			double tmpScore = seAlign_scoreVec_Rcm[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)	
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlign_Rcm_tmp.push_back(tmp);
				finalAlign_Rcm_tmp_score.push_back(tmpScore);
			}
		}

		#ifdef DETECT_CIRCULAR_RNA
		this->filterOutBackSpliceAlignmentsIfNormalReadsExist_selectLeastBackSpliceDistanceOne(
			finalAlign_Nor_tmp, finalAlign_Rcm_tmp, finalAlign_Nor_tmp_score, finalAlign_Rcm_tmp_score,
			true);
		#endif

		//cout << "start to do filterOverlapFinalAlignPair ..." << endl;
		this->filterOverlapFinalAlign_SE(finalAlign_Nor_tmp, finalAlign_Rcm_tmp,
			finalAlign_Nor_tmp_score, finalAlign_Rcm_tmp_score, true);
		return;
	}

	void chooseBestAlignment_final_SE()
	{
		//cout << "start to chooseBestAlignment_final_SE ..." << endl;
		this->getMismatchNumForEveryAlignment_SE();
		this->getEndMatchPosForEveryAlignment_SE();		
		this->generateAlignConfidenceMetric_SE();
		this->getScoreForEachAlignment_SE();
		//cout << "seAlign_scoreVec_Nor.size(): " << seAlign_scoreVec_Nor.size() << endl;
		//cout << "seAlign_scoreVec_Rcm.size(): " << seAlign_scoreVec_Rcm.size() << endl;
		if(seAlign_scoreVec_Nor.size() + seAlign_scoreVec_Rcm.size() == 0)
			return;

		double DifferenceMin = BEST_ALIGNMENT_DIFFERENCE_MAX;

		double tmpBestScore_Nor = 0.0;
		int selectedBestAlignmentNO_Nor = 0;// = 0;
		double tmpBestScore_Rcm = 0.0;
		int selectedBestAlignmentNO_Rcm = 0;// = 0;		

		for(int tmp = 0; tmp < seAlign_scoreVec_Nor.size(); tmp++)
		{
			double tmpScore = seAlign_scoreVec_Nor[tmp];
			if(tmpScore > tmpBestScore_Nor)
			{
				tmpBestScore_Nor = tmpScore;
				selectedBestAlignmentNO_Nor = tmp;
			}
		}
		for(int tmp = 0; tmp < seAlign_scoreVec_Rcm.size(); tmp++)
		{
			double tmpScore = seAlign_scoreVec_Rcm[tmp];
			if(tmpScore > tmpBestScore_Rcm)
			{
				tmpBestScore_Rcm = tmpScore;
				selectedBestAlignmentNO_Rcm = tmp;
			}
		}

		bool bestAlignmentInNor_bool = false;
		bool bestAlignmentInRcm_bool = false;
		int bestAlignemnt_toCheck_index_start_Nor = 0;
		int bestAlignemnt_toCheck_index_start_Rcm = 0;

		if(fabs(tmpBestScore_Nor - tmpBestScore_Rcm) < DifferenceMin)
		{
			highestPairAlignmentScore = tmpBestScore_Nor;

			bestAlignmentInNor_bool = true;
			bestAlignemnt_toCheck_index_start_Nor = selectedBestAlignmentNO_Nor;
			bestAlignmentInRcm_bool = true;
			bestAlignemnt_toCheck_index_start_Rcm = selectedBestAlignmentNO_Rcm;
		}
		else if(tmpBestScore_Nor > tmpBestScore_Rcm)
		{
			highestPairAlignmentScore = tmpBestScore_Nor;

			bestAlignmentInNor_bool = true;
			bestAlignemnt_toCheck_index_start_Nor = selectedBestAlignmentNO_Nor;			
		}
		else
		{
			highestPairAlignmentScore = tmpBestScore_Rcm;

			bestAlignmentInRcm_bool = true;
			bestAlignemnt_toCheck_index_start_Rcm = selectedBestAlignmentNO_Rcm;			
		}

		vector<int> finalAlign_Nor_tmp; 
		vector<int> finalAlign_Rcm_tmp;		
		vector< double > finalAlign_Nor_tmp_score;  
		vector< double > finalAlign_Rcm_tmp_score;	

		if(bestAlignmentInNor_bool)
		{
			//cout << "norAlignmentInfo_PE_1.size(): " << norAlignmentInfo_PE_1.size() << endl;
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor1Rcm2*/; tmp < norAlignmentInfo_PE_1.size(); tmp ++)
			{
				double tmpScore = seAlign_scoreVec_Nor[tmp]; 
				//cout << "tmp: " << tmp << "  tmpScore: " << tmpScore << endl;
				if(fabs(tmpBestScore_Nor - tmpScore) < DifferenceMin)
				{
					//cout << "interBest index in Nor1: " << tmp << endl;
					finalAlign_Nor_tmp.push_back(tmp);
					finalAlign_Nor_tmp_score.push_back(tmpScore);
				}
			}
		}
		if(bestAlignmentInRcm_bool)
		{
			//cout << "rcmAlignmentInfo_PE_1.size(): " << rcmAlignmentInfo_PE_1.size() << endl;
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor2Rcm1*/; tmp < rcmAlignmentInfo_PE_1.size(); tmp ++)
			{
				double tmpScore = seAlign_scoreVec_Rcm[tmp];
				//cout << "tmp: " << tmp << "  tmpScore: " << tmpScore << endl;
				if(fabs(tmpBestScore_Rcm - tmpScore) < DifferenceMin)
				{
					//cout << "interBest index in Rcm1: " << tmp << endl;
					finalAlign_Rcm_tmp.push_back(tmp);
					finalAlign_Rcm_tmp_score.push_back(tmpScore);
				}
			}
		}

		// cout << "before filtering for cirRNA" << endl;
		// cout << "finalAlign_Nor_tmp.size(): " << finalAlign_Nor_tmp.size() << endl;
		// cout << "finalAlign_Rcm_tmp.size(): " << finalAlign_Rcm_tmp.size() << endl;

		#ifdef DETECT_CIRCULAR_RNA
		this->filterOutBackSpliceAlignmentsIfNormalReadsExist_selectLeastBackSpliceDistanceOne(
			finalAlign_Nor_tmp, finalAlign_Rcm_tmp, finalAlign_Nor_tmp_score, finalAlign_Rcm_tmp_score,
			false);
		#endif

		// cout << "after filtering for cirRNA" << endl;
		// cout << "finalAlign_Nor_tmp.size(): " << finalAlign_Nor_tmp.size() << endl;
		// cout << "finalAlign_Rcm_tmp.size(): " << finalAlign_Rcm_tmp.size() << endl;

		this->filterOverlapFinalAlign_SE(finalAlign_Nor_tmp, finalAlign_Rcm_tmp,
			finalAlign_Nor_tmp_score, finalAlign_Rcm_tmp_score, false);
	}

	void chooseBestAlignment_final_PEasSE_end1()
	{
		this->chooseBestAlignment_final_SE();
	}	

	void chooseBestAlignment_final_PEasSE_end2()
	{
		//cout << "chooseBestAlignment_final_PEasSE_end2 starts ......" << endl;
		//cout << "end2_nor.size(): " << norAlignmentInfo_PE_2.size() << endl;
		//cout << "end2_rcm.size(): " << rcmAlignmentInfo_PE_2.size() << endl;
		this->getMismatchNumForEveryAlignment_end2();
		this->getEndMatchPosForEveryAlignment_end2();		
		this->generateAlignConfidenceMetric_end2();
		this->getScoreForEachAlignment_end2();	
		//cout << "end of generating metrics ......" << endl;
		if(seAlign_scoreVec_Nor.size() + seAlign_scoreVec_Rcm.size() == 0)
			return;

		double DifferenceMin = BEST_ALIGNMENT_DIFFERENCE_MAX;

		double tmpBestScore_Nor = 0.0;
		int selectedBestAlignmentNO_Nor = 0;// = 0;
		double tmpBestScore_Rcm = 0.0;
		int selectedBestAlignmentNO_Rcm = 0;// = 0;		

		for(int tmp = 0; tmp < seAlign_scoreVec_Nor.size(); tmp++)
		{
			double tmpScore = seAlign_scoreVec_Nor[tmp];
			if(tmpScore > tmpBestScore_Nor)
			{
				tmpBestScore_Nor = tmpScore;
				selectedBestAlignmentNO_Nor = tmp;
			}
		}
		for(int tmp = 0; tmp < seAlign_scoreVec_Rcm.size(); tmp++)
		{
			double tmpScore = seAlign_scoreVec_Rcm[tmp];
			if(tmpScore > tmpBestScore_Rcm)
			{
				tmpBestScore_Rcm = tmpScore;
				selectedBestAlignmentNO_Rcm = tmp;
			}
		}
		//cout << "end of choosing best Alignment score ...." << endl;
		bool bestAlignmentInNor_bool = false;
		bool bestAlignmentInRcm_bool = false;
		int bestAlignemnt_toCheck_index_start_Nor = 0;
		int bestAlignemnt_toCheck_index_start_Rcm = 0;

		if(fabs(tmpBestScore_Nor - tmpBestScore_Rcm) < DifferenceMin)
		{
			highestPairAlignmentScore = tmpBestScore_Nor;

			bestAlignmentInNor_bool = true;
			bestAlignemnt_toCheck_index_start_Nor = selectedBestAlignmentNO_Nor;
			bestAlignmentInRcm_bool = true;
			bestAlignemnt_toCheck_index_start_Rcm = selectedBestAlignmentNO_Rcm;
		}
		else if(tmpBestScore_Nor > tmpBestScore_Rcm)
		{
			highestPairAlignmentScore = tmpBestScore_Nor;

			bestAlignmentInNor_bool = true;
			bestAlignemnt_toCheck_index_start_Nor = selectedBestAlignmentNO_Nor;			
		}
		else
		{
			highestPairAlignmentScore = tmpBestScore_Rcm;

			bestAlignmentInRcm_bool = true;
			bestAlignemnt_toCheck_index_start_Rcm = selectedBestAlignmentNO_Rcm;			
		}
		//cout << "end of getting best Alignment in Nor or Rcm " << endl;
		vector<int> finalAlign_Nor_tmp; 	
		vector<int> finalAlign_Rcm_tmp;		
		vector< double > finalAlign_Nor_tmp_score;  
		vector< double > finalAlign_Rcm_tmp_score;	

		if(bestAlignmentInNor_bool)
		{
			for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp ++)
			{
				double tmpScore = seAlign_scoreVec_Nor[tmp]; 
				if(fabs(tmpBestScore_Nor - tmpScore) < DifferenceMin)
				{
					finalAlign_Nor_tmp.push_back(tmp);
					finalAlign_Nor_tmp_score.push_back(tmpScore);
				}
			}
		}
		if(bestAlignmentInRcm_bool)
		{
			for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp ++)
			{
				double tmpScore = seAlign_scoreVec_Rcm[tmp];
				if(fabs(tmpBestScore_Rcm - tmpScore) < DifferenceMin)
				{
					finalAlign_Rcm_tmp.push_back(tmp);
					finalAlign_Rcm_tmp_score.push_back(tmpScore);
				}
			}
		}
		//cout << "start to filterOverlapFinalAlign_PEasSE_end2 ......" << endl;
		this->filterOverlapFinalAlign_PEasSE_end2(finalAlign_Nor_tmp, finalAlign_Rcm_tmp,
			finalAlign_Nor_tmp_score, finalAlign_Rcm_tmp_score, false);
	}		

	void chooseBestAlignment_final_PEasSE()
	{
		//cout << "start to chooseBestAlignment_final_PEasSE_end1 ..." << endl;
		this->chooseBestAlignment_final_PEasSE_end1();		
		//cout << "seAlignVec_final_Nor_end1.size(): " << seAlignVec_final_Nor.size() << endl;
		//cout << "seAlignVec_final_Rcm_end1.size(): " << seAlignVec_final_Rcm.size() << endl;
		for(int tmp = 0; tmp < seAlignVec_final_Nor.size(); tmp++)
		{	
			unpairedSEalignVec_final_Nor_end1.push_back(
				seAlignVec_final_Nor[tmp]);
		}
		for(int tmp = 0; tmp < seAlignVec_final_Rcm.size(); tmp++)
		{
			unpairedSEalignVec_final_Rcm_end1.push_back(
				seAlignVec_final_Rcm[tmp]);
		}

		//cout << "seAlignVec_final_Nor.size(): " << seAlignVec_final_Nor.size() << endl;
		//cout << "seAlignVec_final_Rcm.size(): " << seAlignVec_final_Rcm.size() << endl;
		//cout << "clear existing vectors ......" << endl;
		seAlignVec_final_Nor.clear();
		//vector<int>().swap(seAlignVec_final_Nor);
		seAlignVec_final_Rcm.clear();
		//vector<int>().swap(seAlignVec_final_Rcm);
		//cout << "seAlignVec_final_Nor.size(): " << seAlignVec_final_Nor.size() << endl;
		//cout << "seAlignVec_final_Rcm.size(): " << seAlignVec_final_Rcm.size() << endl;
		//cout << "clear existing vectors ......" << endl;
		seAlign_chrNameVec_Nor.clear();
		//vector<string>().swap(seAlign_chrNameVec_Nor);
		seAlign_startMapPosVec_Nor.clear();
		//vector<int>().swap(seAlign_startMapPosVec_Nor);
		seAlign_endMapPosVec_Nor.clear();
		//vector<int>().swap(seAlign_endMapPosVec_Nor);
		seAlign_mappedLengthVec_Nor.clear();
		//vector<int>().swap(seAlign_mappedLengthVec_Nor);
		seAlign_mismatchNumVec_Nor.clear();
		//vector<int>().swap(seAlign_mismatchNumVec_Nor);
		seAlign_insertionLengthVec_Nor.clear();
		//vector<int>().swap(seAlign_insertionLengthVec_Nor);
		seAlign_deletionLengthVec_Nor.clear();
		//vector<int>().swap(seAlign_deletionLengthVec_Nor);
		seAlign_SJconfidenceLevel_Nor.clear();
		//vector<int>().swap(seAlign_SJconfidenceLevel_Nor);
		seAlign_scoreVec_Nor.clear();
		//vector<double>().swap(seAlign_scoreVec_Nor);

		seAlign_chrNameVec_Rcm.clear();
		//vector<string>().swap(seAlign_chrNameVec_Rcm);
		seAlign_startMapPosVec_Rcm.clear();
		//vector<int>().swap(seAlign_startMapPosVec_Rcm);
		seAlign_endMapPosVec_Rcm.clear();
		//vector<int>().swap(seAlign_endMapPosVec_Rcm);
		seAlign_mappedLengthVec_Rcm.clear();
		//vector<int>().swap(seAlign_mappedLengthVec_Rcm);
		seAlign_mismatchNumVec_Rcm.clear();
		//vector<int>().swap(seAlign_mismatchNumVec_Rcm);
		seAlign_insertionLengthVec_Rcm.clear();
		//vector<int>().swap(seAlign_insertionLengthVec_Rcm);
		seAlign_deletionLengthVec_Rcm.clear();
		//vector<int>().swap(seAlign_deletionLengthVec_Rcm);
		seAlign_SJconfidenceLevel_Rcm.clear();
		//vector<int>().swap(seAlign_SJconfidenceLevel_Rcm);
		seAlign_scoreVec_Rcm.clear();
		//vector<double>().swap(seAlign_scoreVec_Rcm);
		//cout << "start to chooseBestAlignment_final_PEasSE_end2 ..." << endl;
		this->chooseBestAlignment_final_PEasSE_end2();
		//cout << "seAlignVec_final_Nor_end2.size(): " << seAlignVec_final_Nor.size() << endl;
		//cout << "seAlignVec_final_Rcm_end2.size(): " << seAlignVec_final_Rcm.size() << endl;
		for(int tmp = 0; tmp < seAlignVec_final_Nor.size(); tmp++)
		{	
			unpairedSEalignVec_final_Nor_end2.push_back(
				seAlignVec_final_Nor[tmp]);
		}
		for(int tmp = 0; tmp < seAlignVec_final_Rcm.size(); tmp++)
		{
			unpairedSEalignVec_final_Rcm_end2.push_back(
				seAlignVec_final_Rcm[tmp]);
		}
	}

	void filterOutBackSpliceAlignmentsIfNormalReadsExist_selectLeastBackSpliceDistanceOne(
		vector<int>& interAlign_vec_Nor, vector<int>& interAlign_vec_Rcm,
		vector< double >& interAlign_vec_Nor_score, vector< double >& interAlign_vec_Rcm_score,
		bool needToCheckAllAlignComplete_or_not) // if all complete, will be reported at phase1
	{
		if(needToCheckAllAlignComplete_or_not)
		{
			bool allTmpAlignmentCompleteBool = this->allTmpAlignment_complete_SE_bool(
				interAlign_vec_Nor, interAlign_vec_Rcm);
			if(!allTmpAlignmentCompleteBool)
				return;
		}
		//cout << "start to do filterOutBackSpliceAlignmentsIfNormalReadsExist_selectLeastBackSpliceDistanceOne ..." << endl;
		int interAlign_vec_size = interAlign_vec_Nor.size() + interAlign_vec_Rcm.size();
		//cout << "interAlign_vec_size: " << interAlign_vec_size << endl;
		if(interAlign_vec_size <= 1)
			return;

		// copy info in interAlign_Nor/Rcm_vec & interAlign_vec_Nor/Rcm_score
		vector<int> interAlign_vec_Nor_tmp;
		vector<int> interAlign_vec_Rcm_tmp;
		vector<double> interAlign_vec_Nor_score_tmp;
		vector<double> interAlign_vec_Rcm_score_tmp;
		for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp++)
		{
			interAlign_vec_Nor_tmp.push_back(interAlign_vec_Nor[tmp]);
			interAlign_vec_Nor_score_tmp.push_back(interAlign_vec_Nor_score[tmp]);
		}
		for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp++)
		{
			interAlign_vec_Rcm_tmp.push_back(interAlign_vec_Rcm[tmp]);
			interAlign_vec_Rcm_score_tmp.push_back(interAlign_vec_Rcm_score[tmp]);
		}

		// check nonBackSplice alignments
		vector<int> interAlign_vec_Nor_tmp_nonBackSplice;
		vector<int> interAlign_vec_Rcm_tmp_nonBackSplice;
		vector<double> interAlign_vec_Nor_score_tmp_nonBackSplice;
		vector<double> interAlign_vec_Rcm_score_tmp_nonBackSplice;		
		this->checkNonBackSpliceAlignments_SE(interAlign_vec_Nor_tmp, interAlign_vec_Rcm_tmp,
			interAlign_vec_Nor_score_tmp, interAlign_vec_Rcm_score_tmp,
			interAlign_vec_Nor_tmp_nonBackSplice, interAlign_vec_Rcm_tmp_nonBackSplice,
			interAlign_vec_Nor_score_tmp_nonBackSplice, interAlign_vec_Rcm_score_tmp_nonBackSplice);
		// cout << "interAlign_vec_Nor_tmp_nonBackSplice.size(): " << interAlign_vec_Nor.size() << endl;
		// cout << "interAlign_vec_Rcm_tmp_nonBackSplice.size(): " << interAlign_vec_Rcm.size() << endl;

		// cout << "before clearing" << endl;
		// cout << "interAlign_vec_Nor.size(): " << interAlign_vec_Nor.size() << endl;
		// cout << "interAlign_vec_Rcm.size(): " << interAlign_vec_Rcm.size() << endl;
		interAlign_vec_Nor.clear();
		interAlign_vec_Rcm.clear();
		interAlign_vec_Nor_score.clear();
		interAlign_vec_Rcm_score.clear();
		// cout << "after clearing" << endl;
		// cout << "interAlign_vec_Nor.size(): " << interAlign_vec_Nor.size() << endl;
		// cout << "interAlign_vec_Rcm.size(): " << interAlign_vec_Rcm.size() << endl;
		if(interAlign_vec_Nor_tmp_nonBackSplice.size() + 
			interAlign_vec_Rcm_tmp_nonBackSplice.size() > 0) // if nonBackSplice alignments exist, copy nonBackSplice alignments
		{
			//cout << "nonBackSplice alignments exist ......" << endl;
			for(int tmp = 0; tmp < interAlign_vec_Nor_tmp_nonBackSplice.size(); tmp++)
			{
				interAlign_vec_Nor.push_back(interAlign_vec_Nor_tmp_nonBackSplice[tmp]);
				interAlign_vec_Nor_score.push_back(interAlign_vec_Nor_score_tmp_nonBackSplice[tmp]);
			}
			for(int tmp = 0; tmp < interAlign_vec_Rcm_tmp_nonBackSplice.size(); tmp++)
			{
				interAlign_vec_Rcm.push_back(interAlign_vec_Rcm_tmp_nonBackSplice[tmp]);
				interAlign_vec_Rcm_score.push_back(interAlign_vec_Rcm_score_tmp_nonBackSplice[tmp]);
			}
		}
		else // if nonBackSplice alignments do not exist, select the best backSplice alignments.
		{
			//cout << "nonBackSplice alignments do not exist ......" << endl;
			this->selectTheBestBackSpliceAlignments_SE_withLeastBackSpliceTotalDistance(
				interAlign_vec_Nor_tmp, interAlign_vec_Rcm_tmp,
				interAlign_vec_Nor_score_tmp, interAlign_vec_Rcm_score_tmp,
				interAlign_vec_Nor, interAlign_vec_Rcm,
				interAlign_vec_Nor_score, interAlign_vec_Rcm_score);
		}
	}

	void checkNonBackSpliceAlignments_SE(vector<int>& interAlign_vec_Nor_tmp, vector<int>& interAlign_vec_Rcm_tmp,
			vector< double >& interAlign_vec_Nor_score_tmp, vector< double >& interAlign_vec_Rcm_score_tmp,
			vector<int>& interAlign_vec_Nor_tmp_nonBackSplice, vector<int>& interAlign_vec_Rcm_tmp_nonBackSplice,
			vector< double >& interAlign_vec_Nor_score_tmp_nonBackSplice, vector< double >& interAlign_vec_Rcm_score_tmp_nonBackSplice)
	{
		for(int tmp = 0; tmp < interAlign_vec_Nor_tmp.size(); tmp++)
		{
			int index_alignInfo_Nor1 = interAlign_vec_Nor_tmp[tmp];
			int score_alignInfo_Nor1 = interAlign_vec_Nor_score_tmp[tmp];
			bool tmpAlignInfo_backSpliceExistsOrNotBool 
				= norAlignmentInfo_PE_1[index_alignInfo_Nor1]->backSpliceExists_bool();
			if(!tmpAlignInfo_backSpliceExistsOrNotBool)
			{
				interAlign_vec_Nor_tmp_nonBackSplice.push_back(index_alignInfo_Nor1);
				interAlign_vec_Nor_score_tmp_nonBackSplice.push_back(score_alignInfo_Nor1);
			}
		}

		for(int tmp = 0; tmp < interAlign_vec_Rcm_tmp.size(); tmp++)
		{
			int index_alignInfo_Rcm1 = interAlign_vec_Rcm_tmp[tmp];
			int score_alignInfo_Rcm1 = interAlign_vec_Rcm_score_tmp[tmp];
			bool tmpAlignInfo_backSpliceExistsOrNotBool 
				= rcmAlignmentInfo_PE_1[index_alignInfo_Rcm1]->backSpliceExists_bool();
			if(!tmpAlignInfo_backSpliceExistsOrNotBool)
			{
				interAlign_vec_Rcm_tmp_nonBackSplice.push_back(index_alignInfo_Rcm1);
				interAlign_vec_Rcm_score_tmp_nonBackSplice.push_back(score_alignInfo_Rcm1);
			}
		}
	}

	void selectTheBestBackSpliceAlignments_SE_withLeastBackSpliceTotalDistance(
		vector<int>& interAlign_vec_Nor_tmp, vector<int>& interAlign_vec_Rcm_tmp,
		vector< double >& interAlign_vec_Nor_score_tmp, vector< double >& interAlign_vec_Rcm_score_tmp,
		vector<int>& interAlign_vec_Nor, vector<int>& interAlign_vec_Rcm,
		vector< double >& interAlign_vec_Nor_score, vector< double >& interAlign_vec_Rcm_score)
	{
		vector<int> backSpliceTotalDistanceVec_Nor1;
		vector<int> backSpliceTotalDistanceVec_Rcm1;

		// generate backSpliceTotalDistanceVec
		for(int tmp = 0; tmp < interAlign_vec_Nor_tmp.size(); tmp++)
		{
			int index_alignInfo_Nor1 = interAlign_vec_Nor_tmp[tmp];
			int tmpAlignInfo_backSpliceDistanceTotal 
				= norAlignmentInfo_PE_1[index_alignInfo_Nor1]->backSplicceDistanceTotal();
			backSpliceTotalDistanceVec_Nor1.push_back(tmpAlignInfo_backSpliceDistanceTotal);
		}

		for(int tmp = 0; tmp < interAlign_vec_Rcm_tmp.size(); tmp++)
		{
			int index_alignInfo_Rcm1 = interAlign_vec_Rcm_tmp[tmp];
			int tmpAlignInfo_backSpliceDistanceTotal
				= rcmAlignmentInfo_PE_1[index_alignInfo_Rcm1]->backSplicceDistanceTotal();
			backSpliceTotalDistanceVec_Rcm1.push_back(tmpAlignInfo_backSpliceDistanceTotal);
		}

		// get the leastBackSpliceTotalDistance among all the alignInfos
		int leastBackSpliceTotalDistance = 1000000;
		// int leastBackSpliceTotalDistance_startMapPos;
		// int leastBackSpliceTotalDistance_endMapPos;
		for(int tmp = 0; tmp < backSpliceTotalDistanceVec_Nor1.size(); tmp++)
		{
			int tmpBackSpliceTotalDistance = backSpliceTotalDistanceVec_Nor1[tmp];
			if(tmpBackSpliceTotalDistance < leastBackSpliceTotalDistance)
				leastBackSpliceTotalDistance = tmpBackSpliceTotalDistance;
		}
		for(int tmp = 0; tmp < backSpliceTotalDistanceVec_Rcm1.size(); tmp++)
		{
			int tmpBackSpliceTotalDistance = backSpliceTotalDistanceVec_Rcm1[tmp];
			if(tmpBackSpliceTotalDistance < leastBackSpliceTotalDistance)
				leastBackSpliceTotalDistance = tmpBackSpliceTotalDistance;
		}

		// select the best backSplice alignments with the least backSpliceDistance
		for(int tmp = 0; tmp < backSpliceTotalDistanceVec_Nor1.size(); tmp++)
		{
			int tmpBackSpliceTotalDistance = backSpliceTotalDistanceVec_Nor1[tmp];
			if(tmpBackSpliceTotalDistance == leastBackSpliceTotalDistance)
			{
				interAlign_vec_Nor.push_back(interAlign_vec_Nor_tmp[tmp]);
				interAlign_vec_Nor_score.push_back(interAlign_vec_Nor_score_tmp[tmp]);
			}
		}
		for(int tmp = 0; tmp < backSpliceTotalDistanceVec_Rcm1.size(); tmp++)
		{
			int tmpBackSpliceTotalDistance = backSpliceTotalDistanceVec_Rcm1[tmp];
			if(tmpBackSpliceTotalDistance == leastBackSpliceTotalDistance)
			{
				interAlign_vec_Rcm.push_back(interAlign_vec_Rcm_tmp[tmp]);
				interAlign_vec_Rcm_score.push_back(interAlign_vec_Rcm_score_tmp[tmp]);
			}
		}		

	}

	void getScoreForEachAlignment_SE()
	{
		int mismatch_weight = MISMATCH_PENALTY;
		int deletion_weight = DELETION_PENALTY;
		int insertion_weight = INSERTION_PENALTY;

		int semiCanonicalSJ_penalty = SEMICANONICALSJ_PENALTY;
		int nonCanonicalSJ_penalty = NONCANONICALSJ_PENALTY;

		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			int tmpMappedLength = seAlign_mappedLengthVec_Nor[tmp];
			int tmpMismatchNum = seAlign_mismatchNumVec_Nor[tmp];
			int tmpInsertionLength = seAlign_insertionLengthVec_Nor[tmp];
			int tmpDeletionLength =  seAlign_deletionLengthVec_Nor[tmp];
			
			int tmpSJpenalty = 0;
			int tmpSJconfidenceLevel = seAlign_SJconfidenceLevel_Nor[tmp];

			if(tmpSJconfidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
			{}
			else if(tmpSJconfidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
				tmpSJpenalty = semiCanonicalSJ_penalty;
			else  // noncanonical
				tmpSJpenalty = nonCanonicalSJ_penalty;

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpSJpenalty;
		
			seAlign_scoreVec_Nor.push_back(tmpScore);
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			int tmpMappedLength = seAlign_mappedLengthVec_Rcm[tmp];
			int tmpMismatchNum = seAlign_mismatchNumVec_Rcm[tmp];
			int tmpInsertionLength = seAlign_insertionLengthVec_Rcm[tmp];
			int tmpDeletionLength =  seAlign_deletionLengthVec_Rcm[tmp];
			
			int tmpSJpenalty = 0;
			int tmpSJconfidenceLevel = seAlign_SJconfidenceLevel_Rcm[tmp];

			if(tmpSJconfidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
			{}
			else if(tmpSJconfidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
				tmpSJpenalty = semiCanonicalSJ_penalty;
			else  // noncanonical
				tmpSJpenalty = nonCanonicalSJ_penalty;

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpSJpenalty;
		
			seAlign_scoreVec_Rcm.push_back(tmpScore);
		}		
	}

	void getScoreForEachAlignment_end1()
	{
		this->getScoreForEachAlignment_SE();
	}

	void getScoreForEachAlignment_end2()
	{
		int mismatch_weight = MISMATCH_PENALTY;
		int deletion_weight = DELETION_PENALTY;
		int insertion_weight = INSERTION_PENALTY;

		int semiCanonicalSJ_penalty = SEMICANONICALSJ_PENALTY;
		int nonCanonicalSJ_penalty = NONCANONICALSJ_PENALTY;

		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			int tmpMappedLength = seAlign_mappedLengthVec_Nor[tmp];
			int tmpMismatchNum = seAlign_mismatchNumVec_Nor[tmp];
			int tmpInsertionLength = seAlign_insertionLengthVec_Nor[tmp];
			int tmpDeletionLength =  seAlign_deletionLengthVec_Nor[tmp];
			
			int tmpSJpenalty = 0;
			int tmpSJconfidenceLevel = seAlign_SJconfidenceLevel_Nor[tmp];

			if(tmpSJconfidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
			{}
			else if(tmpSJconfidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
				tmpSJpenalty = semiCanonicalSJ_penalty;
			else  // noncanonical
				tmpSJpenalty = nonCanonicalSJ_penalty;

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpSJpenalty;
		
			seAlign_scoreVec_Nor.push_back(tmpScore);
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			int tmpMappedLength = seAlign_mappedLengthVec_Rcm[tmp];
			int tmpMismatchNum = seAlign_mismatchNumVec_Rcm[tmp];
			int tmpInsertionLength = seAlign_insertionLengthVec_Rcm[tmp];
			int tmpDeletionLength =  seAlign_deletionLengthVec_Rcm[tmp];
			
			int tmpSJpenalty = 0;
			int tmpSJconfidenceLevel = seAlign_SJconfidenceLevel_Rcm[tmp];

			if(tmpSJconfidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
			{}
			else if(tmpSJconfidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
				tmpSJpenalty = semiCanonicalSJ_penalty;
			else  // noncanonical
				tmpSJpenalty = nonCanonicalSJ_penalty;

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpSJpenalty;
		
			seAlign_scoreVec_Rcm.push_back(tmpScore);
		}		
	}	

	void generateAlignConfidenceMetric_SE()
	{
		for(int tmp_nor = 0; tmp_nor < norAlignmentInfo_PE_1.size(); tmp_nor++)
		{
			string tmpChrName = (norAlignmentInfo_PE_1[tmp_nor])->returnAlignChromName();
			seAlign_chrNameVec_Nor.push_back(tmpChrName);
			int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_1[tmp_nor])->returnAlignChromPos();
			seAlign_startMapPosVec_Nor.push_back(tmpChrMapPosStart_nor);
			int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_1[tmp_nor])->returnEndMatchedPosInChr();
			seAlign_endMapPosVec_Nor.push_back(tmpChrMapPosEnd_nor);
			int tmpMappedLength_nor = (norAlignmentInfo_PE_1[tmp_nor])->mappedLength();
			seAlign_mappedLengthVec_Nor.push_back(tmpMappedLength_nor);
			int tmpMismatchNum_nor = (norAlignmentInfo_PE_1[tmp_nor])->returnMismatchNum();
			seAlign_mismatchNumVec_Nor.push_back(tmpMismatchNum_nor);
			norAlignmentInfo_PE_1[tmp_nor]->generateIndelLength();
			int tmpInsertionLength_nor = norAlignmentInfo_PE_1[tmp_nor]->returnInsertionLength();
			seAlign_insertionLengthVec_Nor.push_back(tmpInsertionLength_nor);
			int tmpDeletionLength_nor = norAlignmentInfo_PE_1[tmp_nor]->returnDeletionLength();
			seAlign_deletionLengthVec_Nor.push_back(tmpDeletionLength_nor);
			int SJconfidenceLevel_nor = (norAlignmentInfo_PE_1[tmp_nor])->spliceJunctionConfidenceLevel();
			seAlign_SJconfidenceLevel_Nor.push_back(SJconfidenceLevel_nor);
		}
		for(int tmp_rcm = 0; tmp_rcm < rcmAlignmentInfo_PE_1.size(); tmp_rcm++)
		{
			string tmpChrName = (rcmAlignmentInfo_PE_1[tmp_rcm])->returnAlignChromName();
			seAlign_chrNameVec_Rcm.push_back(tmpChrName);
			int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_1[tmp_rcm])->returnAlignChromPos();
			seAlign_startMapPosVec_Rcm.push_back(tmpChrMapPosStart_rcm);
			int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_1[tmp_rcm])->returnEndMatchedPosInChr();
			seAlign_endMapPosVec_Rcm.push_back(tmpChrMapPosEnd_rcm);
			int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_1[tmp_rcm])->mappedLength();
			seAlign_mappedLengthVec_Rcm.push_back(tmpMappedLength_rcm);
			int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_1[tmp_rcm])->returnMismatchNum();
			seAlign_mismatchNumVec_Rcm.push_back(tmpMismatchNum_rcm);
			rcmAlignmentInfo_PE_1[tmp_rcm]->generateIndelLength();
			int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_1[tmp_rcm]->returnInsertionLength();
			seAlign_insertionLengthVec_Rcm.push_back(tmpInsertionLength_rcm);
			int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_1[tmp_rcm]->returnDeletionLength();
			seAlign_deletionLengthVec_Rcm.push_back(tmpDeletionLength_rcm);
			int SJconfidenceLevel_rcm = (rcmAlignmentInfo_PE_1[tmp_rcm])->spliceJunctionConfidenceLevel();
			seAlign_SJconfidenceLevel_Rcm.push_back(SJconfidenceLevel_rcm);
		}
	}

	void generateAlignConfidenceMetric_end1()
	{
		this->generateAlignConfidenceMetric_SE();
	}

	void generateAlignConfidenceMetric_end2()
	{
		for(int tmp_nor = 0; tmp_nor < norAlignmentInfo_PE_2.size(); tmp_nor++)
		{
			string tmpChrName = (norAlignmentInfo_PE_2[tmp_nor])->returnAlignChromName();
			seAlign_chrNameVec_Nor.push_back(tmpChrName);
			int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_2[tmp_nor])->returnAlignChromPos();
			seAlign_startMapPosVec_Nor.push_back(tmpChrMapPosStart_nor);
			int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_2[tmp_nor])->returnEndMatchedPosInChr();
			seAlign_endMapPosVec_Nor.push_back(tmpChrMapPosEnd_nor);
			int tmpMappedLength_nor = (norAlignmentInfo_PE_2[tmp_nor])->mappedLength();
			seAlign_mappedLengthVec_Nor.push_back(tmpMappedLength_nor);
			int tmpMismatchNum_nor = (norAlignmentInfo_PE_2[tmp_nor])->returnMismatchNum();
			seAlign_mismatchNumVec_Nor.push_back(tmpMismatchNum_nor);
			norAlignmentInfo_PE_2[tmp_nor]->generateIndelLength();
			int tmpInsertionLength_nor = norAlignmentInfo_PE_2[tmp_nor]->returnInsertionLength();
			seAlign_insertionLengthVec_Nor.push_back(tmpInsertionLength_nor);
			int tmpDeletionLength_nor = norAlignmentInfo_PE_2[tmp_nor]->returnDeletionLength();
			seAlign_deletionLengthVec_Nor.push_back(tmpDeletionLength_nor);
			int SJconfidenceLevel_nor = (norAlignmentInfo_PE_2[tmp_nor])->spliceJunctionConfidenceLevel();
			seAlign_SJconfidenceLevel_Nor.push_back(SJconfidenceLevel_nor);
		}
		for(int tmp_rcm = 0; tmp_rcm < rcmAlignmentInfo_PE_2.size(); tmp_rcm++)
		{
			string tmpChrName = (rcmAlignmentInfo_PE_2[tmp_rcm])->returnAlignChromName();
			seAlign_chrNameVec_Rcm.push_back(tmpChrName);
			int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_2[tmp_rcm])->returnAlignChromPos();
			seAlign_startMapPosVec_Rcm.push_back(tmpChrMapPosStart_rcm);
			int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_2[tmp_rcm])->returnEndMatchedPosInChr();
			seAlign_endMapPosVec_Rcm.push_back(tmpChrMapPosEnd_rcm);
			int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_2[tmp_rcm])->mappedLength();
			seAlign_mappedLengthVec_Rcm.push_back(tmpMappedLength_rcm);
			int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_2[tmp_rcm])->returnMismatchNum();
			seAlign_mismatchNumVec_Rcm.push_back(tmpMismatchNum_rcm);
			rcmAlignmentInfo_PE_2[tmp_rcm]->generateIndelLength();
			int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_2[tmp_rcm]->returnInsertionLength();
			seAlign_insertionLengthVec_Rcm.push_back(tmpInsertionLength_rcm);
			int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_2[tmp_rcm]->returnDeletionLength();
			seAlign_deletionLengthVec_Rcm.push_back(tmpDeletionLength_rcm);
			int SJconfidenceLevel_rcm = (rcmAlignmentInfo_PE_2[tmp_rcm])->spliceJunctionConfidenceLevel();
			seAlign_SJconfidenceLevel_Rcm.push_back(SJconfidenceLevel_rcm);
		}
	}

	// void generateAlignConfidenceMetric_bothEnds()
	// {
	// 	this->generateAlignConfidenceMetric_end1();
	// 	this->generateAlignConfidenceMetric_end2();
	// }

	bool twoAlignOverlapOrNot_SE(int index_1, int index_2, bool Nor_or_Rcm_bool)
	{
		int leftMapPos_1, rightMapPos_1, leftMapPos_2, rightMapPos_2;
		int startMapPos_1, endMapPos_1, startMapPos_2, endMapPos_2;		
		if(Nor_or_Rcm_bool)
		{
			if((norAlignmentInfo_PE_1[index_1]->returnAlignChromName()) 
				!= (norAlignmentInfo_PE_1[index_2]->returnAlignChromName()))
			{
				return false;
			}
			startMapPos_1 = norAlignmentInfo_PE_1[index_1]->returnAlignChromPos();
			endMapPos_1 = norAlignmentInfo_PE_1[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = norAlignmentInfo_PE_1[index_2]->returnAlignChromPos();
			endMapPos_2 = norAlignmentInfo_PE_1[index_2]->returnEndMatchedPosInChr();			
		}
		else
		{
			if((rcmAlignmentInfo_PE_1[index_1]->returnAlignChromName()) 
				!= (rcmAlignmentInfo_PE_1[index_2]->returnAlignChromName()))
			{
				return false;
			}
			startMapPos_1 = rcmAlignmentInfo_PE_1[index_1]->returnAlignChromPos();
			endMapPos_1 = rcmAlignmentInfo_PE_1[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = rcmAlignmentInfo_PE_1[index_2]->returnAlignChromPos();
			endMapPos_2 = rcmAlignmentInfo_PE_1[index_2]->returnEndMatchedPosInChr();
		}
		// generate leftMapPos rightMapPos, just to get rid of cirRNA backSplice reads.
		leftMapPos_1 = startMapPos_1;
		rightMapPos_1 = endMapPos_1;
		leftMapPos_2 = startMapPos_2;
		rightMapPos_2 = endMapPos_2;
		if(startMapPos_1 > endMapPos_1)
		{
			leftMapPos_1 = endMapPos_1;
			rightMapPos_1 = startMapPos_1;
		}
		if(startMapPos_2 > endMapPos_2)
		{
			leftMapPos_2 = endMapPos_2;
			rightMapPos_2 = startMapPos_2;
		}
		if((rightMapPos_1 < leftMapPos_2)||(rightMapPos_2 < leftMapPos_1))
			return false;
		else
			return true;
	}

	bool twoAlignOverlapOrNot_PEasSE_end2(int index_1, int index_2, bool Nor_or_Rcm_bool)
	{
		//cout << "twoAlignOverlapOrNot_PEasSE_end2 starts ..." << endl;
		//cout << "index_1: " << index_1 << endl;
		//cout << "index_2: " << index_2 << endl;
		//cout << "Nor_or_Rcm_bool: " << Nor_or_Rcm_bool << endl;
		int leftMapPos_1, rightMapPos_1, leftMapPos_2, rightMapPos_2;
		int startMapPos_1, endMapPos_1, startMapPos_2, endMapPos_2;		
		if(Nor_or_Rcm_bool)
		{
			if((norAlignmentInfo_PE_2[index_1]->returnAlignChromName()) 
				!= (norAlignmentInfo_PE_2[index_2]->returnAlignChromName()))
			{
				return false;
			}
			startMapPos_1 = norAlignmentInfo_PE_2[index_1]->returnAlignChromPos();
			endMapPos_1 = norAlignmentInfo_PE_2[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = norAlignmentInfo_PE_2[index_2]->returnAlignChromPos();
			endMapPos_2 = norAlignmentInfo_PE_2[index_2]->returnEndMatchedPosInChr();			
		}
		else
		{
			//cout << "chromName_1: " << rcmAlignmentInfo_PE_2[index_1]->returnAlignChromName() << endl;
			//cout << "chromName_2: " << rcmAlignmentInfo_PE_2[index_2]->returnAlignChromName() << endl;			
			if((rcmAlignmentInfo_PE_2[index_1]->returnAlignChromName()) 
				!= (rcmAlignmentInfo_PE_2[index_2]->returnAlignChromName()))
			{
				return false;
			}
			startMapPos_1 = rcmAlignmentInfo_PE_2[index_1]->returnAlignChromPos();
			//cout << "startMapPos_1: " << startMapPos_1 << endl;
			endMapPos_1 = rcmAlignmentInfo_PE_2[index_1]->returnEndMatchedPosInChr();
			//cout << "endMapPos_1: " << endMapPos_1 << endl;
			startMapPos_2 = rcmAlignmentInfo_PE_2[index_2]->returnAlignChromPos();
			//cout << "startMapPos_2: " << startMapPos_2 << endl; 
			endMapPos_2 = rcmAlignmentInfo_PE_2[index_2]->returnEndMatchedPosInChr();
			//cout << "endMapPos_2: " << endMapPos_2 << endl;
		}
		// generate leftMapPos rightMapPos, just to get rid of cirRNA backSplice reads.
		leftMapPos_1 = startMapPos_1;
		rightMapPos_1 = endMapPos_1;
		leftMapPos_2 = startMapPos_2;
		rightMapPos_2 = endMapPos_2;
		if(startMapPos_1 > endMapPos_1)
		{
			leftMapPos_1 = endMapPos_1;
			rightMapPos_1 = startMapPos_1;
		}
		if(startMapPos_2 > endMapPos_2)
		{
			leftMapPos_2 = endMapPos_2;
			rightMapPos_2 = startMapPos_2;
		}
		//cout << "leftMapPos_1: " << leftMapPos_1 << endl;
		//cout << "leftMapPos_2: " << leftMapPos_2 << endl;
		//cout << "rightMapPos_1: " << rightMapPos_1 << endl;
		//cout << "rightMapPos_2: " << rightMapPos_2 << endl; 
		if((rightMapPos_1 < leftMapPos_2)||(rightMapPos_2 < leftMapPos_1))
			return false;
		else
			return true;
	}

	void filterOverlapFinalAlign_SE(vector<int>& interAlign_vec_Nor,
		vector<int>& interAlign_vec_Rcm,
		vector< double >& interAlign_vec_Nor_score,
		vector< double >& interAlign_vec_Rcm_score,
		bool needToCheckAllAlignComplete_or_not)
	{
		set<int> invalidAlign_Nor; // val -- index in interAlign_vec_Nor
		set<int> invalidAlign_Rcm; // val -- index in interAlign_vec_Rcm

		// cout << "start to deal with interAlignPair_vec_Nor" << endl;
		for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlign_vec_Nor.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignOverlapOrNot_SE(interAlign_vec_Nor[tmp],
					interAlign_vec_Nor[tmp2], true);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlign_vec_Nor_score[tmp];
					double tmpScore_2 = interAlign_vec_Nor_score[tmp2];
					bool compare2Align_bool;
					if(tmpScore_1 > tmpScore_2)
					{
						compare2Align_bool = true;
					}
					else if(tmpScore_1 < tmpScore_2)
					{
						compare2Align_bool = false;
					}
					else
					{
						compare2Align_bool = this->selectSmallerMapRegionAlign_SE_bool(interAlign_vec_Nor[tmp],
							interAlign_vec_Nor[tmp2], true); 
					}

					if(compare2Align_bool) // tmp is better
					{
						invalidAlign_Nor.insert(tmp2);
					}
					else
					{	
						invalidAlign_Nor.insert(tmp);
					}
				}
			}
		}

		// cout << "start to deal with interAlign_vec_Rcm" << endl;
		for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlign_vec_Rcm.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignOverlapOrNot_SE(interAlign_vec_Rcm[tmp],
					interAlign_vec_Rcm[tmp2], false);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlign_vec_Rcm_score[tmp];
					double tmpScore_2 = interAlign_vec_Rcm_score[tmp2];
					bool compare2Align_bool;
					if(tmpScore_1 > tmpScore_2)
					{
						compare2Align_bool = true;
					}
					else if(tmpScore_1 < tmpScore_2)
					{
						compare2Align_bool = false;
					}
					else
					{
						compare2Align_bool = this->selectSmallerMapRegionAlign_SE_bool(interAlign_vec_Rcm[tmp],
							interAlign_vec_Rcm[tmp2], false); 
					}

					if(compare2Align_bool) // tmp is better
					{
						invalidAlign_Rcm.insert(tmp2);
					}
					else
					{	
						invalidAlign_Rcm.insert(tmp);
					}
				}
			}
		}		

		if(needToCheckAllAlignComplete_or_not)
		{
			vector< int > finalAlign_Nor_tmp;
			vector< int > finalAlign_Rcm_tmp;
			vector< double > finalAlign_Nor_tmp_score;
			vector< double > finalAlign_Rcm_tmp_score;

			for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp++)
			{
				if(invalidAlign_Nor.find(tmp) != invalidAlign_Nor.end())
				{}
				else
				{
					finalAlign_Nor_tmp.push_back(interAlign_vec_Nor[tmp]);
					finalAlign_Nor_tmp_score.push_back(interAlign_vec_Nor_score[tmp]);
				}
			}
			for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp++)
			{
				if(invalidAlign_Rcm.find(tmp) != invalidAlign_Rcm.end())
				{}
				else
				{
					finalAlign_Rcm_tmp.push_back(interAlign_vec_Rcm[tmp]);
					finalAlign_Rcm_tmp_score.push_back(interAlign_vec_Rcm_score[tmp]);
				}
			}

			bool allTmpFinalAlign_complete = this->allTmpAlignment_complete_SE_bool(
				finalAlign_Nor_tmp, finalAlign_Rcm_tmp);
			if(allTmpFinalAlign_complete)
			{
				for(int tmp = 0; tmp < finalAlign_Nor_tmp.size(); tmp++)
				{
					double tmpScore = finalAlign_Nor_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Nor.push_back(finalAlign_Nor_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlign_Rcm_tmp.size(); tmp++)
				{
					double tmpScore = finalAlign_Rcm_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Rcm.push_back(finalAlign_Rcm_tmp[tmp]);
				}
			}
			else
			{
				for(int tmp = 0; tmp < finalAlign_Nor_tmp.size(); tmp++)
				{
					//double tmpScore = finalAlign_Nor_tmp_score[tmp];
					//if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Nor.push_back(finalAlign_Nor_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlign_Rcm_tmp.size(); tmp++)
				{
					//double tmpScore = finalAlign_Rcm_tmp_score[tmp];
					//if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Rcm.push_back(finalAlign_Rcm_tmp[tmp]);
				}
			}
		}
		else
		{
			for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp++)
			{
				if(invalidAlign_Nor.find(tmp) != invalidAlign_Nor.end())
				{}
				else
				{
					seAlignVec_final_Nor.push_back(interAlign_vec_Nor[tmp]);
				}
			}
			for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp++)
			{
				if(invalidAlign_Rcm.find(tmp) != invalidAlign_Rcm.end())
				{}
				else
				{
					seAlignVec_final_Rcm.push_back(interAlign_vec_Rcm[tmp]);
				}
			}
		}
	}

	void filterOverlapFinalAlign_PEasSE_end2(vector<int>& interAlign_vec_Nor,
		vector<int>& interAlign_vec_Rcm,
		vector< double >& interAlign_vec_Nor_score,
		vector< double >& interAlign_vec_Rcm_score,
		bool needToCheckAllAlignComplete_or_not)
	{
		set<int> invalidAlign_Nor; // val -- index in interAlign_vec_Nor
		set<int> invalidAlign_Rcm; // val -- index in interAlign_vec_Rcm

		//cout << "start to deal with interAlignPair_vec_Nor" << endl;
		for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlign_vec_Nor.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignOverlapOrNot_PEasSE_end2(interAlign_vec_Nor[tmp],
					interAlign_vec_Nor[tmp2], true);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlign_vec_Nor_score[tmp];
					double tmpScore_2 = interAlign_vec_Nor_score[tmp2];
					bool compare2Align_bool;
					if(tmpScore_1 > tmpScore_2)
					{
						compare2Align_bool = true;
					}
					else if(tmpScore_1 < tmpScore_2)
					{
						compare2Align_bool = false;
					}
					else
					{
						compare2Align_bool = this->selectSmallerMapRegionAlign_PEasSE_end2_bool(interAlign_vec_Nor[tmp],
							interAlign_vec_Nor[tmp2], true); 
					}

					if(compare2Align_bool) // tmp is better
					{
						invalidAlign_Nor.insert(tmp2);
					}
					else
					{	
						invalidAlign_Nor.insert(tmp);
					}
				}
			}
		}

		//cout << "start to deal with interAlign_vec_Rcm" << endl;
		//cout << "interAlign_vec_Rcm.size(): " << interAlign_vec_Rcm.size() << endl;
		for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp ++)
		{
			//cout << "tmp in interAlign_vec_Rcm " << tmp << endl;
			for(int tmp2 = tmp+1; tmp2 < interAlign_vec_Rcm.size(); tmp2 ++)
			{
				//cout << "tmp2 in interAlign_vec_Rcm: " << tmp2 << endl;
				//cout << "start to do twoAlignOverlapOrNot_PEasSE_end2 ...." << endl;
				bool overlap_orNot_bool = this->twoAlignOverlapOrNot_PEasSE_end2(interAlign_vec_Rcm[tmp],
					interAlign_vec_Rcm[tmp2], false);
				//cout << "overlap_orNot_bool: " << overlap_orNot_bool << endl;
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlign_vec_Rcm_score[tmp];
					double tmpScore_2 = interAlign_vec_Rcm_score[tmp2];
					bool compare2Align_bool;
					if(tmpScore_1 > tmpScore_2)
					{
						compare2Align_bool = true;
					}
					else if(tmpScore_1 < tmpScore_2)
					{
						compare2Align_bool = false;
					}
					else
					{
						//cout << "start to do selectSmallerMapRegionAlign_PEasSE_end2_bool: " << compare2Align_bool << endl;
						compare2Align_bool = this->selectSmallerMapRegionAlign_PEasSE_end2_bool(interAlign_vec_Rcm[tmp],
							interAlign_vec_Rcm[tmp2], false); 
						//cout << "compare2Align_bool: " << compare2Align_bool << endl;
					}

					if(compare2Align_bool) // tmp is better
					{
						invalidAlign_Rcm.insert(tmp2);
					}
					else
					{	
						invalidAlign_Rcm.insert(tmp);
					}
				}
			}
		}		
		//cout << "needToCheckAllAlignComplete_or_not: " << needToCheckAllAlignComplete_or_not << endl;
		if(needToCheckAllAlignComplete_or_not)
		{
			vector< int > finalAlign_Nor_tmp;
			vector< int > finalAlign_Rcm_tmp;
			vector< double > finalAlign_Nor_tmp_score;
			vector< double > finalAlign_Rcm_tmp_score;

			for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp++)
			{
				if(invalidAlign_Nor.find(tmp) != invalidAlign_Nor.end())
				{}
				else
				{
					finalAlign_Nor_tmp.push_back(interAlign_vec_Nor[tmp]);
					finalAlign_Nor_tmp_score.push_back(interAlign_vec_Nor_score[tmp]);
				}
			}
			for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp++)
			{
				if(invalidAlign_Rcm.find(tmp) != invalidAlign_Rcm.end())
				{}
				else
				{
					finalAlign_Rcm_tmp.push_back(interAlign_vec_Rcm[tmp]);
					finalAlign_Rcm_tmp_score.push_back(interAlign_vec_Rcm_score[tmp]);
				}
			}

			bool allTmpFinalAlign_complete = this->allTmpAlignment_complete_PEasSE_end2_bool(
				finalAlign_Nor_tmp, finalAlign_Rcm_tmp);
			if(allTmpFinalAlign_complete)
			{
				for(int tmp = 0; tmp < finalAlign_Nor_tmp.size(); tmp++)
				{
					double tmpScore = finalAlign_Nor_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Nor.push_back(finalAlign_Nor_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlign_Rcm_tmp.size(); tmp++)
				{
					double tmpScore = finalAlign_Rcm_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Rcm.push_back(finalAlign_Rcm_tmp[tmp]);
				}
			}
			else
			{
				for(int tmp = 0; tmp < finalAlign_Nor_tmp.size(); tmp++)
				{
					//double tmpScore = finalAlign_Nor_tmp_score[tmp];
					//if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Nor.push_back(finalAlign_Nor_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlign_Rcm_tmp.size(); tmp++)
				{
					//double tmpScore = finalAlign_Rcm_tmp_score[tmp];
					//if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						seAlignVec_final_Rcm.push_back(finalAlign_Rcm_tmp[tmp]);
				}
			}
		}
		else
		{
			for(int tmp = 0; tmp < interAlign_vec_Nor.size(); tmp++)
			{
				if(invalidAlign_Nor.find(tmp) != invalidAlign_Nor.end())
				{}
				else
				{
					seAlignVec_final_Nor.push_back(interAlign_vec_Nor[tmp]);
				}
			}
			for(int tmp = 0; tmp < interAlign_vec_Rcm.size(); tmp++)
			{
				if(invalidAlign_Rcm.find(tmp) != invalidAlign_Rcm.end())
				{}
				else
				{
					seAlignVec_final_Rcm.push_back(interAlign_vec_Rcm[tmp]);
				}
			}
		}
	}

	bool allTmpAlignment_complete_SE_bool(vector<int>& nor_align_vec, vector<int>& rcm_align_vec)
	{
		for(int tmp = 0; tmp < nor_align_vec.size(); tmp++)
		{
			int tmpIndex_Nor1 = nor_align_vec[tmp];
			if(!(norAlignmentInfo_PE_1[tmpIndex_Nor1]->noUnfixedHeadTailBool()))
				return false;
		}
		for(int tmp = 0; tmp < rcm_align_vec.size(); tmp++)
		{
			int tmpIndex_Rcm1 = rcm_align_vec[tmp];
			if(!(rcmAlignmentInfo_PE_1[tmpIndex_Rcm1]->noUnfixedHeadTailBool()))
				return false;
		}
		return true;
	}

	bool allTmpAlignment_complete_PEasSE_end2_bool(vector<int>& nor_align_vec, vector<int>& rcm_align_vec)
	{
		for(int tmp = 0; tmp < nor_align_vec.size(); tmp++)
		{
			int tmpIndex_Nor2 = nor_align_vec[tmp];
			if(!(norAlignmentInfo_PE_2[tmpIndex_Nor2]->noUnfixedHeadTailBool()))
				return false;
		}
		for(int tmp = 0; tmp < rcm_align_vec.size(); tmp++)
		{
			int tmpIndex_Rcm2 = rcm_align_vec[tmp];
			if(!(rcmAlignmentInfo_PE_2[tmpIndex_Rcm2]->noUnfixedHeadTailBool()))
				return false;
		}
		return true;
	}

	bool selectSmallerMapRegionAlign_SE_bool(int index_1, int index_2, bool Nor_or_Rcm_bool)
	{
		int startMapPos_1, endMapPos_1, startMapPos_2, endMapPos_2;
		if(Nor_or_Rcm_bool)
		{
			startMapPos_1 = norAlignmentInfo_PE_1[index_1]->returnAlignChromPos();
			endMapPos_1 = norAlignmentInfo_PE_1[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = norAlignmentInfo_PE_1[index_2]->returnAlignChromPos();
			endMapPos_2 = norAlignmentInfo_PE_1[index_2]->returnEndMatchedPosInChr();
		}
		else
		{
			startMapPos_1 = rcmAlignmentInfo_PE_1[index_1]->returnAlignChromPos();
			endMapPos_1 = rcmAlignmentInfo_PE_1[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = rcmAlignmentInfo_PE_1[index_2]->returnAlignChromPos();
			endMapPos_2 = rcmAlignmentInfo_PE_1[index_2]->returnEndMatchedPosInChr();
		}
		int mappingRegionSize_1 = endMapPos_1 - startMapPos_1;
		int mappingRegionSize_2 = endMapPos_2 - startMapPos_2;

		#ifdef DETECT_CIRCULAR_RNA
		int mappingRegionSize_1_relative;
		int mappingRegionSize_2_relative;
		if(mappingRegionSize_1 < 0)
			mappingRegionSize_1_relative = 0 - mappingRegionSize_1;
		else
			mappingRegionSize_1_relative = mappingRegionSize_1;
		if(mappingRegionSize_2 < 0)
			mappingRegionSize_2_relative = 0 - mappingRegionSize_2;
		else
			mappingRegionSize_2_relative = mappingRegionSize_2;
		if(mappingRegionSize_1_relative <= mappingRegionSize_2_relative)
			return true;
		else
			return false;
		#else
		if(mappingRegionSize_1 <= mappingRegionSize_2)
			return true;
		else
			return false;
		#endif
	}

	bool selectSmallerMapRegionAlign_PEasSE_end2_bool(int index_1, int index_2, bool Nor_or_Rcm_bool)
	{
		int startMapPos_1, endMapPos_1, startMapPos_2, endMapPos_2;
		if(Nor_or_Rcm_bool)
		{
			startMapPos_1 = norAlignmentInfo_PE_2[index_1]->returnAlignChromPos();
			endMapPos_1 = norAlignmentInfo_PE_2[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = norAlignmentInfo_PE_2[index_2]->returnAlignChromPos();
			endMapPos_2 = norAlignmentInfo_PE_2[index_2]->returnEndMatchedPosInChr();
		}
		else
		{
			startMapPos_1 = rcmAlignmentInfo_PE_2[index_1]->returnAlignChromPos();
			endMapPos_1 = rcmAlignmentInfo_PE_2[index_1]->returnEndMatchedPosInChr();
			startMapPos_2 = rcmAlignmentInfo_PE_2[index_2]->returnAlignChromPos();
			endMapPos_2 = rcmAlignmentInfo_PE_2[index_2]->returnEndMatchedPosInChr();
		}
		int mappingRegionSize_1 = endMapPos_1 - startMapPos_1;
		int mappingRegionSize_2 = endMapPos_2 - startMapPos_2;

		#ifdef DETECT_CIRCULAR_RNA
		int mappingRegionSize_1_relative;
		int mappingRegionSize_2_relative;
		if(mappingRegionSize_1 < 0)
			mappingRegionSize_1_relative = 0 - mappingRegionSize_1;
		else
			mappingRegionSize_1_relative = mappingRegionSize_1;
		if(mappingRegionSize_2 < 0)
			mappingRegionSize_2_relative = 0 - mappingRegionSize_2;
		else
			mappingRegionSize_2_relative = mappingRegionSize_2;
		if(mappingRegionSize_1_relative <= mappingRegionSize_2_relative)
			return true;
		else
			return false;
		#else
		if(mappingRegionSize_1 <= mappingRegionSize_2)
			return true;
		else
			return false;
		#endif
	}

	void removeDuplicateMismatch(bool SE_or_PE_bool)
	{
		int alignType_max = 4;
		if(SE_or_PE_bool)
			alignType_max = 2;	
		for(int tmp = 1; tmp <= alignType_max; tmp++)
		{
			int tmpAlignInfoVecSize = this->getAlignInfoVecSize(tmp);
			for(int tmp2 = 0; tmp2 < tmpAlignInfoVecSize; tmp2++)
			{
				(this->returnAlignInfoInPeAlignInfo(tmp, tmp2))->removeDuplicateMismatch();				
			}
		}
	}

	int returnAlignInfoMismatchPosVecSize_type_index(int alignInfoKind, 
		int index_alignInfo)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index_alignInfo]->returnMismatchPosVecSize();
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index_alignInfo]->returnMismatchPosVecSize();
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index_alignInfo]->returnMismatchPosVecSize();
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index_alignInfo]->returnMismatchPosVecSize();
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}			
	}

	int returnAlignInfoMismatchPosVecValue_type_index(int alignInfoKind, 
		int index_alignInfo, int index_element)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index_alignInfo]->returnMismatchPosVecValue(
				index_element);
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index_alignInfo]->returnMismatchPosVecValue(
				index_element);
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index_alignInfo]->returnMismatchPosVecValue(
				index_element);
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index_alignInfo]->returnMismatchPosVecValue(
				index_element);
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}			
	}

	char returnAlignInfoMismatchCharVecValue_type_index(int alignInfoKind, 
		int index_alignInfo, int index_element)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index_alignInfo]->returnMismatchCharVecValue(
				index_element);
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index_alignInfo]->returnMismatchCharVecValue(
				index_element);
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index_alignInfo]->returnMismatchCharVecValue(
				index_element);
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index_alignInfo]->returnMismatchCharVecValue(
				index_element);
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}		
	}

	int returnAlignInfoJumpCodeLen_type_index(int alignInfoKind, 
		int index_alignInfo, int index_element)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index_alignInfo]->returnCigarStringJumpCodeLen(
				index_element);
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index_alignInfo]->returnCigarStringJumpCodeLen(
				index_element);
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index_alignInfo]->returnCigarStringJumpCodeLen(
				index_element);
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index_alignInfo]->returnCigarStringJumpCodeLen(
				index_element);
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}			
	}

	Jump_Code returnAlignInfoJumpCodeElement_type_index(int alignInfoKind, 
		int index_alignInfo, int index_element)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index_alignInfo]->returnCigarStringJumpCodeElement(index_element);
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index_alignInfo]->returnCigarStringJumpCodeElement(index_element);
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index_alignInfo]->returnCigarStringJumpCodeElement(index_element);
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index_alignInfo]->returnCigarStringJumpCodeElement(index_element);
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}			
	}

	int returnAlignInfoJumpCodeSize_type_index(int alignInfoKind, int index)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index]->returnCigarStringJumpCodeSize();
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index]->returnCigarStringJumpCodeSize();
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index]->returnCigarStringJumpCodeSize();
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index]->returnCigarStringJumpCodeSize();
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}			
	}

	string returnAlignInfoDirection_type_index(int alignInfoKind, int index)
	{
		if(alignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[index]->returnAlignDirection();
		}
		else if(alignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[index]->returnAlignDirection();
		}
		else if(alignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[index]->returnAlignDirection();
		}
		else if(alignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[index]->returnAlignDirection();
		}
		else
		{
			cout << "error in returnAlignInfoDirection_type_index(), peAlignInfo.h" << endl;
		}	
	}

	void memoryFree()
	{
		//vector< pair< int, vector<int> > >().swap(oriAlignPair_Nor1Rcm2);

		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->memoryFree();
			delete norAlignmentInfo_PE_1[tmp];
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->memoryFree();
			delete rcmAlignmentInfo_PE_1[tmp];
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->memoryFree();
			delete norAlignmentInfo_PE_2[tmp];
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->memoryFree();
			delete rcmAlignmentInfo_PE_2[tmp];
		}
	}

	void setAllPointersNULL()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp] = NULL;
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp] = NULL;
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp] = NULL;
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp] = NULL;
		}		
	}

	PE_Read_Alignment_Info(Path_Info& pathInfo_nor1, Path_Info& pathInfo_rcm1, 
			Path_Info& pathInfo_nor2, Path_Info& pathInfo_rcm2, Index_Info* indexInfo)
	{
		highestPairAlignmentScore = 0;
		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Rcm1 = -1;
		repeatRegion_index_Nor2 = -1;
		repeatRegion_index_Rcm2 = -1;

		//cout << "start PE_Read_Alignment_Info " << endl;
		repeatRegion_index_Nor1 = pathInfo_nor1.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_nor1.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
				= pathInfo_nor1.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).second;
				= pathInfo_nor1.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= (pathInfo_nor1->fixedPathMismatchVec)[tmpPath];
				= pathInfo_nor1.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			//Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
			//	mapChromPosInt, (pathInfo_nor1->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			//	tmpMismatch, indexInfo);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_nor1.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_nor1.fixedPathMismatchPosVec)[tmpPath], (pathInfo_nor1.fixedPathMismatchCharVec)[tmpPath]);
			norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Nor1 ends ..." << endl;
		repeatRegion_index_Rcm1 = pathInfo_rcm1.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_rcm1.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).first;
				= pathInfo_rcm1.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= //(((pathInfo_rcm1->finalPathVec)[tmpPath]).first).second;
				= pathInfo_rcm1.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= //(pathInfo_rcm1->fixedPathMismatchVec)[tmpPath];
				= pathInfo_rcm1.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_rcm1.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_rcm1.fixedPathMismatchPosVec)[tmpPath], (pathInfo_rcm1.fixedPathMismatchCharVec)[tmpPath]);
			rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Rcm1 ends ..." << endl;
		repeatRegion_index_Nor2 = pathInfo_nor2.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_nor2.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= //(((pathInfo_nor2->finalPathVec)[tmpPath]).first).first;
				= pathInfo_nor2.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= //(((pathInfo_nor2->finalPathVec)[tmpPath]).first).second;
				= pathInfo_nor2.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= //(pathInfo_nor2->fixedPathMismatchVec)[tmpPath];
				= pathInfo_nor2.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_nor2.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_nor2.fixedPathMismatchPosVec)[tmpPath], (pathInfo_nor2.fixedPathMismatchCharVec)[tmpPath]);
			norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Nor2 ends ..." << endl;
		repeatRegion_index_Rcm2 = pathInfo_rcm2.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_rcm2.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).first;
				= pathInfo_rcm2.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).second;
				= pathInfo_rcm2.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= //(pathInfo_rcm2->fixedPathMismatchVec)[tmpPath];
				= pathInfo_rcm2.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_rcm2.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_rcm2.fixedPathMismatchPosVec)[tmpPath], (pathInfo_rcm2.fixedPathMismatchCharVec)[tmpPath]);
			rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Rcm2 ends ..." << endl;
	}

	void assignHighestPairAlignmentScore(double highestScore)
	{
		highestPairAlignmentScore = highestScore;
	}
	// void initiatePeAlignInfo_withGroupSegInfoBothDir(
	// 	GroupSeg_Info_BothDir& groupSegInfoBothDir,
	// 	Index_Info* indexInfo)
	// {
	// 	highestPairAlignmentScore = 0;

	// 	repeatRegion_index_Nor1 = -1;
	// 	repeatRegion_index_Rcm1 = -1;		
		
	// }

	void initiatePeAlignInfo_Nor1Rcm2Only(
		Path_Info& pathInfo_nor1, Path_Info& pathInfo_rcm2, Index_Info* indexInfo)
	{
		highestPairAlignmentScore = 0;
		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Nor1 = pathInfo_nor1.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_nor1.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
				= pathInfo_nor1.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).second;
				= pathInfo_nor1.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= (pathInfo_nor1->fixedPathMismatchVec)[tmpPath];
				= pathInfo_nor1.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			//Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
			//	mapChromPosInt, (pathInfo_nor1->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			//	tmpMismatch, indexInfo);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_nor1.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_nor1.fixedPathMismatchPosVec)[tmpPath], (pathInfo_nor1.fixedPathMismatchCharVec)[tmpPath]);
			norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}				
		repeatRegion_index_Rcm2 = -1;
		repeatRegion_index_Rcm2 = pathInfo_rcm2.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_rcm2.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).first;
				= pathInfo_rcm2.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).second;
				= pathInfo_rcm2.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= //(pathInfo_rcm2->fixedPathMismatchVec)[tmpPath];
				= pathInfo_rcm2.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_rcm2.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_rcm2.fixedPathMismatchPosVec)[tmpPath], (pathInfo_rcm2.fixedPathMismatchCharVec)[tmpPath]);
			rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
	}

	void initiatePeAlignInfo(Path_Info& pathInfo_nor1, Path_Info& pathInfo_rcm1, 
			Path_Info& pathInfo_nor2, Path_Info& pathInfo_rcm2, Index_Info* indexInfo, bool SE_or_PE_bool)
	{
		highestPairAlignmentScore = 0;

		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Rcm1 = -1;
		//cout << "start PE_Read_Alignment_Info " << endl;
		repeatRegion_index_Nor1 = pathInfo_nor1.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_nor1.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
				= pathInfo_nor1.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= (((pathInfo_nor1->finalPathVec)[tmpPath]).first).second;
				= pathInfo_nor1.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= (pathInfo_nor1->fixedPathMismatchVec)[tmpPath];
				= pathInfo_nor1.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			//Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
			//	mapChromPosInt, (pathInfo_nor1->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			//	tmpMismatch, indexInfo);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_nor1.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_nor1.fixedPathMismatchPosVec)[tmpPath], (pathInfo_nor1.fixedPathMismatchCharVec)[tmpPath]);
			norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Nor1 ends ..." << endl;
		repeatRegion_index_Rcm1 = pathInfo_rcm1.returnRepeatRegion_index();
		for(int tmpPath = 0; tmpPath < pathInfo_rcm1.returnFinalPathVecSize(); tmpPath++)
		{
			int mapChromNameInt //= (((pathInfo_rcm1->finalPathVec)[tmpPath]).first).first;
				= pathInfo_rcm1.returnMapChrNameIntInFinalPathVec(tmpPath);
			int mapChromPosInt //= //(((pathInfo_rcm1->finalPathVec)[tmpPath]).first).second;
				= pathInfo_rcm1.returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= //(pathInfo_rcm1->fixedPathMismatchVec)[tmpPath];
				= pathInfo_rcm1.returnMismatchNumInFixedPathMismatchVec(tmpPath);
			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
				mapChromPosInt, (pathInfo_rcm1.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
				tmpMismatch, indexInfo, (pathInfo_rcm1.fixedPathMismatchPosVec)[tmpPath], (pathInfo_rcm1.fixedPathMismatchCharVec)[tmpPath]);
			rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			//delete(tmpAlignmentInfo);
		}
		//cout << "pe alignInfo Rcm1 ends ..." << endl;

		if(!SE_or_PE_bool)
		{	
			repeatRegion_index_Nor2 = -1;
			repeatRegion_index_Rcm2 = -1;
			
			repeatRegion_index_Nor2 = pathInfo_nor2.returnRepeatRegion_index();
			for(int tmpPath = 0; tmpPath < pathInfo_nor2.returnFinalPathVecSize(); tmpPath++)
			{
				int mapChromNameInt //= //(((pathInfo_nor2->finalPathVec)[tmpPath]).first).first;
					= pathInfo_nor2.returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= //(((pathInfo_nor2->finalPathVec)[tmpPath]).first).second;
					= pathInfo_nor2.returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= //(pathInfo_nor2->fixedPathMismatchVec)[tmpPath];
					= pathInfo_nor2.returnMismatchNumInFixedPathMismatchVec(tmpPath);
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt), 
					mapChromPosInt, (pathInfo_nor2.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, (pathInfo_nor2.fixedPathMismatchPosVec)[tmpPath], (pathInfo_nor2.fixedPathMismatchCharVec)[tmpPath]);
				norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
				//delete(tmpAlignmentInfo);
			}
			//cout << "pe alignInfo Nor2 ends ..." << endl;
			repeatRegion_index_Rcm2 = pathInfo_rcm2.returnRepeatRegion_index();
			for(int tmpPath = 0; tmpPath < pathInfo_rcm2.returnFinalPathVecSize(); tmpPath++)
			{
				int mapChromNameInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).first;
					= pathInfo_rcm2.returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= //(((pathInfo_rcm2->finalPathVec)[tmpPath]).first).second;
					= pathInfo_rcm2.returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= //(pathInfo_rcm2->fixedPathMismatchVec)[tmpPath];
					= pathInfo_rcm2.returnMismatchNumInFixedPathMismatchVec(tmpPath);
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
					mapChromPosInt, (pathInfo_rcm2.returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, (pathInfo_rcm2.fixedPathMismatchPosVec)[tmpPath], (pathInfo_rcm2.fixedPathMismatchCharVec)[tmpPath]);
				rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
				//delete(tmpAlignmentInfo);
			}
			//cout << "pe alignInfo Rcm2 ends ..." << endl;
		}	
	}

	PE_Read_Alignment_Info()
	{
		highestPairAlignmentScore = 0;
		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Rcm1 = -1;
		repeatRegion_index_Nor2 = -1;
		repeatRegion_index_Rcm2 = -1;
	}

	string returnPeAlignInfoStr()
	{
		string tmpStr = "-- Nor1 alignInfo: \n";
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + ". " + norAlignmentInfo_PE_1[tmp]->returnAlignInfoStr()
				 + " mismatchNum:" + int_to_str(norAlignmentInfo_PE_1[tmp]->returnMismatchNum()) 
				 + " mismatchPos:";
			for(int tmp2 = 0; tmp2 < norAlignmentInfo_PE_1[tmp]->returnMismatchPosVecSize(); tmp2++)
			{
				tmpStr = tmpStr + int_to_str(norAlignmentInfo_PE_1[tmp]->returnMismatchPosVecValue(tmp2)) + ",";
			}
			tmpStr += "\n";
		}
		tmpStr += "-- Rcm1 alignInfo:\n";
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + ". " + rcmAlignmentInfo_PE_1[tmp]->returnAlignInfoStr()
				 + " mismatchNum:" + int_to_str(rcmAlignmentInfo_PE_1[tmp]->returnMismatchNum()) 
				 + " mismatchPos:";
			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_1[tmp]->returnMismatchPosVecSize(); tmp2++)
			{
				tmpStr = tmpStr + int_to_str(rcmAlignmentInfo_PE_1[tmp]->returnMismatchPosVecValue(tmp2)) + ",";
			}
			tmpStr += "\n";
		}
		tmpStr += "-- Nor2 alignInfo:\n";
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + ". " + norAlignmentInfo_PE_2[tmp]->returnAlignInfoStr()
				 + " mismatchNum:" + int_to_str(norAlignmentInfo_PE_2[tmp]->returnMismatchNum()) 
				 + " mismatchPos:";
			for(int tmp2 = 0; tmp2 < norAlignmentInfo_PE_2[tmp]->returnMismatchPosVecSize(); tmp2++)
			{
				tmpStr = tmpStr + int_to_str(norAlignmentInfo_PE_2[tmp]->returnMismatchPosVecValue(tmp2)) + ",";
			}
			tmpStr += "\n";
		}
		tmpStr += "-- Rcm2 alignInfo:\n";
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + ". " + rcmAlignmentInfo_PE_2[tmp]->returnAlignInfoStr()
				 + " mismatchNum:" + int_to_str(rcmAlignmentInfo_PE_2[tmp]->returnMismatchNum()) 
				 + " mismatchPos:";
			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_2[tmp]->returnMismatchPosVecSize(); tmp2++)
			{
				tmpStr = tmpStr + int_to_str(rcmAlignmentInfo_PE_2[tmp]->returnMismatchPosVecValue(tmp2)) + ",";
			}
			tmpStr += "\n";
		}
		return tmpStr;
	}


	int returnRepeatRegion_index(bool nor_rcm_bool, bool end1_end2_bool)
	{
		if(nor_rcm_bool && end1_end2_bool)
		{
			return repeatRegion_index_Nor1;
		}
		else if((!nor_rcm_bool) && end1_end2_bool)
		{
			return repeatRegion_index_Rcm1;
		}
		else if(nor_rcm_bool && (!end1_end2_bool))
		{
			return repeatRegion_index_Nor2;
		}
		else 
			return repeatRegion_index_Rcm2;
	}

	void includeSyntheticSNPtransSeqPeAlignInfo(PE_Read_Alignment_Info& peAlignInfo_syntheticSNPtransSeq_end1, 
		PE_Read_Alignment_Info&	peAlignInfo_syntheticSNPtransSeq_end2, 
		Index_Info* syntheticSNPtransSeqIndexInfo, Index_Info* genomeIndexInfo, PE_Read_Info& readInfo)
	{
		repeatRegion_index_Nor1 = peAlignInfo_syntheticSNPtransSeq_end1.returnRepeatRegion_index(true, true);
		repeatRegion_index_Rcm1 = peAlignInfo_syntheticSNPtransSeq_end1.returnRepeatRegion_index(false, true);
		repeatRegion_index_Nor2 = peAlignInfo_syntheticSNPtransSeq_end2.returnRepeatRegion_index(true, true);
		repeatRegion_index_Rcm2 = peAlignInfo_syntheticSNPtransSeq_end2.returnRepeatRegion_index(false, true);
		peAlignInfo_syntheticSNPtransSeq_end1.getEndMatchPosForEveryAlignment();
		peAlignInfo_syntheticSNPtransSeq_end2.getEndMatchPosForEveryAlignment();		
		int ori_nor1_size = peAlignInfo_syntheticSNPtransSeq_end1.returnNorAlignmentInfo_PE_1_size();
		int ori_nor2_size = peAlignInfo_syntheticSNPtransSeq_end2.returnNorAlignmentInfo_PE_1_size();
		int ori_rcm1_size = peAlignInfo_syntheticSNPtransSeq_end1.returnRcmAlignmentInfo_PE_1_size();
		int ori_rcm2_size = peAlignInfo_syntheticSNPtransSeq_end2.returnRcmAlignmentInfo_PE_1_size();	
		for(int tmp = 0; tmp < ori_nor1_size; tmp++)
		{
			//cout << endl << "Nor1: " << tmp << endl;
			if((peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Nor1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Nor1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end1());
			if(!coverSNPorNot_bool)
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Nor1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			norAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm1_size; tmp++)
		{
			//cout << endl << "Rcm1: " << tmp << endl;
			if((peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Rcm1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Rcm1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end1());
			if(!coverSNPorNot_bool)
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(peAlignInfo_syntheticSNPtransSeq_end1.returnAlignInfo_Rcm1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_nor2_size; tmp++)
		{
			//cout << endl << "Nor2: " << tmp << endl;
			if((peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Nor1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Nor1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end2());
			if(!coverSNPorNot_bool)
				continue;			
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Nor1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			norAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm2_size; tmp++)
		{
			//cout << endl << "Rcm2: " << tmp << endl; 
			if((peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Rcm1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Rcm1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end2());
			if(!coverSNPorNot_bool)
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(peAlignInfo_syntheticSNPtransSeq_end2.returnAlignInfo_Rcm1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
	}

	void convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(PE_Read_Alignment_Info& syntheticSNPtransSeqAlignInfo,
		Index_Info* syntheticSNPtransSeqIndexInfo, Index_Info* genomeIndexInfo, PE_Read_Info& readInfo, bool SE_or_PE_bool)
	{
		repeatRegion_index_Nor1 = syntheticSNPtransSeqAlignInfo.returnRepeatRegion_index(true, true);
		repeatRegion_index_Rcm1 = syntheticSNPtransSeqAlignInfo.returnRepeatRegion_index(false, true);
		repeatRegion_index_Nor2 = syntheticSNPtransSeqAlignInfo.returnRepeatRegion_index(true, false);
		repeatRegion_index_Rcm2 = syntheticSNPtransSeqAlignInfo.returnRepeatRegion_index(false, false);
		syntheticSNPtransSeqAlignInfo.getEndMatchPosForEveryAlignment();
		int ori_nor1_size = syntheticSNPtransSeqAlignInfo.returnNorAlignmentInfo_PE_1_size();
		int ori_nor2_size = syntheticSNPtransSeqAlignInfo.returnNorAlignmentInfo_PE_2_size();
		int ori_rcm1_size = syntheticSNPtransSeqAlignInfo.returnRcmAlignmentInfo_PE_1_size();
		int ori_rcm2_size = syntheticSNPtransSeqAlignInfo.returnRcmAlignmentInfo_PE_2_size();			
		
		for(int tmp = 0; tmp < ori_nor1_size; tmp++)
		{
			//cout << endl << "Nor1: " << tmp << endl;
			if((syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end1());
			if(!coverSNPorNot_bool)
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			norAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm1_size; tmp++)
		{
			//cout << endl << "Rcm1: " << tmp << endl;
			if((syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm1(tmp))->containSJ())
				continue;
			bool coverSNPorNot_bool 
				= (syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm1(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
					readInfo.returnReadLength_end1());
			if(!coverSNPorNot_bool)
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
				(syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm1(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		if(!SE_or_PE_bool)
		{	
			for(int tmp = 0; tmp < ori_nor2_size; tmp++)
			{
				//cout << endl << "Nor2: " << tmp << endl;
				if((syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor2(tmp))->containSJ())
					continue;
				bool coverSNPorNot_bool 
					= (syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor2(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
						readInfo.returnReadLength_end2());
				if(!coverSNPorNot_bool)
					continue;			
				Alignment_Info* newAlignInfo = new Alignment_Info();
				newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
					(syntheticSNPtransSeqAlignInfo.returnAlignInfo_Nor2(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
				norAlignmentInfo_PE_2.push_back(newAlignInfo);
			}
			for(int tmp = 0; tmp < ori_rcm2_size; tmp++)
			{
				//cout << endl << "Rcm2: " << tmp << endl; 
				if((syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm2(tmp))->containSJ())
					continue;
				bool coverSNPorNot_bool 
					= (syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm2(tmp))->syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(
						readInfo.returnReadLength_end2());
				if(!coverSNPorNot_bool)
					continue;
				Alignment_Info* newAlignInfo = new Alignment_Info();
				newAlignInfo->convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(
					(syntheticSNPtransSeqAlignInfo.returnAlignInfo_Rcm2(tmp)), syntheticSNPtransSeqIndexInfo, genomeIndexInfo);
				rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
			}
		}
	}

	void convertTranscriptAlignInfo2GenomeAlignInfo(PE_Read_Alignment_Info& transcriptPeAlignInfo, 
		Transcript_Set* transcriptInfo, Index_Info* genomeIndexInfo)
	{
		//cout << "start to convertTranscriptAlignInfo2GenomeAlignInfo ..." << endl;
		repeatRegion_index_Nor1 = transcriptPeAlignInfo.returnRepeatRegion_index(true, true);
		repeatRegion_index_Rcm1 = transcriptPeAlignInfo.returnRepeatRegion_index(false, true);
		repeatRegion_index_Nor2 = transcriptPeAlignInfo.returnRepeatRegion_index(true, false);
		repeatRegion_index_Rcm2 = transcriptPeAlignInfo.returnRepeatRegion_index(false, false);
		transcriptPeAlignInfo.getEndMatchPosForEveryAlignment();
		int ori_nor1_size = transcriptPeAlignInfo.returnNorAlignmentInfo_PE_1_size();
		int ori_nor2_size = transcriptPeAlignInfo.returnNorAlignmentInfo_PE_2_size();
		int ori_rcm1_size = transcriptPeAlignInfo.returnRcmAlignmentInfo_PE_1_size();
		int ori_rcm2_size = transcriptPeAlignInfo.returnRcmAlignmentInfo_PE_2_size();	
		//cout << "ori_nor1_size: " << ori_nor1_size << endl;
		//cout << "ori_rcm1_size: " << ori_rcm1_size << endl;
		//cout << "ori_nor2_size: " << ori_nor2_size << endl;
		//cout << "ori_rcm2_size: " << ori_rcm2_size << endl;			
		for(int tmp = 0; tmp < ori_nor1_size; tmp++)
		{
			//cout << endl << "Nor1: " << tmp << endl;
			if((transcriptPeAlignInfo.returnAlignInfo_Nor1(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo.returnAlignInfo_Nor1(tmp)),
				transcriptInfo, genomeIndexInfo);
			norAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm1_size; tmp++)
		{
			//cout << endl << "Rcm1: " << tmp << endl;
			if((transcriptPeAlignInfo.returnAlignInfo_Rcm1(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo.returnAlignInfo_Rcm1(tmp)),
				transcriptInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_nor2_size; tmp++)
		{
			//cout << endl << "Nor2: " << tmp << endl;
			if((transcriptPeAlignInfo.returnAlignInfo_Nor2(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo.returnAlignInfo_Nor2(tmp)),
				transcriptInfo, genomeIndexInfo);
			norAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm2_size; tmp++)
		{
			//cout << endl << "Rcm2: " << tmp << endl; 
			if((transcriptPeAlignInfo.returnAlignInfo_Rcm2(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo.returnAlignInfo_Rcm2(tmp)),
				transcriptInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
	}

	void convertTranscriptAlignInfo2GenomeAlignInfo(PE_Read_Alignment_Info* transcriptPeAlignInfo, 
		Transcript_Set* transcriptInfo, Index_Info* genomeIndexInfo)
	{
		//cout << "start to convertTranscriptAlignInfo2GenomeAlignInfo ..." << endl;
		repeatRegion_index_Nor1 = transcriptPeAlignInfo->returnRepeatRegion_index(true, true);
		repeatRegion_index_Rcm1 = transcriptPeAlignInfo->returnRepeatRegion_index(false, true);
		repeatRegion_index_Nor2 = transcriptPeAlignInfo->returnRepeatRegion_index(true, false);
		repeatRegion_index_Rcm2 = transcriptPeAlignInfo->returnRepeatRegion_index(false, false);
		transcriptPeAlignInfo->getEndMatchPosForEveryAlignment();
		int ori_nor1_size = transcriptPeAlignInfo->returnNorAlignmentInfo_PE_1_size();
		int ori_nor2_size = transcriptPeAlignInfo->returnNorAlignmentInfo_PE_2_size();
		int ori_rcm1_size = transcriptPeAlignInfo->returnRcmAlignmentInfo_PE_1_size();
		int ori_rcm2_size = transcriptPeAlignInfo->returnRcmAlignmentInfo_PE_2_size();	
		//cout << "ori_nor1_size: " << ori_nor1_size << endl;
		//cout << "ori_rcm1_size: " << ori_rcm1_size << endl;
		//cout << "ori_nor2_size: " << ori_nor2_size << endl;
		//cout << "ori_rcm2_size: " << ori_rcm2_size << endl;			
		for(int tmp = 0; tmp < ori_nor1_size; tmp++)
		{
			//cout << endl << "Nor1: " << tmp << endl;
			if((transcriptPeAlignInfo->returnAlignInfo_Nor1(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo->returnAlignInfo_Nor1(tmp)),
				transcriptInfo, genomeIndexInfo);
			norAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm1_size; tmp++)
		{
			//cout << endl << "Rcm1: " << tmp << endl;
			if((transcriptPeAlignInfo->returnAlignInfo_Rcm1(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo->returnAlignInfo_Rcm1(tmp)),
				transcriptInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_nor2_size; tmp++)
		{
			//cout << endl << "Nor2: " << tmp << endl;
			if((transcriptPeAlignInfo->returnAlignInfo_Nor2(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo->returnAlignInfo_Nor2(tmp)),
				transcriptInfo, genomeIndexInfo);
			norAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
		for(int tmp = 0; tmp < ori_rcm2_size; tmp++)
		{
			//cout << endl << "Rcm2: " << tmp << endl; 
			if((transcriptPeAlignInfo->returnAlignInfo_Rcm2(tmp))->containSJ())
				continue;
			Alignment_Info* newAlignInfo = new Alignment_Info();
			newAlignInfo->convertTranscriptAlignInfo2GenomeAlignInfo(
				(transcriptPeAlignInfo->returnAlignInfo_Rcm2(tmp)),
				transcriptInfo, genomeIndexInfo);
			rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
	}

	void pairingAlignment2OriPair()
	{
		//cout << "start to pair ... " << endl;
		this->getMismatchNumForEveryAlignment();
		this->getEndMatchPosForEveryAlignment();

		Alignment_Info* tmpAlignInfo_1;
		Alignment_Info* tmpAlignInfo_2;
		bool newEntity = false;
		//cout << "start to pair Nor1Rcm2... " << endl;

		set < pair<int, pair<string, string> > >::iterator setIter;

		set < pair<int, pair<string, string> > > alignInfoSet_nor_1;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++) 
		//pair norAlignmentInfo_PE_1 & rcmAlignmentInfo_PE_2
		{
			//cout << "tmp in norAlignmentInfo_PE_1: " << tmp << endl; 
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmp];

			int tmpAlignInfo_1_chrom_pos = tmpAlignInfo_1->returnAlignChromPos();
			string tmpAlignInfo_1_chrom_name = tmpAlignInfo_1->returnAlignChromName();
			string tmpAlignInfo_1_cigar = tmpAlignInfo_1->returnCigarString();

			setIter = alignInfoSet_nor_1.find(pair<int, pair<string, string> > (
				tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar) ) );
			if(setIter == alignInfoSet_nor_1.end())
				alignInfoSet_nor_1.insert(pair<int, pair<string, string> > (
					tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar)));
			else
				continue;

			set < pair<int, pair<string, string> > > alignInfoSet_rcm_2;
			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_2.size(); tmp2++)
			{
				//cout << "tmp2 in rcmAlignmentInfo_PE_2: " << tmp2 << endl;
				tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmp2];
				int tmpAlignInfo_2_chrom_pos = tmpAlignInfo_2->returnAlignChromPos();
				string tmpAlignInfo_2_chrom_name = tmpAlignInfo_2->returnAlignChromName();
				string tmpAlignInfo_2_cigar = tmpAlignInfo_2->returnCigarString();

				setIter = alignInfoSet_rcm_2.find(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar) ) );
				if(setIter == alignInfoSet_rcm_2.end())
					alignInfoSet_rcm_2.insert(pair<int, pair<string, string> > (
						tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar)));
				else
					continue;

				if((tmpAlignInfo_1->returnAlignChromName()) == (tmpAlignInfo_2->returnAlignChromName()))
				{
					//cout << "chromName_1 == chromName_2 " << endl;
					#ifdef DETECT_CIRCULAR_RNA
					//cout << "DETECT_CIRCULAR_RNA defined " << endl;
					int tmpStartMapPos_1 = tmpAlignInfo_1->returnAlignChromPos();
					int tmpStartMapPos_2 = tmpAlignInfo_2->returnAlignChromPos();
					int tmpEndMapPos_1 = tmpAlignInfo_1->returnEndMatchedPosInChr();
					int tmpEndMapPos_2 = tmpAlignInfo_2->returnEndMatchedPosInChr();
					// cout << "tmpStartMapPos_1: " << tmpStartMapPos_1 << endl;
					// cout << "tmpStartMapPos_2: " << tmpStartMapPos_2 << endl;
					// cout << "tmpEndMapPos_1: " << tmpEndMapPos_1 << endl;
					// cout << "tmpEndMapPos_2: " << tmpEndMapPos_2 << endl;
					int tmpMapArea_startPos = selectTheSmallestAmong4values(
						tmpStartMapPos_1, tmpStartMapPos_2, tmpEndMapPos_1, tmpEndMapPos_2);
					//cout << "tmpMapArea_startPos: " << tmpMapArea_startPos << endl;
					int tmpMapArea_endPos = selectTheLargestAmong4values(
						tmpStartMapPos_1, tmpStartMapPos_2, tmpEndMapPos_1, tmpEndMapPos_2);
					//cout << "tmpMapArea_endPos: " << tmpMapArea_endPos << endl;
					//cout << "tmpMapAreaSize: " << tmpMapArea_endPos - tmpMapArea_startPos + 1 << endl;
					if(tmpMapArea_endPos - tmpMapArea_startPos + 1 <= READ_ALIGN_AREA_LENGTH)
					{	
						if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
						{
							vector<int> newTmpVec;
							newTmpVec.push_back(tmp2);
							oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
							newEntity = false;
						}
						else //tmp has already been in oriAlignPair_Nor1Rcm2
							(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
					}
					else // too far away
					{}
					#else
					if(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnAlignChromPos())
					{
						if( ((tmpAlignInfo_2->returnAlignChromPos()) - (tmpAlignInfo_1->returnEndMatchedPosInChr()))
							< PAIR_READ_DISTANCE_MAX)
						{
							//cout << "type 1" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
						}
						else
						{
							//cout << "too far away ..." << endl;
							//two far away 
						}
					}
					else if((tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr()))
					{
						//cout << "parts of read overlap ..." << endl;
						//cout << "tmp: " << tmp << " tmp 2: " << tmp2 << endl;
						//Note: In addition should check whether they cross the same SJs or not
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct!" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//cout << "overlap error !" << endl;
						}
					}
					else if(
						(tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							 )
					{
						//cout << "type 3" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{}
					}
					else if (
						(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							)
					{
						//cout << "type 4" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{}
					}
					else if (
							(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							&&(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							//&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&( (tmpAlignInfo_1->returnEndMatchedPosInChr())-(tmpAlignInfo_2->returnAlignChromPos()) < PAIR_READ_DISTANCE_MAX)							
							)
					{
						//cout << "type 5" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//cout << "overlap error " << endl;
						}
					}
					else
					{}
					#endif
				}
				else
				{}
			}
		}

		set < pair<int, pair<string, string> > > alignInfoSet_nor_2;
		//cout << "start to pair Nor2Rcm1... " << endl;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			//cout << "tmp in norAlignmentInfo_PE_2: " << tmp << endl; 
			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmp];

			int tmpAlignInfo_1_chrom_pos = tmpAlignInfo_1->returnAlignChromPos();
			string tmpAlignInfo_1_chrom_name = tmpAlignInfo_1->returnAlignChromName();
			string tmpAlignInfo_1_cigar = tmpAlignInfo_1->returnCigarString();

			setIter = alignInfoSet_nor_2.find(pair<int, pair<string, string> > (tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar) ) );
			if(setIter == alignInfoSet_nor_2.end())
				alignInfoSet_nor_2.insert(pair<int, pair<string, string> > (
					tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar)));
			else
				continue;			
			/*if(tmpAlignInfo_1 -> SJstrand == "X")   // FIX ME: do local optimization, shift to get consistant SJ strand
			{
				//cout << "SJstrand_1 = X" << endl;
				continue;
			}*/
			set < pair<int, pair<string, string> > > alignInfoSet_rcm_1;
			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_1.size(); tmp2++)
			{
				//cout << "tmp2 in rcmAlignmentInfo_PE_2: " << tmp2 << endl;
				tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmp2];
				int tmpAlignInfo_2_chrom_pos = tmpAlignInfo_2->returnAlignChromPos();
				string tmpAlignInfo_2_chrom_name = tmpAlignInfo_2->returnAlignChromName();
				string tmpAlignInfo_2_cigar = tmpAlignInfo_2->returnCigarString();

				setIter = alignInfoSet_rcm_1.find(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar) ) );
				if(setIter == alignInfoSet_rcm_1.end())
					alignInfoSet_rcm_1.insert(pair<int, pair<string, string> > (
						tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar)));
				else
					continue;
				
				if((tmpAlignInfo_1->returnAlignChromName()) == (tmpAlignInfo_2->returnAlignChromName()))
				{
					//cout << "chromName_1 == chromName_2 " << endl;
					#ifdef DETECT_CIRCULAR_RNA
					//cout << "DETECT_CIRCULAR_RNA defined " << endl;
					int tmpStartMapPos_1 = tmpAlignInfo_1->returnAlignChromPos();
					int tmpStartMapPos_2 = tmpAlignInfo_2->returnAlignChromPos();
					int tmpEndMapPos_1 = tmpAlignInfo_1->returnEndMatchedPosInChr();
					int tmpEndMapPos_2 = tmpAlignInfo_2->returnEndMatchedPosInChr();
					//cout << "tmpStartMapPos_1: " << tmpStartMapPos_1 << endl;
					//cout << "tmpStartMapPos_2: " << tmpStartMapPos_2 << endl;
					//cout << "tmpEndMapPos_1: " << tmpEndMapPos_1 << endl;
					//cout << "tmpEndMapPos_2: " << tmpEndMapPos_2 << endl;
					int tmpMapArea_startPos = selectTheSmallestAmong4values(
						tmpStartMapPos_1, tmpStartMapPos_2, tmpEndMapPos_1, tmpEndMapPos_2);
					//cout << "tmpMapArea_startPos: " << tmpMapArea_startPos << endl;
					int tmpMapArea_endPos = selectTheLargestAmong4values(
						tmpStartMapPos_1, tmpStartMapPos_2, tmpEndMapPos_1, tmpEndMapPos_2);
					//cout << "tmpMapArea_endPos: " << tmpMapArea_endPos << endl;
					//cout << "tmpMapAreaSize: " << tmpMapArea_endPos - tmpMapArea_startPos + 1 << endl;
					if(tmpMapArea_endPos - tmpMapArea_startPos + 1 <= READ_ALIGN_AREA_LENGTH)
					{	
						if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
						{
							vector<int> newTmpVec;
							newTmpVec.push_back(tmp2);
							oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
							newEntity = false;
						}
						else //tmp has already been in oriAlignPair_Nor1Rcm2
						{
							(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
							newEntity = false;
						}
					}
					else
					{}
					#else
					if(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnAlignChromPos())
					{
						//cout << "type 1" << endl;
						if(((tmpAlignInfo_2->returnAlignChromPos()) - (tmpAlignInfo_1->returnEndMatchedPosInChr())) < PAIR_READ_DISTANCE_MAX)
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{
							//two far away 
						}
					}
					else if((tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr()))
					{
						//cout << "parts of read overlap ..." << endl;
						//cout << "tmp: " << tmp << " tmp 2: " << tmp2 << endl;
						//Note: In addition should check whether they cross the same SJs or not
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct " << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor2Rcm1
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{
							//cout << "overlap error !" << endl;
						}
					}
					else if(
						(tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							)
					{
						//cout << "type 3" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor2Rcm1
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{}
					}
					else if (
						(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							)
					{
						//cout << "type 4" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor2Rcm1
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{}
					}
					else if (
							(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							&&(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							//&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&( (tmpAlignInfo_1->returnEndMatchedPosInChr())-(tmpAlignInfo_2->returnAlignChromPos()) < PAIR_READ_DISTANCE_MAX)							
							)
					{
						//cout << "type 5" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor2Rcm1
							{
								(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
								newEntity = false;
							}
						}
						else
						{}
					}
					else
					{}
					#endif
				}
				else
				{}
			}
		}
		///////////////// 1. only one end read mapped, another one unmapped ///////////////////////

		///////////////// 2. reads can be mapped under both directions //////////////////////////

		///////////////// 3. reads can be mapped to different places in one direction //////////////////
	}

	void pairingAlignment2OriPair_Nor1Rcm2Only(//int mappedAlignPairLen_min
		)
	{
		//cout << "start to pair ... " << endl;
		this->getMismatchNumForEveryAlignment_Nor1Rcm2();
		this->getEndMatchPosForEveryAlignment_Nor1Rcm2();

		Alignment_Info* tmpAlignInfo_1;
		Alignment_Info* tmpAlignInfo_2;
		bool newEntity = false;
		//cout << "start to pair Nor1Rcm2... " << endl;

		set < pair<int, pair<string, string> > >::iterator setIter;

		set < pair<int, pair<string, string> > > alignInfoSet_nor_1;

		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++) 
		//pair norAlignmentInfo_PE_1 & rcmAlignmentInfo_PE_2
		{

			newEntity = true;
			tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmp];

			int tmpAlignInfo_1_chrom_pos = tmpAlignInfo_1->returnAlignChromPos();
			string tmpAlignInfo_1_chrom_name = tmpAlignInfo_1->returnAlignChromName();
			string tmpAlignInfo_1_cigar = tmpAlignInfo_1->returnCigarString();

			setIter = alignInfoSet_nor_1.find(pair<int, pair<string, string> > (tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar) ) );
			if(setIter == alignInfoSet_nor_1.end())
			{
				alignInfoSet_nor_1.insert(pair<int, pair<string, string> > (tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar)));
			}
			else
			{
				continue;
			}

			/*if(tmpAlignInfo_1 -> SJstrand == "X")  // FIX ME: do local optimization, shift to get consistant SJ strand
			{
				//cout << "SJstrand_1 = X" << endl;
				continue;
			}*/

			set < pair<int, pair<string, string> > > alignInfoSet_rcm_2;

			for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_2.size(); tmp2++)
			{

				tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmp2];

				int tmpAlignInfo_2_chrom_pos = tmpAlignInfo_2->returnAlignChromPos();
				string tmpAlignInfo_2_chrom_name = tmpAlignInfo_2->returnAlignChromName();
				string tmpAlignInfo_2_cigar = tmpAlignInfo_2->returnCigarString();

				setIter = alignInfoSet_rcm_2.find(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar) ) );
				if(setIter == alignInfoSet_rcm_2.end())
				{
					alignInfoSet_rcm_2.insert(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar)));
				}
				else
				{
					continue;
				}
				
				/*if((tmpAlignInfo_2 -> SJstrand == "X")   // FIX ME: do local optimization, shift to get consistant SJ strand
					||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
					||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
				{
					//cout << "SJstrand_2 = X or conflict strand: " << tmpAlignInfo_1->SJstrand << ", " << tmpAlignInfo_2->SJstrand << endl;
					continue;
				}*/		

				if((tmpAlignInfo_1->returnAlignChromName()) == (tmpAlignInfo_2->returnAlignChromName()))
				{
					// if((tmpAlignInfo_1->mappedLength() + tmpAlignInfo_2->mappedLength()) < mappedAlignPairLen_min)
					// {
					// 	continue;
					// }

					if(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnAlignChromPos())
					{
						if( ((tmpAlignInfo_2->returnAlignChromPos()) - (tmpAlignInfo_1->returnEndMatchedPosInChr()))< PAIR_READ_DISTANCE_MAX)
						{
							//cout << "type 1" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//cout << "too far away ..." << endl;
							//two far away 
						}
					}
					else if((tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr()))
					{
						//cout << "parts of read overlap ..." << endl;
						//cout << "tmp: " << tmp << " tmp 2: " << tmp2 << endl;
						//Note: In addition should check whether they cross the same SJs or not
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct!" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//cout << "overlap error !" << endl;
						}
					}
					else if(
						(tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							 )
					{
						//cout << "type 3" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{}
					}
					else if (
						(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr())
							&&(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							)
					{
						//cout << "type 4" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{}
					}
					else if (
							(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
							&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
							&&(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
							&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
							//&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
							&&( (tmpAlignInfo_1->returnEndMatchedPosInChr())-(tmpAlignInfo_2->returnAlignChromPos()) < PAIR_READ_DISTANCE_MAX)							
							)
					{
						//cout << "type 5" << endl;
						if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
						{
							//cout << "overlap correct" << endl;
							if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
							{
								vector<int> newTmpVec;
								newTmpVec.push_back(tmp2);
								oriAlignPair_Nor1Rcm2.push_back(pair<int, vector<int> > (tmp, newTmpVec));
								newEntity = false;
							}
							else //tmp has already been in oriAlignPair_Nor1Rcm2
							{
								(oriAlignPair_Nor1Rcm2[oriAlignPair_Nor1Rcm2.size()-1].second).push_back(tmp2);
							}
						}
						else
						{
							//cout << "overlap error " << endl;
						}
					}
					else
					{

					}
				}
				else
				{

				}
			}
		}

		// set < pair<int, pair<string, string> > > alignInfoSet_nor_2;
		// //cout << "start to pair Nor2Rcm1... " << endl;
		// for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		// {
		// 	newEntity = true;
		// 	tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmp];

		// 	int tmpAlignInfo_1_chrom_pos = tmpAlignInfo_1->returnAlignChromPos();
		// 	string tmpAlignInfo_1_chrom_name = tmpAlignInfo_1->returnAlignChromName();
		// 	string tmpAlignInfo_1_cigar = tmpAlignInfo_1->returnCigarString();

		// 	setIter = alignInfoSet_nor_2.find(pair<int, pair<string, string> > (tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar) ) );
		// 	if(setIter == alignInfoSet_nor_2.end())
		// 	{
		// 		alignInfoSet_nor_2.insert(pair<int, pair<string, string> > (tmpAlignInfo_1_chrom_pos, pair<string, string>(tmpAlignInfo_1_chrom_name, tmpAlignInfo_1_cigar)));
		// 	}
		// 	else
		// 	{
		// 		continue;
		// 	}			
		// 	/*if(tmpAlignInfo_1 -> SJstrand == "X")   // FIX ME: do local optimization, shift to get consistant SJ strand
		// 	{
		// 		//cout << "SJstrand_1 = X" << endl;
		// 		continue;
		// 	}*/

		// 	set < pair<int, pair<string, string> > > alignInfoSet_rcm_1;

		// 	for(int tmp2 = 0; tmp2 < rcmAlignmentInfo_PE_1.size(); tmp2++)
		// 	{
		// 		tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmp2];

		// 		int tmpAlignInfo_2_chrom_pos = tmpAlignInfo_2->returnAlignChromPos();
		// 		string tmpAlignInfo_2_chrom_name = tmpAlignInfo_2->returnAlignChromName();
		// 		string tmpAlignInfo_2_cigar = tmpAlignInfo_2->returnCigarString();

		// 		setIter = alignInfoSet_rcm_1.find(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar) ) );
		// 		if(setIter == alignInfoSet_rcm_1.end())
		// 		{
		// 			alignInfoSet_rcm_1.insert(pair<int, pair<string, string> > (tmpAlignInfo_2_chrom_pos, pair<string, string>(tmpAlignInfo_2_chrom_name, tmpAlignInfo_2_cigar)));
		// 		}
		// 		else
		// 		{
		// 			continue;
		// 		}

		// 		/*if((tmpAlignInfo_2 -> SJstrand == "X")   // FIX ME: do local optimization, shift to get consistant SJ strand
		// 			||((tmpAlignInfo_1 -> SJstrand == "+")&&(tmpAlignInfo_2 -> SJstrand == "-"))
		// 			||((tmpAlignInfo_1 -> SJstrand == "-")&&(tmpAlignInfo_2 -> SJstrand == "+")))
		// 		{
		// 			//cout << "SJstrand_2 = X or conflict strand: " << tmpAlignInfo_1->SJstrand << ", " << tmpAlignInfo_2->SJstrand << endl;
		// 			continue;
		// 		}*/
				
		// 		if((tmpAlignInfo_1->returnAlignChromName()) == (tmpAlignInfo_2->returnAlignChromName()))
		// 		{
		// 			// if((tmpAlignInfo_1->mappedLength() + tmpAlignInfo_2->mappedLength()) < mappedAlignPairLen_min)
		// 			// {
		// 			// 	continue;
		// 			// }

		// 			if(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnAlignChromPos())
		// 			{
		// 				//cout << "type 1" << endl;
		// 				if(((tmpAlignInfo_2->returnAlignChromPos()) - (tmpAlignInfo_1->returnEndMatchedPosInChr())) < PAIR_READ_DISTANCE_MAX)
		// 				{
		// 					if(newEntity) //tmp is not in oriAlignPair_Nor1Rcm2
		// 					{
		// 						vector<int> newTmpVec;
		// 						newTmpVec.push_back(tmp2);
		// 						oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
		// 						newEntity = false;
		// 					}
		// 					else //tmp has already been in oriAlignPair_Nor1Rcm2
		// 					{
		// 						(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
		// 						newEntity = false;
		// 					}
		// 				}
		// 				else
		// 				{
		// 					//two far away 
		// 				}
		// 			}
		// 			else if((tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
		// 					&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr()))
		// 			{
		// 				//cout << "parts of read overlap ..." << endl;
		// 				//cout << "tmp: " << tmp << " tmp 2: " << tmp2 << endl;
		// 				//Note: In addition should check whether they cross the same SJs or not
		// 				if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
		// 				{
		// 					//cout << "overlap correct " << endl;
		// 					if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
		// 					{
		// 						vector<int> newTmpVec;
		// 						newTmpVec.push_back(tmp2);
		// 						oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
		// 						newEntity = false;
		// 					}
		// 					else //tmp has already been in oriAlignPair_Nor2Rcm1
		// 					{
		// 						(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
		// 						newEntity = false;
		// 					}
		// 				}
		// 				else
		// 				{
		// 					//cout << "overlap error !" << endl;
		// 				}
		// 			}
		// 			else if(
		// 				(tmpAlignInfo_1->returnAlignChromPos() <= tmpAlignInfo_2->returnAlignChromPos())
		// 					&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
		// 					&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
		// 					)
		// 			{
		// 				//cout << "type 3" << endl;
		// 				if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
		// 				{
		// 					if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
		// 					{
		// 						vector<int> newTmpVec;
		// 						newTmpVec.push_back(tmp2);
		// 						oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
		// 						newEntity = false;
		// 					}
		// 					else //tmp has already been in oriAlignPair_Nor2Rcm1
		// 					{
		// 						(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
		// 						newEntity = false;
		// 					}
		// 				}
		// 				else
		// 				{}
		// 			}
		// 			else if (
		// 				(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
		// 					&&(tmpAlignInfo_1->returnEndMatchedPosInChr() <= tmpAlignInfo_2->returnEndMatchedPosInChr())
		// 					&&(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
		// 					)
		// 			{
		// 				//cout << "type 4" << endl;
		// 				if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
		// 				{
		// 					if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
		// 					{
		// 						vector<int> newTmpVec;
		// 						newTmpVec.push_back(tmp2);
		// 						oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
		// 						newEntity = false;
		// 					}
		// 					else //tmp has already been in oriAlignPair_Nor2Rcm1
		// 					{
		// 						(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
		// 						newEntity = false;
		// 					}
		// 				}
		// 				else
		// 				{}
		// 			}
		// 			else if (
		// 					(tmpAlignInfo_1->unfixedHeadExistsBool()) // pe_1 read has unfixed head
		// 					&&(tmpAlignInfo_2->unfixedTailExistsBool()) // pe_2 read has unfixed tail
		// 					&&(tmpAlignInfo_1->returnAlignChromPos() > tmpAlignInfo_2->returnAlignChromPos())
		// 					&&(tmpAlignInfo_1->returnEndMatchedPosInChr() > tmpAlignInfo_2->returnEndMatchedPosInChr())
		// 					//&&(tmpAlignInfo_1->endMatchedPosInChr <= tmpAlignInfo_2->endMatchedPosInChr)
		// 					&&( (tmpAlignInfo_1->returnEndMatchedPosInChr())-(tmpAlignInfo_2->returnAlignChromPos()) < PAIR_READ_DISTANCE_MAX)							
		// 					)
		// 			{
		// 				//cout << "type 5" << endl;
		// 				if(tmpAlignInfo_1->checkOverlapPairAlignment_new(tmpAlignInfo_2))
		// 				{
		// 					if(newEntity) //tmp is not in oriAlignPair_Nor2Rcm1
		// 					{
		// 						vector<int> newTmpVec;
		// 						newTmpVec.push_back(tmp2);
		// 						oriAlignPair_Nor2Rcm1.push_back(pair<int, vector<int> > (tmp, newTmpVec));
		// 						newEntity = false;
		// 					}
		// 					else //tmp has already been in oriAlignPair_Nor2Rcm1
		// 					{
		// 						(oriAlignPair_Nor2Rcm1[oriAlignPair_Nor2Rcm1.size()-1].second).push_back(tmp2);
		// 						newEntity = false;
		// 					}
		// 				}
		// 				else
		// 				{}
		// 			}
		// 			else
		// 			{

		// 			}
		// 		}
		// 		else
		// 		{

		// 		}
		// 	}
		// }
		///////////////// 1. only one end read mapped, another one unmapped ///////////////////////

		///////////////// 2. reads can be mapped under both directions //////////////////////////

		///////////////// 3. reads can be mapped to different places in one direction //////////////////
	}

	bool check_shortPairEndDistance_orNot(int index, bool Nor1Rcm2_orNot)
	{
		if(Nor1Rcm2_orNot)
		{
			if(oriAlignPairNew_pairDistance_Nor1Rcm2[index] <= PAIR_READ_DISTANCE_SHORT_MAX)
			{
				return true;
			}
			else
				return false;
		}
		else
		{
			if(oriAlignPairNew_pairDistance_Nor2Rcm1[index] <= PAIR_READ_DISTANCE_SHORT_MAX)
			{
				return true;
			}
			else
				return false;
		}
	}


	void output_phase1_SE(
		bool Do_Phase1_Only,
		vector<string>& SeAlignSamStrVec_complete,
		vector<string>& SeAlignInfoStrVec_inComplete,
		vector<string>& SeAlignSamStrVec_inComplete,
		vector<string>& SeAlignSamStrVec_unmapped,
		vector<string>& SeAlignSamStrVec_unmapped_lowScore,
		vector<string>& SeAlignSamStrVec_unmapped_mappedToRepeatRegion,	
		vector<RepeatRegion_Info*>& repeatRegionInfoVec,
		PE_Read_Info& readInfo, Index_Info* indexInfo, int tmpOpenMP, int threadNO,
		Stats_Info* statsInfo, bool FastaOrFastq
		)
	{
		bool finalAlignExistsBool = this->finalAlignExistsBool_SE();
		bool mappedToRepeatRegionBool = this->mappedToRepeatRegionBool_SE();
	
		if(finalAlignExistsBool)
		{
			bool allFinalAlignComplete_bool = this->allFinalAlignmentComplete_SE();
			bool unique_bool = this->checkUniqueOrMulti_SE();
		    //statsInfo->increExistingAlignNum_phase1_SE(threadNO, allFinalAlignComplete_bool, unique_bool)
		    if(allFinalAlignComplete_bool)
		    {
		    	// int alignment_score_min_output = readInfo.returnAlignmentScoreMinOutput_withComplement_perHundredBases_SE(
							// 	ALIGNMENT_SCORE_MIN_OUTPUT_COMPLEMENT_PerHundredBases);
		    	// bool completeAlignmentScore_tooLow = this->alignScoreTooLow_bool_SE(alignment_score_min_output);
		    	bool completeAlignmentScore_tooLow = this->alignScoreTooLow_bool_SE(readInfo);
		    	statsInfo->increCompleteOrTooLowScore_phase1_SE(threadNO, completeAlignmentScore_tooLow, unique_bool);
		    	if(completeAlignmentScore_tooLow)
		    	{
		    		//statsInfo->increLowScoreComplete_phase1_SE(threadNO, unique_bool);
		    		SeAlignSamStrVec_unmapped_lowScore[tmpOpenMP] = this->getSAMformatForUnmapped_SE(readInfo, FastaOrFastq);
		    	}
		    	else
		    	{
		    		SeAlignSamStrVec_complete[tmpOpenMP] = this->getSAMformatForFinalAlignment_SE(readInfo, FastaOrFastq);
		    	}
		    }
		    else
		    {
		    	statsInfo->increIncomplete_phase1_SE(threadNO, unique_bool);
		    	if(!Do_Phase1_Only)
		    	{
		    		SeAlignInfoStrVec_inComplete[tmpOpenMP] = this->getTmpAlignInfoForFinalAlign_SE(readInfo.returnReadName_SE(),
		    			readInfo.returnReadSeq_SE(), readInfo.returnReadQual_SE(), FastaOrFastq);
		    	}
		    	SeAlignSamStrVec_inComplete[tmpOpenMP] = this->getSAMformatForFinalAlignment_SE(readInfo, FastaOrFastq);
		    }
		}
		else
		{	
			statsInfo->increUnmap_phase1_SE(threadNO, mappedToRepeatRegionBool);
			if(mappedToRepeatRegionBool) // one strand of read mapped to repeatRegion
			{
				SeAlignSamStrVec_unmapped_mappedToRepeatRegion[tmpOpenMP] = 
					this->getSamFormatMappedToRepeatRegion_SE(readInfo, repeatRegionInfoVec[threadNO], indexInfo, FastaOrFastq);
			}
			else // SE read unmapped
			{
				SeAlignSamStrVec_unmapped[tmpOpenMP] = 
					this->getSAMformatForUnmapped_SE(readInfo, FastaOrFastq);
			}
		}
	}

	string getTmpAlignInfoForFinalAlign_SE(const string& readName, 
		const string& readSeq, 
		const string& readQual, bool FastaOrFastq)
	{
		string tmpAlignInfoStr = "\n" + readName + "\t" 
			+ int_to_str(seAlignVec_final_Nor.size()) + "\t"
			+ int_to_str(seAlignVec_final_Rcm.size()) + "\n"
			+ readSeq + "\n" + readQual + "\n" + "\n" + "\n";

		Alignment_Info* tmpAlignInfo;
		string tmpNorAlignInfo;
		for(int tmpNor= 0; tmpNor < seAlignVec_final_Nor.size(); tmpNor ++)
		{
			if(seAlignVec_final_Nor.size() >= 40)
				{break;}

			int tmp = seAlignVec_final_Nor[tmpNor];
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo = tmpNorAlignInfo 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo = tmpNorAlignInfo + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo = tmpNorAlignInfo + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo = tmpNorAlignInfo + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo = tmpNorAlignInfo + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo += "\t";				
		}
		string tmpRcmAlignInfo;
		for(int tmpRcm = 0; tmpRcm < seAlignVec_final_Rcm.size(); tmpRcm++)
		{
			if(seAlignVec_final_Rcm.size() >= 40)
				{break;}
			int tmp = seAlignVec_final_Rcm[tmpRcm];
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo = tmpRcmAlignInfo
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";			
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo = tmpRcmAlignInfo + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo = tmpRcmAlignInfo + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo = tmpRcmAlignInfo + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo = tmpRcmAlignInfo + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo += "\t";
		}	
		//string tmpNorAlignInfo_2;
		//string tmpRcmAlignInfo_2;

		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo + "\n"
			+ "Nor_2:\t" + "\n"
			+ "Rcm_2:\t";
		return tmpAlignInfoStr;
	}

	string getSamFormatMappedToRepeatRegion_SE(PE_Read_Info& seReadInfo, RepeatRegion_Info* repeatRegionInfo,
		Index_Info* indexInfo, bool FastaOrFastq)
	{
		bool mappedToRepeatRegionBool_Nor1 = (repeatRegion_index_Nor1 >= 0);
		bool mappedToRepeatRegionBool_Rcm1 = (repeatRegion_index_Rcm1 >= 0);
	
		// fix me:
		int repeatIndex;
		string repeatInfoStr;
		if(mappedToRepeatRegionBool_Nor1)
			repeatIndex = repeatRegion_index_Nor1;
		else
			repeatIndex = repeatRegion_index_Rcm1;
		repeatInfoStr = "Rr:i:" + int_to_str(repeatIndex);
		string tmpQualSeq;
		if(FastaOrFastq)
		{
			tmpQualSeq = "*";
		}
		else
		{
			tmpQualSeq = seReadInfo.returnReadQual_SE();
		}
		string tmpReadName = seReadInfo.returnReadName_SE();

		string seAlignSamStr;
		seAlignSamStr = tmpReadName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + seReadInfo.returnReadSeq_SE()
			+ "\t" + tmpQualSeq + "\tIH:i:0\tHI:i:0\t" + repeatInfoStr; 
		return seAlignSamStr;		
	}

	string getSAMformatForFinalAlignment_SE(PE_Read_Info& seReadInfo, bool FastaOrFastq)
	{
		string tmpReadName = seReadInfo.returnReadName_SE();

		string tmpSamStr;

		int IH_Nor = seAlignVec_final_Nor.size();
		int IH_Rcm = seAlignVec_final_Rcm.size();

		int IH_all = IH_Nor + IH_Rcm;
		for(int tmp = 0; tmp < seAlignVec_final_Nor.size(); tmp ++)
		{
			int HI_Nor_tmp = tmp + 1;
			int tmpNor_index = seAlignVec_final_Nor[tmp];
			Alignment_Info* tmpAlignInfo = norAlignmentInfo_PE_1[tmpNor_index];
			tmpSamStr = tmpSamStr + tmpAlignInfo->getSamFormatString_SE(tmpReadName, 
				seReadInfo.returnReadSeq_SE(), seReadInfo.returnQualitySeq_SE(), IH_all,
				HI_Nor_tmp, (tmp != 0), FastaOrFastq);
			tmpSamStr += "\n";
		}
		for(int tmp = 0; tmp < seAlignVec_final_Rcm.size(); tmp ++)
		{
			int HI_Rcm_tmp = tmp + 1 + IH_Nor;
			int tmpRcm_index = seAlignVec_final_Rcm[tmp];
			Alignment_Info* tmpAlignInfo = rcmAlignmentInfo_PE_1[tmpRcm_index];
			tmpSamStr = tmpSamStr + tmpAlignInfo->getSamFormatString_SE(tmpReadName, 
				seReadInfo.returnRcmReadSeq_SE(), seReadInfo.returnRcmQualitySeq_SE(), IH_all,
				HI_Rcm_tmp, ((tmp != 0)||(IH_Nor > 0)), FastaOrFastq);
			tmpSamStr += "\n";
		}
		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;
	}
	string getSAMformatForUnmapped_SE(PE_Read_Info& seReadInfo, bool FastaOrFastq)
	{
		string tmpQualSeq;
		if(FastaOrFastq)
		{
			tmpQualSeq = "*";
		}
		else
		{
			tmpQualSeq = seReadInfo.returnReadQual_SE();
		}
		string tmpReadName = seReadInfo.returnReadName_SE();

		string seAlignSamStr;
		seAlignSamStr = tmpReadName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + seReadInfo.returnReadSeq_SE()
			+ "\t" + tmpQualSeq + "\tIH:i:0\tHI:i:0";
		return seAlignSamStr;
	}

	int returnAlignmentScoreMinOutput_PE(PE_Read_Info& peReadInfo)
	{
		int peReadLength = peReadInfo.returnReadLength_end1() + peReadInfo.returnReadLength_end2();
		if(peReadLength < 100)
			return (int)(0.5 * peReadLength);
		else
			return 50;
	}

	int returnAlignmentScoreMinOutput_SE(PE_Read_Info& seReadInfo)
	{
		int seReadLength = seReadInfo.returnReadLength_SE();
		if(seReadLength < 50)
			return (int)(0.5 * seReadLength);
		else
			return 25;
	}

	bool alignPairScoreTooLow_bool(PE_Read_Info& peReadInfo)
	{
		//return false;
		int alignment_score_min_output = this->returnAlignmentScoreMinOutput_PE(peReadInfo);
		if(highestPairAlignmentScore < alignment_score_min_output)
			return true;
		else
			return false;
	}

	bool alignScoreTooLow_bool_SE(PE_Read_Info& seReadInfo)
	{
		//return false;
		int alignment_score_min_output = this->returnAlignmentScoreMinOutput_SE(seReadInfo);
		if(highestPairAlignmentScore < alignment_score_min_output)
			return true;
		else
			return false;
	}

	bool checkUniqueOrMulti_SE()
	{
		if(seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size() > 1)
			return false;
		else
			return true;
	}

	bool allFinalAlignmentComplete_SE()
	{
		for(int tmp_1 = 0; tmp_1 < seAlignVec_final_Nor.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				norAlignmentInfo_PE_1[seAlignVec_final_Nor[tmp_1]]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;
		}

		for(int tmp_1 = 0; tmp_1 < seAlignVec_final_Rcm.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				rcmAlignmentInfo_PE_1[seAlignVec_final_Rcm[tmp_1]]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;
		}
		return true;		
	}

	bool finalAlignExistsBool_SE()
	{
		if(seAlignVec_final_Nor.size() + seAlignVec_final_Rcm.size() > 0)
			return true;
		else
			return false;
	}

	int returnTranscriptIndex_finalPeAlign(
		bool Nor1Rcm2_or_Nor2Rcm1, int tmpIndex_peAlignInfo, 
		Transcript_Set* transcriptSetInfo)
	{
		//cout << "returnTranscriptIndex_finalPeAlign starts ...." << endl;
		string tmpTranscriptID;
		if(Nor1Rcm2_or_Nor2Rcm1)
		{
			int tmpNor1NO 
				= finalAlignPair_Nor1Rcm2[tmpIndex_peAlignInfo].first;
			tmpTranscriptID = norAlignmentInfo_PE_1[tmpNor1NO]->returnAlignChromName();
		}
		else
		{
			int tmpNor2NO
				= finalAlignPair_Nor2Rcm1[tmpIndex_peAlignInfo].first;
			tmpTranscriptID = norAlignmentInfo_PE_2[tmpNor2NO]->returnAlignChromName();
		}
		int index_transcriptSet
			= transcriptSetInfo->returnTranscriptSetIndex(tmpTranscriptID);
		return index_transcriptSet;
	}


	void addTranscriptReadCount(Transcript_Count_Vec* transcriptCountVecInfo,
		int index_countVec, Transcript_Set* transcriptSetInfo)
	{
		//cout << "addTranscriptReadCount starts ......" << endl;
		//cout << "index_countVec: " << index_countVec << endl;
		int finalPeAlignNum_Nor1Rcm2 
			= this->returnFinalAlignPairSize_Nor1Rcm2();
		int finalPeAlignNum_Nor2Rcm1
			= this->returnFinalAlignPairSize_Nor2Rcm1();
		int finalPeAlignNum = finalPeAlignNum_Nor1Rcm2
			+ finalPeAlignNum_Nor2Rcm1;
		//cout << "finalPeAlignNum_Nor1Rcm2: " << finalPeAlignNum_Nor1Rcm2 << endl;
		//cout << "finalPeAlignNum_Nor2Rcm1: " << finalPeAlignNum_Nor2Rcm1 << endl;
		double eachCountValue = 1.0;
		if(finalPeAlignNum > 1)
		{
			eachCountValue = 1.0/((double)finalPeAlignNum);
		}
		//cout << "start to process Nor1Rcm2" << endl;
		for(int tmpIndex_peAlignInfo = 0; 
			tmpIndex_peAlignInfo < finalPeAlignNum_Nor1Rcm2;
			tmpIndex_peAlignInfo ++)
		{
			int tmpIndex_transcriptVec = this->returnTranscriptIndex_finalPeAlign(
				true, tmpIndex_peAlignInfo, transcriptSetInfo);
			//cout << "tmpIndex_transcriptVec: " << tmpIndex_transcriptVec << endl;
			transcriptCountVecInfo->addNewReadCount_inOneCountInfo(
				tmpIndex_transcriptVec, eachCountValue, index_countVec);
		}
		//cout << "start to process Nor2Rcm1" << endl;
		for(int tmpIndex_peAlignInfo = 0; 
			tmpIndex_peAlignInfo < finalPeAlignNum_Nor2Rcm1;
			tmpIndex_peAlignInfo ++)
		{
			int tmpIndex_transcriptVec = this->returnTranscriptIndex_finalPeAlign(
				false, tmpIndex_peAlignInfo, transcriptSetInfo);
			//cout << "tmpIndex_transcriptVec: " << tmpIndex_transcriptVec << endl;
			transcriptCountVecInfo->addNewReadCount_inOneCountInfo(
				tmpIndex_transcriptVec, eachCountValue, index_countVec);
		}
		//cout << "end of addTranscriptReadCount" << endl;
	}

	void doTranscriptCount(Transcript_Count_Vec* transcriptCountVecInfo,
		int index_countVec, Transcript_Set* transcriptSetInfo)
	{
		//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		//cout << "doTranscriptCount starts ......" << endl;
		bool pairExistsBool = this->finalPairExistsBool();
		//cout << "pairExistsBool: " << pairExistsBool << endl;
		if(pairExistsBool)
		{
			bool unique_bool = this->checkUniqueOrMulti();
			//cout << "unique_bool: " << unique_bool << endl;
			if(unique_bool)
			{
				transcriptCountVecInfo->uniqueMappedReadCountIncrement_inOneCountInfo(
					index_countVec);
				this->addTranscriptReadCount(transcriptCountVecInfo,
					index_countVec, transcriptSetInfo);
			}
			else
			{
				transcriptCountVecInfo->multiMappedReadCountIncrement_inOneCountInfo(
					index_countVec);
				this->addTranscriptReadCount(transcriptCountVecInfo,
					index_countVec, transcriptSetInfo);
			}
		}
		else
		{
			transcriptCountVecInfo->unmappedReadCountIncrement_inOneCountInfo(
				index_countVec);
		}
	}

	void addGeneReadCount_uniqueMapped_multiMapped_readCountIncrement_inOneCountInfo(
		Gene_Count_Vec* geneCountVecInfo, int index_countVec, 
		Transcript_Set* transcriptSetInfo, Transcript2geneMap_Info* transcript2geneMapInfo)
	{
		int finalPeAlignNum_Nor1Rcm2 = this->returnFinalAlignPairSize_Nor1Rcm2();
		int finalPeAlignNum_Nor2Rcm1 = this->returnFinalAlignPairSize_Nor2Rcm1();
		int finalPeAlignNum = finalPeAlignNum_Nor1Rcm2 + finalPeAlignNum_Nor2Rcm1;

		set<int> geneIndexSet;
	
		for(int tmpIndex_peAlignInfo = 0; tmpIndex_peAlignInfo < finalPeAlignNum_Nor1Rcm2; 
			tmpIndex_peAlignInfo ++)
		{
			int tmpIndex_transcriptVec = this->returnTranscriptIndex_finalPeAlign(
				true, tmpIndex_peAlignInfo, transcriptSetInfo);		
			int tmpIndex_geneVec = transcript2geneMapInfo->getGeneIndexFromTranscriptIndex(
				tmpIndex_transcriptVec);
			geneIndexSet.insert(tmpIndex_geneVec);
		}
		for(int tmpIndex_peAlignInfo = 0; tmpIndex_peAlignInfo < finalPeAlignNum_Nor2Rcm1; 
			tmpIndex_peAlignInfo ++)
		{
			int tmpIndex_transcriptVec = this->returnTranscriptIndex_finalPeAlign(
				false, tmpIndex_peAlignInfo, transcriptSetInfo);
			int tmpIndex_geneVec = transcript2geneMapInfo->getGeneIndexFromTranscriptIndex(
				tmpIndex_transcriptVec);
			geneIndexSet.insert(tmpIndex_geneVec);
		}	

		int geneIndexNum = geneIndexSet.size();
		if(geneIndexNum == 0)
		{
			cout << "error in addGeneReadCount_uniqueMapped_multiMapped_readCountIncrement_inOneCountInfo( " << endl;
			exit(1);
		}
		else if(geneIndexNum == 1)
		{
			geneCountVecInfo->uniqueMappedReadCountIncrement_inOneCountInfo(index_countVec);
			int geneIndex = *(geneIndexSet.begin());
			geneCountVecInfo->addNewReadCount_inOneCountInfo(geneIndex, 1.0, index_countVec);
		}
		else
			geneCountVecInfo->multiMappedReadCountIncrement_inOneCountInfo(index_countVec);
	}

	void doGeneCount(Gene_Count_Vec* geneCountVecInfo, int index_countVec, 
		Transcript_Set* transcriptSetInfo, Transcript2geneMap_Info* transcript2geneMapInfo)
	{
		bool pairExistsBool = this->finalPairExistsBool();
		if(pairExistsBool)
			this->addGeneReadCount_uniqueMapped_multiMapped_readCountIncrement_inOneCountInfo(
				geneCountVecInfo, index_countVec, transcriptSetInfo, transcript2geneMapInfo);
		else
			geneCountVecInfo->unmappedReadCountIncrement_inOneCountInfo(index_countVec);
	}

	void output_phase1(bool outputDirectlyBool_Phase1Only,
		bool Do_Phase1_Only,
		vector<string>& PeAlignSamStrVec_complete,
		vector<string>& PeAlignInfoStrVec_inCompletePair,
		vector<string>& PeAlignInfoStrVec_oneEndUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_lowScore,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
		vector<string>& PeAlignSamStrVec_inCompletePair,
		//vector<string>& PeAlignInfoStrVec_completePaired,
		vector<RepeatRegion_Info*>& repeatRegionInfoVec,
		PE_Read_Info& readInfo, Index_Info* indexInfo, int tmpOpenMP, int threadNO,
		Stats_Info* statsInfo, bool FastaOrFastq, int multiMapSeg_maxLength
		)
	{
		//cout << "start output_phase1 ..." << endl;
		bool pairExistsBool = this->finalPairExistsBool();
		bool mappedToRepeatRegionBool = this->mappedToRepeatRegionBool();

		if(pairExistsBool) // some pair exists
		{
			bool allAlignmentCompleteBool = this->allAlignmentInFinalPairCompleted();
			bool unique_bool = this->checkUniqueOrMulti();

			statsInfo->increPairedNum_phase1(threadNO, allAlignmentCompleteBool, unique_bool);
				
			if(allAlignmentCompleteBool)
			{
				bool completeAlignmentPairScore_tooLow
				 	= this->alignPairScoreTooLow_bool(readInfo); 						

				if(completeAlignmentPairScore_tooLow)
				{
					statsInfo->increLowScoreComplete_phase1(threadNO, unique_bool);
					PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmpOpenMP] = 
						this->getSAMformatForBothEndsUnmapped(readInfo, FastaOrFastq);	
				}
				else
				{
					PeAlignSamStrVec_complete[tmpOpenMP] = 
						this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq, multiMapSeg_maxLength);
				}
			}
			else
			{
				if(!Do_Phase1_Only)
				{
					PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = 
						this->getTmpAlignInfoForFinalPair(
				 		readInfo.returnReadName_1(), readInfo.returnReadName_2(), 
						readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
						readInfo.returnReadQual_1(), readInfo.returnReadQual_2(),
						FastaOrFastq, multiMapSeg_maxLength);
				}
				PeAlignSamStrVec_inCompletePair[tmpOpenMP] 
					= this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq, multiMapSeg_maxLength);
			}
		}
		else // no pair exists: 1.one end unmapped; 2. both ends unmapped
		{
			bool alignmentExistsBool = this->alignInfoExistsBool();
			statsInfo->increUnpairedNum_phase1(threadNO, alignmentExistsBool, mappedToRepeatRegionBool);
			if(alignmentExistsBool) // one end unmapped
			{	
				if(Do_Phase1_Only)
				{
					PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = 
						this->getSAMformatForUnpairedAlignments_secondaryOrNot(
							readInfo, FastaOrFastq);
				}
				else
				{
					PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = 
						this->getTmpAlignInfo(
				 		readInfo.returnReadName_1(), readInfo.returnReadName_2(), 
						readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
						readInfo.returnReadQual_1(), readInfo.returnReadQual_2(),
						FastaOrFastq, multiMapSeg_maxLength);					
				}
			}
			else if(mappedToRepeatRegionBool)
			{
				PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmpOpenMP] = 
					this->getSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
						readInfo, repeatRegionInfoVec[threadNO], indexInfo,
						FastaOrFastq);
			}
			else // both ends unmapped
			{
				PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = 
					this->getSAMformatForBothEndsUnmapped(
						readInfo, FastaOrFastq);						
			}
		}
	}

	void output_phase1(bool outputDirectlyBool_Phase1Only,
		bool Do_Phase1_Only,
		vector<string>& PeAlignSamStrVec_complete,
		vector<string>& PeAlignInfoStrVec_inCompletePair,
		vector<string>& PeAlignInfoStrVec_oneEndUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_lowScore,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
		vector<string>& PeAlignSamStrVec_inCompletePair,
		//vector<string>& PeAlignInfoStrVec_completePaired,
		vector<RepeatRegion_Info*>& repeatRegionInfoVec,
		PE_Read_Info& readInfo, Index_Info* indexInfo, int tmpOpenMP, int threadNO,
		Stats_Info* statsInfo, bool FastaOrFastq)
	{
		int tmp_multiMapSeg_maxLength = 0;
		this->output_phase1(
			outputDirectlyBool_Phase1Only, Do_Phase1_Only,
			PeAlignSamStrVec_complete,
			PeAlignInfoStrVec_inCompletePair,
			PeAlignInfoStrVec_oneEndUnmapped,
			PeAlignSamStrVec_bothEndsUnmapped,
			PeAlignSamStrVec_bothEndsUnmapped_lowScore,
			PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
			PeAlignSamStrVec_inCompletePair,
			repeatRegionInfoVec, readInfo, indexInfo, tmpOpenMP, threadNO,
			statsInfo, FastaOrFastq, tmp_multiMapSeg_maxLength);
	}

	void output_phase1(bool outputDirectlyBool_Phase1Only, bool Do_Phase1_Only, Result_Array* resultArray, 
		vector<RepeatRegion_Info*>& repeatRegionInfoVec, PE_Read_Info& readInfo, Index_Info* indexInfo, 
		int tmpOpenMP, int threadNO, Stats_Info* statsInfo, bool FastaOrFastq)
	{
		this->output_phase1(outputDirectlyBool_Phase1Only, Do_Phase1_Only, resultArray, repeatRegionInfoVec, 
			readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo, FastaOrFastq, 0);
	}

	void output_phase1(bool outputDirectlyBool_Phase1Only, bool Do_Phase1_Only, Result_Array* resultArray, 
		vector<RepeatRegion_Info*>& repeatRegionInfoVec, PE_Read_Info& readInfo, Index_Info* indexInfo, 
		int tmpOpenMP, int threadNO, Stats_Info* statsInfo, bool FastaOrFastq, int multiMapSeg_maxLength)
	{
		//cout << "start output_phase1 ..." << endl;
		bool pairExistsBool = this->finalPairExistsBool();
		bool mappedToRepeatRegionBool = this->mappedToRepeatRegionBool();
		if(pairExistsBool) // some pair exists
		{
			bool allAlignmentCompleteBool = this->allAlignmentInFinalPairCompleted();
			bool unique_bool = this->checkUniqueOrMulti();		
			//cout << "allAlignmentCompleteBool: " << allAlignmentCompleteBool << endl;
			//cout << "unique_bool: " << unique_bool << endl; 
			statsInfo->increPairedNum_phase1(threadNO, allAlignmentCompleteBool, unique_bool);	
			if(allAlignmentCompleteBool)
			{
				bool completeAlignmentPairScore_tooLow = this->alignPairScoreTooLow_bool(readInfo);
				if(completeAlignmentPairScore_tooLow)
				{
					statsInfo->increLowScoreComplete_phase1(threadNO, unique_bool);
					resultArray->insert_PeAlignSamStrVec_bothEndsUnmapped_lowScore(
					this->getSAMformatForBothEndsUnmapped(readInfo, FastaOrFastq), tmpOpenMP);
				}
				else
					resultArray->insert_PeAlignSamStrVec_complete(this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq), tmpOpenMP);
			}
			else
			{
				resultArray->insert_PeAlignInfoStrVec_inCompletePair(this->getTmpAlignInfoForFinalPair(
				 		readInfo.returnReadName_1(), readInfo.returnReadName_2(), readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
						readInfo.returnReadQual_1(), readInfo.returnReadQual_2(), FastaOrFastq), tmpOpenMP);
				resultArray->insert_PeAlignSamStrVec_inCompletePair( 
							this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq), tmpOpenMP);
			}
		}
		else // no pair exists: 1.one end unmapped; 2. both ends unmapped
		{
			bool alignmentExistsBool = this->alignInfoExistsBool();
			statsInfo->increUnpairedNum_phase1(threadNO, alignmentExistsBool, mappedToRepeatRegionBool);
			if(alignmentExistsBool) // one end unmapped
				resultArray->insert_PeAlignInfoStrVec_oneEndUnmapped(this->getTmpAlignInfo(readInfo.returnReadName_1(), 
					readInfo.returnReadName_2(), readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(), 
					readInfo.returnReadQual_1(), readInfo.returnReadQual_2(), FastaOrFastq), tmpOpenMP);
			else if(mappedToRepeatRegionBool)
				resultArray->insert_PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion( this->getSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
							readInfo, repeatRegionInfoVec[threadNO], indexInfo, FastaOrFastq), tmpOpenMP);
			else // both ends unmapped
				resultArray->insert_PeAlignSamStrVec_bothEndsUnmapped(this->getSAMformatForBothEndsUnmapped(readInfo, FastaOrFastq), tmpOpenMP);						
		}
	}

	void output_phase2_fixOneEndUnpaired()
	{}
	
	void output_phase2_fixHeadTail()
	{}

	int checkUniqueOrMulti()
	{
		int pairSize = (finalAlignPair_Nor1Rcm2.size() + finalAlignPair_Nor2Rcm1.size());
		if(pairSize > 1)
		{
			return false;
		}
		else
			return true;
	}

	Alignment_Info* returnAlignInfo_Nor1(int tmp)
	{
		return norAlignmentInfo_PE_1[tmp];
	}
	Alignment_Info* returnAlignInfo_Rcm1(int tmp)
	{
		return rcmAlignmentInfo_PE_1[tmp];
	}
	Alignment_Info* returnAlignInfo_Nor2(int tmp)
	{
		return norAlignmentInfo_PE_2[tmp];
	}
	Alignment_Info* returnAlignInfo_Rcm2(int tmp)
	{
		return rcmAlignmentInfo_PE_2[tmp];
	}

	int returnOriAlignPair_Nor1Rcm2_size()
	{
		return oriAlignPair_Nor1Rcm2.size();
	}
	int returnOriAlignPair_Nor2Rcm1_size()
	{
		return oriAlignPair_Nor2Rcm1.size();
	}

	int returnNorAlignmentInfo_PE_1_size()
	{
		return norAlignmentInfo_PE_1.size();
	}
	int returnRcmAlignmentInfo_PE_1_size()
	{
		return rcmAlignmentInfo_PE_1.size();
	}
	int returnNorAlignmentInfo_PE_2_size()
	{
		return norAlignmentInfo_PE_2.size();
	}
	int returnRcmAlignmentInfo_PE_2_size()
	{
		return rcmAlignmentInfo_PE_2.size();
	}	
	int getAlignInfoVecSize(int alignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		//cout << "getAlignInfoVecSize starts ... " << endl;
		//cout << "alignInfoKind: " << alignInfoKind << endl;
		int size = 100;
		if(alignInfoKind == 1)
		{
			size = norAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 2)
		{
			size = rcmAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 3)
		{
			size = norAlignmentInfo_PE_2.size();
		}
		else if(alignInfoKind == 4)
		{
			size = rcmAlignmentInfo_PE_2.size();
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfoVecSize()" << endl;
		}		
		//cout << "size: " << size << endl;
		return size;
	}

	int returnAllAlignmentInfoSize()
	{
		return (norAlignmentInfo_PE_1.size()
			+ rcmAlignmentInfo_PE_1.size()
			+ norAlignmentInfo_PE_2.size()
			+ rcmAlignmentInfo_PE_2.size());
	}

	bool mappedToRepeatRegionBool()
	{
		if( (repeatRegion_index_Nor1 > 0) || (repeatRegion_index_Rcm1 > 0) 
			|| (repeatRegion_index_Nor2 > 0) || (repeatRegion_index_Rcm2 > 0) )
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool mappedToRepeatRegionBool_SE()
	{
		if( (repeatRegion_index_Nor1 > 0) || (repeatRegion_index_Rcm1 > 0) )
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	Alignment_Info* returnAlignInfoInPeAlignInfo(bool Nor1Rcm2OrNor2Rcm1, bool NorOrRcm, int index)
	{
		if(Nor1Rcm2OrNor2Rcm1)
		{
			if(NorOrRcm)
			{
				return norAlignmentInfo_PE_1[index];
			}
			else
			{
				return rcmAlignmentInfo_PE_2[index];
			}
		}
		else
		{	
			if(NorOrRcm)
			{
				return norAlignmentInfo_PE_2[index];
			}
			else
			{
				return rcmAlignmentInfo_PE_1[index];
			}
		}
	}

	Alignment_Info* returnAlignInfoInPeAlignInfo(int tmpAlignInfoType, int tmpIndex_peAlignInfo)
	{
		if(tmpAlignInfoType == 1)
			return norAlignmentInfo_PE_1[tmpIndex_peAlignInfo];
		else if(tmpAlignInfoType == 2)
			return rcmAlignmentInfo_PE_1[tmpIndex_peAlignInfo];
		else if(tmpAlignInfoType == 3)
			return norAlignmentInfo_PE_2[tmpIndex_peAlignInfo];
		else
			return rcmAlignmentInfo_PE_2[tmpIndex_peAlignInfo];
	}

	void oriAlignPair2oriAlignPairNew()
	{
		for(int tmp1 = 0; tmp1 < oriAlignPair_Nor1Rcm2.size(); tmp1++)
		{
			int tmpPair_NO_1 =  (oriAlignPair_Nor1Rcm2[tmp1].first);
			int tmpSize = (oriAlignPair_Nor1Rcm2[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
			{
				int tmpPair_NO_2 = (oriAlignPair_Nor1Rcm2[tmp1].second)[tmp2];

				string tmpChrName = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromName();
				int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromPos();
				int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnEndMatchedPosInChr();
				int tmpMappedLength_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->mappedLength();
				int tmpMismatchNum_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnMismatchNum();
				norAlignmentInfo_PE_1[tmpPair_NO_1]->generateIndelLength();
				int tmpInsertionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnInsertionLength();
				int tmpDeletionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnDeletionLength();

				int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnAlignChromPos();
				int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnEndMatchedPosInChr();
				int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->mappedLength();
				int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnMismatchNum();
				rcmAlignmentInfo_PE_2[tmpPair_NO_2]->generateIndelLength();
				int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnInsertionLength();
				int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnDeletionLength();

				int tmpChrMapPosStart, tmpChrMapPosEnd;
				#ifdef DETECT_CIRCULAR_RNA
				tmpChrMapPosStart = selectTheSmallestAmong4values(
						tmpChrMapPosStart_nor, tmpChrMapPosStart_rcm, tmpChrMapPosEnd_nor, tmpChrMapPosEnd_rcm);
				tmpChrMapPosEnd = selectTheLargestAmong4values(
						tmpChrMapPosStart_nor, tmpChrMapPosStart_rcm, tmpChrMapPosEnd_nor, tmpChrMapPosEnd_rcm);
				#else
				if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
					tmpChrMapPosStart = tmpChrMapPosStart_nor;
				else
					tmpChrMapPosStart = tmpChrMapPosStart_rcm;

				if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
					tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
				else
					tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
				#endif

				int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
				int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
				int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
				int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
				int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

				oriAlignPair_Nor1Rcm2_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
				oriAlignPairNew_chrNameVec_Nor1Rcm2.push_back(tmpChrName);
				oriAlignPairNew_startMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosStart);
				oriAlignPairNew_endMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosEnd);
				oriAlignPairNew_mappedLength_Nor1Rcm2.push_back(tmpMappedLength);
				oriAlignPairNew_mismatchNum_Nor1Rcm2.push_back(tmpMismatchNum);
				oriAlignPairNew_pairDistance_Nor1Rcm2.push_back(tmpPairDistance);
				oriAlignPairNew_insertionLength_Nor1Rcm2.push_back(tmpInsertionLength);
				oriAlignPairNew_deletionLength_Nor1Rcm2.push_back(tmpDeletionLength);

				int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
					(norAlignmentInfo_PE_1[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
					(rcmAlignmentInfo_PE_2[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
				oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2.push_back(SJconfidenceLevel);
			}
		}

		for(int tmp1 = 0; tmp1 < oriAlignPair_Nor2Rcm1.size(); tmp1++)
		{
			int tmpPair_NO_1 = (oriAlignPair_Nor2Rcm1[tmp1]).first;
			int tmpSize = (oriAlignPair_Nor2Rcm1[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
			{
				int tmpPair_NO_2 = (oriAlignPair_Nor2Rcm1[tmp1].second)[tmp2];

				string tmpChrName = (norAlignmentInfo_PE_2[tmpPair_NO_1])->returnAlignChromName();
				int tmpChrMapPosStart_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnAlignChromPos();
				int tmpChrMapPosEnd_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnEndMatchedPosInChr();
				int tmpMappedLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->mappedLength();
				int tmpMismatchNum_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnMismatchNum();
				norAlignmentInfo_PE_2[tmpPair_NO_1]->generateIndelLength();
				int tmpInsertionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnInsertionLength();
				int tmpDeletionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnDeletionLength();

				int tmpChrMapPosStart_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnAlignChromPos();
				int tmpChrMapPosEnd_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnEndMatchedPosInChr();
				int tmpMappedLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->mappedLength();
				int tmpMismatchNum_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnMismatchNum();
				rcmAlignmentInfo_PE_1[tmpPair_NO_2]->generateIndelLength();
				int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnInsertionLength();
				int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnDeletionLength();

				int tmpChrMapPosStart, tmpChrMapPosEnd;
				#ifdef DETECT_CIRCULAR_RNA
				tmpChrMapPosStart = selectTheSmallestAmong4values(
						tmpChrMapPosStart_nor, tmpChrMapPosStart_rcm, tmpChrMapPosEnd_nor, tmpChrMapPosEnd_rcm);
				tmpChrMapPosEnd = selectTheLargestAmong4values(
						tmpChrMapPosStart_nor, tmpChrMapPosStart_rcm, tmpChrMapPosEnd_nor, tmpChrMapPosEnd_rcm);
				#else
				if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
					tmpChrMapPosStart = tmpChrMapPosStart_nor;
				else
					tmpChrMapPosStart = tmpChrMapPosStart_rcm;

				if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
					tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
				else
					tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
				#endif

				int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
				int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
				int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
				int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
				int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

				oriAlignPair_Nor2Rcm1_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
				oriAlignPairNew_chrNameVec_Nor2Rcm1.push_back(tmpChrName);
				oriAlignPairNew_startMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosStart);
				oriAlignPairNew_endMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosEnd);
				oriAlignPairNew_mappedLength_Nor2Rcm1.push_back(tmpMappedLength);
				oriAlignPairNew_mismatchNum_Nor2Rcm1.push_back(tmpMismatchNum);
				oriAlignPairNew_pairDistance_Nor2Rcm1.push_back(tmpPairDistance);
				oriAlignPairNew_insertionLength_Nor2Rcm1.push_back(tmpInsertionLength);
				oriAlignPairNew_deletionLength_Nor2Rcm1.push_back(tmpDeletionLength);

				int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
					(norAlignmentInfo_PE_2[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
					(rcmAlignmentInfo_PE_1[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
				oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1.push_back(SJconfidenceLevel);
			}
		}
	}

	void oriAlignPair2oriAlignPairNew_Nor1Rcm2Only()
	{
		for(int tmp1 = 0; tmp1 < oriAlignPair_Nor1Rcm2.size(); tmp1++)
		{
			int tmpPair_NO_1 =  (oriAlignPair_Nor1Rcm2[tmp1].first);
			int tmpSize = (oriAlignPair_Nor1Rcm2[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
			{
				int tmpPair_NO_2 = (oriAlignPair_Nor1Rcm2[tmp1].second)[tmp2];

				string tmpChrName = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromName();
				int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromPos();
				int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnEndMatchedPosInChr();
				int tmpMappedLength_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->mappedLength();
				int tmpMismatchNum_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnMismatchNum();
				norAlignmentInfo_PE_1[tmpPair_NO_1]->generateIndelLength();
				int tmpInsertionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnInsertionLength();
				int tmpDeletionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnDeletionLength();

				int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnAlignChromPos();
				int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnEndMatchedPosInChr();
				int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->mappedLength();
				int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnMismatchNum();
				rcmAlignmentInfo_PE_2[tmpPair_NO_2]->generateIndelLength();
				int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnInsertionLength();
				int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnDeletionLength();

				int tmpChrMapPosStart, tmpChrMapPosEnd;
				if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
				{
					tmpChrMapPosStart = tmpChrMapPosStart_nor;
				}
				else
				{
					tmpChrMapPosStart = tmpChrMapPosStart_rcm;
				}

				if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
				}
				else
				{
					tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
				}
				int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
				int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
				int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
				int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
				int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

				oriAlignPair_Nor1Rcm2_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
				oriAlignPairNew_chrNameVec_Nor1Rcm2.push_back(tmpChrName);
				oriAlignPairNew_startMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosStart);
				oriAlignPairNew_endMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosEnd);
				oriAlignPairNew_mappedLength_Nor1Rcm2.push_back(tmpMappedLength);
				oriAlignPairNew_mismatchNum_Nor1Rcm2.push_back(tmpMismatchNum);
				oriAlignPairNew_pairDistance_Nor1Rcm2.push_back(tmpPairDistance);
				oriAlignPairNew_insertionLength_Nor1Rcm2.push_back(tmpInsertionLength);
				oriAlignPairNew_deletionLength_Nor1Rcm2.push_back(tmpDeletionLength);

				int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
					(norAlignmentInfo_PE_1[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
					(rcmAlignmentInfo_PE_2[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
				oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2.push_back(SJconfidenceLevel);
			}
		}

		// for(int tmp1 = 0; tmp1 < oriAlignPair_Nor2Rcm1.size(); tmp1++)
		// {
		// 	int tmpPair_NO_1 = (oriAlignPair_Nor2Rcm1[tmp1]).first;
		// 	int tmpSize = (oriAlignPair_Nor2Rcm1[tmp1].second).size();
		// 	for(int tmp2 = 0; tmp2 < tmpSize; tmp2++)
		// 	{
		// 		int tmpPair_NO_2 = (oriAlignPair_Nor2Rcm1[tmp1].second)[tmp2];

		// 		string tmpChrName = (norAlignmentInfo_PE_2[tmpPair_NO_1])->returnAlignChromName();
		// 		int tmpChrMapPosStart_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnAlignChromPos();
		// 		int tmpChrMapPosEnd_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnEndMatchedPosInChr();
		// 		int tmpMappedLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->mappedLength();
		// 		int tmpMismatchNum_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnMismatchNum();
		// 		norAlignmentInfo_PE_2[tmpPair_NO_1]->generateIndelLength();
		// 		int tmpInsertionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnInsertionLength();
		// 		int tmpDeletionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnDeletionLength();

		// 		int tmpChrMapPosStart_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnAlignChromPos();
		// 		int tmpChrMapPosEnd_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnEndMatchedPosInChr();
		// 		int tmpMappedLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->mappedLength();
		// 		int tmpMismatchNum_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnMismatchNum();
		// 		rcmAlignmentInfo_PE_1[tmpPair_NO_2]->generateIndelLength();
		// 		int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnInsertionLength();
		// 		int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnDeletionLength();

		// 		int tmpChrMapPosStart, tmpChrMapPosEnd;
		// 		if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
		// 		{
		// 			tmpChrMapPosStart = tmpChrMapPosStart_nor;
		// 		}
		// 		else
		// 		{
		// 			tmpChrMapPosStart = tmpChrMapPosStart_rcm;
		// 		}

		// 		if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
		// 		{
		// 			tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
		// 		}
		// 		else
		// 		{
		// 			tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
		// 		}
		// 		int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
		// 		int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
		// 		int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
		// 		int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
		// 		int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

		// 		oriAlignPair_Nor2Rcm1_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
		// 		oriAlignPairNew_chrNameVec_Nor2Rcm1.push_back(tmpChrName);
		// 		oriAlignPairNew_startMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosStart);
		// 		oriAlignPairNew_endMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosEnd);
		// 		oriAlignPairNew_mappedLength_Nor2Rcm1.push_back(tmpMappedLength);
		// 		oriAlignPairNew_mismatchNum_Nor2Rcm1.push_back(tmpMismatchNum);
		// 		oriAlignPairNew_pairDistance_Nor2Rcm1.push_back(tmpPairDistance);
		// 		oriAlignPairNew_insertionLength_Nor2Rcm1.push_back(tmpInsertionLength);
		// 		oriAlignPairNew_deletionLength_Nor2Rcm1.push_back(tmpDeletionLength);

		// 		int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
		// 			(norAlignmentInfo_PE_2[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
		// 			(rcmAlignmentInfo_PE_1[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
		// 		oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1.push_back(SJconfidenceLevel);
		// 	}
		// }
	}

	int getPairedAlignmentSJconfidenceLevel(int confidenceLevel_1, int confidenceLevel_2)
	{
		if(confidenceLevel_1 > confidenceLevel_2)
		{
			return confidenceLevel_1;
		}
		else 
			return confidenceLevel_2;
	}

	/*
	void oriAlignPairGroupedByRegion()
	{
		vector < string > oriPairRegionChrVec_Nor1Rcm2;
		vector< pair <int, int> > oriPairRegionPosVec_Nor1Rcm2;

		vector < string > oriPairRegionChrVec_Nor2Rcm1;
		vector< pair <int, int> > oriPairRegionPosVec_Nor2Rcm1;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			int tmpPairStartMapPos = oriAlignPairNew_startMapPosVec_Nor1Rcm2[tmp];
			int tmpPairEndMapPos = oriAlignPairNew_endMapPosVec_Nor1Rcm2[tmp];
			string tmpPairChrName = oriAlignPairNew_chrNameVec_Nor1Rcm2[tmp];

			int regionVecSize = oriPairRegionPosVec_Nor1Rcm2.size();
			bool overlapRegionFound = false;
			for(int tmp2 = 0; tmp2 < regionVecSize; tmp2++)
			{
				string tmpRegionChr = oriPairRegionChrVec_Nor1Rcm2[tmp2];
				int tmpRegionStartPos = oriPairRegionPosVec_Nor1Rcm2[tmp2].first;
				int tmpRegionEndPos = oriPairRegionPosVec_Nor1Rcm2[tmp2].second;
				if ( ( tmpPairChrName == tmpRegionChr ) &&
					( !( (tmpPairStartMapPos >= tmpRegionEndPos) || ( tmpPairEndMapPos <= tmpRegionStartPos ) ) ) )
				{
					overlapRegionFound = true;					
					if(tmpPairStartMapPos < tmpRegionStartPos)
						oriPairRegionPosVec_Nor1Rcm2[tmp2].first = tmpPairStartMapPos;
					if(tmpPairEndMapPos > tmpRegionEndPos)
						oriPairRegionPosVec_Nor1Rcm2[tmp2].second = tmpPairEndMapPos;
					oriAlignPairGroupedByRegion_Nor1Rcm2[tmp2].push_back(tmp);
					break;
				}
			}

			if(!overlapRegionFound)
			{
				vector<int> newIntVec;
				newIntVec.push_back(tmp);
				oriAlignPairGroupedByRegion_Nor1Rcm2.push_back(newIntVec);
				oriPairRegionChrVec_Nor1Rcm2.push_back(tmpPairChrName);
				oriPairRegionPosVec_Nor1Rcm2.push_back(pair<int,int> (tmpPairStartMapPos, tmpPairEndMapPos));
			}
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			int tmpPairStartMapPos = oriAlignPairNew_startMapPosVec_Nor2Rcm1[tmp];
			int tmpPairEndMapPos = oriAlignPairNew_endMapPosVec_Nor2Rcm1[tmp];
			string tmpPairChrName = oriAlignPairNew_chrNameVec_Nor2Rcm1[tmp];

			int regionVecSize = oriPairRegionPosVec_Nor2Rcm1.size();
			bool overlapRegionFound = false;
			for(int tmp2 = 0; tmp2 < regionVecSize; tmp2++)
			{
				string tmpRegionChr = oriPairRegionChrVec_Nor2Rcm1[tmp2];
				int tmpRegionStartPos = oriPairRegionPosVec_Nor2Rcm1[tmp2].first;
				int tmpRegionEndPos = oriPairRegionPosVec_Nor2Rcm1[tmp2].second;
				if ( ( tmpPairChrName == tmpRegionChr ) &&
					( !( (tmpPairStartMapPos >= tmpRegionEndPos) || ( tmpPairEndMapPos <= tmpRegionStartPos ) ) ) )
				{
					overlapRegionFound = true;					
					if(tmpPairStartMapPos < tmpRegionStartPos)
						oriPairRegionPosVec_Nor2Rcm1[tmp2].first = tmpPairStartMapPos;
					if(tmpPairEndMapPos > tmpRegionEndPos)
						oriPairRegionPosVec_Nor2Rcm1[tmp2].second = tmpPairEndMapPos;
					oriAlignPairGroupedByRegion_Nor2Rcm1[tmp2].push_back(tmp);
					break;
				}
			}

			if(!overlapRegionFound)
			{
				vector<int> newIntVec;
				newIntVec.push_back(tmp);
				oriAlignPairGroupedByRegion_Nor2Rcm1.push_back(newIntVec);
				oriPairRegionChrVec_Nor2Rcm1.push_back(tmpPairChrName);
				oriPairRegionPosVec_Nor2Rcm1.push_back(pair<int,int> (tmpPairStartMapPos, tmpPairEndMapPos));
			}
		}
	}*/

	int getSJpenalty(bool Nor1Rcm2_or_Nor2Rcm1_bool, int indexInOriAlignPair, 
		int semiCanonicalSJ_penalty, int nonCanonicalSJ_penalty)
	{
		int confidenceLevel;
		int SJ_penalty = 0;
		if(Nor1Rcm2_or_Nor2Rcm1_bool)
		{
			confidenceLevel = oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2[indexInOriAlignPair];
		}
		else
		{
			confidenceLevel = oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1[indexInOriAlignPair]; 
		}

		if(confidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
		{}
		else if(confidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
			SJ_penalty = semiCanonicalSJ_penalty;
		else  // noncanonical
			SJ_penalty = nonCanonicalSJ_penalty;

		return SJ_penalty;
	}

	int getMappingAreaSizePenalty(bool Nor1Rcm2_or_Nor2Rcm1_bool, int indexInOriAlignPair, 
		int distant_insert_size_min, int distant_insert_size_penalty)
	{
		int map_start_pos;
		int map_end_pos;
		if(Nor1Rcm2_or_Nor2Rcm1_bool)
		{	
			map_start_pos = oriAlignPairNew_startMapPosVec_Nor1Rcm2[indexInOriAlignPair];
			map_end_pos = oriAlignPairNew_endMapPosVec_Nor1Rcm2[indexInOriAlignPair];
		}
		else
		{	
			map_start_pos = oriAlignPairNew_startMapPosVec_Nor2Rcm1[indexInOriAlignPair];
			map_end_pos = oriAlignPairNew_endMapPosVec_Nor2Rcm1[indexInOriAlignPair];
		}
		if((map_end_pos - map_start_pos) >= distant_insert_size_min)
			return distant_insert_size_penalty;
		else
			return 0;
	}

	void getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty()
	{
		int mismatch_weight = MISMATCH_PENALTY;
		int deletion_weight = DELETION_PENALTY;
		int insertion_weight = INSERTION_PENALTY;

		int semiCanonicalSJ_penalty = SEMICANONICALSJ_PENALTY;
		int nonCanonicalSJ_penalty = NONCANONICALSJ_PENALTY;
	
		int distant_insert_size_min = DISTANT_INSERT_SIZE_MIN;
		int distant_insert_size_penalty = DISTANT_INSERT_SIZE_PENALTY;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor1Rcm2[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor1Rcm2[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor1Rcm2[tmp];

			double tmpPairDistance_score;
			#ifdef DETECT_CIRCULAR_RNA
			if((tmpPairDistance >= 500)||(tmpPairDistance <= -500))
				tmpPairDistance_score = 0.1;
			else
				tmpPairDistance_score = 0;
			#else
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			#endif

			int SJ_penalty = this->getSJpenalty(true, tmp, semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);

			int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
				true, tmp, distant_insert_size_min, distant_insert_size_penalty);

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score - SJ_penalty
				- mappingAreaSizePenalty;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor2Rcm1[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor2Rcm1[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor2Rcm1[tmp];

			double tmpPairDistance_score;
			#ifdef DETECT_CIRCULAR_RNA
			if((tmpPairDistance >= 500)||(tmpPairDistance <= -500))
				tmpPairDistance_score = 0.1;
			else
				tmpPairDistance_score = 0;
			#else
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			#endif

			int SJ_penalty = this->getSJpenalty(false, tmp, semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);

			int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
				false, tmp, distant_insert_size_min, distant_insert_size_penalty);

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score - SJ_penalty
				- mappingAreaSizePenalty;

			oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		}
	}

	void getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty_Nor1Rcm2()
	{
		int mismatch_weight = MISMATCH_PENALTY;
		int deletion_weight = DELETION_PENALTY;
		int insertion_weight = INSERTION_PENALTY;

		int semiCanonicalSJ_penalty = SEMICANONICALSJ_PENALTY;
		int nonCanonicalSJ_penalty = NONCANONICALSJ_PENALTY;
	
		int distant_insert_size_min = DISTANT_INSERT_SIZE_MIN;
		int distant_insert_size_penalty = DISTANT_INSERT_SIZE_PENALTY;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor1Rcm2[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor1Rcm2[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor1Rcm2[tmp];

			double tmpPairDistance_score;
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			
			int SJ_penalty = this->getSJpenalty(true, tmp, semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);

			int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
				true, tmp, distant_insert_size_min, distant_insert_size_penalty);

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score - SJ_penalty
				- mappingAreaSizePenalty;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		// for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		// {
		// 	int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
		// 	int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
		// 	int tmpPairDistance = oriAlignPairNew_pairDistance_Nor2Rcm1[tmp];
		// 	int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor2Rcm1[tmp];
		// 	int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor2Rcm1[tmp];

		// 	double tmpPairDistance_score;
		// 	if(tmpPairDistance > 500)
		// 		tmpPairDistance_score = 0.1;
		// 	else if(tmpPairDistance <= 0)
		// 		tmpPairDistance_score = 0;
		// 	else
		// 		tmpPairDistance_score = (double)tmpPairDistance/5000;

		// 	int SJ_penalty = this->getSJpenalty(false, tmp, semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);

		// 	int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
		// 		false, tmp, distant_insert_size_min, distant_insert_size_penalty);

		// 	double tmpScore = tmpMappedLength 
		//         - tmpMismatchNum * mismatch_weight 
		// 		- tmpInsertionLength * insertion_weight 
		// 		- tmpDeletionLength * deletion_weight
		// 		- tmpPairDistance_score - SJ_penalty
		// 		- mappingAreaSizePenalty;

		// 	oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		// }
	}
	/*
	void getScoreForEachPair_mapLen_mis_peDis()
	{
		int mismatch_weight = 2;
		int deletion_weight = 3;
		int insertion_weight = 4;
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor1Rcm2[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor1Rcm2[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor1Rcm2[tmp];

			double tmpPairDistance_score;
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			
			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
			int tmpPairDistance = oriAlignPairNew_pairDistance_Nor2Rcm1[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor2Rcm1[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor2Rcm1[tmp];

			double tmpPairDistance_score;
			if(tmpPairDistance > 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score;

			oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		}
	}

	void getScoreForEachPair_mapLen_mis()
	{
		int mismatch_weight = 2;
		int deletion_weight = 3;
		int insertion_weight = 4;
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor1Rcm2[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor1Rcm2[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor1Rcm2[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor1Rcm2[tmp];
			//int tmpPairDistance = oriAlignPairNew_pairDistance_Nor1Rcm2[tmp];
			//if(tmpPairDistance > 500)
			//	tmpPairDistance_score = 0.1;
			//else
			//	tmpPairDistance_score = 0;
			double tmpScore = tmpMappedLength 		        
			    - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight;

			oriAlignPairNew_score_Nor1Rcm2.push_back(tmpScore);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			int tmpMappedLength = oriAlignPairNew_mappedLength_Nor2Rcm1[tmp];
			int tmpMismatchNum = oriAlignPairNew_mismatchNum_Nor2Rcm1[tmp];
			int tmpInsertionLength = oriAlignPairNew_insertionLength_Nor2Rcm1[tmp];
			int tmpDeletionLength = oriAlignPairNew_deletionLength_Nor2Rcm1[tmp];
			//int tmpPairDistance = oriAlignPairNew_pairDistance_Nor2Rcm1[tmp];
			//if(tmpPairDistance > 500)
			//	tmpPairDistance_score = 0.1;
			//else
			//	tmpPairDistance_score = 0;
			double tmpScore = tmpMappedLength 		        
			    - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight;

			oriAlignPairNew_score_Nor2Rcm1.push_back(tmpScore);
		}
	}*/

	void chooseAllPairedAlignment()
	{
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			finalAlignPair_Nor1Rcm2.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
		}

		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			finalAlignPair_Nor2Rcm1.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
		}		
	}

	void alignmentFilter_fixPhase1_SJpenalty_Nor1Rcm2Only(int readLength_1, int readLength_2)
	{
		double min_map_score_keptAlignmentPair = ((readLength_1 + readLength_2)/100) * Min_Map_Score_keptAlignmentPair_PerHundredBases;
			//Min_Map_Score_keptAlignmentPair;
		//cout << "start to do pairing ..." << endl;
		this->pairingAlignment2OriPair_Nor1Rcm2Only();
		//cout << "start to do oriAlignPair2oriAlignPairNew ..." << endl;
		this->oriAlignPair2oriAlignPairNew_Nor1Rcm2Only();
		if(oriAlignPair_Nor1Rcm2_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
		{
			return;
			//oriAlignPair_Nor1Rcm2_new.clear();
			//oriAlignPair_Nor2Rcm1_new.clear();
		}
		// if(oriAlignPair_Nor2Rcm1_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
		// {
		// 	return;
		// 	//oriAlignPair_Nor1Rcm2_new.clear();
		// 	//oriAlignPair_Nor2Rcm1_new.clear();
		// }
		//cout << "start to do getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty ..." << endl;
		this->getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty_Nor1Rcm2();
		
		if(oriAlignPair_Nor1Rcm2_new.size()// + oriAlignPair_Nor2Rcm1_new.size() 
			== 0)
		{
			return;
		}
		//cout << "finish to do getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty ..." << endl;
		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp;
		//vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp;
		vector< double > finalAlignPair_Nor1Rcm2_tmp_score;
		//vector< double > finalAlignPair_Nor2Rcm1_tmp_score;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{		
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlignPair_Nor1Rcm2_tmp.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
				finalAlignPair_Nor1Rcm2_tmp_score.push_back(tmpScore);
			}
		}	
		// for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		// {		
		// 	double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
		// 	if(tmpScore > min_map_score_keptAlignmentPair)	
		// 	{
		// 		if(tmpScore > highestPairAlignmentScore)
		// 			highestPairAlignmentScore = tmpScore;
		// 		finalAlignPair_Nor2Rcm1_tmp.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
		// 		finalAlignPair_Nor2Rcm1_tmp_score.push_back(tmpScore);
		// 	}
		// }
		//cout << "start to do filterOverlapFinalAlignPair ..." << endl;
		this->filterOverlapFinalAlignPair_Nor1Rcm2Only(
			finalAlignPair_Nor1Rcm2_tmp, //finalAlignPair_Nor2Rcm1_tmp,
			finalAlignPair_Nor1Rcm2_tmp_score, //finalAlignPair_Nor2Rcm1_tmp_score,
			true);
		return;
	}

	void alignmentFilter_fixPhase1_SJpenalty(int readLength_1, int readLength_2)
	{
		double min_map_score_keptAlignmentPair 
			= //((readLength_1 + readLength_2)/100) * Min_Map_Score_keptAlignmentPair_PerHundredBases;
			Min_Map_Score_keptAlignmentPair;
		//cout << "start to do pairing ..." << endl;
		this->pairingAlignment2OriPair();
		//cout << "oriAlignPair_Nor1Rcm2.size(): " << oriAlignPair_Nor1Rcm2.size() << endl;
		//cout << "(oriAlignPair_Nor1Rcm2[0].second).size(): " << (oriAlignPair_Nor1Rcm2[0].second).size() << endl;
		//cout << "oriAlignPair_Nor2Rcm1.size(): " << oriAlignPair_Nor2Rcm1.size() << endl;
		//cout << "start to do oriAlignPair2oriAlignPairNew ..." << endl;
		this->oriAlignPair2oriAlignPairNew();
		if(oriAlignPair_Nor1Rcm2_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
			return;
		if(oriAlignPair_Nor2Rcm1_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
			return;
		//cout << "start to do getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty ..." << endl;
		this->getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty();
		
		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
			return;
		//cout << "finish to do getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty ..." << endl;
		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp;
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp;
		vector< double > finalAlignPair_Nor1Rcm2_tmp_score;
		vector< double > finalAlignPair_Nor2Rcm1_tmp_score;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{		
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlignPair_Nor1Rcm2_tmp.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
				finalAlignPair_Nor1Rcm2_tmp_score.push_back(tmpScore);
			}
		}	
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{		
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)	
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlignPair_Nor2Rcm1_tmp.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
				finalAlignPair_Nor2Rcm1_tmp_score.push_back(tmpScore);
			}
		}
		//cout << "finalAlignPair_Nor1Rcm2_tmp.size(): " << finalAlignPair_Nor1Rcm2_tmp.size() << endl;
		//cout << "finalAlignPair_Nor2Rcm1_tmp.size(): " << finalAlignPair_Nor2Rcm1_tmp.size() << endl;
		//cout << "start to do filterOverlapFinalAlignPair ..." << endl;
		this->filterOverlapFinalAlignPair(finalAlignPair_Nor1Rcm2_tmp, finalAlignPair_Nor2Rcm1_tmp,
			finalAlignPair_Nor1Rcm2_tmp_score, finalAlignPair_Nor2Rcm1_tmp_score,
			true);
		return;
	}

	void alignmentFilter_fixOneEndUnmapped_SJpenalty(int readLength_1, int readLength_2)
	{
		double min_map_score_keptAlignmentPair = //((readLength_1 + readLength_2)/100) * Min_Map_Score_keptAlignmentPair_PerHundredBases;
			Min_Map_Score_keptAlignmentPair;
		this->pairingAlignment2OriPair();
		this->oriAlignPair2oriAlignPairNew();
		this->getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty();
		
		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
			return;
		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp;
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp;
		vector< double > finalAlignPair_Nor1Rcm2_tmp_score;
		vector< double > finalAlignPair_Nor2Rcm1_tmp_score;

		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{		
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;				
				finalAlignPair_Nor1Rcm2_tmp.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
				finalAlignPair_Nor1Rcm2_tmp_score.push_back(tmpScore);
			}
		}	
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{		
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)	
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				finalAlignPair_Nor2Rcm1_tmp.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
				finalAlignPair_Nor2Rcm1_tmp_score.push_back(tmpScore);
			}
		}
		this->filterOverlapFinalAlignPair(finalAlignPair_Nor1Rcm2_tmp, finalAlignPair_Nor2Rcm1_tmp,
			finalAlignPair_Nor1Rcm2_tmp_score, finalAlignPair_Nor2Rcm1_tmp_score,
			true);
		return;
	}

	void chooseBestAlignment_phase1_subRegion()
	{
		this->getMismatchNumForEveryAlignment();
		this->getEndMatchPosForEveryAlignment();
		
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			int tmpPair_NO_1 = tmp;
			int tmpPair_NO_2 = tmp;
		
			string tmpChrName = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromName();
			int tmpChrMapPosStart_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnAlignChromPos();
			int tmpChrMapPosEnd_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnEndMatchedPosInChr();
			int tmpMappedLength_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->mappedLength();
			int tmpMismatchNum_nor = (norAlignmentInfo_PE_1[tmpPair_NO_1])->returnMismatchNum();
			norAlignmentInfo_PE_1[tmpPair_NO_1]->generateIndelLength();
			int tmpInsertionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnInsertionLength();
			int tmpDeletionLength_nor = norAlignmentInfo_PE_1[tmpPair_NO_1]->returnDeletionLength();

			int tmpChrMapPosStart_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnAlignChromPos();
			int tmpChrMapPosEnd_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnEndMatchedPosInChr();
			int tmpMappedLength_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->mappedLength();
			int tmpMismatchNum_rcm = (rcmAlignmentInfo_PE_2[tmpPair_NO_2])->returnMismatchNum();
			rcmAlignmentInfo_PE_2[tmpPair_NO_2]->generateIndelLength();
			int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnInsertionLength();
			int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_2[tmpPair_NO_2]->returnDeletionLength();

			int tmpChrMapPosStart, tmpChrMapPosEnd;
			if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
			{
				tmpChrMapPosStart = tmpChrMapPosStart_nor;
			}
			else
			{
				tmpChrMapPosStart = tmpChrMapPosStart_rcm;
			}

			if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
			{
				tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
			}
			else
			{
				tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
			}
			int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
			int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
			int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
			int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
			int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

			oriAlignPair_Nor1Rcm2_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
			oriAlignPairNew_chrNameVec_Nor1Rcm2.push_back(tmpChrName);
			oriAlignPairNew_startMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosStart);
			oriAlignPairNew_endMapPosVec_Nor1Rcm2.push_back(tmpChrMapPosEnd);
			oriAlignPairNew_mappedLength_Nor1Rcm2.push_back(tmpMappedLength);
			oriAlignPairNew_mismatchNum_Nor1Rcm2.push_back(tmpMismatchNum);
			oriAlignPairNew_pairDistance_Nor1Rcm2.push_back(tmpPairDistance);
			oriAlignPairNew_insertionLength_Nor1Rcm2.push_back(tmpInsertionLength);
			oriAlignPairNew_deletionLength_Nor1Rcm2.push_back(tmpDeletionLength);

			int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
				(norAlignmentInfo_PE_1[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
				(rcmAlignmentInfo_PE_2[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
			oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2.push_back(SJconfidenceLevel);
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			int tmpPair_NO_1 = tmp;
			int tmpPair_NO_2 = tmp;

			string tmpChrName = (norAlignmentInfo_PE_2[tmpPair_NO_1])->returnAlignChromName();
			int tmpChrMapPosStart_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnAlignChromPos();
			int tmpChrMapPosEnd_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnEndMatchedPosInChr();
			int tmpMappedLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->mappedLength();
			int tmpMismatchNum_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnMismatchNum();
			norAlignmentInfo_PE_2[tmpPair_NO_1]->generateIndelLength();
			int tmpInsertionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnInsertionLength();
			int tmpDeletionLength_nor = norAlignmentInfo_PE_2[tmpPair_NO_1]->returnDeletionLength();

			int tmpChrMapPosStart_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnAlignChromPos();
			int tmpChrMapPosEnd_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnEndMatchedPosInChr();
			int tmpMappedLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->mappedLength();
			int tmpMismatchNum_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnMismatchNum();
			rcmAlignmentInfo_PE_1[tmpPair_NO_2]->generateIndelLength();
			int tmpInsertionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnInsertionLength();
			int tmpDeletionLength_rcm = rcmAlignmentInfo_PE_1[tmpPair_NO_2]->returnDeletionLength();

			int tmpChrMapPosStart, tmpChrMapPosEnd;
			if(tmpChrMapPosStart_nor <= tmpChrMapPosStart_rcm)
			{
				tmpChrMapPosStart = tmpChrMapPosStart_nor;
			}
			else
			{
				tmpChrMapPosStart = tmpChrMapPosStart_rcm;
			}

			if(tmpChrMapPosEnd_nor <= tmpChrMapPosEnd_rcm)
			{
				tmpChrMapPosEnd = tmpChrMapPosEnd_rcm;
			}
			else
			{
				tmpChrMapPosEnd = tmpChrMapPosEnd_nor;
			}
			int tmpMappedLength = tmpMappedLength_nor + tmpMappedLength_rcm;
			int tmpMismatchNum = tmpMismatchNum_nor + tmpMismatchNum_rcm;
			int tmpPairDistance = tmpChrMapPosStart_rcm - tmpChrMapPosEnd_nor;
			int tmpInsertionLength = tmpInsertionLength_nor + tmpInsertionLength_rcm;
			int tmpDeletionLength = tmpDeletionLength_nor + tmpDeletionLength_rcm;

			oriAlignPair_Nor2Rcm1_new.push_back(pair<int, int> (tmpPair_NO_1, tmpPair_NO_2) );
			oriAlignPairNew_chrNameVec_Nor2Rcm1.push_back(tmpChrName);
			oriAlignPairNew_startMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosStart);
			oriAlignPairNew_endMapPosVec_Nor2Rcm1.push_back(tmpChrMapPosEnd);
			oriAlignPairNew_mappedLength_Nor2Rcm1.push_back(tmpMappedLength);
			oriAlignPairNew_mismatchNum_Nor2Rcm1.push_back(tmpMismatchNum);
			oriAlignPairNew_pairDistance_Nor2Rcm1.push_back(tmpPairDistance);
			oriAlignPairNew_insertionLength_Nor2Rcm1.push_back(tmpInsertionLength);
			oriAlignPairNew_deletionLength_Nor2Rcm1.push_back(tmpDeletionLength);

			int SJconfidenceLevel = this->getPairedAlignmentSJconfidenceLevel(
				(norAlignmentInfo_PE_2[tmpPair_NO_1])->spliceJunctionConfidenceLevel(),
				(rcmAlignmentInfo_PE_1[tmpPair_NO_2])->spliceJunctionConfidenceLevel());
			oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1.push_back(SJconfidenceLevel);			
		}
		
		if(oriAlignPair_Nor1Rcm2_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
		{
			return;
		}
		if(oriAlignPair_Nor2Rcm1_new.size() > MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES)
		{
			return;
		}
		
		this->getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty();
		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
		{
			return;
		}

		double DifferenceMin = BEST_ALIGNMENT_DIFFERENCE_MAX;

		double tmpBestScore_Nor1Rcm2 = 0.0;
		int selectedBestAlignmentNO_Nor1Rcm2 = 0;// = 0;
		double tmpBestScore_Nor2Rcm1 = 0.0;
		int selectedBestAlignmentNO_Nor2Rcm1 = 0;// = 0;
		//bool selectedBestAlignment_Nor1Rcm2 = true;// = true;
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > tmpBestScore_Nor1Rcm2)
			{
				tmpBestScore_Nor1Rcm2 = tmpScore;
				selectedBestAlignmentNO_Nor1Rcm2 = tmp;
			}
		}
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > tmpBestScore_Nor2Rcm1)
			{
				tmpBestScore_Nor2Rcm1 = tmpScore;
				selectedBestAlignmentNO_Nor2Rcm1 = tmp;
			}
		}
		
		bool bestAlignmentInNor1Rcm2_bool = false;
		bool bestAlignmentInNor2Rcm1_bool = false;
		int bestAlignemnt_toCheck_index_start_Nor1Rcm2 = 0;
		int bestAlignemnt_toCheck_index_start_Nor2Rcm1 = 0;		
		
		if(fabs(tmpBestScore_Nor1Rcm2 - tmpBestScore_Nor2Rcm1) < DifferenceMin)
		{
			highestPairAlignmentScore = tmpBestScore_Nor1Rcm2;

			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;
			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;
		}
		else if(tmpBestScore_Nor1Rcm2 > tmpBestScore_Nor2Rcm1)
		{
			highestPairAlignmentScore = tmpBestScore_Nor1Rcm2;

			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;			
		}
		else
		{
			highestPairAlignmentScore = tmpBestScore_Nor2Rcm1;

			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;			
		}

		if(bestAlignmentInNor1Rcm2_bool)
		{
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor1Rcm2*/; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp]; 
				if(fabs(tmpBestScore_Nor1Rcm2 - tmpScore) < DifferenceMin)
				{
					finalAlignPair_Nor1Rcm2.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					//finalAlignPair_Nor1Rcm2_tmp_score.push_back(tmpScore);
				}
			}
		}
		if(bestAlignmentInNor2Rcm1_bool)
		{
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor2Rcm1*/; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
				if(fabs(tmpBestScore_Nor2Rcm1 - tmpScore) < DifferenceMin)
				{
					finalAlignPair_Nor2Rcm1.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					//finalAlignPair_Nor2Rcm1_tmp_score.push_back(tmpScore);
				}
			}
		}		
	}

	void chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty()
	{
		this->pairingAlignment2OriPair();
		this->oriAlignPair2oriAlignPairNew();
		this->getScoreForEachPair_mapLen_mis_peDis_SJconfidenceLevel_mappingAreaSizePenalty();
		
		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
			return;
		double DifferenceMin = BEST_ALIGNMENT_DIFFERENCE_MAX;

		double tmpBestScore_Nor1Rcm2 = 0.0;
		int selectedBestAlignmentNO_Nor1Rcm2 = 0;// = 0;
		double tmpBestScore_Nor2Rcm1 = 0.0;
		int selectedBestAlignmentNO_Nor2Rcm1 = 0;// = 0;
		//bool selectedBestAlignment_Nor1Rcm2 = true;// = true;
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > tmpBestScore_Nor1Rcm2)
			{
				//if(tmpScore > highestPairAlignmentScore)
				//	highestPairAlignmentScore = tmpScore;
				tmpBestScore_Nor1Rcm2 = tmpScore;
				selectedBestAlignmentNO_Nor1Rcm2 = tmp;
				//selectedBestAlignment_Nor1Rcm2 = true;
			}
		}
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > tmpBestScore_Nor2Rcm1)
			{
				//if(tmpScore > highestPairAlignmentScore)
				//	highestPairAlignmentScore = tmpScore;
				tmpBestScore_Nor2Rcm1 = tmpScore;
				selectedBestAlignmentNO_Nor2Rcm1 = tmp;
				//selectedBestAlignment_Nor1Rcm2 = false;
			}
		}

		bool bestAlignmentInNor1Rcm2_bool = false;
		bool bestAlignmentInNor2Rcm1_bool = false;
		int bestAlignemnt_toCheck_index_start_Nor1Rcm2 = 0;
		int bestAlignemnt_toCheck_index_start_Nor2Rcm1 = 0;		
		if(fabs(tmpBestScore_Nor1Rcm2 - tmpBestScore_Nor2Rcm1) < DifferenceMin)
		{
			highestPairAlignmentScore = tmpBestScore_Nor1Rcm2;

			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;
			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;
		}
		else if(tmpBestScore_Nor1Rcm2 > tmpBestScore_Nor2Rcm1)
		{
			highestPairAlignmentScore = tmpBestScore_Nor1Rcm2;

			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;			
		}
		else
		{
			highestPairAlignmentScore = tmpBestScore_Nor2Rcm1;

			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;			
		}

		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp; 
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp;		

		vector< double > finalAlignPair_Nor1Rcm2_tmp_score;  
		vector< double > finalAlignPair_Nor2Rcm1_tmp_score;	

		if(bestAlignmentInNor1Rcm2_bool)
		{
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor1Rcm2*/; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp]; 
				if(fabs(tmpBestScore_Nor1Rcm2 - tmpScore) < DifferenceMin)
				{
					finalAlignPair_Nor1Rcm2_tmp.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					finalAlignPair_Nor1Rcm2_tmp_score.push_back(tmpScore);
				}
			}
		}
		if(bestAlignmentInNor2Rcm1_bool)
		{
			for(int tmp = 0/*bestAlignemnt_toCheck_index_start_Nor2Rcm1*/; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
				if(fabs(tmpBestScore_Nor2Rcm1 - tmpScore) < DifferenceMin)
				{
					finalAlignPair_Nor2Rcm1_tmp.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					finalAlignPair_Nor2Rcm1_tmp_score.push_back(tmpScore);
				}
			}
		}

		this->filterOverlapFinalAlignPair(finalAlignPair_Nor1Rcm2_tmp, finalAlignPair_Nor2Rcm1_tmp,
			finalAlignPair_Nor1Rcm2_tmp_score, finalAlignPair_Nor2Rcm1_tmp_score, false);

		return;
	}

	/*
	void chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical()
	{
		this->pairingAlignment2OriPair();
		//this->pairingAlignment();

		//cout << "***************" << endl;
		//cout << this->printOutPairingResults() << endl;
		//cout << "***************" << endl;

		this->oriAlignPair2oriAlignPairNew();
		this->getScoreForEachPair_mapLen_mis_peDis();
		
		if(oriAlignPair_Nor1Rcm2_new.size() + oriAlignPair_Nor2Rcm1_new.size() == 0)
		{
			return;
		}
		double DifferenceMin = BEST_ALIGNMENT_DIFFERENCE_MAX;

		double tmpBestScore_Nor1Rcm2 = 0.0;
		int selectedBestAlignmentNO_Nor1Rcm2 = 0;// = 0;
		double tmpBestScore_Nor2Rcm1 = 0.0;
		int selectedBestAlignmentNO_Nor2Rcm1 = 0;// = 0;
		//bool selectedBestAlignment_Nor1Rcm2 = true;// = true;
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp];
			if(tmpScore > tmpBestScore_Nor1Rcm2)
			{
				tmpBestScore_Nor1Rcm2 = tmpScore;
				selectedBestAlignmentNO_Nor1Rcm2 = tmp;
				//selectedBestAlignment_Nor1Rcm2 = true;
			}
		}
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
		{
			double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
			if(tmpScore > tmpBestScore_Nor2Rcm1)
			{
				tmpBestScore_Nor2Rcm1 = tmpScore;
				selectedBestAlignmentNO_Nor2Rcm1 = tmp;
				//selectedBestAlignment_Nor1Rcm2 = false;
			}
		}

		bool bestAlignmentInNor1Rcm2_bool = false;
		bool bestAlignmentInNor2Rcm1_bool = false;
		int bestAlignemnt_toCheck_index_start_Nor1Rcm2 = 0;
		int bestAlignemnt_toCheck_index_start_Nor2Rcm1 = 0;		
		if(fabs(tmpBestScore_Nor1Rcm2 - tmpBestScore_Nor2Rcm1) < DifferenceMin)
		{
			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;
			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;
		}
		else if(tmpBestScore_Nor1Rcm2 > tmpBestScore_Nor2Rcm1)
		{
			bestAlignmentInNor1Rcm2_bool = true;
			bestAlignemnt_toCheck_index_start_Nor1Rcm2 = selectedBestAlignmentNO_Nor1Rcm2;			
		}
		else
		{
			bestAlignmentInNor2Rcm1_bool = true;
			bestAlignemnt_toCheck_index_start_Nor2Rcm1 = selectedBestAlignmentNO_Nor2Rcm1;			
		}
		//if(selectedBestAlignment_Nor1Rcm2)
		//{
		//vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp;
		//vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp;

		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp_SJcanonical;  // all SJ canonical or no SJ
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp_SJcanonical;		

		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp_SJsemicanonical;
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp_SJsemicanonical;		

		vector<pair<int, int> > finalAlignPair_Nor1Rcm2_tmp_SJnoncanonical;
		vector<pair<int, int> > finalAlignPair_Nor2Rcm1_tmp_SJnoncanonical;

		if(bestAlignmentInNor1Rcm2_bool)
		{
			for(int tmp = bestAlignemnt_toCheck_index_start_Nor1Rcm2; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor1Rcm2[tmp]; 
				if(fabs(tmpBestScore_Nor1Rcm2 - tmpScore) < DifferenceMin)
				{
					//finalAlignPair_Nor1Rcm2_tmp.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					if(oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2[tmp] <= SPLICE_JUNCTION_CANONICAL_ONLY ) // all SJ canonical or no SJ
					{
						finalAlignPair_Nor1Rcm2_tmp_SJcanonical.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					}					
					else if(oriAlignPairNew_SJconfidenceLevel_Nor1Rcm2[tmp] == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL)
					{				
						finalAlignPair_Nor1Rcm2_tmp_SJsemicanonical.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					}
					else
					{						
						finalAlignPair_Nor1Rcm2_tmp_SJnoncanonical.push_back(oriAlignPair_Nor1Rcm2_new[tmp]);
					}
				}
			}
		}
		if(bestAlignmentInNor2Rcm1_bool)
		{
			for(int tmp = bestAlignemnt_toCheck_index_start_Nor2Rcm1; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp ++)
			{
				double tmpScore = oriAlignPairNew_score_Nor2Rcm1[tmp];
				if(fabs(tmpBestScore_Nor2Rcm1 - tmpScore) < DifferenceMin)
				{
					//finalAlignPair_Nor2Rcm1_tmp.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					if(oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1[tmp] <= SPLICE_JUNCTION_CANONICAL_ONLY ) // all SJ canonical or no SJ
					{
						finalAlignPair_Nor2Rcm1_tmp_SJcanonical.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					}	
					else if(oriAlignPairNew_SJconfidenceLevel_Nor2Rcm1[tmp] == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL)
					{
						finalAlignPair_Nor2Rcm1_tmp_SJsemicanonical.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					}
					else
					{
						finalAlignPair_Nor2Rcm1_tmp_SJnoncanonical.push_back(oriAlignPair_Nor2Rcm1_new[tmp]);
					}  
				}
			}
		}
		vector < pair<int,int> > finalAlignPair_Nor1Rcm2_inter_vec;
		vector < pair<int,int> > finalAlignPair_Nor2Rcm1_inter_vec;	

		if((finalAlignPair_Nor1Rcm2_tmp_SJcanonical.size() + finalAlignPair_Nor2Rcm1_tmp_SJcanonical.size()) > 0)
		{
			for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp_SJcanonical.size(); tmp++)
			{
				finalAlignPair_Nor1Rcm2_inter_vec.push_back(finalAlignPair_Nor1Rcm2_tmp_SJcanonical[tmp]);
				//return;
			}
			for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp_SJcanonical.size(); tmp++)
			{
				finalAlignPair_Nor2Rcm1_inter_vec.push_back(finalAlignPair_Nor2Rcm1_tmp_SJcanonical[tmp]);
				//return;
			}		
		}
		else if((finalAlignPair_Nor1Rcm2_tmp_SJsemicanonical.size() 
			+ finalAlignPair_Nor2Rcm1_tmp_SJsemicanonical.size()) > 0)
		{
			for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp_SJsemicanonical.size(); tmp++)
			{
				finalAlignPair_Nor1Rcm2_inter_vec.push_back(finalAlignPair_Nor1Rcm2_tmp_SJsemicanonical[tmp]);
				//return;
			}
			for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp_SJsemicanonical.size(); tmp++)
			{
				finalAlignPair_Nor2Rcm1_inter_vec.push_back(finalAlignPair_Nor2Rcm1_tmp_SJsemicanonical[tmp]);
				//return;
			}				
		}
		else
		{
			for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp_SJnoncanonical.size(); tmp++)
			{
				finalAlignPair_Nor1Rcm2_inter_vec.push_back(finalAlignPair_Nor1Rcm2_tmp_SJnoncanonical[tmp]);
				//return;
			}
			for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp_SJnoncanonical.size(); tmp++)
			{
				finalAlignPair_Nor2Rcm1_inter_vec.push_back(finalAlignPair_Nor2Rcm1_tmp_SJnoncanonical[tmp]);
				//return;
			}	
		}

		this->filterOverlapFinalAlignPair(finalAlignPair_Nor1Rcm2_inter_vec, finalAlignPair_Nor2Rcm1_inter_vec);

		return;
	}*/

	void filterOverlapFinalAlignPair(vector< pair<int,int> >& interAlignPair_vec_Nor1Rcm2,
		vector< pair<int,int> >& interAlignPair_vec_Nor2Rcm1,
		vector< double >& interAlignPair_vec_Nor1Rcm2_score,
		vector< double >& interAlignPair_vec_Nor2Rcm1_score,
		bool needToCheckAllAlignComplete_or_not) 
		// filter overlapped alignments, always choose 1st one, no matter if others are better than the 1st one.
	{
		set<int> invalidAlignPair_Nor1Rcm2;
		set<int> invalidAlignPair_Nor2Rcm1;
		// cout << "start to deal with interAlignPair_vec_Nor1Rcm2" << endl;
		for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlignPair_vec_Nor1Rcm2.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignPairOverlapOrNot(interAlignPair_vec_Nor1Rcm2[tmp].first,
					interAlignPair_vec_Nor1Rcm2[tmp2].first, interAlignPair_vec_Nor1Rcm2[tmp].second,
					interAlignPair_vec_Nor1Rcm2[tmp2].second, true);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlignPair_vec_Nor1Rcm2_score[tmp];
					double tmpScore_2 = interAlignPair_vec_Nor1Rcm2_score[tmp2];
					bool compare2AlignPair_bool;
					if(tmpScore_1 > tmpScore_2)
						compare2AlignPair_bool = true;
					else if(tmpScore_1 < tmpScore_2)
						compare2AlignPair_bool = false;
					else
						compare2AlignPair_bool = this->compare2AlignPair(interAlignPair_vec_Nor1Rcm2[tmp].first,
							interAlignPair_vec_Nor1Rcm2[tmp2].first, interAlignPair_vec_Nor1Rcm2[tmp].second,
							interAlignPair_vec_Nor1Rcm2[tmp2].second, true); 

					if(compare2AlignPair_bool) // tmp is better
						invalidAlignPair_Nor1Rcm2.insert(tmp2);
					else	
						invalidAlignPair_Nor1Rcm2.insert(tmp);
				}
			}
		}
		//cout << "start to deal with interAlignPair_vec_Nor2Rcm1" << endl;
		for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlignPair_vec_Nor2Rcm1.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignPairOverlapOrNot(interAlignPair_vec_Nor2Rcm1[tmp].first,
					interAlignPair_vec_Nor2Rcm1[tmp2].first, interAlignPair_vec_Nor2Rcm1[tmp].second,
					interAlignPair_vec_Nor2Rcm1[tmp2].second, false);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlignPair_vec_Nor2Rcm1_score[tmp];
					double tmpScore_2 = interAlignPair_vec_Nor2Rcm1_score[tmp2];
					bool compare2AlignPair_bool;
					if(tmpScore_1 > tmpScore_2)
						compare2AlignPair_bool = true;
					else if(tmpScore_1 < tmpScore_2)
						compare2AlignPair_bool = false;
					else
						compare2AlignPair_bool = this->compare2AlignPair(interAlignPair_vec_Nor2Rcm1[tmp].first,
							interAlignPair_vec_Nor2Rcm1[tmp2].first, interAlignPair_vec_Nor2Rcm1[tmp].second,
							interAlignPair_vec_Nor2Rcm1[tmp2].second, false); 
					if(compare2AlignPair_bool) // tmp is better
						invalidAlignPair_Nor2Rcm1.insert(tmp2);
					else	
						invalidAlignPair_Nor2Rcm1.insert(tmp);
				}
			}
		}

		if(needToCheckAllAlignComplete_or_not)
		{	
			vector< pair<int,int> > finalAlignPair_Nor1Rcm2_tmp;
			vector< pair<int,int> > finalAlignPair_Nor2Rcm1_tmp;
			vector< double > finalAlignPair_Nor1Rcm2_tmp_score;
			vector< double > finalAlignPair_Nor2Rcm1_tmp_score;

			for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp++)
			{
				if(invalidAlignPair_Nor1Rcm2.find(tmp) != invalidAlignPair_Nor1Rcm2.end())
					continue;
				else
				{
					finalAlignPair_Nor1Rcm2_tmp.push_back(interAlignPair_vec_Nor1Rcm2[tmp]);
					finalAlignPair_Nor1Rcm2_tmp_score.push_back(interAlignPair_vec_Nor1Rcm2_score[tmp]);
				}
			}
			
			for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp++)
			{
				if(invalidAlignPair_Nor2Rcm1.find(tmp) != invalidAlignPair_Nor2Rcm1.end())
					continue;
				else
				{
					finalAlignPair_Nor2Rcm1_tmp.push_back(interAlignPair_vec_Nor2Rcm1[tmp]);
					finalAlignPair_Nor2Rcm1_tmp_score.push_back(interAlignPair_vec_Nor2Rcm1_score[tmp]);
				}
			}

			bool allTmpFinalPair_complete = this->allTmpAlignmentPair_complete_bool(
				finalAlignPair_Nor1Rcm2_tmp, finalAlignPair_Nor2Rcm1_tmp);
			if(allTmpFinalPair_complete)
			{
				for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor1Rcm2_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						finalAlignPair_Nor1Rcm2.push_back(finalAlignPair_Nor1Rcm2_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor2Rcm1_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						finalAlignPair_Nor2Rcm1.push_back(finalAlignPair_Nor2Rcm1_tmp[tmp]);
				}			
			}
			else
			{
				for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor1Rcm2_tmp_score[tmp];
					finalAlignPair_Nor1Rcm2.push_back(finalAlignPair_Nor1Rcm2_tmp[tmp]);
				}
				for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor2Rcm1_tmp_score[tmp];
					finalAlignPair_Nor2Rcm1.push_back(finalAlignPair_Nor2Rcm1_tmp[tmp]);
				}	
			}
		}
		else
		{
			for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp++)
			{
				if(invalidAlignPair_Nor1Rcm2.find(tmp) != invalidAlignPair_Nor1Rcm2.end())
					continue;
				else
					finalAlignPair_Nor1Rcm2.push_back(interAlignPair_vec_Nor1Rcm2[tmp]);
			}
			
			for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp++)
			{
				if(invalidAlignPair_Nor2Rcm1.find(tmp) != invalidAlignPair_Nor2Rcm1.end())
					continue;
				else
					finalAlignPair_Nor2Rcm1.push_back(interAlignPair_vec_Nor2Rcm1[tmp]);
			}			
		}
	}

	void filterOverlapFinalAlignPair_Nor1Rcm2Only(
		vector< pair<int,int> >& interAlignPair_vec_Nor1Rcm2,
		//vector< pair<int,int> >& interAlignPair_vec_Nor2Rcm1,
		vector< double >& interAlignPair_vec_Nor1Rcm2_score,
		//vector< double >& interAlignPair_vec_Nor2Rcm1_score,
		bool needToCheckAllAlignComplete_or_not) 
		// filter overlapped alignments, always choose 1st one, no matter if others are better than the 1st one.
	{
		set<int> invalidAlignPair_Nor1Rcm2;
		//set<int> invalidAlignPair_Nor2Rcm1;

		// cout << "start to generate invalidAlignpair_" << endl;
		// cout << "interAlignPair_vec_Nor1Rcm2.size(): " << interAlignPair_vec_Nor1Rcm2.size() << endl;
		// cout << "interAlignPair_vec_Nor2Rcm1.size(): " << interAlignPair_vec_Nor2Rcm1.size() << endl;		

		// cout << "start to deal with interAlignPair_vec_Nor1Rcm2" << endl;
		for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp ++)
		{
			for(int tmp2 = tmp+1; tmp2 < interAlignPair_vec_Nor1Rcm2.size(); tmp2 ++)
			{
				bool overlap_orNot_bool = this->twoAlignPairOverlapOrNot(interAlignPair_vec_Nor1Rcm2[tmp].first,
					interAlignPair_vec_Nor1Rcm2[tmp2].first, interAlignPair_vec_Nor1Rcm2[tmp].second,
					interAlignPair_vec_Nor1Rcm2[tmp2].second, true);
				if(overlap_orNot_bool)
				{
					double tmpScore_1 = interAlignPair_vec_Nor1Rcm2_score[tmp];
					double tmpScore_2 = interAlignPair_vec_Nor1Rcm2_score[tmp2];
					bool compare2AlignPair_bool;
					if(tmpScore_1 > tmpScore_2)
					{
						compare2AlignPair_bool = true;
					}
					else if(tmpScore_1 < tmpScore_2)
					{
						compare2AlignPair_bool = false;
					}
					else
					{
						compare2AlignPair_bool = this->compare2AlignPair(interAlignPair_vec_Nor1Rcm2[tmp].first,
							interAlignPair_vec_Nor1Rcm2[tmp2].first, interAlignPair_vec_Nor1Rcm2[tmp].second,
							interAlignPair_vec_Nor1Rcm2[tmp2].second, true); 
					}

					if(compare2AlignPair_bool) // tmp is better
					{
						invalidAlignPair_Nor1Rcm2.insert(tmp2);
					}
					else
					{	
						invalidAlignPair_Nor1Rcm2.insert(tmp);
					}
				}
			}
		}
		//cout << "start to deal with interAlignPair_vec_Nor2Rcm1" << endl;
		// for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp ++)
		// {
		// 	for(int tmp2 = tmp+1; tmp2 < interAlignPair_vec_Nor2Rcm1.size(); tmp2 ++)
		// 	{
		// 		bool overlap_orNot_bool = this->twoAlignPairOverlapOrNot(interAlignPair_vec_Nor2Rcm1[tmp].first,
		// 			interAlignPair_vec_Nor2Rcm1[tmp2].first, interAlignPair_vec_Nor2Rcm1[tmp].second,
		// 			interAlignPair_vec_Nor2Rcm1[tmp2].second, false);
		// 		if(overlap_orNot_bool)
		// 		{
		// 			double tmpScore_1 = interAlignPair_vec_Nor2Rcm1_score[tmp];
		// 			double tmpScore_2 = interAlignPair_vec_Nor2Rcm1_score[tmp2];
		// 			bool compare2AlignPair_bool;
		// 			if(tmpScore_1 > tmpScore_2)
		// 			{
		// 				compare2AlignPair_bool = true;
		// 			}
		// 			else if(tmpScore_1 < tmpScore_2)
		// 			{
		// 				compare2AlignPair_bool = false;
		// 			}
		// 			else
		// 			{	
		// 				compare2AlignPair_bool = this->compare2AlignPair(interAlignPair_vec_Nor2Rcm1[tmp].first,
		// 					interAlignPair_vec_Nor2Rcm1[tmp2].first, interAlignPair_vec_Nor2Rcm1[tmp].second,
		// 					interAlignPair_vec_Nor2Rcm1[tmp2].second, false); 
		// 			}
		// 			if(compare2AlignPair_bool) // tmp is better
		// 			{
		// 				invalidAlignPair_Nor2Rcm1.insert(tmp2);
		// 			}
		// 			else
		// 			{	
		// 				invalidAlignPair_Nor2Rcm1.insert(tmp);
		// 			}
		// 		}
		// 	}
		// }

		//cout << "finish generating invalidAlignPair_" << endl;

		if(needToCheckAllAlignComplete_or_not)
		{	
			vector< pair<int,int> > finalAlignPair_Nor1Rcm2_tmp;
			//vector< pair<int,int> > finalAlignPair_Nor2Rcm1_tmp;
			vector< double > finalAlignPair_Nor1Rcm2_tmp_score;
			//vector< double > finalAlignPair_Nor2Rcm1_tmp_score;


			for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp++)
			{
				if(invalidAlignPair_Nor1Rcm2.find(tmp) != invalidAlignPair_Nor1Rcm2.end())
				{
					continue;
				}
				else
				{
					finalAlignPair_Nor1Rcm2_tmp.push_back(interAlignPair_vec_Nor1Rcm2[tmp]);
					finalAlignPair_Nor1Rcm2_tmp_score.push_back(interAlignPair_vec_Nor1Rcm2_score[tmp]);
				}
			}
			
			// for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp++)
			// {
			// 	if(invalidAlignPair_Nor2Rcm1.find(tmp) != invalidAlignPair_Nor2Rcm1.end())
			// 	{
			// 		continue;
			// 	}
			// 	else
			// 	{
			// 		finalAlignPair_Nor2Rcm1_tmp.push_back(interAlignPair_vec_Nor2Rcm1[tmp]);
			// 		finalAlignPair_Nor2Rcm1_tmp_score.push_back(interAlignPair_vec_Nor2Rcm1_score[tmp]);
			// 	}
			// }

			bool allTmpFinalPair_complete = this->allTmpAlignmentPair_complete_bool_Nor1Rcm2Only(
				finalAlignPair_Nor1Rcm2_tmp);//, finalAlignPair_Nor2Rcm1_tmp);
			if(allTmpFinalPair_complete)
			{
				for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor1Rcm2_tmp_score[tmp];
					if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
						finalAlignPair_Nor1Rcm2.push_back(finalAlignPair_Nor1Rcm2_tmp[tmp]);
				}
				// for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
				// {
				// 	double tmpScore = finalAlignPair_Nor2Rcm1_tmp_score[tmp];
				// 	if(fabs(highestPairAlignmentScore - tmpScore) < BEST_ALIGNMENT_DIFFERENCE_MAX)
				// 		finalAlignPair_Nor2Rcm1.push_back(finalAlignPair_Nor2Rcm1_tmp[tmp]);
				// }			
			}
			else
			{
				for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
				{
					double tmpScore = finalAlignPair_Nor1Rcm2_tmp_score[tmp];
					finalAlignPair_Nor1Rcm2.push_back(finalAlignPair_Nor1Rcm2_tmp[tmp]);
				}
				// for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
				// {
				// 	double tmpScore = finalAlignPair_Nor2Rcm1_tmp_score[tmp];
				// 	finalAlignPair_Nor2Rcm1.push_back(finalAlignPair_Nor2Rcm1_tmp[tmp]);
				// }	
			}
		}
		else
		{
			for(int tmp = 0; tmp < interAlignPair_vec_Nor1Rcm2.size(); tmp++)
			{
				if(invalidAlignPair_Nor1Rcm2.find(tmp) != invalidAlignPair_Nor1Rcm2.end())
					continue;
				else
				{
					finalAlignPair_Nor1Rcm2.push_back(interAlignPair_vec_Nor1Rcm2[tmp]);
				}
			}
			
			// for(int tmp = 0; tmp < interAlignPair_vec_Nor2Rcm1.size(); tmp++)
			// {
			// 	if(invalidAlignPair_Nor2Rcm1.find(tmp) != invalidAlignPair_Nor2Rcm1.end())
			// 		continue;
			// 	else
			// 	{
			// 		finalAlignPair_Nor2Rcm1.push_back(interAlignPair_vec_Nor2Rcm1[tmp]);
			// 	}
			// }			
		}
	}


	bool allTmpAlignmentPair_complete_bool(
		vector< pair<int,int> >& finalAlignPair_Nor1Rcm2_tmp,
		vector< pair<int,int> >& finalAlignPair_Nor2Rcm1_tmp)
	{
		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
		{
			int Nor1_alignInfo_index = finalAlignPair_Nor1Rcm2_tmp[tmp].first;
			int Rcm2_alignInfo_index = finalAlignPair_Nor1Rcm2_tmp[tmp].second;
			if(!((norAlignmentInfo_PE_1[Nor1_alignInfo_index]->noUnfixedHeadTailBool()) 
				&&(rcmAlignmentInfo_PE_2[Rcm2_alignInfo_index]->noUnfixedHeadTailBool())))
			{
				return false;
			}
		}
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
		{
			int Nor2_alignInfo_index = finalAlignPair_Nor2Rcm1_tmp[tmp].first;
			int Rcm1_alignInfo_index = finalAlignPair_Nor2Rcm1_tmp[tmp].second;
			if(!((norAlignmentInfo_PE_2[Nor2_alignInfo_index]->noUnfixedHeadTailBool()) 
				&&(rcmAlignmentInfo_PE_1[Rcm1_alignInfo_index]->noUnfixedHeadTailBool())))
			{
				return false;
			}
		}
		return true;;
	}

	bool allTmpAlignmentPair_complete_bool_Nor1Rcm2Only(
		vector< pair<int,int> >& finalAlignPair_Nor1Rcm2_tmp)
		//vector< pair<int,int> >& finalAlignPair_Nor2Rcm1_tmp)
	{
		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2_tmp.size(); tmp++)
		{
			int Nor1_alignInfo_index = finalAlignPair_Nor1Rcm2_tmp[tmp].first;
			int Rcm2_alignInfo_index = finalAlignPair_Nor1Rcm2_tmp[tmp].second;
			if(!((norAlignmentInfo_PE_1[Nor1_alignInfo_index]->noUnfixedHeadTailBool()) 
				&&(rcmAlignmentInfo_PE_2[Rcm2_alignInfo_index]->noUnfixedHeadTailBool())))
			{
				return false;
			}
		}
		// for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1_tmp.size(); tmp++)
		// {
		// 	int Nor2_alignInfo_index = finalAlignPair_Nor2Rcm1_tmp[tmp].first;
		// 	int Rcm1_alignInfo_index = finalAlignPair_Nor2Rcm1_tmp[tmp].second;
		// 	if(!((norAlignmentInfo_PE_2[Nor2_alignInfo_index]->noUnfixedHeadTailBool()) 
		// 		&&(rcmAlignmentInfo_PE_1[Rcm1_alignInfo_index]->noUnfixedHeadTailBool())))
		// 	{
		// 		return false;
		// 	}
		// }
		return true;;
	}	

	bool compare2AlignPair_sub(Alignment_Info* alignInfo_nor, Alignment_Info* alignInfo_rcm,
		Alignment_Info* alignInfo_nor_toCompare, Alignment_Info* alignInfo_rcm_toCompare)
	{
		int pair1_end1_startPos = alignInfo_nor->returnAlignChromPos();
		int pair1_end1_endPos = alignInfo_nor->returnEndMatchedPosInChr();
		int pair1_end2_startPos	= alignInfo_rcm->returnAlignChromPos();
		int pair1_end2_endPos = alignInfo_rcm->returnEndMatchedPosInChr();
		int pair2_end1_startPos = alignInfo_nor_toCompare->returnAlignChromPos(); 
		int pair2_end1_endPos = alignInfo_nor_toCompare->returnEndMatchedPosInChr();
		int pair2_end2_startPos = alignInfo_rcm_toCompare->returnAlignChromPos();
		int pair2_end2_endPos = alignInfo_rcm_toCompare->returnEndMatchedPosInChr();		

		int pair1_startPos = pair1_end1_startPos;
		int pair1_endPos = pair1_end2_endPos;
		int pair2_startPos = pair2_end1_startPos;
		int pair2_endPos = pair2_end2_endPos;
		#ifdef DETECT_CIRCULAR_RNA
			pair1_startPos = selectTheSmallestAmong4values(
				pair1_end1_startPos, pair1_end1_endPos, pair1_end2_startPos, pair1_end2_endPos);
			pair2_startPos = selectTheSmallestAmong4values(
				pair2_end1_startPos, pair2_end1_endPos, pair2_end2_startPos, pair2_end2_endPos);
			pair1_endPos = selectTheLargestAmong4values(
				pair1_end1_startPos, pair1_end1_endPos, pair1_end2_startPos, pair1_end2_endPos);
			pair2_endPos = selectTheLargestAmong4values(
				pair2_end1_startPos, pair2_end1_endPos, pair2_end2_startPos, pair2_end2_endPos);
		#else
		if(pair1_end2_startPos < pair1_startPos)
			pair1_startPos = pair1_end2_startPos;
		if(pair1_end1_endPos > pair1_endPos)
			pair1_endPos = pair1_end1_endPos;
		if(pair2_end2_startPos < pair2_startPos)
			pair2_startPos = pair2_end2_startPos;
		if(pair2_end1_endPos > pair2_endPos)
			pair2_endPos = pair2_end1_endPos;
		#endif

		int pair1_mappingAreaSize = pair1_endPos - pair1_startPos;
		int pair2_mappingAreaSize = pair2_endPos - pair2_startPos;

		if(pair1_mappingAreaSize <= pair2_mappingAreaSize)
			return true;
		else
			return false;
	}

	bool compare2AlignPair(int norAlignInfo_index, int norAlignInfo_index_toCompare, 
		int rcmAlignInfo_index, int rcmAlignInfo_index_toCompare,
		bool Nor1Rcm2_or_Nor2Rcm1_bool)
	{
		if(Nor1Rcm2_or_Nor2Rcm1_bool)
		{
			int norAlignInfo_1_index = norAlignInfo_index;
			int rcmAlignInfo_2_index = rcmAlignInfo_index;

			int norAlignInfo_1_index_toCompare = norAlignInfo_index_toCompare;
			int rcmAlignInfo_2_index_toCompare = rcmAlignInfo_index_toCompare; 
			return (this->compare2AlignPair_sub(
				norAlignmentInfo_PE_1[norAlignInfo_1_index],
				rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index],
				norAlignmentInfo_PE_1[norAlignInfo_1_index_toCompare],
				rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index_toCompare]));
		}
		else
		{
			int norAlignInfo_2_index = norAlignInfo_index;
			int rcmAlignInfo_1_index = rcmAlignInfo_index;

			int norAlignInfo_2_index_toCompare = norAlignInfo_index_toCompare;
			int rcmAlignInfo_1_index_toCompare = rcmAlignInfo_index_toCompare; 
			return (this->compare2AlignPair_sub(
				norAlignmentInfo_PE_2[norAlignInfo_2_index],
				rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index],
				norAlignmentInfo_PE_2[norAlignInfo_2_index_toCompare],
				rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index_toCompare]));
		}			
	}

	bool twoAlignPairOverlapOrNot(int norAlignInfo_index, int norAlignInfo_index_toCompare, 
		int rcmAlignInfo_index, int rcmAlignInfo_index_toCompare,
		bool Nor1Rcm2_or_Nor2Rcm1_bool)
	{
		if(Nor1Rcm2_or_Nor2Rcm1_bool)
		{
			int norAlignInfo_1_index = norAlignInfo_index;
			int rcmAlignInfo_2_index = rcmAlignInfo_index;

			int norAlignInfo_1_index_toCompare = norAlignInfo_index_toCompare;
			int rcmAlignInfo_2_index_toCompare = rcmAlignInfo_index_toCompare; 

			if(norAlignmentInfo_PE_1[norAlignInfo_1_index]->returnAlignChromName() 
				!= norAlignmentInfo_PE_1[norAlignInfo_1_index_toCompare]->returnAlignChromName())
			{
				return false;
			}
			else
			{	
				return (this->twoAlignPairOverlapOrNot_sub(
					norAlignmentInfo_PE_1[norAlignInfo_1_index]->returnAlignChromPos(), 
					norAlignmentInfo_PE_1[norAlignInfo_1_index]->returnEndMatchedPosInChr(),
					rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index]->returnAlignChromPos(),
					rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index]->returnEndMatchedPosInChr(),
					norAlignmentInfo_PE_1[norAlignInfo_1_index_toCompare]->returnAlignChromPos(), 
					norAlignmentInfo_PE_1[norAlignInfo_1_index_toCompare]->returnEndMatchedPosInChr(),
					rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index_toCompare]->returnAlignChromPos(),
					rcmAlignmentInfo_PE_2[rcmAlignInfo_2_index_toCompare]->returnEndMatchedPosInChr()));
			}
		}
		else
		{
			int norAlignInfo_2_index = norAlignInfo_index;
			int rcmAlignInfo_1_index = rcmAlignInfo_index;

			int norAlignInfo_2_index_toCompare = norAlignInfo_index_toCompare;
			int rcmAlignInfo_1_index_toCompare = rcmAlignInfo_index_toCompare; 

			if(norAlignmentInfo_PE_2[norAlignInfo_2_index]->returnAlignChromName() 
				!= norAlignmentInfo_PE_2[norAlignInfo_2_index_toCompare]->returnAlignChromName())
			{
				return false;
			}
			else
			{
				return (this->twoAlignPairOverlapOrNot_sub(
					norAlignmentInfo_PE_2[norAlignInfo_2_index]->returnAlignChromPos(), 
					norAlignmentInfo_PE_2[norAlignInfo_2_index]->returnEndMatchedPosInChr(),
					rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index]->returnAlignChromPos(),
					rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index]->returnEndMatchedPosInChr(),
					norAlignmentInfo_PE_2[norAlignInfo_2_index_toCompare]->returnAlignChromPos(), 
					norAlignmentInfo_PE_2[norAlignInfo_2_index_toCompare]->returnEndMatchedPosInChr(),
					rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index_toCompare]->returnAlignChromPos(),
					rcmAlignmentInfo_PE_1[rcmAlignInfo_1_index_toCompare]->returnEndMatchedPosInChr()));
			}
		}		
	}

	bool twoAlignPairOverlapOrNot_sub(int pair1_end1_startPos, int pair1_end1_endPos,
		int pair1_end2_startPos, int pair1_end2_endPos,
		int pair2_end1_startPos, int pair2_end1_endPos,
		int pair2_end2_startPos, int pair2_end2_endPos)
	{
		int pair1_startPos = pair1_end1_startPos;
		int pair1_endPos = pair1_end2_endPos;
		int pair2_startPos = pair2_end1_startPos;
		int pair2_endPos = pair2_end2_endPos;
		#ifdef DETECT_CIRCULAR_RNA
			pair1_startPos = selectTheSmallestAmong4values(
				pair1_end1_startPos, pair1_end1_endPos, pair1_end2_startPos, pair1_end2_endPos);
			pair2_startPos = selectTheSmallestAmong4values(
				pair2_end1_startPos, pair2_end1_endPos, pair2_end2_startPos, pair2_end2_endPos);
			pair1_endPos = selectTheLargestAmong4values(
				pair1_end1_startPos, pair1_end1_endPos, pair1_end2_startPos, pair1_end2_endPos);
			pair2_endPos = selectTheLargestAmong4values(
				pair2_end1_startPos, pair2_end1_endPos, pair2_end2_startPos, pair2_end2_endPos);
		#else
		if(pair1_end2_startPos < pair1_startPos)
			pair1_startPos = pair1_end2_startPos;
		if(pair1_end1_endPos > pair1_endPos)
			pair1_endPos = pair1_end1_endPos;
		if(pair2_end2_startPos < pair2_startPos)
			pair2_startPos = pair2_end2_startPos;
		if(pair2_end1_endPos > pair2_endPos)
			pair2_endPos = pair2_end1_endPos;
		#endif
		if( (pair1_endPos < pair2_startPos) || (pair1_startPos > pair2_endPos) )
			return false;
		else
			return true;
	}
	
	string printOutPairingResults()
	{
		string pairingResults;
		pairingResults = "pairing results information:\n";
		pairingResults += "Nor1Rcm2 pair:\n";
		for(int tmp = 0; tmp < oriAlignPair_Nor1Rcm2_new.size(); tmp++)
		{
			pairingResults = "Pair " + int_to_str(tmp + 1) + ": ";
			pairingResults = pairingResults + int_to_str(oriAlignPair_Nor1Rcm2_new[tmp].first)
				+ ", " + int_to_str(oriAlignPair_Nor1Rcm2_new[tmp].second) + "\n";
		}  
		pairingResults += "Nor2Rcm1 pair:\n";
		for(int tmp = 0; tmp < oriAlignPair_Nor2Rcm1_new.size(); tmp++)
		{
			pairingResults = "Pair " + int_to_str(tmp + 1) + ": ";
			pairingResults = pairingResults + int_to_str(oriAlignPair_Nor2Rcm1_new[tmp].first)
				+ ", " + int_to_str(oriAlignPair_Nor2Rcm1_new[tmp].second) + "\n";
		}  	
		return pairingResults;
	}

	void generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_forLocalIndexMap2filterFusion(
		string& tmpOriPeSamStr_left, string& tmpOriPeSamStr_right, string& tmpDetectedFusionOtherEndSamStr, 
		PE_Read_Info& peReadInfo, Index_Info* indexInfo, 
		bool fasta_or_fastq_bool, int& tmp_multiMapSeg_maxLength)
	{
		bool SE_or_PE_bool = false;
		string line1, line2, line3, line4, line5, line6, line7, line8, line9, line10;

		int startLoc = 0;

		vector<string> samFieldVec_left;
		vector<string> samFieldVec_right;
		vector<string> samFieldVec_fusionOtherEnd;
		//cout << "start to extract SAM_1" << endl;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpOriPeSamStr_left.find("\t", startLoc);
			string tmpSamField = tmpOriPeSamStr_left.substr(startLoc, tabLoc-startLoc);
			samFieldVec_left.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		startLoc = 0;
		//cout << "start to extract SAM_2" << endl;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpOriPeSamStr_right.find("\t", startLoc);
			string tmpSamField = tmpOriPeSamStr_right.substr(startLoc, tabLoc-startLoc);
			samFieldVec_right.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		startLoc = 0;
		//cout << "start to extract SAM_fusionOtherEnd" << endl;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpDetectedFusionOtherEndSamStr.find("\t", startLoc);
			string tmpSamField = tmpDetectedFusionOtherEndSamStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpSamField: " << tmpSamField << endl;
			samFieldVec_fusionOtherEnd.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string tmpReadSeq_fusionOtherEnd = samFieldVec_fusionOtherEnd[9];
		int tmpReadSeqLength_fusionOtherEnd = tmpReadSeq_fusionOtherEnd.length();
		//cout << "tmpReadSeqLength_fusionOtherEnd: " << tmpReadSeqLength_fusionOtherEnd << endl;
		bool leftRead_or_rightRead_read1_bool;
		int flag_left = atoi(samFieldVec_left[1].c_str());
		if(flag_left & 0x40)
			leftRead_or_rightRead_read1_bool = true;
		else
			leftRead_or_rightRead_read1_bool = false;
		//cout << "flag_left: " << flag_left << endl;
		int loc_XM = tmpOriPeSamStr_left.find("XM:i:");
		int loc_nextTab = tmpOriPeSamStr_left.find("\t", loc_XM+1);
		string XM_fieldStr = tmpOriPeSamStr_left.substr(loc_XM+5, loc_nextTab-1-(loc_XM+5)+1);
		//cout << "XM_fieldStr: " << XM_fieldStr << endl;
		bool leftHead_or_rightTail_bool;
		bool for_or_rcm_bool;
		int loc_ZF = tmpDetectedFusionOtherEndSamStr.find("ZF:Z:FUS");
		string ZF_fieldStr = tmpDetectedFusionOtherEndSamStr.substr(loc_ZF);
		//cout << "ZF_fieldStr: " << ZF_fieldStr << endl;
		if(ZF_fieldStr == "ZF:Z:FUS_leftHead_for")
		{
			leftHead_or_rightTail_bool = true;
			for_or_rcm_bool = true;
		}
		else if(ZF_fieldStr == "ZF:Z:FUS_leftHead_rev")
		{
			leftHead_or_rightTail_bool = true;
			for_or_rcm_bool = false;
		}
		else if(ZF_fieldStr == "ZF:Z:FUS_rightTail_for")
		{
			leftHead_or_rightTail_bool = false;
			for_or_rcm_bool = true;
		}
		else if(ZF_fieldStr == "ZF:Z:FUS_rightTail_rev")
		{
			leftHead_or_rightTail_bool = false;
			for_or_rcm_bool = false;
		}
		else
		{
			cout << "incorrect ZF:Z: field" << endl;
			exit(1);
		}

		string tmpReadName_left = samFieldVec_left[0];
		string tmpReadName_right = samFieldVec_right[0];
		string tmpReadSeq_left = samFieldVec_left[9];
		string tmpReadSeq_right = samFieldVec_right[9];
		int tmpReadSeqLength_left = tmpReadSeq_left.length();
		int tmpReadSeqLength_right = tmpReadSeq_right.length();
		string tmpReadSeq_left_ori = tmpReadSeq_left;
		string tmpReadSeq_right_ori = convertStringToReverseComplement(tmpReadSeq_right);
		string tmpQualSeq_left = samFieldVec_left[10];
		string tmpQualSeq_right = samFieldVec_right[10];
		string tmpQualSeq_left_ori;
		string tmpQualSeq_right_ori;
		string tmpCigarString_fusionOtherEnd = samFieldVec_fusionOtherEnd[5];
		string tmpChrName_fusionOtherEnd = samFieldVec_fusionOtherEnd[2];
		string tmpChrPos_fusionOtherEnd = samFieldVec_fusionOtherEnd[3];
		if(fasta_or_fastq_bool)
		{
			tmpQualSeq_left_ori = "*";
			tmpQualSeq_right_ori = "*";
		}
		else
		{
			tmpQualSeq_left_ori = tmpQualSeq_left;
			tmpQualSeq_right_ori = convertQualityScoreString2Reverse(tmpQualSeq_right);
		}
		
		line2 = tmpReadSeq_left_ori;
		line3 = tmpQualSeq_left_ori;
		line5 = tmpReadSeq_right_ori;
		line6 = tmpQualSeq_right_ori;
		line7 = "Nor1:";
		line8 = "Rcm1:";
		line9 = "Nor2:";
		line10 = "Rcm2:";		
		if(leftHead_or_rightTail_bool)
		{
			int tmpSoftClipLength = tmpReadSeqLength_left - tmpReadSeqLength_fusionOtherEnd;
			if(for_or_rcm_bool)
			{
				line1 = tmpReadName_left + "/1\t1\t0\t" + XM_fieldStr;
				line4 = tmpReadName_right + "/2\t0\t0";
				string tmpSoftClipCigarString = tmpCigarString_fusionOtherEnd + int_to_str(tmpSoftClipLength) + "S";
				line7 = line7 + "\t" + tmpChrName_fusionOtherEnd + "," + tmpChrPos_fusionOtherEnd + "," 
					+ tmpSoftClipCigarString + ",0,,,";
			}
			else
			{
				line1 = tmpReadName_left + "/1\t0\t1\t" + XM_fieldStr;
				line4 = tmpReadName_right + "/2\t0\t0";
				string tmpSoftClipCigarString = int_to_str(tmpSoftClipLength) + "S" + tmpCigarString_fusionOtherEnd;
				line8 = line8 + "\t" + tmpChrName_fusionOtherEnd + "," + tmpChrPos_fusionOtherEnd + "," 
					+ tmpSoftClipCigarString + ",0,,,";
			}
		}
		else
		{
			int tmpSoftClipLength = tmpReadSeqLength_right - tmpReadSeqLength_fusionOtherEnd;
			if(for_or_rcm_bool)
			{
				line1 = tmpReadName_left + "/1\t0\t0\t" + XM_fieldStr;
				line4 = tmpReadName_right + "/2\t1\t0";
				string tmpSoftClipCigarString = tmpCigarString_fusionOtherEnd + int_to_str(tmpSoftClipLength) + "S";
				line9 = line9 + "\t" + tmpChrName_fusionOtherEnd + "," + tmpChrPos_fusionOtherEnd + "," 
					+ tmpSoftClipCigarString + ",0,,,";
			}
			else
			{
				line1 = tmpReadName_left + "/1\t0\t0\t" + XM_fieldStr;
				line4 = tmpReadName_right + "/2\t0\t1";
				string tmpSoftClipCigarString = int_to_str(tmpSoftClipLength) + "S" + tmpCigarString_fusionOtherEnd;
				line10 = line10 + "\t" + tmpChrName_fusionOtherEnd + "," + tmpChrPos_fusionOtherEnd + "," 
					+ tmpSoftClipCigarString + ",0,,,";
			}		
		}
		// cout << "********************************************\npeAlignInfo:" << endl;
		// cout << line1 << endl << line2 << endl << line3 << endl << line4 << endl << line5 << endl << line6 << endl
		// 	 << line7 << endl << line8 << endl << line9 << endl << line10 << endl;

		this->generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_getline(
			line1, line2, line3, line4, line5, line6, line7, line8, line9, line10,
			peReadInfo, indexInfo, fasta_or_fastq_bool, SE_or_PE_bool, tmp_multiMapSeg_maxLength);
	}

	void generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_getline(
		const string& line1, const string& line2, const string& line3,
		const string& line4, const string& line5, const string& line6,
		const string& line7, const string& line8,
		const string& line9, const string& line10, PE_Read_Info& peReadInfo, 
		Index_Info* indexInfo, bool fasta_or_fastq_bool, bool SE_or_PE_bool, 
		int& multiMapSeg_maxLength)
	{
		highestPairAlignmentScore = 0;

		int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
		string readNameStr_1, readNameStr_2;

		int startSearchPos = 0, foundSearchPos;	foundSearchPos = line1.find("\t", startSearchPos);
		readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
			
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);
		Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
			
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		multiMapSeg_maxLength = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());


		startSearchPos = 0; foundSearchPos = line4.find("\t", startSearchPos);
		readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);
		Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				
		startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);		
		Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

		int readLength_1 = line2.length() - 1;
		int readLength_2 = line5.length() - 1;

		//peReadInfo.getFastaFormatReadInfo_new(readNameStr_1, readNameStr_2,
		//	line2, line5);
		peReadInfo.initiateReadInfo(readNameStr_1, readNameStr_2,
			line2, line5, line3, line6, fasta_or_fastq_bool, SE_or_PE_bool);
		this->getPeReadAlignmentInfo(line7, line8, line9, line10,
			Nor1Num, Rcm1Num, Nor2Num, Rcm2Num, indexInfo);
	}		

	void generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(
		const string& line1, const string& line2, const string& line3,
		const string& line4, const string& line5, const string& line6,
		const string& line7, const string& line8,
		const string& line9, const string& line10, PE_Read_Info& peReadInfo, 
		Index_Info* indexInfo, bool fasta_or_fastq_bool, bool SE_or_PE_bool,
		int& multiMapSeg_maxLength)
	{
		//cout << "generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline( starts ..." << endl;
		highestPairAlignmentScore = 0;

		int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
		string readNameStr_1, readNameStr_2;

		int startSearchPos = 0, foundSearchPos;	foundSearchPos = line1.find("\t", startSearchPos);
		readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
		//cout << "readNameStr_1: " << readNameStr_1 << endl;
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);
		Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
		
		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

		startSearchPos = foundSearchPos + 1; foundSearchPos = line1.find("\t", startSearchPos);		
		multiMapSeg_maxLength = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
		// cout << "Nor1Num: " << Nor1Num << endl;
		// cout << "Rcm1Num: " << Rcm1Num << endl;

		if(!SE_or_PE_bool)
		{	
			startSearchPos = 0; foundSearchPos = line4.find("\t", startSearchPos);
			readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
					
			startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);
			Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
					
			startSearchPos = foundSearchPos + 1; foundSearchPos = line4.find("\t", startSearchPos);		
			Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	
		}
		//peReadInfo.get_PE_Read_Info(readNameStr_1, readNameStr_2,
		//	line2, line5);
		//cout << "start to initiateReadInfo  ...." << endl;
		peReadInfo.initiateReadInfo(readNameStr_1, readNameStr_2,
			line2, line5, line3, line6, fasta_or_fastq_bool, SE_or_PE_bool);
		//cout << "end of initiating peReadInfo ..." << endl;
		this->getPeReadAlignmentInfo(line7, line8, line9, line10,
			Nor1Num, Rcm1Num, Nor2Num, Rcm2Num, indexInfo, SE_or_PE_bool);
	}

	bool finalPairExistsBool()
	{
		int finalPairNum = finalAlignPair_Nor1Rcm2.size() + finalAlignPair_Nor2Rcm1.size();
		if(finalPairNum == 0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	bool oriPairExistsBool()
	{
		int oriPairNum = oriAlignPair_Nor1Rcm2.size() + oriAlignPair_Nor2Rcm1.size();
		if(oriPairNum == 0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	bool alignInfoExistsBool()
	{
		int alignInfoNum = norAlignmentInfo_PE_1.size() + rcmAlignmentInfo_PE_1.size()
			+ norAlignmentInfo_PE_2.size() + rcmAlignmentInfo_PE_2.size();
		if(alignInfoNum == 0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	bool oneEndUnmapBool()
	{
		bool oriPairBool = this->oriPairExistsBool();
		bool alignInfoBool = this->alignInfoExistsBool();
		if((!oriPairBool)&&(alignInfoBool))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool allAlignmentInOriPairCompleted()
	{
		bool incomleteAlignmentExists = true;

		for(int tmp_1 = 0; tmp_1 < oriAlignPair_Nor1Rcm2.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				norAlignmentInfo_PE_1[(oriAlignPair_Nor1Rcm2[tmp_1].first)]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;
			for(int tmp_2 = 0; tmp_2 < 
				(oriAlignPair_Nor1Rcm2[tmp_1].second).size(); tmp_2++)
			{
				bool tmpBool_2 = 
					rcmAlignmentInfo_PE_2[((oriAlignPair_Nor1Rcm2[tmp_1].second)[tmp_2])]->noUnfixedHeadTailBool();
				if(!tmpBool_2)
					return false;
			}
		}

		for(int tmp_1 = 0; tmp_1 < oriAlignPair_Nor2Rcm1.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				norAlignmentInfo_PE_2[(oriAlignPair_Nor2Rcm1[tmp_1].first)]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;
			for(int tmp_2 = 0; tmp_2 < 
				(oriAlignPair_Nor2Rcm1[tmp_1].second).size(); tmp_2++)
			{
				bool tmpBool_2 = 
					rcmAlignmentInfo_PE_1[((oriAlignPair_Nor2Rcm1[tmp_1].second)[tmp_2])]->noUnfixedHeadTailBool();
				if(!tmpBool_2)
					return false;
			}
		}
		return true;
	}

	bool allAlignmentInFinalPairCompleted()
	{
		bool incomleteAlignmentExists = true;

		for(int tmp_1 = 0; tmp_1 < finalAlignPair_Nor1Rcm2.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				norAlignmentInfo_PE_1[(finalAlignPair_Nor1Rcm2[tmp_1].first)]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;

			bool tmpBool_2 = 
				rcmAlignmentInfo_PE_2[(finalAlignPair_Nor1Rcm2[tmp_1].second)]->noUnfixedHeadTailBool();
			if(!tmpBool_2)
				return false;
		}

		for(int tmp_1 = 0; tmp_1 < finalAlignPair_Nor2Rcm1.size(); tmp_1++)
		{
			bool tmpBool_1 = 
				norAlignmentInfo_PE_2[(finalAlignPair_Nor2Rcm1[tmp_1].first)]->noUnfixedHeadTailBool();
			if(!tmpBool_1)
				return false;

			bool tmpBool_2 = 
				rcmAlignmentInfo_PE_1[(finalAlignPair_Nor2Rcm1[tmp_1].second)]->noUnfixedHeadTailBool();
			if(!tmpBool_2)
				return false;
		}

		return true;
	}

	bool allFinalUnpairedAlignmentCompleted()
	{
		for(int tmp = 0; tmp < unpairedSEalignVec_final_Nor_end1.size(); tmp++)
		{
			int tmpIndex = unpairedSEalignVec_final_Nor_end1[tmp];
			bool tmpBool = 
				norAlignmentInfo_PE_1[tmpIndex]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}
		for(int tmp = 0; tmp < unpairedSEalignVec_final_Rcm_end1.size(); tmp++)
		{
			int tmpIndex = unpairedSEalignVec_final_Rcm_end1[tmp];
			bool tmpBool = 
				rcmAlignmentInfo_PE_1[tmpIndex]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}
		for(int tmp = 0; tmp < unpairedSEalignVec_final_Nor_end2.size(); tmp++)
		{
			int tmpIndex = unpairedSEalignVec_final_Nor_end2[tmp];
			bool tmpBool = 
				norAlignmentInfo_PE_2[tmpIndex]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}
		for(int tmp = 0; tmp < unpairedSEalignVec_final_Rcm_end2.size(); tmp++)
		{
			int tmpIndex = unpairedSEalignVec_final_Rcm_end2[tmp];
			bool tmpBool = 
				rcmAlignmentInfo_PE_2[tmpIndex]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}
		return true;
	}

	bool allUnpairedAlignmentCompleted()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			bool tmpBool = 
				norAlignmentInfo_PE_1[tmp]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}

		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			bool tmpBool = 
				norAlignmentInfo_PE_2[tmp]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;			
		}

		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			bool tmpBool = 
				rcmAlignmentInfo_PE_1[tmp]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;
		}

		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			bool tmpBool = 
				rcmAlignmentInfo_PE_2[tmp]->noUnfixedHeadTailBool();
			if(!tmpBool)
				return false;			
		}

		return true;
	}
	
	void getPeReadAlignmentInfo(const string& line5, const string& line6, const string& line7,
		const string& line8, int nor1Num, int rcm1Num, int nor2Num, int rcm2Num, Index_Info* indexInfo)
	{
		this->getPeReadAlignmentInfo(line5, line6, line7, line8, nor1Num, rcm1Num, nor2Num, rcm2Num, indexInfo, false);
	}

	void getPeReadAlignmentInfo(const string& line5, const string& line6, const string& line7,
		const string& line8, int nor1Num, int rcm1Num, int nor2Num, int rcm2Num, Index_Info* indexInfo, bool SE_or_PE_bool)
	{
		//Alignment_Info* tmpAlignmentInfo;
		//int tmp 
		highestPairAlignmentScore = 0;

		int tmpAlignInfoStrStartPos;
		int tmpAlignInfoStrEndPos;
		string tmpAlignInfoStr;

		norAlignmentInfo_PE_1.clear();
		rcmAlignmentInfo_PE_1.clear();
		norAlignmentInfo_PE_2.clear();
		rcmAlignmentInfo_PE_2.clear();		
		//cout << 
		if(nor1Num < 40)
		{
			tmpAlignInfoStrStartPos = line5.find("\t", 0) + 1;
			for(int tmp = 0; tmp < nor1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line5.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line5.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+", indexInfo);
				norAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		
		if(rcm1Num < 40)
		{	
			tmpAlignInfoStrStartPos = line6.find("\t", 0) + 1;
			for(int tmp = 0; tmp < rcm1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line6.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line6.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-", indexInfo);
				rcmAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}	

		if(!SE_or_PE_bool)
		{	
			//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
			//cout << "norNum_2: " << nor2Num << endl;
			if(nor2Num < 40)
			{
				//cout << "line7: " << line7 << endl;
				tmpAlignInfoStrStartPos = line7.find("\t", 0) + 1;
				for(int tmp = 0; tmp < nor2Num; tmp++)
				{
					//cout << "start at:" << tmpAlignInfoStrStartPos << " end at:" 
					//	<< tmpAlignInfoStrEndPos << endl;
					tmpAlignInfoStrEndPos = line7.find("\t", tmpAlignInfoStrStartPos) - 1;
					tmpAlignInfoStr = line7.substr(tmpAlignInfoStrStartPos, 
						tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
					//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
					Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+", indexInfo);
					norAlignmentInfo_PE_2.push_back(tmpAlignInfo);

					tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
				}
			}
			//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
			//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
			//cout << "rcm2Num: " << rcm2Num << endl;
			if(rcm2Num < 40)
			{
				tmpAlignInfoStrStartPos = line8.find("\t", 0) + 1;
				for(int tmp = 0; tmp < rcm2Num; tmp++)
				{
					tmpAlignInfoStrEndPos = line8.find("\t", tmpAlignInfoStrStartPos) - 1;
					tmpAlignInfoStr = line8.substr(tmpAlignInfoStrStartPos, 
						tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
					//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
					Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-", indexInfo);
					rcmAlignmentInfo_PE_2.push_back(tmpAlignInfo);

					tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
				}
			}
			//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
		}
	}


	PE_Read_Alignment_Info(const string& line5, const string& line6, const string& line7,
		const string& line8, int nor1Num, int rcm1Num, int nor2Num, int rcm2Num, Index_Info* indexInfo)
	{
		highestPairAlignmentScore = 0;
		//Alignment_Info* tmpAlignmentInfo;
		//int tmp 
		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Rcm1 = -1;
		repeatRegion_index_Nor2 = -1;
		repeatRegion_index_Rcm2 = -1;


		int tmpAlignInfoStrStartPos;
		int tmpAlignInfoStrEndPos;
		string tmpAlignInfoStr;

		norAlignmentInfo_PE_1.clear();
		rcmAlignmentInfo_PE_1.clear();
		norAlignmentInfo_PE_2.clear();
		rcmAlignmentInfo_PE_2.clear();		
		//cout << 
		if(nor1Num < 40)
		{
			tmpAlignInfoStrStartPos = line5.find("\t", 0) + 1;
			for(int tmp = 0; tmp < nor1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line5.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line5.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+", indexInfo);
				norAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		
		if(rcm1Num < 40)
		{	
			tmpAlignInfoStrStartPos = line6.find("\t", 0) + 1;
			for(int tmp = 0; tmp < rcm1Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line6.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line6.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-", indexInfo);
				rcmAlignmentInfo_PE_1.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}	
		//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
		//cout << "norNum_2: " << nor2Num << endl;
		if(nor2Num < 40)
		{
			//cout << "line7: " << line7 << endl;
			tmpAlignInfoStrStartPos = line7.find("\t", 0) + 1;
			for(int tmp = 0; tmp < nor2Num; tmp++)
			{
				//cout << "start at:" << tmpAlignInfoStrStartPos << " end at:" 
				//	<< tmpAlignInfoStrEndPos << endl;
				tmpAlignInfoStrEndPos = line7.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line7.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "+", indexInfo);
				norAlignmentInfo_PE_2.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		//cout << "norAlignmentInfo_PE_2.size():  " << norAlignmentInfo_PE_2.size() << endl;
		//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
		//cout << "rcm2Num: " << rcm2Num << endl;
		if(rcm2Num < 40)
		{
			tmpAlignInfoStrStartPos = line8.find("\t", 0) + 1;
			for(int tmp = 0; tmp < rcm2Num; tmp++)
			{
				tmpAlignInfoStrEndPos = line8.find("\t", tmpAlignInfoStrStartPos) - 1;
				tmpAlignInfoStr = line8.substr(tmpAlignInfoStrStartPos, 
					tmpAlignInfoStrEndPos - tmpAlignInfoStrStartPos + 1);
				//cout << "tmpAlignInfoStr: " << tmpAlignInfoStr << endl;
				Alignment_Info* tmpAlignInfo = new Alignment_Info(tmpAlignInfoStr, "-", indexInfo);
				rcmAlignmentInfo_PE_2.push_back(tmpAlignInfo);

				tmpAlignInfoStrStartPos = tmpAlignInfoStrEndPos + 2;
			}
		}
		//cout << "rcmAlignmentInfo_PE_2.size():  " << rcmAlignmentInfo_PE_2.size() << endl;
	}

	void pushBackPathInfo2PeAlignInfo(Path_Info* newPathInfo, bool End1OrEnd2, 
		bool NorOrRcm, Index_Info* indexInfo) // Note: End1OrEnd2, NorOrRcm are both describing the other end alignment
	{
		//cout << "End1OrEnd2: " << End1OrEnd2 << endl;
		//cout << "NorOrRcm: " << NorOrRcm << endl;
		//cout << "PathSize: " << (newPathInfo->finalPathVec).size() << endl;
		if(End1OrEnd2 && NorOrRcm)
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->returnFinalPathVecSize(); tmpPath++)
			{
				int mapChromNameInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).first;
					= newPathInfo->returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).second;
					= newPathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
					= newPathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
					mapChromPosInt, (newPathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, 
					(newPathInfo->fixedPathMismatchPosVec)[tmpPath], 
					(newPathInfo->fixedPathMismatchCharVec)[tmpPath]);
				rcmAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			}			
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->returnFinalPathVecSize(); tmpPath++)
			{
				//int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				//int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				//int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
				int mapChromNameInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).first;
					= newPathInfo->returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).second;
					= newPathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
					= newPathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt),
					mapChromPosInt, (newPathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, 
					(newPathInfo->fixedPathMismatchPosVec)[tmpPath], 
					(newPathInfo->fixedPathMismatchCharVec)[tmpPath]);
				norAlignmentInfo_PE_2.push_back(tmpAlignmentInfo);
			}		
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->returnFinalPathVecSize(); tmpPath++)
			{
				//int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				//int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				//int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
				int mapChromNameInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).first;
					= newPathInfo->returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).second;
					= newPathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
					= newPathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->returnChrNameStr(mapChromNameInt), 
					mapChromPosInt, (newPathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, 
					(newPathInfo->fixedPathMismatchPosVec)[tmpPath], 
					(newPathInfo->fixedPathMismatchCharVec)[tmpPath]);
				rcmAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}		
		}
		else
		{
			for(int tmpPath = 0; tmpPath < newPathInfo->returnFinalPathVecSize(); tmpPath++)
			{
				//int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				//int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				//int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
				int mapChromNameInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).first;
					= newPathInfo->returnMapChrNameIntInFinalPathVec(tmpPath);
				int mapChromPosInt //= (((newPathInfo->finalPathVec)[tmpPath]).first).second;
					= newPathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
				int tmpMismatch //= (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;
					= newPathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);				
				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->returnChrNameStr(mapChromNameInt),  
					mapChromPosInt, (newPathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
					tmpMismatch, indexInfo, 
					(newPathInfo->fixedPathMismatchPosVec)[tmpPath], 
					(newPathInfo->fixedPathMismatchCharVec)[tmpPath]);
				norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}		
		}

	}

	Alignment_Info* fixHeadTail_getAlignInfo(
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else //if(replacedAlignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
	}

	bool fixHeadTail_alignInfo_headSoftClippingOrNot(
		int replacedAlignInfoKind, int replacedAlignInfoNO)
	// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		bool headSoftClippingOrNot_bool = false;
		if(replacedAlignInfoKind == 1)
		{
			headSoftClippingOrNot_bool 
				= (((norAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode)[0].type == "S");
		}
		else if(replacedAlignInfoKind == 2)
		{
			headSoftClippingOrNot_bool 
				= (((rcmAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode)[0].type == "S");  
		}
		else if(replacedAlignInfoKind == 3)
		{
			headSoftClippingOrNot_bool 
				= (((norAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode)[0].type == "S");
		}
		else //if(replacedAlignInfoKind == 4)
		{
			headSoftClippingOrNot_bool 
				= (((rcmAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode)[0].type == "S");
		}
		return headSoftClippingOrNot_bool;
	}

	bool fixHeadTail_alignInfo_tailSoftClippingOrNot(
		int replacedAlignInfoKind, int replacedAlignInfoNO)
	// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		bool tailSoftClippingOrNot_bool = false;
		if(replacedAlignInfoKind == 1)
		{
			int cigarStringJumpCodeSize 
				= ((norAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode).size();
			tailSoftClippingOrNot_bool 
				= (((norAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type 
					== "S");
		}
		else if(replacedAlignInfoKind == 2)
		{
			int cigarStringJumpCodeSize 
				= ((rcmAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode).size();			
			tailSoftClippingOrNot_bool 
				= (((rcmAlignmentInfo_PE_1[replacedAlignInfoNO])->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type 
					== "S");  
		}
		else if(replacedAlignInfoKind == 3)
		{
			int cigarStringJumpCodeSize 
				= ((norAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode).size();			
			tailSoftClippingOrNot_bool 
				= (((norAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type 
					== "S");
		}
		else //if(replacedAlignInfoKind == 4)
		{
			int cigarStringJumpCodeSize 
				= ((rcmAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode).size();
			tailSoftClippingOrNot_bool 
				= (((rcmAlignmentInfo_PE_2[replacedAlignInfoNO])->cigarStringJumpCode)[cigarStringJumpCodeSize-1].type 
					== "S");
		}
		return tailSoftClippingOrNot_bool;
	}

	Alignment_Info* getAlignInfo(
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			return norAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 2)
		{
			return rcmAlignmentInfo_PE_1[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 3)
		{
			return norAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else if(replacedAlignInfoKind == 4)
		{
			return rcmAlignmentInfo_PE_2[replacedAlignInfoNO];
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfo()" << endl;
		}
	}

	void replaceWithNewAlignInfo(int tmpAlignInfoType, 
		int tmpIndex_peAlignInfo, Alignment_Info* newAlignInfo)
	{
		/*
		if(tmpAlignInfoType == 1)
		{
			//delete norAlignmentInfo_PE_1[tmpIndex_peAlignInfo];
			norAlignmentInfo_PE_1[tmpIndex_peAlignInfo]->copyAlignInfo(newAlignInfo); //= newAlignInfo;
		}
		else if(tmpAlignInfoType == 2)
		{
			//delete rcmAlignmentInfo_PE_1[tmpIndex_peAlignInfo];
			rcmAlignmentInfo_PE_1[tmpIndex_peAlignInfo]->copyAlignInfo(newAlignInfo); //= newAlignInfo;
		}
		else if(tmpAlignInfoType == 3)
		{		
			//delete norAlignmentInfo_PE_2[tmpIndex_peAlignInfo];	
			norAlignmentInfo_PE_2[tmpIndex_peAlignInfo]->copyAlignInfo(newAlignInfo); //= newAlignInfo;
		}
		else
		{		
			//delete rcmAlignmentInfo_PE_2[tmpIndex_peAlignInfo];
			rcmAlignmentInfo_PE_2[tmpIndex_peAlignInfo]->copyAlignInfo(newAlignInfo); //= newAlignInfo;
		}*/
		this->pushBackNewAlignInfo(tmpAlignInfoType, newAlignInfo);
	}

	void pushBackNewAlignInfo(int tmpAlignInfoType, Alignment_Info* newAlignInfo)
	{
		if(tmpAlignInfoType == 1)
		{
			norAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		else if(tmpAlignInfoType == 2)
		{
			rcmAlignmentInfo_PE_1.push_back(newAlignInfo);
		}
		else if(tmpAlignInfoType == 3)
		{
			norAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
		else
		{
			rcmAlignmentInfo_PE_2.push_back(newAlignInfo);
		}
	}

	/*
	void incompleteHead_replaceAndPushBackPathInfo2PeAlignInfo(
		Path_Info* pathInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo,
		int midPartSegLength, Index_Info* indexInfo, Alignment_Info* oriAlignInfo)
	{
		//Alignment_Info* oriAlignInfo = (this->returnAlignInfoInPeAlignInfo(
		//		tmpAlignInfoType, tmpIndex_peAlignInfo));
		int oldMismatch = oriAlignInfo->returnMismatchNum();
		for(int tmpPath = 0; tmpPath < pathInfo->returnFinalPathVecSize(); tmpPath++ )
		{
			//int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt //= (((pathInfo->finalPathVec)[tmpPath]).first).second;
				= pathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= (pathInfo->fixedPathMismatchVec)[tmpPath];
				= pathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
			int newMismatch = tmpMismatch + oldMismatch;

			//vector<int> newMismatchPosVec;
			//vector<char> newMisamtchCharVec;
			//pathInfo->

			Alignment_Info* tmpNewAlignInfo = new Alignment_Info();
			tmpNewAlignInfo = oriAlignInfo->newAlignInfo_addIncompleteHead(
				pathInfo, tmpPath, mapChromPosInt, newMismatch, indexInfo);

			if(tmpPath == 0)
			{
				this->replaceWithNewAlignInfo(
					tmpAlignInfoType, tmpIndex_peAlignInfo,
					tmpNewAlignInfo);
			}
			else
			{
				this->pushBackNewAlignInfo(
					tmpAlignInfoType, tmpNewAlignInfo);
			}
		}
	}

	void incompleteTail_replaceAndPushBackPathInfo2PeAlignInfo(
		Path_Info* pathInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo,
		int midPartSegLength, Index_Info* indexInfo, Alignment_Info* oriAlignInfo, int midPartStartLocInRead)
	{
		//Alignment_Info* oriAlignInfo = (this->returnAlignInfoInPeAlignInfo(
		//		tmpAlignInfoType, tmpIndex_peAlignInfo));
		int oldMismatch = oriAlignInfo->returnMismatchNum();
		int mapChromPosInt = oriAlignInfo->returnAlignChromPos();
		for(int tmpPath = 0; tmpPath < pathInfo->returnFinalPathVecSize(); tmpPath++ )
		{
			//int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			//int mapChromPosInt = (((pathInfo->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch //= (pathInfo->fixedPathMismatchVec)[tmpPath];
				= pathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
			int newMismatch = tmpMismatch + oldMismatch;

			Alignment_Info* tmpNewAlignInfo = new Alignment_Info();
			//Alignment_Info* 
				tmpNewAlignInfo = oriAlignInfo->newAlignInfo_addIncompleteTail(
					pathInfo, tmpPath, mapChromPosInt, newMismatch, indexInfo, midPartStartLocInRead);

			// stop 5
			if(tmpPath == 0)
			{
				this->replaceWithNewAlignInfo(
					tmpAlignInfoType, tmpIndex_peAlignInfo,
					tmpNewAlignInfo);
				//delete tmpNewAlignInfo;
			}//stop 6
			else
			{
				this->pushBackNewAlignInfo(
					tmpAlignInfoType, tmpNewAlignInfo);
			}
		}
	}*/


	void incompleteHead_replaceAndPushBackPathInfo2PeAlignInfo_new(
		Path_Info* pathInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo,
		int midPartSegLength, Index_Info* indexInfo, Alignment_Info* oriAlignInfo)
	{
		//Alignment_Info* oriAlignInfo = (this->returnAlignInfoInPeAlignInfo(
		//		tmpAlignInfoType, tmpIndex_peAlignInfo));
		int oldMismatch = oriAlignInfo->returnMismatchNum();
		for(int tmpPath = 0; tmpPath < pathInfo->returnFinalPathVecSize(); tmpPath++ )
		{
			//int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt //= (((pathInfo->finalPathVec)[tmpPath]).first).second;
				= pathInfo->returnMapChrPosIntInFinalPathVec(tmpPath);
			int tmpMismatch //= (pathInfo->fixedPathMismatchVec)[tmpPath];
				= pathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
			int newMismatch = tmpMismatch + oldMismatch;

			Alignment_Info* tmpNewAlignInfo = new Alignment_Info();
			tmpNewAlignInfo->newAlignInfo_addIncompleteHead_new(
				pathInfo, tmpPath, mapChromPosInt, newMismatch, indexInfo, oriAlignInfo);

			if(tmpPath == 0)
			{
				this->replaceWithNewAlignInfo(
					tmpAlignInfoType, tmpIndex_peAlignInfo,
					tmpNewAlignInfo);
			}
			else
			{
				this->pushBackNewAlignInfo(
					tmpAlignInfoType, tmpNewAlignInfo);
			}
		}
	}

	void incompleteTail_replaceAndPushBackPathInfo2PeAlignInfo_new(
		Path_Info* pathInfo, int tmpAlignInfoType, int tmpIndex_peAlignInfo,
		int midPartSegLength, Index_Info* indexInfo, Alignment_Info* oriAlignInfo, int midPartStartLocInRead)
	{
		//Alignment_Info* oriAlignInfo = (this->returnAlignInfoInPeAlignInfo(
		//		tmpAlignInfoType, tmpIndex_peAlignInfo));
		int oldMismatch = oriAlignInfo->returnMismatchNum();
		int mapChromPosInt = oriAlignInfo->returnAlignChromPos();
		for(int tmpPath = 0; tmpPath < pathInfo->returnFinalPathVecSize(); tmpPath++ )
		{
			//int mapChromNameInt = (((pathInfo_nor1->finalPathVec)[tmpPath]).first).first;
			//int mapChromPosInt = (((pathInfo->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch //= (pathInfo->fixedPathMismatchVec)[tmpPath];
				= pathInfo->returnMismatchNumInFixedPathMismatchVec(tmpPath);
			int newMismatch = tmpMismatch + oldMismatch;

			Alignment_Info* tmpNewAlignInfo = new Alignment_Info();
			tmpNewAlignInfo->newAlignInfo_addIncompleteTail_new(
				pathInfo, tmpPath, mapChromPosInt, newMismatch, indexInfo, midPartStartLocInRead, oriAlignInfo);

			// stop 5
			if(tmpPath == 0)
			{
				this->replaceWithNewAlignInfo(
					tmpAlignInfoType, tmpIndex_peAlignInfo,
					tmpNewAlignInfo);
				//delete tmpNewAlignInfo;
			}//stop 6
			else
			{
				this->pushBackNewAlignInfo(
					tmpAlignInfoType, tmpNewAlignInfo);
			}
		}
	}



	void fixHeadTail_replaceWithNewAlignInfo(Alignment_Info* newAlignmentInfo, 
		int replacedAlignInfoKind, int replacedAlignInfoNO)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(replacedAlignInfoKind == 1)
		{
			delete norAlignmentInfo_PE_1[replacedAlignInfoNO];
			norAlignmentInfo_PE_1[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 2)
		{
			delete rcmAlignmentInfo_PE_1[replacedAlignInfoNO];
			rcmAlignmentInfo_PE_1[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 3)
		{
			delete norAlignmentInfo_PE_2[replacedAlignInfoNO];
			norAlignmentInfo_PE_2[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else if(replacedAlignInfoKind == 4)
		{
			delete rcmAlignmentInfo_PE_2[replacedAlignInfoNO];
			rcmAlignmentInfo_PE_2[replacedAlignInfoNO] = newAlignmentInfo;
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_replaceWithNewAlignInfo()" << endl;
		}
		//this->fixHeadTail_addNewAlignInfo(newAlignmentInfo, replacedAlignInfoKind);
	}

	void fixHeadTail_addNewAlignInfo(Alignment_Info* newAlignmentInfo, 
		int addedAlignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		if(addedAlignInfoKind == 1)
		{
			norAlignmentInfo_PE_1.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 2)
		{
			rcmAlignmentInfo_PE_1.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 3)
		{
			norAlignmentInfo_PE_2.push_back(newAlignmentInfo);
		}
		else if(addedAlignInfoKind == 4)
		{
			rcmAlignmentInfo_PE_2.push_back(newAlignmentInfo);
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_addNewAlignInfo()" << endl;
		}
		//delete newAlignmentInfo;
	}

	int fixHeadTail_getAlignInfoVecSize(int alignInfoKind)// 1. Nor_1; 2. Rcm_1; 3. Nor_2; 4. Rcm_2
	{
		int size = 100;
		if(alignInfoKind == 1)
		{
			size = norAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 2)
		{
			size = rcmAlignmentInfo_PE_1.size();
		}
		else if(alignInfoKind == 3)
		{
			size = norAlignmentInfo_PE_2.size();
		}
		else if(alignInfoKind == 4)
		{
			size = rcmAlignmentInfo_PE_2.size();
		}
		else
		{
			cout << "error in PE_Align_Info.fixHeadTail_getAlignInfoVecSize()" << endl;
		}		
		return size;
	}


	void addHeadInfo2OriAligmentInfo(vector<Alignment_Info*>& alignInfoVec, int alignInfoVecNO,
		int firstMatchLength, int spliceJunctionDistance)
	{

	}

	void getMismatchNumForEveryAlignment_SE()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		}		
	}

	void getMismatchNumForEveryAlignment_end1()
	{
		this->getMismatchNumForEveryAlignment_SE();
	}

	void getMismatchNumForEveryAlignment_end2()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		}
	}

	void getMismatchNumForEveryAlignment()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		}
	}

	void getMismatchNumForEveryAlignment_Nor1Rcm2()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		}
		// for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		// {
		// 	rcmAlignmentInfo_PE_1[tmp]->generateMismatchNum();
		// }
		// for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		// {
		// 	norAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		// }
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->generateMismatchNum();
		}
	}

	void getEndMatchPosForEveryAlignment_SE()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
	}

	void getEndMatchPosForEveryAlignment_end1()
	{
		this->getEndMatchPosForEveryAlignment_SE();
	}

	void getEndMatchPosForEveryAlignment_end2()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
	}

	void getEndMatchPosForEveryAlignment()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			rcmAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			norAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
	}	

	void getEndMatchPosForEveryAlignment_Nor1Rcm2()
	{
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			norAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		}
		// for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		// {
		// 	rcmAlignmentInfo_PE_1[tmp]->getEndMatchedPosInChr();
		// }
		// for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		// {
		// 	norAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		// }
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			rcmAlignmentInfo_PE_2[tmp]->getEndMatchedPosInChr();
		}
	}	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output AlignInfo functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output AlignInfo functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////


	string getTmpAlignInfo(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq,
		int multiMapSeg_maxLength)
	{

		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(norAlignmentInfo_PE_1.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_1.size()) + "\t" + int_to_str(multiMapSeg_maxLength)
			+ "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(norAlignmentInfo_PE_2.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;
		string tmpNorAlignInfo_1;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			if(norAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ","; //+ "\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";
		}
		string tmpRcmAlignInfo_1;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			if(norAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";				
		}
		string tmpRcmAlignInfo_2;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}	

		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}

	string getTmpAlignInfo(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq)
	{
		int multiMapSeg_maxLength = 0;
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(norAlignmentInfo_PE_1.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_1.size()) + "\t" + int_to_str(multiMapSeg_maxLength)
			+ "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(norAlignmentInfo_PE_2.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;
		string tmpNorAlignInfo_1;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			if(norAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ","; //+ "\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";
		}
		string tmpRcmAlignInfo_1;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			if(norAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";				
		}
		string tmpRcmAlignInfo_2;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}	

		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}

	string getTmpAlignInfoForFinalPair(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq, 
		int multiMapSeg_maxLength)
	{
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t" + int_to_str(multiMapSeg_maxLength)
			+ "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;

		string tmpNorAlignInfo_1;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}

			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].first;
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";				
		}
		string tmpRcmAlignInfo_1;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";			
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}		
		string tmpNorAlignInfo_2;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].first;
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";		
		}
		string tmpRcmAlignInfo_2;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}		
		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}

	string getTmpAlignInfoForFinalPair(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq)
	{
		int multiMapSeg_maxLength = 0;
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t" + int_to_str(multiMapSeg_maxLength)
			+ "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;

		string tmpNorAlignInfo_1;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}

			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].first;
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";				
		}
		string tmpRcmAlignInfo_1;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";			
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}		
		string tmpNorAlignInfo_2;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].first;
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";		
		}
		string tmpRcmAlignInfo_2;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}		
		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}
	/////////////////////  output directly functions ///////////////////////////

	void outputTmpAlignInfo(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, ofstream* output, bool FastaOrFastq)
	{
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(norAlignmentInfo_PE_1.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(norAlignmentInfo_PE_2.size()) + "\t"
			+ int_to_str(rcmAlignmentInfo_PE_2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;
		string tmpNorAlignInfo_1;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_1.size(); tmp++)
		{
			if(norAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";		
		}
		string tmpRcmAlignInfo_1;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_1.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_1.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmp = 0; tmp < norAlignmentInfo_PE_2.size(); tmp++)
		{
			if(norAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";		
		}
		string tmpRcmAlignInfo_2;
		for(int tmp = 0; tmp < rcmAlignmentInfo_PE_2.size(); tmp++)
		{
			if(rcmAlignmentInfo_PE_2.size() >= 40)
			{
				break;
			}
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}	

		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		*output << tmpAlignInfoStr << endl;
		//return tmpAlignInfoStr;
	}

	void outputTmpAlignInfoForFinalPair(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, ofstream* output, bool FastaOrFastq)
	{
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(finalAlignPair_Nor2Rcm1.size()) + "\t"
			+ int_to_str(finalAlignPair_Nor1Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		Alignment_Info* tmpAlignInfo;

		string tmpNorAlignInfo_1;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}

			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].first;
			tmpAlignInfo = norAlignmentInfo_PE_1[tmp];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_1 += "\t";		
		}
		string tmpRcmAlignInfo_1;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_1[tmp];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";			
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_1 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_1 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_1 += "\t";
		}		
		string tmpNorAlignInfo_2;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < finalAlignPair_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(finalAlignPair_Nor2Rcm1.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor2Rcm1[tmpNor2Rcm1].first;
			tmpAlignInfo = norAlignmentInfo_PE_2[tmp];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpNorAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpNorAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpNorAlignInfo_2 += "\t";		
		}
		string tmpRcmAlignInfo_2;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < finalAlignPair_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(finalAlignPair_Nor1Rcm2.size() >= 40)
				{break;}
			int tmp = finalAlignPair_Nor1Rcm2[tmpNor1Rcm2].second;
			tmpAlignInfo = rcmAlignmentInfo_PE_2[tmp];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ tmpAlignInfo->returnAlignChromName() + ","
				+ int_to_str(tmpAlignInfo->returnAlignChromPos()) + ","
				+ tmpAlignInfo->jumpCodeVec2Str() + ","
				+ int_to_str(tmpAlignInfo->returnMismatchNum()) + ",";//\t";
			if(STORE_MISMATCH_POS)
			{
				int tmpMismatchPosVecSize = tmpAlignInfo->returnMismatchPosVecSize();
				if(tmpMismatchPosVecSize >= 1)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(0));
				}
				for(int tmp2 = 1; tmp2 < tmpMismatchPosVecSize; tmp2++)
				{
					tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(tmpAlignInfo->returnMismatchPosVecValue(tmp2));
				}
				tmpRcmAlignInfo_2 += ",";
				if(STORE_MISMATCH_CHA)
				{
					int tmpMismatchCharVecSize = tmpAlignInfo->returnMismatchCharVecSize();
					if(tmpMismatchCharVecSize >= 1)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (tmpAlignInfo->returnMismatchCharVecValue(0));
					}
					for(int tmp3 = 1; tmp3 < tmpMismatchCharVecSize; tmp3++)
					{
						tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (tmpAlignInfo->returnMismatchCharVecValue(tmp3));
					}
					tmpRcmAlignInfo_2 += ",";
				}
			}
			else
			{}
			tmpRcmAlignInfo_2 += "\t";	
		}		
		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		*output << tmpAlignInfoStr << endl;
		//return tmpAlignInfoStr;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output SAM functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// Output SAM functions For Fasta Reads ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	string getSAMformatForFinalPair_secondaryOrNot(PE_Read_Info& peReadInfo,
		bool FastaOrFastq, int multiMapSeg_maxLength)
		/* after pairing candidate alignments and choosing the best pairs, according to
		 finalAlignPair_Nor1Rcm2 and finalAlignPair_Nor2Rcm1;*/
	{
		//return "paired";

 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();
		string tmpSamStr;

		int IH_Nor1Rcm2 = finalAlignPair_Nor1Rcm2.size();
		int IH_Nor2Rcm1 = finalAlignPair_Nor2Rcm1.size();

		int IH_allPair = IH_Nor1Rcm2 + IH_Nor2Rcm1;

		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp++)
		{

			int HI_Nor1Rcm2_tmp = tmp + 1;
			int tmpNor1NO = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcm2NO = finalAlignPair_Nor1Rcm2[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmpNor1NO];
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmpRcm2NO];

			int template_start = tmpAlignInfo_1->returnAlignChromPos();
			int template_end = tmpAlignInfo_2->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_1, readSeq_1);
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(), 
					tmpAlignInfo_2, true, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length, FastaOrFastq,
					multiMapSeg_maxLength);

			//outputFile << tmpSamStr << endl; 
			tmpSamStr += "\n";

			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
					tmpAlignInfo_1, false, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length, FastaOrFastq,
					multiMapSeg_maxLength);

			tmpSamStr += "\n";
		}
		
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp++)
		{

			int HI_Nor2Rcm1_tmp = tmp + 1
				+ IH_Nor1Rcm2;

			int tmpNor2NO = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcm1NO = finalAlignPair_Nor2Rcm1[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmpNor2NO];		
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmpRcm1NO];
			
			int template_start = tmpAlignInfo_1->returnAlignChromPos();
			int template_end = tmpAlignInfo_2->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;		

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_2, readSeq_2);
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
					tmpAlignInfo_2, false, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
					template_length, FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\n";
			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
					tmpAlignInfo_1, true, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
					template_length, FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\n";
		}
		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;
	}	

	string getSAMformatForFinalPair_secondaryOrNot(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
		/* after pairing candidate alignments and choosing the best pairs, according to
		 finalAlignPair_Nor1Rcm2 and finalAlignPair_Nor2Rcm1;*/
	{
		//return "paired";
		int multiMapSeg_maxLength = 0;
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();
		string tmpSamStr;

		int IH_Nor1Rcm2 = finalAlignPair_Nor1Rcm2.size();
		int IH_Nor2Rcm1 = finalAlignPair_Nor2Rcm1.size();

		int IH_allPair = IH_Nor1Rcm2 + IH_Nor2Rcm1;

		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp++)
		{

			int HI_Nor1Rcm2_tmp = tmp + 1;
			int tmpNor1NO = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcm2NO = finalAlignPair_Nor1Rcm2[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_1[tmpNor1NO];
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_2[tmpRcm2NO];

			int template_start = tmpAlignInfo_1->returnAlignChromPos();
			int template_end = tmpAlignInfo_2->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_1, readSeq_1);
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(), 
					tmpAlignInfo_2, true, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length, FastaOrFastq,
					multiMapSeg_maxLength);

			//outputFile << tmpSamStr << endl; 
			tmpSamStr += "\n";

			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
					tmpAlignInfo_1, false, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length, FastaOrFastq,
					multiMapSeg_maxLength);

			tmpSamStr += "\n";
		}
		
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp++)
		{

			int HI_Nor2Rcm1_tmp = tmp + 1
				+ IH_Nor1Rcm2;

			int tmpNor2NO = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcm1NO = finalAlignPair_Nor2Rcm1[tmp].second;
			Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmpNor2NO];		
			Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmpRcm1NO];
			
			int template_start = tmpAlignInfo_1->returnAlignChromPos();
			int template_end = tmpAlignInfo_2->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;		

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_2, readSeq_2);
				= tmpSamStr + tmpAlignInfo_1->getSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
					tmpAlignInfo_2, false, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
					template_length, FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\n";
			tmpSamStr 
				= tmpSamStr + tmpAlignInfo_2->getSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
					tmpAlignInfo_1, true, IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
					template_length, FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\n";
		}
		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asUnmapped(
			peReadInfo, FastaOrFastq);
		//return this->getSAMformatForBothEndsUnmapped(peReadInfo, FastaOrFastq);
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot_asUnmapped(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		return this->getSAMformatForBothEndsUnmapped(peReadInfo, FastaOrFastq);
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
		PE_Read_Info& peReadInfo, bool FastaOrFastq, int multiMapSeg_maxLength)
	{
		int IH_Nor1 = unpairedSEalignVec_final_Nor_end1.size();
		int IH_Rcm1 = unpairedSEalignVec_final_Rcm_end1.size();
		int IH_Nor2 = unpairedSEalignVec_final_Nor_end2.size();
		int IH_Rcm2 = unpairedSEalignVec_final_Rcm_end2.size();
		if((IH_Nor1 + IH_Rcm1 == 1)&&(IH_Nor2 + IH_Rcm2 == 1))
			return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asMapped(
						peReadInfo, FastaOrFastq, multiMapSeg_maxLength);
		else
			return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asUnmapped(
						peReadInfo, FastaOrFastq);
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
		PE_Read_Info& peReadInfo, bool FastaOrFastq)
	{
		int multiMapSeg_maxLength = 0;
		int IH_Nor1 = unpairedSEalignVec_final_Nor_end1.size();
		int IH_Rcm1 = unpairedSEalignVec_final_Rcm_end1.size();
		int IH_Nor2 = unpairedSEalignVec_final_Nor_end2.size();
		int IH_Rcm2 = unpairedSEalignVec_final_Rcm_end2.size();
		if((IH_Nor1 + IH_Rcm1 == 1)&&(IH_Nor2 + IH_Rcm2 == 1))
			return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asMapped(
						peReadInfo, FastaOrFastq);
		else
			return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asUnmapped(
						peReadInfo, FastaOrFastq);
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
		PE_Read_Info& peReadInfo, bool FastaOrFastq, int multiMapSeg_maxLength)
	{
		return this->getSAMformatForUnpairedAlignments_secondaryOrNot_asMapped(
						peReadInfo, FastaOrFastq, multiMapSeg_maxLength);
	}


	string getSAMformatForUnpairedAlignments_secondaryOrNot_asMapped(PE_Read_Info& peReadInfo, 
		bool FastaOrFastq, int multiMapSeg_maxLength)
	{
		string tmpSamStr;
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		int IH_Nor1 = unpairedSEalignVec_final_Nor_end1.size();
		int IH_Rcm1 = unpairedSEalignVec_final_Rcm_end1.size();
		int IH_Nor2 = unpairedSEalignVec_final_Nor_end2.size();
		int IH_Rcm2 = unpairedSEalignVec_final_Rcm_end2.size();
		int IH_end1 = IH_Nor1 + IH_Rcm1;
		int IH_end2 = IH_Nor2 + IH_Rcm2;

		for(int tmp = 0; tmp < IH_Nor1; tmp++)
		{
			int tmpHI = tmp + 1;
			int tmpNor1NO = unpairedSEalignVec_final_Nor_end1[tmp];
			tmpSamStr = tmpSamStr 
				+ norAlignmentInfo_PE_1[tmpNor1NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(), 
					true, IH_end1, tmpHI, (IH_end2 > 0), (tmp != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";
		}
		for(int tmp = 0; tmp < IH_Rcm1; tmp++)
		{
			int tmpHI = IH_Nor1 + tmp + 1;
			int tmpRcm1NO = unpairedSEalignVec_final_Rcm_end1[tmp];
			tmpSamStr = tmpSamStr
				+ rcmAlignmentInfo_PE_1[tmpRcm1NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
					true, IH_end1, tmpHI, (IH_end2 > 0), ((tmp + IH_Nor1) != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";
		}
		if(IH_end1 == 0) // 1st end unmapped
		{
			tmpSamStr = tmpSamStr + readName_1 + "\t69\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1()
				+ "\t";
			if(FastaOrFastq)
				tmpSamStr = tmpSamStr + "*\tIH:i:0\tHI:i:0";
			else 
				tmpSamStr = tmpSamStr + peReadInfo.returnQualitySeq_1() + "\tIH:i:0\tHI:i:0";
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";			
		}

		for(int tmp = 0; tmp < IH_Nor2; tmp++)
		{
			int tmpHI = tmp + 1;
			int tmpNor2NO = unpairedSEalignVec_final_Nor_end2[tmp];
			tmpSamStr = tmpSamStr
				+ norAlignmentInfo_PE_2[tmpNor2NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
					false, IH_end2, tmpHI, (IH_end1 > 0), (tmp != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}
		for(int tmp = 0; tmp < IH_Rcm2; tmp++)
		{
			int tmpHI = IH_Nor2 + tmp + 1;
			int tmpRcm2NO = unpairedSEalignVec_final_Rcm_end2[tmp];
			tmpSamStr = tmpSamStr
				+ rcmAlignmentInfo_PE_2[tmpRcm2NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
					false, IH_end2, tmpHI, (IH_end1 > 0), ((tmp + IH_Nor2) != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}
		if(IH_end2 == 0) // 2nd end unmapped
		{
			tmpSamStr = tmpSamStr + readName_2 + "\t133\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1()
				+ "\t";
			if(FastaOrFastq)
				tmpSamStr = tmpSamStr + "*\tIH:i:0\tHI:i:0";
			else
				tmpSamStr = tmpSamStr + peReadInfo.returnQualitySeq_2() + "\tIH:i:0\tHI:i:0";
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}

		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;		
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot_asMapped(PE_Read_Info& peReadInfo, 
		bool FastaOrFastq)
	{
		int multiMapSeg_maxLength = 0;
		string tmpSamStr;
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		int IH_Nor1 = unpairedSEalignVec_final_Nor_end1.size();
		int IH_Rcm1 = unpairedSEalignVec_final_Rcm_end1.size();
		int IH_Nor2 = unpairedSEalignVec_final_Nor_end2.size();
		int IH_Rcm2 = unpairedSEalignVec_final_Rcm_end2.size();
		int IH_end1 = IH_Nor1 + IH_Rcm1;
		int IH_end2 = IH_Nor2 + IH_Rcm2;

		for(int tmp = 0; tmp < IH_Nor1; tmp++)
		{
			int tmpHI = tmp + 1;
			int tmpNor1NO = unpairedSEalignVec_final_Nor_end1[tmp];
			tmpSamStr = tmpSamStr 
				+ norAlignmentInfo_PE_1[tmpNor1NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(), 
					true, IH_end1, tmpHI, (IH_end2 > 0), (tmp != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";
		}
		for(int tmp = 0; tmp < IH_Rcm1; tmp++)
		{
			int tmpHI = IH_Nor1 + tmp + 1;
			int tmpRcm1NO = unpairedSEalignVec_final_Rcm_end1[tmp];
			tmpSamStr = tmpSamStr
				+ rcmAlignmentInfo_PE_1[tmpRcm1NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
					true, IH_end1, tmpHI, (IH_end2 > 0), ((tmp + IH_Nor1)!= 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";
		}
		if(IH_end1 == 0) // 1st end unmapped
		{
			tmpSamStr = tmpSamStr + readName_1 + "\t69\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1()
				+ "\t";
			if(FastaOrFastq)
				tmpSamStr = tmpSamStr + "*\tIH:i:0\tHI:i:0";
			else 
				tmpSamStr = tmpSamStr + peReadInfo.returnQualitySeq_1() + "\tIH:i:0\tHI:i:0";
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end2);
			tmpSamStr += "\n";			
		}

		for(int tmp = 0; tmp < IH_Nor2; tmp++)
		{
			int tmpHI = tmp + 1;
			int tmpNor2NO = unpairedSEalignVec_final_Nor_end2[tmp];
			tmpSamStr = tmpSamStr
				+ norAlignmentInfo_PE_2[tmpNor2NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
					false, IH_end2, tmpHI, (IH_end1 > 0), (tmp != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}
		for(int tmp = 0; tmp < IH_Rcm2; tmp++)
		{
			int tmpHI = IH_Nor2 + tmp + 1;
			int tmpRcm2NO = unpairedSEalignVec_final_Rcm_end2[tmp];
			tmpSamStr = tmpSamStr
				+ rcmAlignmentInfo_PE_2[tmpRcm2NO]->getSamFormatString_unpaired_secondaryOrNot_PEasSE(
					readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
					false, IH_end2, tmpHI, (IH_end1 > 0), ((tmp + IH_Nor2) != 0), FastaOrFastq, multiMapSeg_maxLength);
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}
		if(IH_end2 == 0) // 2nd end unmapped
		{
			tmpSamStr = tmpSamStr + readName_2 + "\t133\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1()
				+ "\t";
			if(FastaOrFastq)
				tmpSamStr = tmpSamStr + "*\tIH:i:0\tHI:i:0";
			else
				tmpSamStr = tmpSamStr + peReadInfo.returnQualitySeq_2() + "\tIH:i:0\tHI:i:0";
			tmpSamStr += "\tXI:i:";
			tmpSamStr += int_to_str(IH_end1);
			tmpSamStr += "\n";
		}

		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		return returnStr;		
	}

	string getSAMformatForBothEndsUnmapped(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();		

		string peAlignSamStr;
		peAlignSamStr = readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1() 
			+ "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0\n" +
			readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_2() 
			+ "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0";
		return peAlignSamStr;
	}

	string getSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
		PE_Read_Info& peReadInfo, RepeatRegion_Info* repeatRegionInfo, 
		Index_Info* indexInfo, bool FastaOrFastq)
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}
		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		bool mappedToRepeatRegionBool_Nor1 = (repeatRegion_index_Nor1 >= 0);
		bool mappedToRepeatRegionBool_Rcm1 = (repeatRegion_index_Rcm1 >= 0);
		bool mappedToRepeatRegionBool_Nor2 = (repeatRegion_index_Nor2 >= 0);
		bool mappedToRepeatRegionBool_Rcm2 = (repeatRegion_index_Rcm2 >= 0);

		string peAlignSamStr;

		if((!mappedToRepeatRegionBool_Nor1)&&(!mappedToRepeatRegionBool_Rcm1))
		{
			peAlignSamStr = peAlignSamStr + readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_1() + "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0\n";			
		}
		else
		{
			peAlignSamStr = peAlignSamStr + readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_1() + "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0";
			if(repeatRegion_index_Nor1 >= 0)
			{
				peAlignSamStr =  peAlignSamStr + "\tRf:i:" + int_to_str(repeatRegion_index_Nor1);
			}
			if(repeatRegion_index_Rcm1 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRr:i:" + int_to_str(repeatRegion_index_Rcm1);
			}
			peAlignSamStr += "\n";
		}

		if((repeatRegion_index_Nor2 < 0)&&(repeatRegion_index_Rcm2 < 0))
		{
			peAlignSamStr = peAlignSamStr + readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_2() + "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0\n";		
		}
		else
		{
			peAlignSamStr = peAlignSamStr + readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_2() + "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0";

			if(repeatRegion_index_Nor2 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRf:i:" + int_to_str(repeatRegion_index_Nor2);
			}
			if(repeatRegion_index_Rcm2 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRr:i:" + int_to_str(repeatRegion_index_Rcm2);
			}
			peAlignSamStr += "\n";
		}		

		return peAlignSamStr.substr(0,peAlignSamStr.length()-1);
	}

	//////////////////// output directly functions ////////////////////////
	/*
	void outputSAMformatForFinalPair_secondaryOrNot(
		PE_Read_Info& peReadInfo, ofstream* outputFile, bool FastaOrFastq)
		// after pairing candidate alignments and choosing the best pairs, according to
		// finalAlignPair_Nor1Rcm2 and finalAlignPair_Nor2Rcm1;
	{
		//return "paired";
  		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		int IH_Nor1Rcm2 = finalAlignPair_Nor1Rcm2.size();
		int IH_Nor2Rcm1 = finalAlignPair_Nor2Rcm1.size();

		int IH_allPair = IH_Nor1Rcm2 + IH_Nor2Rcm1;

		for(int tmp = 0; tmp < finalAlignPair_Nor1Rcm2.size(); tmp++)
		{
			int HI_Nor1Rcm2_tmp = tmp + 1;
			int tmpNor1NO = finalAlignPair_Nor1Rcm2[tmp].first;
			int tmpRcm2NO = finalAlignPair_Nor1Rcm2[tmp].second;

			int template_start = norAlignmentInfo_PE_1[tmpNor1NO]->returnAlignChromPos();
			int template_end = rcmAlignmentInfo_PE_2[tmpRcm2NO]->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;
			//(*outputFile) << "clear" << endl;
			//cout << "a" << endl;
			norAlignmentInfo_PE_1[tmpNor1NO]->outputSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(),
					rcmAlignmentInfo_PE_2[tmpRcm2NO], true, IH_allPair, 
					HI_Nor1Rcm2_tmp, (tmp != 0), template_length, outputFile, FastaOrFastq);
			//cout << "b" << endl;
			rcmAlignmentInfo_PE_2[tmpRcm2NO]->outputSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
					norAlignmentInfo_PE_1[tmpNor1NO], false, IH_allPair, 
					HI_Nor1Rcm2_tmp, (tmp != 0), template_length, outputFile, FastaOrFastq);
			//tmpSamStr += "\n";
		}
		
		for(int tmp = 0; tmp < finalAlignPair_Nor2Rcm1.size(); tmp++)
		{
			int HI_Nor2Rcm1_tmp = tmp + 1
				+ IH_Nor1Rcm2;

			int tmpNor2NO = finalAlignPair_Nor2Rcm1[tmp].first;
			int tmpRcm1NO = finalAlignPair_Nor2Rcm1[tmp].second;
			//Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmpNor2NO];		
			//Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmpRcm1NO];
			
			int template_start = norAlignmentInfo_PE_2[tmpNor2NO]->returnAlignChromPos();
			int template_end = rcmAlignmentInfo_PE_1[tmpRcm1NO]->returnEndMatchedPosInChr();

			int template_length = template_end - template_start + 1;		

			//tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_2, readSeq_2);
			//	= tmpSamStr + 
			//cout << "c" << endl;
			norAlignmentInfo_PE_2[tmpNor2NO]->outputSamFormatString_paired_secondaryOrNot(
					readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
					rcmAlignmentInfo_PE_1[tmpRcm1NO], false, 
					IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), template_length, 
					outputFile, FastaOrFastq);
			//tmpSamStr += "\n";
			//tmpSamStr 
			//	= tmpSamStr + 
			//cout << "d" << endl;
			rcmAlignmentInfo_PE_1[tmpRcm1NO]->outputSamFormatString_paired_secondaryOrNot(
					readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
					norAlignmentInfo_PE_2[tmpNor2NO], true, 
					IH_allPair, HI_Nor2Rcm1_tmp, ((tmp != 0)||(IH_Nor1Rcm2 > 0)), template_length, outputFile,
					FastaOrFastq);
			//tmpSamStr += "\n";
		}
		//string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		//return returnStr;
	}	

	void outputSAMformatForUnpairedAlignments_secondaryOrNot(
		PE_Read_Info& peReadInfo, ofstream* outputFile, bool FastaOrFastq)
	{
		this->outputSAMformatForBothEndsUnmapped(peReadInfo, outputFile, FastaOrFastq);
	}

	void outputSAMformatForBothEndsUnmapped(
		PE_Read_Info& peReadInfo, ofstream* outputFile, bool FastaOrFastq)//, 
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}

		*outputFile << peReadInfo.returnReadName_beforeSlash_1() 
			<< "\t77\t*\t0\t0\t*\t*\t0\t0\t" 
			<< peReadInfo.returnReadSeq_1()
			<< "\t" << qualitySeq_1 << "\tIH:i:0\tHI:i:0\n" 
			<< peReadInfo.returnReadName_beforeSlash_2()
			<< "\t141\t*\t0\t0\t*\t*\t0\t0\t" 
			<< peReadInfo.returnReadSeq_2()
			<< "\t" << qualitySeq_2 << "\tIH:i:0\tHI:i:0" << endl;
	}

	void outputSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
		PE_Read_Info& peReadInfo, RepeatRegion_Info* repeatRegionInfo, 
		Index_Info* indexInfo, ofstream* outputFile, bool FastaOrFastq)
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}

  		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		bool mappedToRepeatRegionBool_Nor1 = (repeatRegion_index_Nor1 >= 0);
		bool mappedToRepeatRegionBool_Rcm1 = (repeatRegion_index_Rcm1 >= 0);
		bool mappedToRepeatRegionBool_Nor2 = (repeatRegion_index_Nor2 >= 0);
		bool mappedToRepeatRegionBool_Rcm2 = (repeatRegion_index_Rcm2 >= 0);

		//string peAlignSamStr;

		if((!mappedToRepeatRegionBool_Nor1)&&(!mappedToRepeatRegionBool_Rcm1))
		{
			*outputFile << readName_1 << "\t77\t*\t0\t0\t*\t*\t0\t0\t" << 
				peReadInfo.returnReadSeq_1() << "\t" << qualitySeq_1 << "\tIH:i:0\tHI:i:0" << endl;			
		}
		else
		{
			*outputFile << readName_1 << "\t77\t*\t0\t0\t*\t*\t0\t0\t" << 
				peReadInfo.returnReadSeq_1() << "\t" << qualitySeq_1 << "\tIH:i:0\tHI:i:0";
			if(repeatRegion_index_Nor1 >= 0)
			{
				*outputFile << "\tRf:i:" << int_to_str(repeatRegion_index_Nor1);
			}
			if(repeatRegion_index_Rcm1 >= 0)
			{
				*outputFile << "\tRr:i:" << int_to_str(repeatRegion_index_Rcm1);
			}
			*outputFile << endl;
		}

		if((repeatRegion_index_Nor2 < 0)&&(repeatRegion_index_Rcm2 < 0))
		{
			*outputFile << readName_2 << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << 
				(peReadInfo.returnReadSeq_2()) << "\t" << qualitySeq_2 << "\tIH:i:0\tHI:i:0" << endl;		
		}
		else
		{
			*outputFile << readName_2 << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << 
				(peReadInfo.returnReadSeq_2()) << "\t" << qualitySeq_2 << "\tIH:i:0\tHI:i:0";
			if(repeatRegion_index_Nor2 >= 0)
			{
				*outputFile << "\tRf:i:" << int_to_str(repeatRegion_index_Nor2);
			}
			if(repeatRegion_index_Rcm2 >= 0)
			{
				*outputFile << "\tRr:i:" << int_to_str(repeatRegion_index_Rcm2);
			}
			*outputFile << endl;
		}		
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
};


#endif