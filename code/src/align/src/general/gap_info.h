// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
/* 
To fix:
	1. get # of mismatches: score_string, ...
	2. in fix_insertion: if secondMatchLength == 0, then .......
	3. extend back 2nd segment to avoid pendingSeq.length() > 31;
*/
#ifndef GAP_INFO_H
#define GAP_INFO_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;

class Gap_Info
{
//public:
private:
	int LengthOfSeqPerMismatchAllowed;
	int FixSpliceBuffer;
public:
	Gap_Info()
	{
		LengthOfSeqPerMismatchAllowed = 15;
		FixSpliceBuffer = FIX_SPLICE_BUFFER_MAX;
	}

	bool fixGapInPath_phase2(Path_Info* pathInfo, Seg_Info* segInfo, Index_Info* indexInfo, 
		const string& readSeq_inProcess, int readLength, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		//cout << "start fixGapInPath function..." << endl;
		bool fixGapInPathBool = false;
		int pathVecSize = pathInfo->returnPathVecSegSize();

		for(int tmpPathNO = 0; tmpPathNO < pathVecSize; tmpPathNO ++)
		{

			//cout << "...... start to fix path " << int_to_str(tmpPathNO + 1) << endl;
			if(!(pathInfo->returnPathValidBoolVec(tmpPathNO)))
			{
				//newPathFixed = false;
				pathInfo->pushBackPathFixedBoolVec(false);				
				continue;
			}

			//cout << "...... start to fix path " << int_to_str(tmpPathNO + 1) << endl;

			vector<int> newMismatchPosVec;
			vector<char> newMismatchCharVec;

			int firstSegGroupNO //= (pathInfo->PathVec_seg[tmpPathNO])[0].first;
				= pathInfo->returnSegGroupNOinPathVec_seg(tmpPathNO, 0);
			int firstSegLength = segInfo->returnSegmentLength(firstSegGroupNO);


			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(firstSegLength, "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);

			int newPathMismatchNum = 0;			
			bool newPathFixed = true;			
			int tmpPathSegSize //= pathInfo->PathVec_seg[tmpPathNO].size();
				= pathInfo->returnPathSegSizeInPathVec_seg(tmpPathNO);

			for(int tmpPathSegNO = 0; tmpPathSegNO < tmpPathSegSize-1; tmpPathSegNO++)
			{
				int tmpMismatchNum = 0;
				int tmpSegGroupNO //= (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].first;
					= pathInfo->returnSegGroupNOinPathVec_seg(tmpPathNO, tmpPathSegNO);
				int tmpSegCandiNO //= (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO].second;
					= pathInfo->returnSegCandiNOinPathVec_seg(tmpPathNO, tmpPathSegNO);

				int tmpSegGroupNO_next //= (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].first;
					= pathInfo->returnSegGroupNOinPathVec_seg(tmpPathNO, tmpPathSegNO + 1);
				int tmpSegCandiNO_next //= (pathInfo->PathVec_seg[tmpPathNO])[tmpPathSegNO+1].second;
					= pathInfo->returnSegCandiNOinPathVec_seg(tmpPathNO, tmpPathSegNO + 1);

				int tmpRelation = segInfo->checkSegRelation_phase2(tmpSegGroupNO, tmpSegCandiNO, 
					tmpSegGroupNO_next, tmpSegCandiNO_next, spliceJunctionDistanceMax);

				int tmpSegmentLocInRead_1 //= (segInfo->norSegmentLocInRead)[tmpSegGroupNO];
					= segInfo->returnSegmentLocInRead(tmpSegGroupNO);
				int tmpSegmentLocInRead_2 //= (segInfo->norSegmentLocInRead)[tmpSegGroupNO_next];
					= segInfo->returnSegmentLocInRead(tmpSegGroupNO_next);
				int tmpSegmentLength_1 //= (segInfo->norSegmentLength)[tmpSegGroupNO];
					= segInfo->returnSegmentLength(tmpSegGroupNO);
				int tmpSegmentLength_2 //= (segInfo->norSegmentLength)[tmpSegGroupNO_next];
					= segInfo->returnSegmentLength(tmpSegGroupNO_next);

				unsigned int tmpSegmentMapPosInWholeGenome_1 //= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO * CANDALILOC + tmpSegCandiNO);
					= segInfo->returnSegmentMapPos(tmpSegGroupNO, tmpSegCandiNO);
				unsigned int tmpSegmentMapPosInWholeGenome_2 //= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO_next * CANDALILOC + tmpSegCandiNO_next);
					= segInfo->returnSegmentMapPos(tmpSegGroupNO_next, tmpSegCandiNO_next);

				unsigned int tmpChrNameInt, tmpChrPosInt;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_1, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_1 = tmpChrPosInt;
				//cout << "...... tmpChrMapPos_1: " << tmpSegmentMapPos_1 << endl;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_2, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_2 = tmpChrPosInt;
				//cout << "...... tmpChrMapPos_2: " << tmpSegmentMapPos_2 << endl;
				
				string tmpChrNameStr;
				if(tmpChrNameStr_1 == tmpChrNameStr_2)
				{
					tmpChrNameStr = tmpChrNameStr_1;
				}
				else
				{
					newPathFixed = false;
					pathInfo->pushBackPathFixedBoolVec(newPathFixed);
					delete newPathSpliceInfo;
					break;
				}

				bool tmpDoubleAnchorFixed = fixDoubleAnchor_extendBack_phase2(newPathSpliceInfo, 
					tmpRelation, tmpSegmentLocInRead_1, tmpSegmentLocInRead_2,
					tmpSegmentLength_1, tmpSegmentLength_2, tmpSegmentMapPos_1, 
					tmpSegmentMapPos_2, readSeq_inProcess, indexInfo, tmpChrNameStr, &tmpMismatchNum,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, newMismatchPosVec, newMismatchCharVec);

				newPathMismatchNum = newPathMismatchNum + tmpMismatchNum;
				if(!tmpDoubleAnchorFixed)
				{
					if((tmpRelation == FIX_DELETION_GAP)||(tmpRelation == FIX_DELETION_NEIGHBOUR)
						||(tmpRelation == FIX_INSERTION_GAP)||(tmpRelation == FIX_INSERTION_NEIGHBOUR)
						||(tmpRelation == FIX_MATCH))
					{

					}
					newPathFixed = false;
					pathInfo->pushBackPathFixedBoolVec(newPathFixed);
					delete newPathSpliceInfo;
					break;					
				}
			}

			if(newPathFixed) // if newPathFixed = false, has been out of the loop above; newPathSpliceInfo has been deleted
			{
				newPathSpliceInfo->getFinalJumpCode();	
				bool allJumpCodeValidBool = newPathSpliceInfo->allFinalJumpCodeValid();
				if(allJumpCodeValidBool)
				{
					pathInfo->pushBackPathFixedBoolVec(allJumpCodeValidBool);
					pathInfo->pushBackFixedPathVec(tmpPathNO, newPathSpliceInfo);
					pathInfo->pushBackFixedPathMismatchVec(newPathMismatchNum);
					if(STORE_MISMATCH_POS)
					{
						pathInfo->fixedPathMismatchPosVec.push_back(newMismatchPosVec);	
						if(STORE_MISMATCH_CHA)
						{
							pathInfo->fixedPathMismatchCharVec.push_back(newMismatchCharVec);
						}
					}
				}
				else
				{
					delete newPathSpliceInfo;
					pathInfo->pushBackPathFixedBoolVec(allJumpCodeValidBool);
				}
			}
		}

		pathInfo->getFinalPath_extend2HeadTail_new(indexInfo, segInfo, readLength, readSeq_inProcess, Do_extendHeadTail);

		//cout << "finish getting finalPath" << endl;
		fixGapInPathBool = true;
		return fixGapInPathBool;
	}

	bool fixDoubleAnchor_extendBack_phase2(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, Index_Info* indexInfo, const string& chromName, int* mismatchNum,
		bool annotation_provided_bool, bool Do_annotation_only_bool, 
		Annotation_Info* annotationInfo, vector<int>& mismatchPosVec, vector<char>& mismatchCharVec
		)
	{
		//cout << "fixDoubleAnchor_extendBack starts!" << endl;
		bool fixDoubleAnchorBool = false;

		int chrNameInt = indexInfo->convertStringToInt(chromName);
		int extendBackNumMax = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		if(extendBackNumMax > segmentMapPos_2 - 1)
		{
			extendBackNumMax = segmentMapPos_2 - 1; 
		}
		int extendBackNum = extendBackInChromSeq(segmentLocInRead_2, readSeq_inProcess, 
			segmentMapPos_2, indexInfo->returnChromStr(chrNameInt), extendBackNumMax);

		segmentLocInRead_2 = segmentLocInRead_2 - extendBackNum;
		segmentLength_2 = segmentLength_2 + extendBackNum;
		segmentMapPos_2 = segmentMapPos_2 - extendBackNum;

		if((relation == FIX_TOO_CLOSE) || (relation == FIX_TOO_FAR) || (relation == FIX_NO_RELATIONSHIP))
		{
			//return false;
		}
		else if(relation == FIX_MATCH)
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Match(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum, mismatchPosVec, mismatchCharVec);
			//return fixDoubleAnchorBool;
		}
		else if((relation == FIX_INSERTION_NEIGHBOUR) || (relation == FIX_INSERTION_GAP))
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Insertion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum, mismatchPosVec, mismatchCharVec);
		}
		else if((relation == FIX_DELETION_NEIGHBOUR) || (relation == FIX_DELETION_GAP))
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Deletion(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum, mismatchPosVec, mismatchCharVec);
		}
		else if((relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_SPLICE_GAP))
		{
			fixDoubleAnchorBool = fixDoubleAnchor_Splice_phase2(cigarInfo, relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chromName, mismatchNum, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
				mismatchPosVec, mismatchCharVec);
		}
		else
		{
			cout << "error in fixDoubleAnchor ... " << endl;
		}
		return fixDoubleAnchorBool;
	}

	void insertMismatchPosVec_fixMatchInfo_SmallGap(int startPosInRead, int toProcessSeq_len, vector<int>& mismatchPosVec)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			int tmpMismatchPos = startPosInRead + tmp;
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixMatchInfo_SmallGap(int toProcessSeq_len, 
		int chrNameInt, int startPosInChrom, Index_Info* indexInfo, vector<char>& mismatchCharVec)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			char tmpMismatchChar = indexInfo->returnOneBaseCharInGenome(chrNameInt, startPosInChrom + tmp);
			//int tmpMismatchPos = startPosInRead + tmp;
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}
	void insertMismatchPosVec_fixMatchInfo(FixDoubleAnchor_Match_Info* fixMatchInfo, vector<int>& mismatchPosVec)
	{
		int vecSize = fixMatchInfo->returnMismatchPosInReadVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixMatchInfo->returnMismatchPosInRead(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixMatchInfo(FixDoubleAnchor_Match_Info* fixMatchInfo, vector<char>& mismatchCharVec)
	{
		int vecSize = fixMatchInfo->returnMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixMatchInfo->returnMismatchChar(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}	

	bool fixDoubleAnchor_Match(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum, vector<int>& mismatchPosVec, vector<char>& mismatchCharVec
		)
	{
		//cout << " fixDoubleAnchor_Match starts ... " << endl;
		bool fixDoubleAnchorBool_Match = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int chrNameInt = indexInfo->convertStringToInt(chromName);
		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;

		if(subSeqLengthInProcess < 2)
		{
			Jump_Code matchJumpCode(//segmentLength_1 + 
				subSeqLengthInProcess + segmentLength_2, "M");
			cigarInfo->jump_code.push_back(matchJumpCode);
			fixDoubleAnchorBool_Match = true;
			(*mismatchNum) = subSeqLengthInProcess;
		
			if(STORE_MISMATCH_POS)
			{
				this->insertMismatchPosVec_fixMatchInfo_SmallGap(segmentLocInRead_1 + segmentLength_1,
					subSeqLengthInProcess, mismatchPosVec);
				if(STORE_MISMATCH_CHA)
				{
					this->insertMismatchCharVec_fixMatchInfo_SmallGap(subSeqLengthInProcess, 
						chrNameInt, segmentMapPos_1 + segmentLength_1, indexInfo, mismatchCharVec);
				}
			}
		}
		else
		{
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);
		
			int chrNameInt = indexInfo->convertStringToInt(chromName);
		
			//string chromSubSeqInProcess = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1, subSeqLengthInProcess);
			string chromSubSeqInProcess = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_1 + segmentLength_1, subSeqLengthInProcess);
			/*
			size_t max_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;
			//cout << "readSeqInProcess: " << endl << readSubSeqInProcess << endl;
			//cout << "chromSeqInProcess: " << endl << chromSubSeqInProcess << endl;

			bool scoreStringBool = score_string(readSubSeqInProcess, chromSubSeqInProcess, max_mismatch, mismatch_bits, comb_bits);// need to debug
			*/
			//cout << "scoreStringBool: " << scoreStringBool << endl;
			
			int max_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;
			FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
			bool scoreStringBool = fixMatchInfo->fixMatch(readSubSeqInProcess, chromSubSeqInProcess,
				max_mismatch, segmentLocInRead_1 + segmentLength_1);
			if(scoreStringBool)
			{
				//(*mismatchNum) = mismatch_bits;
				(*mismatchNum) = fixMatchInfo->returnMismatchNum();
				Jump_Code matchJumpCode(//segmentLength_1 + 
					subSeqLengthInProcess + segmentLength_2, "M");
				cigarInfo->jump_code.push_back(matchJumpCode);
				if(STORE_MISMATCH_POS)
				{
					this->insertMismatchPosVec_fixMatchInfo(fixMatchInfo, mismatchPosVec);
					if(STORE_MISMATCH_CHA)
					{
						this->insertMismatchCharVec_fixMatchInfo(fixMatchInfo, mismatchCharVec);
					}
				}
			}
			else // score string failed, insert sudo-match jump code
			{
				//cout << " score_string failed !" << endl;
				//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");	
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);			
			}
			delete fixMatchInfo;
			fixDoubleAnchorBool_Match = scoreStringBool;
		}
		return fixDoubleAnchorBool_Match;
	}
	void insertMismatchPosVec_fixInsInfo(FixDoubleAnchor_Insertion_Info* fixInsInfo, vector<int>& mismatchPosVec)
	{
		int vecSize = fixInsInfo->returnBestInsertionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixInsInfo->returnBestInsertionMismatchPos(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);	
		}	
	}

	void insertMismatchCharVec_fixInsInfo(FixDoubleAnchor_Insertion_Info* fixInsInfo, vector<char>& mismatchCharVec)
	{
		int vecSize = fixInsInfo->returnBestInsertionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixInsInfo->returnBestInsertionMismatchChar(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}
	bool fixDoubleAnchor_Insertion(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum, vector<int>& mismatchPosVec, vector<char>& mismatchCharVec
		)
	{
		bool fixDoubleAnchorBool_Insertion = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int insertionLength = (segmentMapPos_1 - segmentLocInRead_1) - (segmentMapPos_2 - segmentLocInRead_2);
		int chrNameInt = indexInfo->convertStringToInt(chromName);

		if(subSeqLengthInProcess <= insertionLength)
		{
			//(*mismatchNum) = subSeqLengthInProcess;
			//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
			//cigarInfo->jump_code.push_back(firstMatchJumpCode);

			int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - insertionLength;
			if(secondMatchLength > 0)
			{
				Jump_Code midInsertionJumpCode(insertionLength, "I");	
				Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess - insertionLength, "M");							
				cigarInfo->jump_code.push_back(midInsertionJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				fixDoubleAnchorBool_Insertion = true;
				(*mismatchNum) = 0;
			}
			else
			{
				Jump_Code midInsertionJumpCode(insertionLength, "i");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");
				cigarInfo->jump_code.push_back(midInsertionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				fixDoubleAnchorBool_Insertion = false;
			}
			//fixDoubleAnchorBool_Insertion = true;	
		}	  
		else
		{
			/*
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeqLengthInProcess);
			//string chromSubSeqInProcess = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1, segmentMapPos_2 - 1 - (segmentMapPos_1 + segmentLength_1) + 1);
			string chromSubSeqInProcess = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_1 + segmentLength_1, segmentMapPos_2 - 1 - (segmentMapPos_1 + segmentLength_1) + 1);
			size_t prefix_length = 0;
			size_t mismatch_bits = 0; //?
			size_t max_ins_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 2;
			size_t comb_bits_ins = 0;

			GenomeScan* genome_scan = new GenomeScan;
			bool insertion_fixed = (*genome_scan).Double_anchored_score_ins(readSubSeqInProcess, chromSubSeqInProcess, max_ins_mismatch, prefix_length, comb_bits_ins, mismatch_bits); //X: fix insertion
			*/
			FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info();
			int tmp_toFix_insertion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_insertion_read_end = segmentLocInRead_2 - 1;
			int tmp_toFix_insertion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_insertion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;  // how to set the max-mismatch parameter

			bool insertion_fixed = insInfo->detectBestInsertion_lessMismatch(
				tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
				tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
				chrNameInt, tmp_max_allowed_mismatchNum);

			if(insertion_fixed)
			{
				//int firstMatchLength = //segmentLength_1 + 
				//		prefix_length;
				
				int firstMatchLength = insInfo->return_best_insertion_prefix_match_length();
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - firstMatchLength - insertionLength;
				int mismatch_bits = insInfo->return_best_insertion_mismatch();
				if(secondMatchLength > 0)
				{
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");
					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(insertionJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);
					(*mismatchNum) = mismatch_bits;
					if(STORE_MISMATCH_POS)
					{
						//insInfo->generateBestInsertionMismatchVec((segmentLocInRead_1 + segmentLength_1));
						this->insertMismatchPosVec_fixInsInfo(insInfo, mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							this->insertMismatchCharVec_fixInsInfo(insInfo, mismatchCharVec);
						}
					}
				}	
				else
				{
					insertion_fixed = false;
					//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
					Jump_Code insertionJumpCode(insertionLength, "i");
					Jump_Code midMatchJumpCode(subSeqLengthInProcess - insertionLength, "m");
					Jump_Code secondMatchJumpCode(segmentLength_2, "M");
					//cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(insertionJumpCode);
					cigarInfo->jump_code.push_back(midMatchJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);
				}			
			}
			else
			{
				//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
				Jump_Code insertionJumpCode(insertionLength, "i");
				Jump_Code midMatchJumpCode(subSeqLengthInProcess - insertionLength, "m");
				Jump_Code secondMatchJumpCode(segmentLength_2, "M");
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(insertionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);				
			}	
			//delete(genome_scan);
			delete insInfo;
			fixDoubleAnchorBool_Insertion = insertion_fixed;
		}	
		return fixDoubleAnchorBool_Insertion;
	}
	void insertMismatchPosVec_fixDeletionInfo_SmallGap(
		int startPosInRead, int toProcessSeq_len, vector<int>& mismatchPosVec)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			int tmpMismatchPos = startPosInRead + tmp;
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void insertMismatchCharVec_fixDeletionInfo_SmallGap(
		int toProcessSeq_len, int chrNameInt, int startPosInChrom, 
		Index_Info* indexInfo, vector<char>& mismatchCharVec)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			char tmpMismatchChar = indexInfo->returnOneBaseCharInGenome(chrNameInt, startPosInChrom + tmp);
			//int tmpMismatchPos = startPosInRead + tmp;
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}
	void insertMismatchPosVec_fixDelInfo(FixDoubleAnchor_Deletion_Info* fixDelInfo, vector<int>& mismatchPosVec)
	{
		int vecSize = fixDelInfo->returnBestDeletionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixDelInfo->returnBestDeletionMismatchPos(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);	
		}	
	}

	void insertMismatchCharVec_fixDelInfo(FixDoubleAnchor_Deletion_Info* fixDelInfo, vector<char>& mismatchCharVec)
	{
		int vecSize = fixDelInfo->returnBestDeletionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixDelInfo->returnBestDeletionMismatchChar(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}	
	bool fixDoubleAnchor_Deletion(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum, 
		vector<int>& mismatchPosVec, vector<char>& mismatchCharVec)
	{
		//cout << "start to fix deletion ... " << endl;
		bool fixDoubleAnchorBool_Deletion = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
		int deletionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);
		//cout << "deletionLength: " << deletionLength << endl;
		int chrNameInt = indexInfo->convertStringToInt(chromName);

		if(subSeqLengthInProcess < 2)
		{
			//cout << " subSeqLengthInProcess < 2 " << endl;
			//Jump_Code firstMatchJumpCode(segmentLength_1, "M");
			Jump_Code deletionJumpCode(deletionLength, "D");
			Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess, "M");
			//cigarInfo->jump_code.push_back(firstMatchJumpCode);
			cigarInfo->jump_code.push_back(deletionJumpCode);
			cigarInfo->jump_code.push_back(secondMatchJumpCode);
			fixDoubleAnchorBool_Deletion = true;
			(*mismatchNum) = subSeqLengthInProcess;
			if(STORE_MISMATCH_POS)
			{
				this->insertMismatchPosVec_fixDeletionInfo_SmallGap(
					segmentLocInRead_1 + segmentLength_1, subSeqLengthInProcess, mismatchPosVec);					
				if(STORE_MISMATCH_CHA)
				{
					this->insertMismatchCharVec_fixDeletionInfo_SmallGap(
						subSeqLengthInProcess, chrNameInt, 
						segmentMapPos_2 - subSeqLengthInProcess, indexInfo, mismatchCharVec);
				}
			}
		}
		else
		{
			/*
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1, subSeqLengthInProcess);
			int chromSubSeqLengthInProcess= subSeqLengthInProcess + 2;
			//string left_chrom_seq = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_1 + segmentLength_1 - 1, chromSubSeqLengthInProcess);
			string left_chrom_seq = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_1 + segmentLength_1, chromSubSeqLengthInProcess);
			//string right_chrom_seq = indexInfo->chromStr[chrNameInt].substr(segmentMapPos_2 - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);
			string right_chrom_seq = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_2 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);

			bool small_deletion = true;
			size_t prefix_length = 0;
			size_t max_double_splice_mismatch = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 2;
			size_t mismatch_bits = 0;
			size_t comb_bits = 0;
			GenomeScan* genome_scan = new GenomeScan;
			bool deletion_fixed = (*genome_scan).Double_anchored_score_least_mis(
				readSubSeqInProcess, left_chrom_seq, right_chrom_seq, 
				prefix_length, max_double_splice_mismatch, comb_bits, small_deletion, mismatch_bits);*/
	
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
				Jump_Code firstMatchJumpCode(firstMatchLength, "M");
				Jump_Code deletionJumpCode(deletionLength, "D");
				Jump_Code secondMatchJumpCode(secondMatchLength, "M");
				cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(deletionJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);
				(*mismatchNum) = mismatch_bits;
				if(STORE_MISMATCH_POS)
				{
					//delInfo->generateBestDeletionMismatchVec((segmentLocInRead_1 + segmentLength_1));
					this->insertMismatchPosVec_fixDelInfo(delInfo, mismatchPosVec);
					if(STORE_MISMATCH_CHA)
					{
						this->insertMismatchCharVec_fixDelInfo(delInfo, mismatchCharVec);
					}
				}	
			}
			else
			{
				//int firstMatchLength = segmentLength_1;
				int midMatchLength = subSeqLengthInProcess;
				int secondMatchLength = segmentLength_2;
				//Jump_Code firstMatchJumpCode(firstMatchLength, "M");
				Jump_Code deletionJumpCode(deletionLength, "d");
				Jump_Code midMatchJumpCode(midMatchLength, "m");
				Jump_Code secondMatchJumpCode(secondMatchLength, "M");				
				//cigarInfo->jump_code.push_back(firstMatchJumpCode);
				cigarInfo->jump_code.push_back(deletionJumpCode);
				cigarInfo->jump_code.push_back(midMatchJumpCode);
				cigarInfo->jump_code.push_back(secondMatchJumpCode);				
			}
			delete delInfo;
			//delete(genome_scan);
			fixDoubleAnchorBool_Deletion = deletion_fixed;
		}
		return fixDoubleAnchorBool_Deletion;
	}

	void insertMismatchPosVec_fixSpliceInfo(FixDoubleAnchor_Splice_Info* fixSpliceInfo, vector<int>& mismatchPosVec)
	{
		int vecSize = fixSpliceInfo->returnBestSpliceMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixSpliceInfo->returnMismatchPosInRead(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixSpliceInfo(FixDoubleAnchor_Splice_Info* fixSpliceInfo, vector<char>& mismatchCharVec)
	{
		int vecSize = fixSpliceInfo->returnBestSpliceMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixSpliceInfo->returnMismatchChar(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}	

	void insertMismatchPosVec_fixComplicatedSpliceInfo(
		FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo, vector<int>& mismatchPosVec)
	{
		int vecSize = fixComplicateSpliceInfo->returnBestComplicateSpliceMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixComplicateSpliceInfo->returnMismatchPosInRead(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixComplicateSpliceInfo(
		FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo, vector<char>& mismatchCharVec)
	{
		int vecSize = fixComplicateSpliceInfo->returnBestComplicateSpliceMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixComplicateSpliceInfo->returnMismatchChar(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}	

	// DO CANONICAL SPLICE-JUNCTION ONLY 
	bool fixDoubleAnchor_Splice_phase2(Splice_Info* cigarInfo, int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, const string& chromName, int* mismatchNum,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
		vector<int>& mismatchPosVec, vector<char>& mismatchCharVec)
	{
		//cout << "start to fix splice" << endl;
		bool fixDoubleAnchorBool_Splice = false;
		bool final_splice_fixed = false;
		int chrNameInt = indexInfo->convertStringToInt(chromName);
		//cout << "chrNameInt: " << chrNameInt << endl;

		int tmpBuffer_left = FixSpliceBuffer;
		if(tmpBuffer_left > segmentLength_1 - 2) //anchor >= 2
		{
			tmpBuffer_left = segmentLength_1 - 2;
		}
		int tmpBuffer_right = FixSpliceBuffer;
		if(tmpBuffer_right > segmentLength_2 - 2)
		{
			tmpBuffer_right = segmentLength_2 - 2;
		}

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1 + tmpBuffer_left + tmpBuffer_right;
		int spliceJunctionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);
		int chromSubSeqLengthInProcess = subSeqLengthInProcess + 2;

		size_t max_double_splice_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed;
		if((annotation_provided_bool && Do_annotation_only_bool) || !(annotation_provided_bool))
		{	
			FixDoubleAnchor_Splice_Info* fixSpliceInfo = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo->detectBestSpliceSite//_canonicalOnly_lessMismatch
				( //detectBestSpliceSite_prefer_canonical_lessMismatch(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);		

			if(FIX_INDEL_AROUND_SJ_BOOL)
			{
				if(splice_fixed && (fixSpliceInfo->fixSpliceResultConfident()) )
				{
					int prefix_length = fixSpliceInfo->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");

					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(spliceJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);

					(*mismatchNum) = fixSpliceInfo->returnBestSplice_mismatchNum();
					final_splice_fixed = true;
				
					//if(STORE_MISMATCH_POS)
					//{
						//fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, 
						//	(segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, mismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, mismatchCharVec);
						//}
					//}
				}
				else
				{
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

						final_splice_fixed = true;
						(*mismatchNum) = fixComplicateSpliceInfo->return_mismatch_bestComplicatedSplice();
						cigarInfo->jump_code.push_back(prefixMatchJumpCode);	
						cigarInfo->jump_code.push_back(firstJumpCode);
						cigarInfo->jump_code.push_back(midMatchJumpCode);	
						cigarInfo->jump_code.push_back(secondJumpCode);	
						cigarInfo->jump_code.push_back(suffixMatchJumpCode);			
						//if(STORE_MISMATCH_POS)
						//{
							//fixComplicateSpliceInfo->generateBestComplicateSJMismatchVec(
							//	segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
							this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo, mismatchPosVec);
							//if(STORE_MISMATCH_CHA)
							//{
								this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo, mismatchCharVec);
							//}
						//}		
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

						final_splice_fixed = true;
						(*mismatchNum) = fixSpliceInfo->returnBestSplice_mismatchNum();
						cigarInfo->jump_code.push_back(firstMatchJumpCode);
						cigarInfo->jump_code.push_back(spliceJumpCode);	
						cigarInfo->jump_code.push_back(secondMatchJumpCode);
						//if(STORE_MISMATCH_POS)
						//{
							//fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, mismatchPosVec);
							//if(STORE_MISMATCH_CHA)
							//{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, mismatchCharVec);
							//}
						//}
					}
					else
					{
						final_splice_fixed = false;
					}
					delete fixComplicateSpliceInfo;
				}
			}
			else
			{
				if(splice_fixed)
				{
					int prefix_length = fixSpliceInfo->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");

					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(spliceJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);

					(*mismatchNum) = fixSpliceInfo->returnBestSplice_mismatchNum();
					final_splice_fixed = true;
					if(STORE_MISMATCH_POS)
					{
						//fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, mismatchCharVec);
						}
					}
				}
				else
				{
					final_splice_fixed = false;
				}
			}
			delete fixSpliceInfo;
			//delete genome_scan;
			fixDoubleAnchorBool_Splice = final_splice_fixed;
		}
		else
		{
			bool annotated_splice_fixed_bool = false;
			FixDoubleAnchor_Splice_Info* fixSpliceInfo_annotated = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo_annotated->detectBestSpliceSite//_canonicalOnly_lessMismatch
				( //detectBestSpliceSite_prefer_canonical_lessMismatch(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, true, annotationInfo);		
			if(FIX_INDEL_AROUND_SJ_BOOL)
			{
				if(splice_fixed && ((fixSpliceInfo_annotated->returnBestSplice_canonicalOrNot()) 
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

					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(spliceJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);

					(*mismatchNum) = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
					final_splice_fixed = true;
					// fixed with annotated SJ
					annotated_splice_fixed_bool = true;
					
					if(STORE_MISMATCH_POS)
					{
						//fixSpliceInfo_annotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, 
						//	(segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchCharVec);
						}
					}
				}
				else
				{
					FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_annotated = new FixDoubleAnchor_Splice_Complicate_Info();
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

						final_splice_fixed = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;

						(*mismatchNum) = fixComplicateSpliceInfo_annotated->return_mismatch_bestComplicatedSplice();
						cigarInfo->jump_code.push_back(prefixMatchJumpCode);	
						cigarInfo->jump_code.push_back(firstJumpCode);
						cigarInfo->jump_code.push_back(midMatchJumpCode);	
						cigarInfo->jump_code.push_back(secondJumpCode);	
						cigarInfo->jump_code.push_back(suffixMatchJumpCode);			
						if(STORE_MISMATCH_POS)
						{
							//fixComplicateSpliceInfo_annotated->generateBestComplicateSJMismatchVec(
							//	segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
							this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_annotated, mismatchPosVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_annotated, mismatchCharVec);
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

						final_splice_fixed = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;

						(*mismatchNum) = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
						cigarInfo->jump_code.push_back(firstMatchJumpCode);
						cigarInfo->jump_code.push_back(spliceJumpCode);	
						cigarInfo->jump_code.push_back(secondMatchJumpCode);
						if(STORE_MISMATCH_POS)
						{
							//fixSpliceInfo_annotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchPosVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchCharVec);
							}
						}
					}
					else
					{
						final_splice_fixed = false;
						// unfixed with annotated SJ
						annotated_splice_fixed_bool = false;
					}
					delete fixComplicateSpliceInfo_annotated;
				}
			}
			else
			{
				if(splice_fixed)
				{
					int prefix_length = fixSpliceInfo_annotated->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");

					cigarInfo->jump_code.push_back(firstMatchJumpCode);
					cigarInfo->jump_code.push_back(spliceJumpCode);
					cigarInfo->jump_code.push_back(secondMatchJumpCode);

					(*mismatchNum) = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
					final_splice_fixed = true;
					// fixed with annotated SJ
					annotated_splice_fixed_bool = true;

					if(STORE_MISMATCH_POS)
					{
						//fixSpliceInfo_annotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, mismatchCharVec);
						}
					}
				}
				else
				{
					final_splice_fixed = false;
					// unfixed with annotated SJ
					annotated_splice_fixed_bool = false;
				}
			}
			delete fixSpliceInfo_annotated;

			if(annotated_splice_fixed_bool)
			{}
			else
			{
				FixDoubleAnchor_Splice_Info* fixSpliceInfo_unannotated = new FixDoubleAnchor_Splice_Info();
				bool splice_fixed = fixSpliceInfo_unannotated->detectBestSpliceSite//_canonicalOnly_lessMismatch
					( //detectBestSpliceSite_prefer_canonical_lessMismatch(
					segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
					segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
					readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
					false, false, annotationInfo);		
				if(FIX_INDEL_AROUND_SJ_BOOL)
				{
					if(splice_fixed && ((fixSpliceInfo_unannotated->returnBestSplice_canonicalOrNot()) 
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

						cigarInfo->jump_code.push_back(firstMatchJumpCode);
						cigarInfo->jump_code.push_back(spliceJumpCode);
						cigarInfo->jump_code.push_back(secondMatchJumpCode);

						(*mismatchNum) = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
						final_splice_fixed = true;
						// fixed with unannotated SJ
						
						if(STORE_MISMATCH_POS)
						{
							//fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, 
							//	(segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchPosVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchCharVec);
							}
						}
					}
					else
					{
						FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_unannotated = new FixDoubleAnchor_Splice_Complicate_Info();
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

							final_splice_fixed = true;
							// fixed with unannotated SJ

							(*mismatchNum) = fixComplicateSpliceInfo_unannotated->return_mismatch_bestComplicatedSplice();
							cigarInfo->jump_code.push_back(prefixMatchJumpCode);	
							cigarInfo->jump_code.push_back(firstJumpCode);
							cigarInfo->jump_code.push_back(midMatchJumpCode);	
							cigarInfo->jump_code.push_back(secondJumpCode);	
							cigarInfo->jump_code.push_back(suffixMatchJumpCode);			
							if(STORE_MISMATCH_POS)
							{
								//fixComplicateSpliceInfo_unannotated->generateBestComplicateSJMismatchVec(
								//	segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
								this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_unannotated, mismatchPosVec);
								if(STORE_MISMATCH_CHA)
								{
									this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_unannotated, mismatchCharVec);
								}
							}		
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

							final_splice_fixed = true;
							// fixed with unannotated SJ

							(*mismatchNum) = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
							cigarInfo->jump_code.push_back(firstMatchJumpCode);
							cigarInfo->jump_code.push_back(spliceJumpCode);	
							cigarInfo->jump_code.push_back(secondMatchJumpCode);
							if(STORE_MISMATCH_POS)
							{
								//fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
								this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchPosVec);
								if(STORE_MISMATCH_CHA)
								{
									this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchCharVec);
								}
							}
						}
						else
						{
							final_splice_fixed = false;
							// unfixed with unannotated SJ
						}
						delete fixComplicateSpliceInfo_unannotated;
					}
				}
				else
				{
					if(splice_fixed)
					{
						int prefix_length = fixSpliceInfo_unannotated->returnBestSplice_prefixMatchLength();
						int firstMatchLength = //segmentLength_1 
							- tmpBuffer_left + prefix_length;
						int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
						//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

						Jump_Code firstMatchJumpCode(firstMatchLength, "M");
						Jump_Code spliceJumpCode(spliceJunctionLength, "N");
						Jump_Code secondMatchJumpCode(secondMatchLength, "M");

						cigarInfo->jump_code.push_back(firstMatchJumpCode);
						cigarInfo->jump_code.push_back(spliceJumpCode);
						cigarInfo->jump_code.push_back(secondMatchJumpCode);

						(*mismatchNum) = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
						final_splice_fixed = true;
						// fixed with unannotated SJ

						if(STORE_MISMATCH_POS)
						{
							//fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchPosVec);
							if(STORE_MISMATCH_CHA)
							{
								this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, mismatchCharVec);
							}
						}
					}
					else
					{
						final_splice_fixed = false;
						// unfixed with unannotated SJ
					}
				}
				delete fixSpliceInfo_unannotated;
			}

			//delete genome_scan;
			fixDoubleAnchorBool_Splice = final_splice_fixed;						
		}
		return fixDoubleAnchorBool_Splice;
	}		
};

#endif