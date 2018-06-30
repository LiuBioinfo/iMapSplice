// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PATH_INFO_H
#define PATH_INFO_H

#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <map>
#include "index_info.h"
//#include "fixDoubleAnchor_annotation_info.h"

using namespace std;

class Path_Info
{
//public:
private:
	//*********************   Intermediate results   ********************************//
	int LengthOfSeqPerMismatchAllowed;
	int FixSpliceBuffer;
	map< pair< pair<int,int>, pair<int,int> >, int > fixGapMap; // < <seg_group_1, seg_candi_1> <seg_group_2, seg_candi_2> >,  index_in_vector< pair<bool, vector<int, vector<Jump_Code> > > >
	vector< pair<bool, pair<int, vector<Jump_Code> > > > fixGapVec; // < fixed_or_not, vector <mismatch, vector<Jump_Code> >  >
	int repeatRegion_index;
	vector< vector<int> > fixGapMismatchPosVec;
	vector< vector<char> > fixGapMismatchCharVec;	
	
	//************************  penultimate final Path vector  **********************************//
	// PathVec_seg.size() == gapIndexVec.size() == PathValidBoolVec.size() == PathFixedBoolVec.size();
	vector< vector< pair<int, int> > > PathVec_seg; // vector <vector <seg_group, seg_candi> >
	vector< vector<int> > PathVec_seg_gapCase;	
	// (gapIndexVec[A])[B] is the gap_index  in fixGapVec for the gap between (PathVec_seg[A])[B] and (PathVec_seg[A])[B+1]
	vector< vector<int> > PathVec_seg_gapIndex; // gapIndexVec[A].size() = PathVec_seg[A].size() - 1
	vector< bool > PathValidBoolVec; // PathValidBoolVec.size() = PathVec_seg.size()
	vector< bool > PathFixedBoolVec;

	//************************   Final Path vector  **********************************//
	// fixedPathVec.size() == fixedPathMismatchVec.size() == finapPathVec.size() == true # in PathFixedBoolVec
	//  == fixedPathMismatchPosVec.size() == fixedPathMismatchCharVec.size()
	vector< pair<int, Splice_Info*> > fixedPathVec; 
	vector< int > fixedPathMismatchVec; 

	vector < pair< pair<int, int>, Splice_Info*> > finalPathVec; // <chrNameInt, chrPosInt>, cigarStringJumpCodeVec

public:
	vector< vector<int> > fixedPathMismatchPosVec;
	vector< vector<char> > fixedPathMismatchCharVec;	


	void pushBack_CandiPathPairInfo(vector< pair<int, int> >& tmpPathVec_seg, vector<Jump_Code>& tmpJumpCodeVec,
		vector<int>& tmpMismatchPosVec, vector<char>& tmpMismatchCharVec, int tmpPathNumIndex)
	{
		PathVec_seg.push_back(tmpPathVec_seg);
		PathValidBoolVec.push_back(true);
		PathFixedBoolVec.push_back(true);
		Splice_Info* tmpSpliceInfo = new Splice_Info();
		for(int tmp = 0; tmp < tmpJumpCodeVec.size(); tmp++)
		{
			tmpSpliceInfo->jump_code.push_back(tmpJumpCodeVec[tmp]);
		}
		fixedPathVec.push_back(pair<int, Splice_Info*>(tmpPathNumIndex, tmpSpliceInfo));
		int tmpMismatchNum = tmpMismatchPosVec.size();
		fixedPathMismatchVec.push_back(tmpMismatchNum);
		fixedPathMismatchPosVec.push_back(tmpMismatchPosVec);
		fixedPathMismatchCharVec.push_back(tmpMismatchCharVec);
	}	

	void memoryFree()
	{
		vector< vector< pair<int, int> > >().swap(PathVec_seg);
		vector< vector<int> >().swap(PathVec_seg_gapIndex);
		vector< bool >().swap(PathValidBoolVec);
		vector< bool >().swap(PathFixedBoolVec);

		for(int tmp = 0; tmp < fixedPathVec.size(); tmp++)
		{
			delete fixedPathVec[tmp].second;
		}
		vector< pair<int, Splice_Info*> >().swap(fixedPathVec);
		vector<int>().swap(fixedPathMismatchVec);

		for(int tmp = 0; tmp < finalPathVec.size(); tmp++)
		{
			delete finalPathVec[tmp].second;
		}
		vector < pair< pair<int, int>, Splice_Info*> >().swap(finalPathVec);
		map< pair< pair<int,int>, pair<int,int> >, int >().swap(fixGapMap);
		vector< pair<bool, pair<int, vector<Jump_Code> > > >().swap(fixGapVec);
		vector< vector<int> >().swap(fixGapMismatchPosVec);
		vector< vector<char> >().swap(fixGapMismatchCharVec);
		vector< vector<int> >().swap(fixedPathMismatchPosVec);
		vector< vector<char> >().swap(fixedPathMismatchCharVec);
	}

	string getFixGapInfoStr()
	{
		string fixGapInfoStr = "fixGapMapInfo: \n";
		for(map < pair< pair<int,int>, pair<int,int> >, int >::iterator mapIter = fixGapMap.begin(); 
			mapIter != fixGapMap.end(); mapIter++)
		{
			fixGapInfoStr = fixGapInfoStr + "seg_1: " + int_to_str(((mapIter->first).first).first+1) 
				+ "-" + int_to_str(((mapIter->first).first).second+1)
				+ " ~ " + "seg_2: " +  int_to_str(((mapIter->first).second).first+1) 
				+ "-" + int_to_str(((mapIter->first).second).second+1) + " index: " + int_to_str(mapIter->second+1) + "\n"; 
		}
		fixGapInfoStr += "\nfixGapVec:";
		for(int tmp = 0; tmp < fixGapVec.size(); tmp ++)
		{
			fixGapInfoStr = fixGapInfoStr + "\nGap:" + int_to_str(tmp+1);
			if(fixGapVec[tmp].first)
			{
				fixGapInfoStr = fixGapInfoStr + " mismatch#: " + int_to_str((fixGapVec[tmp].second).first) 
					+ " mismatchPos: ";
				for(int tmp2 = 0; tmp2 < fixGapMismatchPosVec[tmp].size(); tmp2++)
				{
					fixGapInfoStr = fixGapInfoStr + int_to_str((fixGapMismatchPosVec[tmp])[tmp2]) + ",";
				}
			}
			else
			{
				fixGapInfoStr += " -- false";
			}
		}
		return fixGapInfoStr;
	}

	Path_Info()
	{
		LengthOfSeqPerMismatchAllowed = 10;//15;
		FixSpliceBuffer = FIX_SPLICE_BUFFER_MAX;
		repeatRegion_index = -1;
	}

	void pushBack2MismatchPosVec(vector<int>& newMismatchPosVec, int index_fixGapVec)
	{
		int vecSize = fixGapMismatchPosVec[index_fixGapVec].size();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			newMismatchPosVec.push_back((fixGapMismatchPosVec[index_fixGapVec])[tmp]);
		}
	}

	void pushBack2MismatchCharVec(vector<char>& newMismatchCharVec, int index_fixGapVec)
	{
		int vecSize = fixGapMismatchCharVec[index_fixGapVec].size();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			newMismatchCharVec.push_back((fixGapMismatchCharVec[index_fixGapVec])[tmp]);
		}
	}

	void pushBack2FinalPathMismatchPosVec(vector<int>& toPushBackMismatchPosVec, int index_fixedPathMismatchPosVec)
	{
		for(int tmp = 0; tmp < toPushBackMismatchPosVec.size(); tmp++)
		{
			fixedPathMismatchPosVec[index_fixedPathMismatchPosVec].push_back(toPushBackMismatchPosVec[tmp]);
		}
	}

	void pushBack2FinalPathMismatchCharVec(vector<char>& toPushBackMismatchCharVec, int index_fixedPathMismatchCharVec)
	{
		for(int tmp = 0; tmp < toPushBackMismatchCharVec.size(); tmp++)
		{
			fixedPathMismatchCharVec[index_fixedPathMismatchCharVec].push_back(toPushBackMismatchCharVec[tmp]);
		}
	}	

	bool fixAllPath_fixGapAlso(Seg_Info* segInfo, Index_Info* indexInfo, const string& readSeq_inProcess, bool Do_extendHeadTail)
	{
		//cout << endl << "*********************************"
		//	<< endl << "start to do fixAllPath_fixGapAlso ..." << endl;
		// #ifdef CAL_TIME
		// fixGap_1_begin = clock();
		// #endif
 		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}
		
		int pathVecSize = this->PathVec_seg.size();
 		//cout << endl << "pathVecSize: " << pathVecSize << endl;
 		for(int tmpPath = 0; tmpPath < pathVecSize; tmpPath++)
 		{		
 			//cout << "tmpPath: " << tmpPath << endl;
			if(!(PathValidBoolVec[tmpPath]))
			{
				//newPathFixed = false;
				this->PathFixedBoolVec.push_back(false);				
				continue;
			}
	 		//cout << "valid tmpPath: " << tmpPath << endl;
			int firstSegGroupNO = (PathVec_seg[tmpPath])[0].first;
			int firstSegLength = segInfo->returnSegmentLength(firstSegGroupNO);//(segInfo->norSegmentLength)[firstSegGroupNO];
			//cout << "firstSegGroupNO: " << firstSegGroupNO << endl;
			//cout << "firstSegLength: " << firstSegLength << endl;
			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(firstSegLength, "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);

			int newPathMismatchNum = 0;
			
			bool newPathFixed = true;
			
			int tmpPathSegSize = PathVec_seg[tmpPath].size();

			vector<int> newMismatchPosVec;
			vector<char> newMismatchCharVec;
 			for(int tmpGap = 0; tmpGap < tmpPathSegSize - 1; tmpGap++)
 			{
 				int tmpGapIndex = (PathVec_seg_gapIndex[tmpPath])[tmpGap];
 				bool tmpGapFixedOrNot = fixGapVec[tmpGapIndex].first;
 				int tmpMismatchNum = (fixGapVec[tmpGapIndex].second).first;
 				int tmpGapJumpCodeVecSize = ((fixGapVec[tmpGapIndex].second).second).size();
 				for(int tmpGapJumpCode = 0; tmpGapJumpCode < tmpGapJumpCodeVecSize;
 					tmpGapJumpCode++)
 				{
 					newPathSpliceInfo->jump_code.push_back(
 						((fixGapVec[tmpGapIndex].second).second)[tmpGapJumpCode]);
 				}
 				newPathMismatchNum += tmpMismatchNum;

 				// if(STORE_MISMATCH_POS)
 				// {
 					this->pushBack2MismatchPosVec(newMismatchPosVec, tmpGapIndex);
 					//if(STORE_MISMATCH_CHA)
 					//{
 						this->pushBack2MismatchCharVec(newMismatchCharVec, tmpGapIndex);		
 					//}
 				//}
 			}
 			//cout << "newPathFixed: " << newPathFixed << endl;
 			if(newPathFixed)
 			{
				newPathSpliceInfo->getFinalJumpCode();	

				bool allJumpCodeValidBool = newPathSpliceInfo->allFinalJumpCodeValid();
				//cout << "allJumpCodeValidBool: " << allJumpCodeValidBool << endl;
				if(allJumpCodeValidBool)
				{
					(PathFixedBoolVec).push_back(allJumpCodeValidBool);
					(fixedPathVec).push_back(pair <int, Splice_Info*> (tmpPath, newPathSpliceInfo) );
					(fixedPathMismatchVec).push_back(newPathMismatchNum);
					//if(STORE_MISMATCH_POS)
					//{
						fixedPathMismatchPosVec.push_back(newMismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							fixedPathMismatchCharVec.push_back(newMismatchCharVec);
						//}
					//}
				}
				else
				{
					(PathFixedBoolVec).push_back(allJumpCodeValidBool);
				} 				
 			} 			
 			//cout << "end of checking path " << tmpPath << endl;
 		}
		
		// #ifdef CAL_TIME
		// fixGap_1_end = clock();
		// fixGap_1_cost = fixGap_1_cost + fixGap_1_end - fixGap_1_begin;
		// fixGap_2_begin = clock();
		// #endif	
		
		int readLength = readSeq_inProcess.length();
		this->getFinalPath_extend2HeadTail_new(indexInfo, segInfo, readLength, readSeq_inProcess, Do_extendHeadTail);
		
		// #ifdef CAL_TIME
		// fixGap_2_end = clock();
		// fixGap_2_cost = fixGap_2_cost + fixGap_2_end - fixGap_2_begin;
		// #endif
		
		return true; 			
	}

	bool getPossiPathFromSeg_fixGapAlso(Seg_Info* segInfo, Index_Info* indexInfo, const string& readSeq_inProcess,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		//cout << endl << endl << "*************************************************" << endl 
		//	<< "new getPossiPathFromSeg_fixGapAlso" << endl; 
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO =  segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);

				vector<int> tmpPathElementVec_gapIndex;
				PathVec_seg_gapIndex.push_back(tmpPathElementVec_gapIndex);

				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
				vector<int>().swap(tmpPathElementVec_gapIndex);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec_fixGapAlso(segInfo, tmpSegNO, indexInfo, readSeq_inProcess,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
		}

		possiPathExists = true;
		return possiPathExists;

		//cout << "end of getPossiPathFromSeg_fixGapAlso ...." << endl
		//	<< "*****************************************************************" << endl << endl;
	}


	void getFinalPath_extend2HeadTail_new(Index_Info* indexInfo, Seg_Info* segInfo, int readLength, 
		const string& readSeq_inProcess, bool Do_extendHeadTail)
	{
		//cout << "start to get Final path" << endl;

		int minMatchBaseToSupport_perMismatch = MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH;
		int minMatchBaseToSupport_twoMismatch = MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH;
		if(!Do_extendHeadTail)
		{
			minMatchBaseToSupport_perMismatch = 100000;
			minMatchBaseToSupport_twoMismatch = 100000;
		}

		if(segInfo->returnSegmentNum() < 1)
		{return;}

		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{
			///////////// try to capture whether path # is a key factor for running time
			//if(tmpPath > 0)
			//	return;
			//cout << "tmpPath: " << tmpPath << endl;
			int mismatchNumToAdd = 0;
			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;

			int tmpPathElementSize = (PathVec_seg[fixedPathNO]).size();
			int tmpPathLastSegGroupNO = (PathVec_seg[fixedPathNO])[tmpPathElementSize - 1].first;

			unsigned int PathMapPos //= *(segInfo->norSegmentAlignLoc + tmpPath1stSegGroupNO * CANDALILOC + tmpPath1stSegCandiNO);
				= segInfo->returnSegmentMapPos(tmpPath1stSegGroupNO, tmpPath1stSegCandiNO);

			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);

			
			int tmpPathFinalMapPos = tmpChrPosInt;
			int tmpPathFinalMapChr = tmpChrNameInt;

			Splice_Info* tmpSpliceInfo = new Splice_Info();
			tmpSpliceInfo->jump_code.clear(); 

			int tmpUnfixedHeadLength 
				= segInfo->returnSegmentLocInRead(tmpPath1stSegGroupNO) - 1;
			bool scoreStringBool_head; 
			vector<int> tmpNewMismatchPosVec_head;
			vector<char> tmpNewMismatchCharVec_head;
			//cout << "tmpUnfixedHeadLength: " << tmpUnfixedHeadLength << endl;
			int tmp_extension_length_head = 0;
			if(tmpUnfixedHeadLength > 0)
			{	
				if((tmpPathFinalMapPos - tmpUnfixedHeadLength - 1 < 1))
				{
					scoreStringBool_head = false;
				}
				else
				{
					tmp_extension_length_head = extensionBackwards_errorTolerated(tmpChrNameInt, 
						tmpChrPosInt-1, indexInfo, readSeq_inProcess, tmpUnfixedHeadLength, 1, 
						tmpNewMismatchPosVec_head, tmpNewMismatchCharVec_head, 
						minMatchBaseToSupport_perMismatch, minMatchBaseToSupport_twoMismatch);	
					if(tmp_extension_length_head == 0)
						scoreStringBool_head = false;
					else
						scoreStringBool_head = true;
				}

				if(scoreStringBool_head)
				{
					if(tmpUnfixedHeadLength - tmp_extension_length_head > 0)
					{
						Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength - tmp_extension_length_head, "S");					
						tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);	
					}
					Jump_Code tmpMatchJumpCode(tmp_extension_length_head, "M");
					tmpSpliceInfo->jump_code.push_back(tmpMatchJumpCode);	
					
					mismatchNumToAdd = mismatchNumToAdd + tmpNewMismatchPosVec_head.size();			
					//if(STORE_MISMATCH_POS)
					//{
						this->pushBack2FinalPathMismatchPosVec(tmpNewMismatchPosVec_head, tmpPath);
						if(STORE_MISMATCH_CHA)
						{
							this->pushBack2FinalPathMismatchCharVec(tmpNewMismatchCharVec_head, tmpPath);
						}
					//}
				}
				else
				{
					Jump_Code tmpHeadJumpCode(tmpUnfixedHeadLength - tmp_extension_length_head, "S");
					//Jump_Code tmpMatchJumpCode(tmp_extension_head, "M");
					tmpSpliceInfo->jump_code.push_back(tmpHeadJumpCode);	
					//tmpSpliceInfo->jump_code.push_back(tmpMatchJumpCode);				
				}
			}
			else
			{}

			tmpSpliceInfo->appendJumpCode(fixedPathVec[tmpPath].second);
			//////////////////////////// add last jump code /////////////////////////////////////////////////////
			int endMappedBaseMapPos = (fixedPathVec[tmpPath].second)->getEndBaseMapPos_jump_code(tmpPathFinalMapPos);
			//int readLength = segInfo->norSegmentLocInRead[segInfo->segmentNum - 1] + segInfo->norSegmentLength[segInfo->segmentNum - 1] - 1;
			int tmpUnfixedTailLength = readLength - (segInfo->returnSegmentLocInRead(tmpPathLastSegGroupNO) 
										+ segInfo->returnSegmentLength(tmpPathLastSegGroupNO) - 1);
			vector<int> tmpNewMismatchPosVec_tail;
			vector<char> tmpNewMismatchCharVec_tail;
			bool scoreStringBool_tail;
			int tmp_extension_length_tail;
			//cout << "tmpUnfixedTailLength: " << tmpUnfixedTailLength << endl;
			if(tmpUnfixedTailLength > 0)
			{
				if((endMappedBaseMapPos + tmpUnfixedTailLength > indexInfo->returnChromLength(tmpChrNameInt)))
				{
					scoreStringBool_tail = false;
				}
				else
				{
					tmp_extension_length_tail = extensionForwards_errorTolerated(tmpChrNameInt, endMappedBaseMapPos + 1,
						indexInfo, readSeq_inProcess, readLength - tmpUnfixedTailLength + 1, readLength,
						tmpNewMismatchPosVec_tail, tmpNewMismatchCharVec_tail, 
						minMatchBaseToSupport_perMismatch, minMatchBaseToSupport_twoMismatch);	
					if(tmp_extension_length_tail == 0)
						scoreStringBool_tail = false;
					else
						scoreStringBool_tail = true;
				}
				
				if(scoreStringBool_tail)
				{
					Jump_Code tmpMatchJumpCode(tmp_extension_length_tail, "M");
					tmpSpliceInfo->jump_code.push_back(tmpMatchJumpCode);
					if(tmpUnfixedTailLength - tmp_extension_length_tail > 0)
					{		
						Jump_Code tmpTailJumpCode(tmpUnfixedTailLength - tmp_extension_length_tail, "S");
						tmpSpliceInfo->jump_code.push_back(tmpTailJumpCode);		
					}
					mismatchNumToAdd = mismatchNumToAdd + tmpNewMismatchPosVec_tail.size();			
					//if(STORE_MISMATCH_POS)
					//{
						this->pushBack2FinalPathMismatchPosVec(tmpNewMismatchPosVec_tail, tmpPath);
						if(STORE_MISMATCH_CHA)
						{
							this->pushBack2FinalPathMismatchCharVec(tmpNewMismatchCharVec_tail, tmpPath);
						}
					//}
				}
				else
				{
					Jump_Code tmpTailJumpCode(tmpUnfixedTailLength, "S");
					tmpSpliceInfo->jump_code.push_back(tmpTailJumpCode);					
				}
			}
			else
			{}
			//cout << "after fixing head and tail ..." << endl;
			tmpSpliceInfo->getFinalJumpCode();			
			
			if(tmpUnfixedHeadLength > 0)
			{
				if(scoreStringBool_head)
				{
					tmpPathFinalMapPos = tmpPathFinalMapPos - tmp_extension_length_head;
				}
				else
				{}
			}

			int oldMismatchNum = fixedPathMismatchVec[tmpPath];
			int newMismatchNum = oldMismatchNum + mismatchNumToAdd;

			fixedPathMismatchVec[tmpPath] = newMismatchNum;

			finalPathVec.push_back(pair< pair<int, int>, Splice_Info*> (pair<int,int> (tmpPathFinalMapChr, tmpPathFinalMapPos), tmpSpliceInfo) );
		}		
		return;
	}


	void insertMismatchPosVec(vector<int>& mismatchPosVec, int fixGapVec_index)
	{
		int vecSize = mismatchPosVec.size();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			fixGapMismatchPosVec[fixGapVec_index].push_back(mismatchPosVec[tmp]);
		}
	}
	void insertMismatchCharVec(vector<char>& mismatchCharVec, int fixGapVec_index)
	{
		int vecSize = mismatchCharVec.size();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			fixGapMismatchCharVec[fixGapVec_index].push_back(mismatchCharVec[tmp]);
		}
	}
	void insertMismatchPosVec_fixMatchInfo(FixDoubleAnchor_Match_Info* fixMatchInfo, int fixGapVec_index)
	{
		int vecSize = fixMatchInfo->returnMismatchPosInReadVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixMatchInfo->returnMismatchPosInRead(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixMatchInfo(FixDoubleAnchor_Match_Info* fixMatchInfo, int fixGapVec_index)
	{
		int vecSize = fixMatchInfo->returnMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixMatchInfo->returnMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}	

	void insertMismatchPosVec_fixNwdpInfo(FixDoubleAnchor_NWDP_Info* fixNwdpInfo, int fixGapVec_index,
		int toProcessSeq_startLoc)
	{
		int vecSize = fixNwdpInfo->returnBestGapMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixNwdpInfo->returnBestGapMismatchPos(tmp) 
				+ toProcessSeq_startLoc - 1;
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixNwdpInfo(FixDoubleAnchor_NWDP_Info* fixNwdpInfo, int fixGapVec_index)
	{
		int vecSize = fixNwdpInfo->returnBestGapMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixNwdpInfo->returnBestGapMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}		
	}

	void insertMismatchPosVec_fixSpliceInfo(FixDoubleAnchor_Splice_Info* fixSpliceInfo, int fixGapVec_index)
	{
		int vecSize = fixSpliceInfo->returnBestSpliceMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixSpliceInfo->returnMismatchPosInRead(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixSpliceInfo(FixDoubleAnchor_Splice_Info* fixSpliceInfo, int fixGapVec_index)
	{
		int vecSize = fixSpliceInfo->returnBestSpliceMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixSpliceInfo->returnMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}	

	void insertMismatchPosVec_fixCirRnaInfo(FixDoubleAnchor_CirRNA_Info* fixCirRnaInfo, int fixGapVec_index)
	{
		int vecSize = fixCirRnaInfo->returnBestSpliceMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixCirRnaInfo->returnMismatchPosInRead(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixCirRnaInfo(FixDoubleAnchor_CirRNA_Info* fixCirRnaInfo, int fixGapVec_index)
	{
		int vecSize = fixCirRnaInfo->returnBestSpliceMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixCirRnaInfo->returnMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}	

	void insertMismatchPosVec_fixComplicatedSpliceInfo(
		FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo, int fixGapVec_index)
	{
		int vecSize = fixComplicateSpliceInfo->returnBestComplicateSpliceMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixComplicateSpliceInfo->returnMismatchPosInRead(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);
		}
	}
	void insertMismatchCharVec_fixComplicateSpliceInfo(FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo, int fixGapVec_index)
	{
		int vecSize = fixComplicateSpliceInfo->returnBestComplicateSpliceMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixComplicateSpliceInfo->returnMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}	

	void insertMismatchPosVec_fixInsInfo(FixDoubleAnchor_Insertion_Info* fixInsInfo, int fixGapVec_index)
	{
		int vecSize = fixInsInfo->returnBestInsertionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixInsInfo->returnBestInsertionMismatchPos(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);	
		}	
	}

	void insertMismatchCharVec_fixInsInfo(FixDoubleAnchor_Insertion_Info* fixInsInfo, int fixGapVec_index)
	{
		int vecSize = fixInsInfo->returnBestInsertionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixInsInfo->returnBestInsertionMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}

	void insertMismatchPosVec_fixDelInfo(FixDoubleAnchor_Deletion_Info* fixDelInfo, int fixGapVec_index)
	{
		int vecSize = fixDelInfo->returnBestDeletionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = fixDelInfo->returnBestDeletionMismatchPos(tmp);
			fixGapMismatchPosVec[fixGapVec_index].push_back(tmpMismatchPos);	
		}	
	}

	void insertMismatchCharVec_fixDelInfo(FixDoubleAnchor_Deletion_Info* fixDelInfo, int fixGapVec_index)
	{
		int vecSize = fixDelInfo->returnBestDeletionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = fixDelInfo->returnBestDeletionMismatchChar(tmp);
			fixGapMismatchCharVec[fixGapVec_index].push_back(tmpMismatchChar);
		}
	}

	int returnRepeatRegion_index()
	{
		return repeatRegion_index;
	}
	int returnMismatchNumInFixedPathMismatchVec(int tmpPath)
	{
		return fixedPathMismatchVec[tmpPath];
	}
	int returnMapChrNameIntInFinalPathVec(int tmpPath)
	{
		return (finalPathVec[tmpPath].first).first;
	}
	int returnMapChrPosIntInFinalPathVec(int tmpPath)
	{
		return (finalPathVec[tmpPath].first).second;
	}
	int returnFinalPathVecSize()
	{
		return finalPathVec.size();
	}
	Splice_Info* returnSpliceInfoInFinalPathVec(int tmpPath)
	{
		return finalPathVec[tmpPath].second;
	}
	void pushBackFixedPathVec(int tmpPathNO, Splice_Info* newPathSpliceInfo)
	{
		fixedPathVec.push_back(pair<int, Splice_Info*> (tmpPathNO, newPathSpliceInfo));
	}
	void pushBackFixedPathMismatchVec(int newPathMismatchNum)
	{
		fixedPathMismatchVec.push_back(newPathMismatchNum);
	}
	int returnPathSegSizeInPathVec_seg(int path_index)
	{
		return PathVec_seg[path_index].size();
	}
	int returnSegGroupNOinPathVec_seg(int path_index, int seg_index)
	{
		return ((PathVec_seg[path_index])[seg_index]).first;
	}
	int returnSegCandiNOinPathVec_seg(int path_index, int seg_index)
	{
		return ((PathVec_seg[path_index])[seg_index]).second;
	}
	int returnPathVecSegSize()
	{
		return PathVec_seg.size();
	}
	bool returnPathValidBoolVec(int index)
	{
		return PathValidBoolVec[index];
	}
	void pushBackPathFixedBoolVec(bool trueOrNot)
	{
		PathFixedBoolVec.push_back(trueOrNot);
	}

	void filterPath_incompleteHead(Seg_Info* segInfo)
	{
		int segmentNum = segInfo->returnSegmentNum();
		for(int tmpPath = 0; tmpPath < PathVec_seg.size(); tmpPath ++)
		{
			if(!PathValidBoolVec[tmpPath])
				continue;
			else
			{
				int tmpVecSize = PathVec_seg[tmpPath].size();
				int tmpLastSegNO = (PathVec_seg[tmpPath])[tmpVecSize - 1].first;
				if( ((segmentNum-1) != tmpLastSegNO) || (tmpVecSize == 1) )
				{
					PathValidBoolVec[tmpPath] = false;
				}
			}
		}
	}

	void filterPath_incompleteTail()
	{
		//int segmentNum = segInfo->segmentNum;
		for(int tmpPath = 0; tmpPath < PathVec_seg.size(); tmpPath ++)
		{
			if(!PathValidBoolVec[tmpPath])
				continue;
			else
			{
				//int tmpVecSize = PathVec_seg[tmpPath].size();
				int tmp1stSegNO = (PathVec_seg[tmpPath])[0].first;
				if( (tmp1stSegNO != 0) || ((PathVec_seg[tmpPath]).size() == 1) )
				{
					PathValidBoolVec[tmpPath] = false;
				}
			}
		}		
	}
 
 	bool foundInFixGapMap(int group_1, int candi_1, int group_2, int candi_2)
	{
		map< pair< pair<int,int>, pair<int,int> >, int >::iterator iter = fixGapMap.find(
			pair< pair<int,int>, pair<int, int> > (pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ) );
		if(iter != fixGapMap.end())
		{
			return true;
		}
		else
		{
			return false;
		}	
	}

 	bool foundInFixGapMap_index(int group_1, int candi_1, int group_2, int candi_2, int& index)
	{
		map< pair< pair<int,int>, pair<int,int> >, int >::iterator iter = fixGapMap.find(
			pair< pair<int,int>, pair<int, int> > (pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ) );
		if(iter != fixGapMap.end())
		{
			index = iter->second;
			return true;
		}
		else
		{
			return false;
		}	
	}

	bool gapCanBefixedOrNotBool(int group_1, int candi_1, int group_2, int candi_2)
	{
		map< pair< pair<int,int>, pair<int,int> >, int >::iterator iter = fixGapMap.find(
			pair< pair<int,int>, pair<int, int> > (pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ) );
		return fixGapVec[iter->second].first;
	}

 	void fixAllGap_cirRNA(Seg_Info* segInfo,
		Index_Info* indexInfo, const string& readSeq_inProcess,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo)
 	{
 		//cout << "start to do fixAllGap_cirRNA" << endl;
 		for(int tmpPath = 0; tmpPath < PathVec_seg.size(); tmpPath++)
 		{
 			for(int tmpGap = 0; tmpGap < PathVec_seg[tmpPath].size() - 1; tmpGap++)
 			{
 				int group_1 = (PathVec_seg[tmpPath])[tmpGap].first;
 				int candi_1 = (PathVec_seg[tmpPath])[tmpGap].second;
  				int group_2 = (PathVec_seg[tmpPath])[tmpGap+1].first;
 				int candi_2 = (PathVec_seg[tmpPath])[tmpGap+1].second;
 				//cout << "group_1: " << group_1 << endl;
 				//cout << "candi_1: " << candi_1 << endl;
  				//cout << "group_2: " << group_2 << endl;
 				//cout << "candi_2: " << candi_2 << endl;				
 				if(this->foundInFixGapMap(group_1, candi_1, group_2, candi_2))
 				{
 					continue;
 				}
 				else
 				{
 					//cout << "start to fix gap ...." << endl;
 					int fixGapVecSize = fixGapVec.size();
 					this->pushBackAndFixGap_cirRNA(group_1, candi_1, group_2, candi_2, fixGapVecSize, segInfo, indexInfo, 
 						readSeq_inProcess, annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
 				}
 			}
 		}
 	}
	
 	bool fixAllPath_cirRNA(Seg_Info* segInfo, Index_Info* indexInfo, const string& readSeq_inProcess, 
 		bool Do_extendHeadTail,
 		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo)
 	{
 		//cout << "start to do fixAllPath_cirRNA" << endl;
 		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

 		this->fixAllGap_cirRNA(segInfo, indexInfo, readSeq_inProcess, 
 			annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
		
		//bool fixGapInPathBool = false;
		int pathVecSize = this->PathVec_seg.size();

 		for(int tmpPath = 0; tmpPath < pathVecSize; tmpPath++)
 		{
			if(!(PathValidBoolVec[tmpPath]))
			{
				//newPathFixed = false;
				this->PathFixedBoolVec.push_back(false);				
				continue;
			}

			int firstSegGroupNO = (PathVec_seg[tmpPath])[0].first;
			int firstSegLength = segInfo->returnSegmentLength(firstSegGroupNO); //(segInfo->norSegmentLength)[firstSegGroupNO];


			Splice_Info* newPathSpliceInfo = new Splice_Info();
			Jump_Code firstJumpCode(firstSegLength, "M");
			newPathSpliceInfo->jump_code.push_back(firstJumpCode);

			int newPathMismatchNum = 0;
			
			bool newPathFixed = true;
			
			int tmpPathSegSize = PathVec_seg[tmpPath].size();

 			for(int tmpGap = 0; tmpGap < tmpPathSegSize - 1; tmpGap++)
 			{
 				int group_1 = (PathVec_seg[tmpPath])[tmpGap].first;
 				int candi_1 = (PathVec_seg[tmpPath])[tmpGap].second;
  				int group_2 = (PathVec_seg[tmpPath])[tmpGap+1].first;
 				int candi_2 = (PathVec_seg[tmpPath])[tmpGap+1].second;

 				bool foundInGapMap = this->foundInFixGapMap(group_1, candi_1, group_2, candi_2);
 				if(!foundInGapMap)
 				{
					newPathFixed = false;
					PathFixedBoolVec.push_back(newPathFixed);
					break;
 				}
 				else
 				{
					map< pair< pair<int,int>, pair<int,int> >, int >::iterator iter = fixGapMap.find(
						pair< pair<int,int>, pair<int, int> > (pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ) ); 					
 					int tmpIndex = iter->second;
 					bool tmpGapFixedOrNot = fixGapVec[tmpIndex].first;
 					if(!tmpGapFixedOrNot)
 					{
 						newPathFixed = false;
						PathFixedBoolVec.push_back(newPathFixed);
						break;					
 					}
 					else
 					{
 						int tmpMismatchNum = (fixGapVec[tmpIndex].second).first;
 						newPathMismatchNum += tmpMismatchNum;
 						int tmpJumpCodeSize = ((fixGapVec[tmpIndex].second).second).size();
 						for(int tmpJumpCodeNO = 0; tmpJumpCodeNO < tmpJumpCodeSize; tmpJumpCodeNO++)
 						{
 							newPathSpliceInfo->jump_code.push_back(((fixGapVec[tmpIndex].second).second)[tmpJumpCodeNO]);
 						} 
 					}

 				}
 			}

 			if(newPathFixed)
 			{
				newPathSpliceInfo->getFinalJumpCode();	
				//cout << "finish getting final jumpcode" << endl;
				bool allJumpCodeValidBool = newPathSpliceInfo->allFinalJumpCodeValid_cirRNA();
				//cout << "allJumpCodeValidBool" << allJumpCodeValidBool << endl;
				if(allJumpCodeValidBool)
				{
					(PathFixedBoolVec).push_back(allJumpCodeValidBool);
					(fixedPathVec).push_back(pair <int, Splice_Info*> (tmpPath, newPathSpliceInfo) );
					(fixedPathMismatchVec).push_back(newPathMismatchNum);
				}
				else
				{
					(PathFixedBoolVec).push_back(allJumpCodeValidBool);
				} 				
 			}
 		} 		

		int readLength = readSeq_inProcess.length();
		this->getFinalPath_extend2HeadTail_new(indexInfo, segInfo, readLength, readSeq_inProcess, Do_extendHeadTail);

		//cout << "finish getting finalPath" << endl;
		//fixGapInPathBool = true;
		return true;

 	}

	void pushBackAndFixGap(int group_1, int candi_1, int group_2, int candi_2, int currentFixGapVecSize, Seg_Info* segInfo,
		Index_Info* indexInfo, const string& readSeq_inProcess, 
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
		int spliceJunctionDistanceMax)
	{
		fixGapMap.insert(pair< pair< pair<int,int>, pair<int,int> >, int > (
			pair< pair<int,int>, pair<int, int> > (
				pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ), 
				currentFixGapVecSize) );

		bool gapFixedOrNotBool = false;
		int tmpMismatchNum = 0;
		vector<Jump_Code> tmpJumpCodeVec;

		int tmpSegGroupNO = group_1;
		int tmpSegCandiNO = candi_1;

		int tmpSegGroupNO_next = group_2;
		int tmpSegCandiNO_next = candi_2;

		int tmpRelation = segInfo->checkSegRelation(tmpSegGroupNO, tmpSegCandiNO, 
			tmpSegGroupNO_next, tmpSegCandiNO_next, spliceJunctionDistanceMax);

		int tmpSegmentLocInRead_1 = segInfo->returnSegmentLocInRead(tmpSegGroupNO); //(segInfo->norSegmentLocInRead)[tmpSegGroupNO];
		int tmpSegmentLocInRead_2 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_next); //(segInfo->norSegmentLocInRead)[tmpSegGroupNO_next];
		int tmpSegmentLength_1 = segInfo->returnSegmentLength(tmpSegGroupNO); //(segInfo->norSegmentLength)[tmpSegGroupNO];
		int tmpSegmentLength_2 = segInfo->returnSegmentLength(tmpSegGroupNO_next); //(segInfo->norSegmentLength)[tmpSegGroupNO_next];

		unsigned int tmpSegmentMapPosInWholeGenome_1 
			//= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO * CANDALILOC + tmpSegCandiNO);
			= segInfo->returnSegmentMapPos(tmpSegGroupNO, tmpSegCandiNO);
		unsigned int tmpSegmentMapPosInWholeGenome_2 
			//= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO_next * CANDALILOC + tmpSegCandiNO_next);
			= segInfo->returnSegmentMapPos(tmpSegGroupNO_next, tmpSegCandiNO_next);

				unsigned int tmpChrNameInt, tmpChrPosInt;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_1, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_1 = tmpChrPosInt;
				//cout << "...... tmpChrMapPos_1: " << tmpSegmentMapPos_1 << endl;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_2, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_2 = tmpChrPosInt;
				
				string tmpChrNameStr;
				if(tmpChrNameStr_1 == tmpChrNameStr_2)
				{
					tmpChrNameStr = tmpChrNameStr_1;
					//cout << "...... tmpChrName: " << tmpChrNameStr << endl;
					gapFixedOrNotBool = false;
					fixGapVec.push_back(pair<bool, pair<int, vector<Jump_Code> > > (false,
						pair<int, vector<Jump_Code> >(tmpMismatchNum, tmpJumpCodeVec)));
					//if(STORE_MISMATCH_POS)
					//{
						vector<int> mismatchPosVec;
						fixGapMismatchPosVec.push_back(mismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							vector<char> mismatchCharVec;
							fixGapMismatchCharVec.push_back(mismatchCharVec);
						//}
					//}
				}
				else
				{
					//cout << "...... different chrName " << endl;
					gapFixedOrNotBool = false;
					fixGapVec.push_back(pair<bool, pair<int, vector<Jump_Code> > > (false,
						pair<int, vector<Jump_Code> >(tmpMismatchNum, tmpJumpCodeVec)));
					//if(STORE_MISMATCH_POS)
					//{
						vector<int> mismatchPosVec;
						fixGapMismatchPosVec.push_back(mismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							vector<char> mismatchCharVec;
							fixGapMismatchCharVec.push_back(mismatchCharVec);
						//}
					//}
					return;
				}


		this->fixDoubleAnchor_extendBack(tmpRelation, tmpSegmentLocInRead_1, tmpSegmentLocInRead_2,
				tmpSegmentLength_1, tmpSegmentLength_2, tmpSegmentMapPos_1, tmpSegmentMapPos_2, 
				readSeq_inProcess, indexInfo, tmpChrNameStr, currentFixGapVecSize, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);

	}

	void pushBackAndFixGap_cirRNA(int group_1, int candi_1, int group_2, int candi_2, int currentFixGapVecSize, Seg_Info* segInfo,
		Index_Info* indexInfo, const string& readSeq_inProcess,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo)
	{
		fixGapMap.insert(pair< pair< pair<int,int>, pair<int,int> >, int > (
			pair< pair<int,int>, pair<int, int> > (
				pair<int, int>(group_1, candi_1), pair<int,int>(group_2, candi_2) ), 
				currentFixGapVecSize) );

		bool gapFixedOrNotBool = false;
		int tmpMismatchNum = 0;
		vector<Jump_Code> tmpJumpCodeVec;

		int tmpSegGroupNO = group_1;
		int tmpSegCandiNO = candi_1;

		int tmpSegGroupNO_next = group_2;
		int tmpSegCandiNO_next = candi_2;

		int tmpRelation = segInfo->checkSegRelation_cirRNA(tmpSegGroupNO, tmpSegCandiNO, tmpSegGroupNO_next, tmpSegCandiNO_next);

		int tmpSegmentLocInRead_1 = segInfo->returnSegmentLocInRead(tmpSegGroupNO); //(segInfo->norSegmentLocInRead)[tmpSegGroupNO];
		int tmpSegmentLocInRead_2 = segInfo->returnSegmentLocInRead(tmpSegGroupNO_next); //(segInfo->norSegmentLocInRead)[tmpSegGroupNO_next];
		int tmpSegmentLength_1 = segInfo->returnSegmentLength(tmpSegGroupNO); //(segInfo->norSegmentLength)[tmpSegGroupNO];
		int tmpSegmentLength_2 = segInfo->returnSegmentLength(tmpSegGroupNO_next); //(segInfo->norSegmentLength)[tmpSegGroupNO_next];

		unsigned int tmpSegmentMapPosInWholeGenome_1 
			//= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO * CANDALILOC + tmpSegCandiNO);
			= segInfo->returnSegmentMapPos(tmpSegGroupNO, tmpSegCandiNO);
		unsigned int tmpSegmentMapPosInWholeGenome_2 
			//= *(segInfo->norSegmentAlignLoc + tmpSegGroupNO_next * CANDALILOC + tmpSegCandiNO_next);
			= segInfo->returnSegmentMapPos(tmpSegGroupNO_next, tmpSegCandiNO_next);

				unsigned int tmpChrNameInt, tmpChrPosInt;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_1, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_1 = tmpChrPosInt;
				//cout << "...... tmpChrMapPos_1: " << tmpSegmentMapPos_1 << endl;
				indexInfo->getChrLocation(tmpSegmentMapPosInWholeGenome_2, &tmpChrNameInt, &tmpChrPosInt);
				string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpSegmentMapPos_2 = tmpChrPosInt;
				
				string tmpChrNameStr;
				if(tmpChrNameStr_1 == tmpChrNameStr_2)
				{
					tmpChrNameStr = tmpChrNameStr_1;
					//cout << "...... tmpChrName: " << tmpChrNameStr << endl;
					gapFixedOrNotBool = false;
					fixGapVec.push_back(pair<bool, pair<int, vector<Jump_Code> > > (false,
						pair<int, vector<Jump_Code> >(tmpMismatchNum, tmpJumpCodeVec)));
					if(STORE_MISMATCH_POS)
					{
						vector<int> mismatchPosVec;
						fixGapMismatchPosVec.push_back(mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							vector<char> mismatchCharVec;
							fixGapMismatchCharVec.push_back(mismatchCharVec);
						}
					}
				}
				else
				{
					//cout << "...... different chrName " << endl;
					gapFixedOrNotBool = false;
					fixGapVec.push_back(pair<bool, pair<int, vector<Jump_Code> > > (false,
						pair<int, vector<Jump_Code> >(tmpMismatchNum, tmpJumpCodeVec)));
					
					if(STORE_MISMATCH_POS)
					{
						vector<int> mismatchPosVec;
						fixGapMismatchPosVec.push_back(mismatchPosVec);
						if(STORE_MISMATCH_CHA)
						{
							vector<char> mismatchCharVec;
							fixGapMismatchCharVec.push_back(mismatchCharVec);
						}
					}
					return;
				}

		this->fixDoubleAnchor_extendBack_cirRNA(tmpRelation, tmpSegmentLocInRead_1, tmpSegmentLocInRead_2,
				tmpSegmentLength_1, tmpSegmentLength_2, tmpSegmentMapPos_1, tmpSegmentMapPos_2, 
				readSeq_inProcess, indexInfo, tmpChrNameStr, currentFixGapVecSize, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);

	}

	void fixDoubleAnchor_extendBack(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, 
		Index_Info* indexInfo, const string& chromName, int index_fixGapVec,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo)
	{

		/*cout << "fixDoubleAnchor_extendBack starts!" << endl << " seg_1_start: " << segmentLocInRead_1 << 
			" seg_2_start: " << segmentLocInRead_2 << " seg_1_mapPos: " << segmentMapPos_1 
			<< " seg_2_mapPos: " << segmentMapPos_2 << endl;*/
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
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Match(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_INSERTION_NEIGHBOUR) || (relation == FIX_INSERTION_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Insertion(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_DELETION_NEIGHBOUR) || (relation == FIX_DELETION_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Deletion(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_SPLICE_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Splice(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
		}
		else
		{
			#ifdef DETECT_CIRCULAR_RNA
			if((relation == FIX_CIRCULAR_RNA_NEIGHBOUR) || (relation == FIX_CIRCULAR_RNA_GAP))
			{
				this->fixDoubleAnchor_Circular_RNA(relation, segmentLocInRead_1, segmentLocInRead_2,
					segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
					indexInfo, chrNameInt, index_fixGapVec);
			}
			else
				cout << "error in fixDoubleAnchor ... " << endl;
			#else
			cout << "error in fixDoubleAnchor ... " << endl;
			#endif
		}

		return;// fixDoubleAnchorBool;

	}

	void fixDoubleAnchor_extendBack_cirRNA(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, 
		Index_Info* indexInfo, const string& chromName, int index_fixGapVec,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo
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
		if(relation == FIX_CIRCULAR_RNA)
		{
			this->fixDoubleAnchor_cirRNA(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);			
		}
		else if(relation == FIX_MATCH)
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Match(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_INSERTION_NEIGHBOUR) || (relation == FIX_INSERTION_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Insertion(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_DELETION_NEIGHBOUR) || (relation == FIX_DELETION_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Deletion(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec);
		}
		else if((relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_SPLICE_GAP))
		{
			//fixDoubleAnchorBool = 
			this->fixDoubleAnchor_Splice(relation, segmentLocInRead_1, segmentLocInRead_2,
				segmentLength_1, segmentLength_2, segmentMapPos_1, segmentMapPos_2, readSeq_inProcess, 
				indexInfo, chrNameInt, index_fixGapVec, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
		}
		else
		{
			cout << "error in fixDoubleAnchor ... " << endl;
		}
		return;// fixDoubleAnchorBool;
	}

	void fixDoubleAnchor_cirRNA(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec
		)
	{
		//cout << "start to fix splice" << endl;

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

		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << " spliceJunctionLength: " << spliceJunctionLength << endl;
		//cout << "segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1: " << segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1 << endl;
		string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left - 1, subSeqLengthInProcess);
		//cout << "segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1: " << segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1 << endl;
		//string left_chrom_seq = (indexInfo->returnChromStr(chrNameInt)).substr(segmentMapPos_1 + segmentLength_1 - tmpBuffer_left - 1, chromSubSeqLengthInProcess);
		string left_chrom_seq = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, chromSubSeqLengthInProcess);
		//cout << "segmentMapPos_2 - 1 - chromSubSeqLengthInProcess: " << segmentMapPos_2 - 1 - chromSubSeqLengthInProcess << endl;
		//string right_chrom_seq = (indexInfo->chromStr[chrNameInt]).substr(segmentMapPos_2 + tmpBuffer_right - 1 - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess); 
		string right_chrom_seq = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_2 + tmpBuffer_right - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);

		size_t prefix_length = 0;
		size_t max_double_splice_mismatch = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 1;
		size_t mismatch_bits = 0;
		size_t comb_bits = 0;
		//bool adjacent_segments = false;
		bool double_anchor_noncanonical = DETECT_NONCANONICAL_SJ; //false;//true;//false;//DO_NONCANONICAL; ////debug
		string flank_seq;
		GenomeScan* genome_scan = new GenomeScan;
		bool splice_fixed = (*genome_scan).Double_anchored_score(readSubSeqInProcess, left_chrom_seq, right_chrom_seq, prefix_length, 
			max_double_splice_mismatch, comb_bits,
		 	//(!adjacent_segments) && 
			double_anchor_noncanonical, flank_seq, mismatch_bits);
		//cout << "splice_fixed: " << endl;
		if(splice_fixed)
		{
			int firstMatchLength = //segmentLength_1 
				- tmpBuffer_left + prefix_length;
			int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
			//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

			Jump_Code firstMatchJumpCode(firstMatchLength, "M");
			Jump_Code spliceJumpCode(spliceJunctionLength, "N");
			Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

			fixGapVec[index_fixGapVec].first = true;
			(fixGapVec[index_fixGapVec].second).first = mismatch_bits;
			((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
			((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
			((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);
			
			/*if(STORE_MISMATCH_POS)
			{
				this->insertMismatchPosVec()
			}*/
		}
		else
		{
			fixGapVec[index_fixGapVec].first = false;
		}
		delete genome_scan;
		//fixDoubleAnchorBool_Splice = splice_fixed;
		return;// fixDoubleAnchorBool_Splice;
	}

	void insertMismatchPosVec_fixMatchInfo_SmallGap(int index_fixGapVec, int startPosInRead, int toProcessSeq_len)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			int tmpMismatchPos = startPosInRead + tmp;
			fixGapMismatchPosVec[index_fixGapVec].push_back(tmpMismatchPos);
		}
	}

	void insertMismatchCharVec_fixMatchInfo_SmallGap(int index_fixGapVec, int toProcessSeq_len, 
		int chrNameInt, int startPosInChrom, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			char tmpMismatchChar = indexInfo->returnOneBaseCharInGenome(chrNameInt, startPosInChrom + tmp);
			//int tmpMismatchPos = startPosInRead + tmp;
			fixGapMismatchCharVec[index_fixGapVec].push_back(tmpMismatchChar);
		}
	}

	void fixDoubleAnchor_Match(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec
		)
	{
		//cout << endl << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl
		//	<< "start to do fixDoubleAnchor_Match" << endl;

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		
		//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl; 
		
		if(subSeqLengthInProcess < 2)
		{
			
			Jump_Code matchJumpCode(//segmentLength_1 + 
				subSeqLengthInProcess + segmentLength_2, "M");
		
			fixGapVec[index_fixGapVec].first = true;
			(fixGapVec[index_fixGapVec].second).first = subSeqLengthInProcess;
			((fixGapVec[index_fixGapVec].second).second).push_back(matchJumpCode);
			
			this->insertMismatchPosVec_fixMatchInfo_SmallGap(index_fixGapVec, segmentLocInRead_1 + segmentLength_1,
					subSeqLengthInProcess);
			this->insertMismatchCharVec_fixMatchInfo_SmallGap(index_fixGapVec, 
						subSeqLengthInProcess, chrNameInt, segmentMapPos_1 + segmentLength_1, indexInfo);
		}
		else
		{
			
			string readSubSeqInProcess = readSeq_inProcess.substr(segmentLocInRead_1 + segmentLength_1 - 1,
				subSeqLengthInProcess);
			string chromSubSeqInProcess = indexInfo->returnChromStrSubstr(chrNameInt, segmentMapPos_1 + segmentLength_1, subSeqLengthInProcess);
			
			//cout << "segStartInRead: " << segmentLocInRead_1 + segmentLength_1 << endl;
			//cout << "scoreStringBool: " << scoreStringBool << endl;
			int max_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;
			FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
			bool scoreStringBool = fixMatchInfo->fixMatch(readSubSeqInProcess, chromSubSeqInProcess,
				max_mismatch, segmentLocInRead_1 + segmentLength_1);
			if(scoreStringBool)
			{
				Jump_Code matchJumpCode(//segmentLength_1 + 
					subSeqLengthInProcess + segmentLength_2, "M");

				fixGapVec[index_fixGapVec].first = true;
				(fixGapVec[index_fixGapVec].second).first = fixMatchInfo->returnMismatchNum();//mismatch_bits;
				((fixGapVec[index_fixGapVec].second).second).push_back(matchJumpCode);	
				this->insertMismatchPosVec_fixMatchInfo(fixMatchInfo, index_fixGapVec);
				this->insertMismatchCharVec_fixMatchInfo(fixMatchInfo, index_fixGapVec);			
			}
			else // score string failed, insert sudo-match jump code
			{
				//cout << "fix-match failed !" << endl;
				//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
				fixGapVec[index_fixGapVec].first = false;	
			}
			delete fixMatchInfo;
		}
		//cout << "end of doing fixDoubleAnchor_Match " << endl << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl << endl;
		//return fixDoubleAnchorBool_Match;
	}

	void fixDoubleAnchor_Insertion(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chromNameInt, int index_fixGapVec
		)
	{
		// cout << "start to fix Insertion" << endl;
		// cout << "segmentLocInRead_1: " << segmentLocInRead_1 << endl;
		// cout << "segmentLocInRead_2: " << segmentLocInRead_2 << endl;
		// cout << "segmentLength_1: " << segmentLength_1 << endl;
		// cout << "segmentLength_2: " << segmentLength_2 << endl;
		// cout << "segmentMapPos_1: " << segmentMapPos_1 << endl;
		// cout << "segmentMapPos_2: " << segmentMapPos_2 << endl; 

		bool fixDoubleAnchorBool_Insertion_bool = false;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
		int insertionLength = (segmentMapPos_1 - segmentLocInRead_1) - (segmentMapPos_2 - segmentLocInRead_2);
		//int chrNameInt = indexInfo->convertStringToInt(chromName);

		// cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
		// cout << "insertionLength: " << insertionLength << endl;

		if(subSeqLengthInProcess <= insertionLength)
		{
			int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - insertionLength;
			if(secondMatchLength > 0)
			{
				Jump_Code midInsertionJumpCode(insertionLength, "I");	
				Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess - insertionLength, "M");							

				fixDoubleAnchorBool_Insertion_bool = true;
				fixGapVec[index_fixGapVec].first = true;
				(fixGapVec[index_fixGapVec].second).first = 0;
				((fixGapVec[index_fixGapVec].second).second).push_back(midInsertionJumpCode);
				((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);			
				// no mismatch pos or char generated 		
			}
			else
			{
				fixDoubleAnchorBool_Insertion_bool = false;
				fixGapVec[index_fixGapVec].first = false;
			}
		}	  
		else
		{
			FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info();
			int tmp_toFix_insertion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_insertion_read_end = segmentLocInRead_2 - 1;
			int tmp_toFix_insertion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_insertion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;  // how to set the max-mismatch parameter

			bool insertion_fixed = insInfo->detectBestInsertion_lessMismatch(
				tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
				tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
				chromNameInt, tmp_max_allowed_mismatchNum);

			if(insertion_fixed)
			{
				int firstMatchLength = insInfo->return_best_insertion_prefix_match_length();
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - firstMatchLength - insertionLength;
				int mismatch_bits = insInfo->return_best_insertion_mismatch();
				if(secondMatchLength > 0)
				{
					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code insertionJumpCode(insertionLength, "I");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M");

					fixDoubleAnchorBool_Insertion_bool = true;
					fixGapVec[index_fixGapVec].first = true;
					(fixGapVec[index_fixGapVec].second).first = mismatch_bits;
					((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
					((fixGapVec[index_fixGapVec].second).second).push_back(insertionJumpCode);		
					((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);	

					this->insertMismatchPosVec_fixInsInfo(insInfo, index_fixGapVec);
					this->insertMismatchCharVec_fixInsInfo(insInfo, index_fixGapVec);	
				}	
				else
				{
					fixDoubleAnchorBool_Insertion_bool = false;
					fixGapVec[index_fixGapVec].first = false;
				}			
			}
			else
			{
				//cout << "fix-insertion failed !" << endl;
				//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
				fixDoubleAnchorBool_Insertion_bool = false;
				fixGapVec[index_fixGapVec].first = false;			
			}
			delete insInfo;
		}

		//cout << "fixDoubleAnchorBool_Insertion_bool: " << fixDoubleAnchorBool_Insertion_bool << endl;
		if(!fixDoubleAnchorBool_Insertion_bool)
		{
			FixDoubleAnchor_NWDP_Info* nwdpInfo = new FixDoubleAnchor_NWDP_Info();
			int tmp_toFix_insertion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_insertion_read_end = segmentLocInRead_2 - 1;
			int tmp_toFix_insertion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_insertion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;  // how to set the max-mismatch parameter
			//cout << "start to do fix_NWDP ..." << endl;
			bool fix_NWDP_bool = nwdpInfo->fix_NWDP(
				tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
				tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
				chromNameInt, tmp_max_allowed_mismatchNum);
			//cout << "end of fixing NWDP ..." << endl;
			//cout << "fix_NWDP_bool: " << fix_NWDP_bool << endl;
			if(fix_NWDP_bool)
			{
				fixGapVec[index_fixGapVec].first = true;
				(fixGapVec[index_fixGapVec].second).first = nwdpInfo->returnBestGapMismatchPosVecSize();
				nwdpInfo->copyJumpCodeVec2TargetVec(((fixGapVec[index_fixGapVec].second).second));
				Jump_Code secondSegMatchJumpCode(segmentLength_2, "M");
				((fixGapVec[index_fixGapVec].second).second).push_back(secondSegMatchJumpCode);
				this->insertMismatchPosVec_fixNwdpInfo(nwdpInfo, index_fixGapVec, tmp_toFix_insertion_read_start);
				this->insertMismatchCharVec_fixNwdpInfo(nwdpInfo, index_fixGapVec);
			}
			else
			{
				fixGapVec[index_fixGapVec].first = false;
			}
			delete nwdpInfo;
		}

		return;// fixDoubleAnchorBool_Insertion;
	}
	
	void insertMismatchPosVec_fixDeletionInfo_SmallGap(int index_fixGapVec, int startPosInRead, int toProcessSeq_len)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			int tmpMismatchPos = startPosInRead + tmp;
			fixGapMismatchPosVec[index_fixGapVec].push_back(tmpMismatchPos);
		}
	}

	void insertMismatchCharVec_fixDeletionInfo_SmallGap(int index_fixGapVec, int toProcessSeq_len, 
		int chrNameInt, int startPosInChrom, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < toProcessSeq_len; tmp++)
		{
			char tmpMismatchChar = indexInfo->returnOneBaseCharInGenome(chrNameInt, startPosInChrom + tmp);
			//int tmpMismatchPos = startPosInRead + tmp;
			fixGapMismatchCharVec[index_fixGapVec].push_back(tmpMismatchChar);
		}
	}

	void fixDoubleAnchor_Deletion(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec
		)
	{
		bool fixDoubleAnchor_deletion_bool = false;
		//int buffer_left = 2, buffer_right = 2;
		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1;
			//+ buffer_left + buffer_right;
		int deletionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);

		if(subSeqLengthInProcess < 2)
		{
			Jump_Code deletionJumpCode(deletionLength, "D");
			Jump_Code secondMatchJumpCode(segmentLength_2 + subSeqLengthInProcess, "M");

			fixDoubleAnchor_deletion_bool = true;
			fixGapVec[index_fixGapVec].first = true;

			(fixGapVec[index_fixGapVec].second).first = subSeqLengthInProcess;
			((fixGapVec[index_fixGapVec].second).second).push_back(deletionJumpCode);
			((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);	

			this->insertMismatchPosVec_fixDeletionInfo_SmallGap(index_fixGapVec, segmentLocInRead_1 + segmentLength_1,
				subSeqLengthInProcess);
			this->insertMismatchCharVec_fixDeletionInfo_SmallGap(index_fixGapVec, subSeqLengthInProcess, 
				chrNameInt, segmentMapPos_2 - subSeqLengthInProcess, indexInfo);
		}
		else
		{	
			//cout << "start to do fixDeletion ...." << endl;
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
			//cout << "deletion_fixed: " << deletion_fixed << endl;
			if(deletion_fixed)
			{
				int firstMatchLength = //segmentLength_1 + 
						//prefix_length;
						delInfo->return_best_deletion_prefix_match_length();
				int secondMatchLength = segmentLength_2 + subSeqLengthInProcess - firstMatchLength;
				int mismatch_bits = delInfo->return_best_deletion_mismatch();
				Jump_Code firstMatchJumpCode(firstMatchLength
					//- buffer_left
					, "M");
				Jump_Code deletionJumpCode(deletionLength, "D");
				Jump_Code secondMatchJumpCode(secondMatchLength
					//- buffer_right
					, "M");

				fixDoubleAnchor_deletion_bool = true;
				fixGapVec[index_fixGapVec].first = true;
				(fixGapVec[index_fixGapVec].second).first = mismatch_bits;
				((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
				((fixGapVec[index_fixGapVec].second).second).push_back(deletionJumpCode);	
				((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);	
				
				this->insertMismatchPosVec_fixDelInfo(delInfo, index_fixGapVec);
				this->insertMismatchCharVec_fixDelInfo(delInfo, index_fixGapVec);	
			}
			else
			{
				//cout << "fix-deletion failed !" << endl;
				//cout << "subSeqLengthInProcess: " << subSeqLengthInProcess << endl;
				fixDoubleAnchor_deletion_bool = false;
				fixGapVec[index_fixGapVec].first = false;			
			}
			delete delInfo;
		}
		if(!fixDoubleAnchor_deletion_bool)
		{
			FixDoubleAnchor_NWDP_Info* nwdpInfo = new FixDoubleAnchor_NWDP_Info();
			int tmp_toFix_deletion_read_start = segmentLocInRead_1 + segmentLength_1;
			int tmp_toFix_deletion_read_end = segmentLocInRead_2 - 1;
			//int subSeqLengthInProcess = tmp_toFix_deletion_read_end - tmp_toFix_deletion_read_start + 1;
			int tmp_toFix_deletion_chrom_start = segmentMapPos_1 + segmentLength_1;
			int tmp_toFix_deletion_chrom_end = segmentMapPos_2 - 1;
			int tmp_max_allowed_mismatchNum = subSeqLengthInProcess/LengthOfSeqPerMismatchAllowed + 1;			
			bool fix_NWDP_bool = nwdpInfo->fix_NWDP(
				tmp_toFix_deletion_read_start, tmp_toFix_deletion_read_end,
				tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, readSeq_inProcess, indexInfo,
				chrNameInt, tmp_max_allowed_mismatchNum);
			if(fix_NWDP_bool)
			{
				fixGapVec[index_fixGapVec].first = true;
				(fixGapVec[index_fixGapVec].second).first = nwdpInfo->returnBestGapMismatchPosVecSize();
				nwdpInfo->copyJumpCodeVec2TargetVec(((fixGapVec[index_fixGapVec].second).second));
				Jump_Code secondSegMatchJumpCode(segmentLength_2, "M");
				((fixGapVec[index_fixGapVec].second).second).push_back(secondSegMatchJumpCode);
				this->insertMismatchPosVec_fixNwdpInfo(nwdpInfo, index_fixGapVec, tmp_toFix_deletion_read_start);
				this->insertMismatchCharVec_fixNwdpInfo(nwdpInfo, index_fixGapVec);
			}
			else
			{
				fixGapVec[index_fixGapVec].first = false;
			}
			delete nwdpInfo;
		}

		return;// fixDoubleAnchorBool_Deletion;
	}

	void fixDoubleAnchor_Splice(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec, 
		bool annotation_provided_bool, 
		bool Do_annotation_only_bool, Annotation_Info* annotationInfo)
	{
		// cout << endl << "**************************************" << endl;
		// cout << "start to fix double anchor splice ...." << endl;
		// cout << "segmentLocInRead_1: " << segmentLocInRead_1 << endl;
		// cout << "segmentLocInRead_2: " << segmentLocInRead_2 << endl;
		// cout << "segmentMapPos_1: " << segmentMapPos_1 << endl;
		// cout << "segmentMapPos_2: " << segmentMapPos_2 << endl;
		int tmpBuffer_left = FixSpliceBuffer;//+1;
		if(tmpBuffer_left > segmentLength_1 - 2) //anchor >= 2
		{
			tmpBuffer_left = segmentLength_1 - 2;
		}
		int tmpBuffer_right = FixSpliceBuffer;//+1;
		if(tmpBuffer_right > segmentLength_2 - 2)
		{
			tmpBuffer_right = segmentLength_2 - 2;
		}

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1 + tmpBuffer_left + tmpBuffer_right;
		int spliceJunctionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);

		//size_t prefix_length = 0;
		size_t max_double_splice_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;

		if((annotation_provided_bool && Do_annotation_only_bool) || !(annotation_provided_bool))
		{
			//cout << "start to do FixDoubleAnchor_Splice_Info" << endl;
			FixDoubleAnchor_Splice_Info* fixSpliceInfo = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo->detectBestSpliceSite(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
				//cout << "splice_fixed: " << splice_fixed << endl;
				if(splice_fixed && fixSpliceInfo->fixSpliceResultConfident())//((fixSpliceInfo->returnBestSplice_canonicalOrNot()) || (fixSpliceInfo->returnBestSplice_mismatchNum() <= 2) ) )
				{
					int prefix_length = fixSpliceInfo->returnBestSplice_prefixMatchLength();
					int firstMatchLength = //segmentLength_1 
						- tmpBuffer_left + prefix_length;
					int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;
					//cout << "prefix_length: " << prefix_length << ", spliceJunctionLength: " << spliceJunctionLength << endl;

					Jump_Code firstMatchJumpCode(firstMatchLength, "M");
					Jump_Code spliceJumpCode(spliceJunctionLength, "N");
					Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

					fixGapVec[index_fixGapVec].first = true;
					(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo->returnBestSplice_mismatchNum();
					((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
					((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
					((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

					//fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
					this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);

					this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
				}
				else
				{
					//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
					FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo = new FixDoubleAnchor_Splice_Complicate_Info();
					bool complicate_splice_fixed = fixComplicateSpliceInfo->detectComplicateSplice(
						segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
						segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
						readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
					//cout << "complicate_splice_fixed: " << complicate_splice_fixed << endl;
					//cout << "finish fixing complicated SJ ....." << endl;	
					bool complicated_splice_fixed_better_bool = false;
					if(complicate_splice_fixed)
					{
						complicated_splice_fixed_better_bool = fixComplicateSpliceInfo->compared2SpliceInfo(splice_fixed,
							fixSpliceInfo);
					}
					else
					{}
					//cout << "complicated_splice_fixed_better_bool: " << complicated_splice_fixed_better_bool << endl;
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

						fixGapVec[index_fixGapVec].first = true;
						(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo->return_mismatch_bestComplicatedSplice();
						((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

						//fixComplicateSpliceInfo->generateBestComplicateSJMismatchVec(
						//	segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1);
						this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo, index_fixGapVec);

						this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo, index_fixGapVec);
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

						fixGapVec[index_fixGapVec].first = true;
						(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo->returnBestSplice_mismatchNum();
						((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						//fixSpliceInfo->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);

						this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo, index_fixGapVec);
					}
					else
					{
						fixGapVec[index_fixGapVec].first = false;
					}
					delete fixComplicateSpliceInfo;
				}

			delete fixSpliceInfo;
		}
		else // annotation_provided_bool == true && Do_annotation_only_bool == false
		{ 
			bool annotated_splice_fixed_bool = false;
			FixDoubleAnchor_Splice_Info* fixSpliceInfo_annotated = new FixDoubleAnchor_Splice_Info();
			bool splice_fixed = fixSpliceInfo_annotated->detectBestSpliceSite//_prefer_canonical_lessMismatch
				(
				segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
				segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
				readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
				annotation_provided_bool, //Do_annotation_only_bool, 
				true, annotationInfo);
			//if(FIX_INDEL_AROUND_SJ_BOOL)
			//{
				if(splice_fixed && 
					((fixSpliceInfo_annotated->returnBestSplice_canonicalOrNot()) 
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

					fixGapVec[index_fixGapVec].first = true;					
					// fixed with annotated SJ
					annotated_splice_fixed_bool = true;

					(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
					((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
					((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
					((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

					this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
					this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
				}
				else
				{
					/////  test fixing complicated SJ  ///////
					//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
					FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_annotated
						= new FixDoubleAnchor_Splice_Complicate_Info();
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

						fixGapVec[index_fixGapVec].first = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;

						(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo_annotated->return_mismatch_bestComplicatedSplice();
						((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

						this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_annotated, index_fixGapVec);
						this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_annotated, index_fixGapVec);
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

						fixGapVec[index_fixGapVec].first = true;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = true;			

						(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_annotated->returnBestSplice_mismatchNum();
						((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
						this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_annotated, index_fixGapVec);
					}
					else
					{
						fixGapVec[index_fixGapVec].first = false;
						// fixed with annotated SJ
						annotated_splice_fixed_bool = false;	
					}
					delete fixComplicateSpliceInfo_annotated;
				}
			delete fixSpliceInfo_annotated;

			if(annotated_splice_fixed_bool) // fixed with annotated SJ
			{}
			else
			{
				FixDoubleAnchor_Splice_Info* fixSpliceInfo_unannotated = new FixDoubleAnchor_Splice_Info();
				bool splice_fixed = fixSpliceInfo_unannotated->detectBestSpliceSite//_prefer_canonical_lessMismatch
					(
					segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
					segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
					readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch, 
					false, false, annotationInfo);

					if(splice_fixed && 
						((fixSpliceInfo_unannotated->returnBestSplice_canonicalOrNot()) 
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

						fixGapVec[index_fixGapVec].first = true;					
						// fixed with unannotated SJ

						(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
						((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
						((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
						((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

						this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
						this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
					}
					else
					{
						/////  test fixing complicated SJ  ///////
						//cout << "\nstart to try fixing complicated SJ ....\n" << endl;
						FixDoubleAnchor_Splice_Complicate_Info* fixComplicateSpliceInfo_unannotated
							= new FixDoubleAnchor_Splice_Complicate_Info();
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

							fixGapVec[index_fixGapVec].first = true;
							// fixed with unannotated SJ

							(fixGapVec[index_fixGapVec].second).first = fixComplicateSpliceInfo_unannotated->return_mismatch_bestComplicatedSplice();
							((fixGapVec[index_fixGapVec].second).second).push_back(prefixMatchJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(firstJumpCode);
							((fixGapVec[index_fixGapVec].second).second).push_back(midMatchJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(secondJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(suffixMatchJumpCode);	

							this->insertMismatchPosVec_fixComplicatedSpliceInfo(fixComplicateSpliceInfo_unannotated, index_fixGapVec);
							this->insertMismatchCharVec_fixComplicateSpliceInfo(fixComplicateSpliceInfo_unannotated, index_fixGapVec);	
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

							fixGapVec[index_fixGapVec].first = true;
							// fixed with unannotated SJ	

							(fixGapVec[index_fixGapVec].second).first = fixSpliceInfo_unannotated->returnBestSplice_mismatchNum();
							((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
							((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
							((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

							//fixSpliceInfo_unannotated->generateBestSpliceMismatchVec(subSeqLengthInProcess, (segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left));
							this->insertMismatchPosVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
							this->insertMismatchCharVec_fixSpliceInfo(fixSpliceInfo_unannotated, index_fixGapVec);
						}
						else
						{
							fixGapVec[index_fixGapVec].first = false;
							// fixed with annotated SJ
							//annotated_splice_fixed_bool = false;	
						}
						delete fixComplicateSpliceInfo_unannotated;
					}

				delete fixSpliceInfo_unannotated;
			}
		}
		return;
	}

	void fixDoubleAnchor_Lariat_RNA(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec)
	{

	}

	void fixDoubleAnchor_Circular_RNA(int relation, int segmentLocInRead_1, int segmentLocInRead_2,
		int segmentLength_1, int segmentLength_2, int segmentMapPos_1, int segmentMapPos_2,
		const string& readSeq_inProcess, //Index_Info* indexInfo, 
		Index_Info* indexInfo, int chrNameInt, int index_fixGapVec
		)
	{
		int tmpBuffer_left = FixSpliceBuffer;//+1;
		if(tmpBuffer_left > segmentLength_1 - 2) //anchor >= 2
		{
			tmpBuffer_left = segmentLength_1 - 2;
		}
		int tmpBuffer_right = FixSpliceBuffer;//+1;
		if(tmpBuffer_right > segmentLength_2 - 2)
		{
			tmpBuffer_right = segmentLength_2 - 2;
		}

		int subSeqLengthInProcess = segmentLocInRead_2 - 1 - (segmentLocInRead_1 + segmentLength_1) + 1 + tmpBuffer_left + tmpBuffer_right;
		int cirRNAspliceJunctionLength = (segmentMapPos_2 - segmentLocInRead_2) - (segmentMapPos_1 - segmentLocInRead_1);

		//size_t prefix_length = 0;
		size_t max_double_splice_mismatch = (subSeqLengthInProcess)/LengthOfSeqPerMismatchAllowed + 1;
		FixDoubleAnchor_CirRNA_Info* fixCirRnaInfo = new FixDoubleAnchor_CirRNA_Info();
		bool cirRNA_fixed_bool = fixCirRnaInfo->detectBestCirRNAspliceSite(
			segmentLocInRead_1 + segmentLength_1 - tmpBuffer_left, segmentLocInRead_2 + tmpBuffer_right - 1,
			segmentMapPos_1 + segmentLength_1 - tmpBuffer_left, segmentMapPos_2 + tmpBuffer_right - 1,
			readSeq_inProcess, indexInfo, chrNameInt, max_double_splice_mismatch);
		if(cirRNA_fixed_bool)
		{
			int prefix_length = fixCirRnaInfo->returnBestCirRNAspliceSite_prefixMatchLength();
			int firstMatchLength = - tmpBuffer_left + prefix_length;
			int secondMatchLength = segmentLength_2 - tmpBuffer_right + subSeqLengthInProcess - prefix_length;

			Jump_Code firstMatchJumpCode(firstMatchLength, "M");
			Jump_Code spliceJumpCode(cirRNAspliceJunctionLength, "N");
			Jump_Code secondMatchJumpCode(secondMatchLength, "M"); 

			fixGapVec[index_fixGapVec].first = true;
			(fixGapVec[index_fixGapVec].second).first = fixCirRnaInfo->returnBestCirRNAspliceSite_mismatchNum();
			((fixGapVec[index_fixGapVec].second).second).push_back(firstMatchJumpCode);
			((fixGapVec[index_fixGapVec].second).second).push_back(spliceJumpCode);	
			((fixGapVec[index_fixGapVec].second).second).push_back(secondMatchJumpCode);

			this->insertMismatchPosVec_fixCirRnaInfo(fixCirRnaInfo, index_fixGapVec);
			this->insertMismatchCharVec_fixCirRnaInfo(fixCirRnaInfo, index_fixGapVec);
		}
		else
		{
			fixGapVec[index_fixGapVec].first = false;
		}
		delete fixCirRnaInfo;
		return;// fixDoubleAnchorBool_Deletion;
	}

	int pathValidNumInt()
	{
		int pathValidNum = 0;
		for(int tmp = 0; tmp < PathValidBoolVec.size(); tmp++)
		{
			if(PathValidBoolVec[tmp])
				pathValidNum++;
		}
		return pathValidNum;
	}

	int checkTwoStringMatchOrNot(const string& string_1, const string& string_2, int maxMismatch)
	{
		/*#ifdef CAL_TIME
		checkTwoStringMatch_begin = clock();
		#endif*/
		bool matchBool = true;
		if(string_1.length() != string_2.length())
		{
			cout << "string_1 != string_2 in terms of length" << endl;
			matchBool = false;
		}
		int mismatchNum = 0;
		for(int tmp = 0; tmp < string_1.length(); tmp++)
		{
			if(string_1[tmp] != string_2[tmp])
			{	mismatchNum++;
				if(mismatchNum > maxMismatch)
				{
					matchBool = false;; // not match
					break;
				}
			}
		}
		/*#ifdef CAL_TIME
		checkTwoStringMatch_end = clock();
		checkTwoStringMatch_cost = checkTwoStringMatch_cost + checkTwoStringMatch_end - checkTwoStringMatch_begin;
		#endif*/
		if(matchBool)
			return mismatchNum;
		else
			return -1;
	}


	int checkTwoStringMatchOrNot(const string& string_1, const string& string_2, int maxMismatch, 
		vector<int>& mismatchPosVec, vector<char>& mismatchCharVec, int startLocInRead)
	{
		/*#ifdef CAL_TIME
		checkTwoStringMatch_begin = clock();
		#endif*/
		bool matchBool = true;
		if(string_1.length() != string_2.length())
		{
			cout << "string_1 != string_2 in terms of length" << endl;
			matchBool = false;
		}
		int mismatchNum = 0;
		for(int tmp = 0; tmp < string_1.length(); tmp++)
		{
			if(string_1[tmp] != string_2[tmp])
			{	
				mismatchNum++;
				if(STORE_MISMATCH_POS)
				{
					mismatchPosVec.push_back(startLocInRead + tmp);
					if(STORE_MISMATCH_CHA)
					{
						mismatchCharVec.push_back(string_2.at(tmp));
					}
				}
				if(mismatchNum > maxMismatch)
				{
					matchBool = false;; // not match
					break;
				}
			}
		}
		/*#ifdef CAL_TIME
		checkTwoStringMatch_end = clock();
		checkTwoStringMatch_cost = checkTwoStringMatch_cost + checkTwoStringMatch_end - checkTwoStringMatch_begin;
		#endif*/
		if(matchBool)
			return mismatchNum;
		else
			return -1;
	}

	int extendBack2Head(const string& string_1, const string& string_2)
	{
		//int extendBackLength = 0; 
		int stringLength = string_1.length();
		for(int tmp = 0; tmp < stringLength; tmp++)
		{
			if(string_1[stringLength - 1 - tmp] != string_2[stringLength - 1 - tmp])
			{	
				return tmp;
			}
		}
		return string_1.length();		
	}

	string finalFixedPathStr(Index_Info* indexInfo)
	{
		string finalFixedPathStr = "final Fixed Path Info: \n";
		for(int tmpPath = 0; tmpPath < finalPathVec.size(); tmpPath++)
		{

			int tmpChrNameInt = (finalPathVec[tmpPath].first).first;
			string tmpChrNameStr = indexInfo->returnChrNameStr(tmpChrNameInt);
			int tmpSegmentMapPos = (finalPathVec[tmpPath].first).second;		
			
			finalFixedPathStr = finalFixedPathStr + int_to_str(tmpPath+1) + " ... fixed Path  " + ": "		
				+ tmpChrNameStr + " " + int_to_str(tmpSegmentMapPos) + " " + (finalPathVec[tmpPath].second)->printFinalJumpCode() 
				+ " mismatch#: " + int_to_str(fixedPathMismatchVec[tmpPath]) + " mismatch_pos_base: ";
			if(STORE_MISMATCH_POS)
			{	
				for(int tmp = 0; tmp < fixedPathMismatchPosVec[tmpPath].size(); tmp++)
				{
					finalFixedPathStr = finalFixedPathStr + int_to_str((fixedPathMismatchPosVec[tmpPath])[tmp]) + ",";
				}
			}
			finalFixedPathStr += "\n";
		}

		return finalFixedPathStr;
	}

	string fixedPathVecStr(Index_Info* indexInfo, Seg_Info* segInfo)
	{
		string fixedPathStr = "fixed Path Info: \n";

		for(int tmpPath = 0; tmpPath < fixedPathVec.size(); tmpPath++)
		{

			int fixedPathNO = fixedPathVec[tmpPath].first;

			int tmpPath1stSegGroupNO = (PathVec_seg[fixedPathNO])[0].first;
			int tmpPath1stSegCandiNO = (PathVec_seg[fixedPathNO])[0].second;
			unsigned int PathMapPos //= *(segInfo->norSegmentAlignLoc + tmpPath1stSegGroupNO * CANDALILOC + tmpPath1stSegCandiNO);
				= segInfo->returnSegmentMapPos(tmpPath1stSegGroupNO, tmpPath1stSegCandiNO);
			unsigned int tmpChrNameInt, tmpChrPosInt;
			indexInfo->getChrLocation(PathMapPos, &tmpChrNameInt, &tmpChrPosInt);
			string tmpChrNameStr = indexInfo->returnChrNameStr(tmpChrNameInt);
			int tmpSegmentMapPos = tmpChrPosInt;			
			
			fixedPathStr = fixedPathStr + int_to_str(tmpPath+1) + " ... fixed Path " + int_to_str(fixedPathNO+1) + ": "			
				+ tmpChrNameStr + " " + int_to_str(tmpSegmentMapPos) + " " + (fixedPathVec[tmpPath].second)->printJumpCode();
			fixedPathStr += "\n";
		}

		return fixedPathStr;
	}
	
	void addNewSegGroupToCurrentPathInfoVec(Seg_Info* segInfo, int segGroupNO, int spliceJunctionDistanceMax)
	{
		bool longSegBool = //true; // 0806 
			segInfo->checkSegLongOrNot(segGroupNO);
		//if(!longSegBool) // 0906
		//	return;	
		//cout << "start to add segGroupNO: " << segGroupNO << endl;
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(segGroupNO); tmpSegCandiLoc++)
		{
			//this->addNewSegCandiToCurrentPathInfoVec_uniquePath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//this->addNewSegCandiToCurrentPathInfoVec_allPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//cout << "start to add segCandiNO: " << tmpSegCandiLoc << endl;
			//this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			this->addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction(
				segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum, spliceJunctionDistanceMax);
		}	
	}

	
	void addNewSegGroupToCurrentPathInfoVec_fixGapAlso(
		Seg_Info* segInfo, int segGroupNO, Index_Info* indexInfo, const string& readSeq_inProcess, 
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		bool longSegBool = segInfo->checkSegLongOrNot(segGroupNO);//0906 true;// 0806 segInfo->checkSegLongOrNot(segGroupNO);
		//cout << "start to add segGroupNO: " << segGroupNO << endl;
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(segGroupNO); tmpSegCandiLoc++)
		{
			//this->addNewSegCandiToCurrentPathInfoVec_uniquePath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//this->addNewSegCandiToCurrentPathInfoVec_allPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//cout << "start to add segCandiNO: " << tmpSegCandiLoc << endl;
			//this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			this->addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction_fixGapAlso(
				segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum, indexInfo, readSeq_inProcess,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
		}	
	}

	void addNewSegGroupToCurrentPathInfoVec_localGreedyMapping(Seg_Info* segInfo, int segGroupNO,
		int spliceJunctionDistanceMax)
	{
		bool longSegBool = //true; // 0806 
			segInfo->checkSegLongOrNot(segGroupNO);
		if(!longSegBool)  // fix me: never introduce short seg to get a path
			return;     
		//cout << "start to add segGroupNO: " << segGroupNO << endl;
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(segGroupNO); tmpSegCandiLoc++)
		{
			//this->addNewSegCandiToCurrentPathInfoVec_uniquePath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//this->addNewSegCandiToCurrentPathInfoVec_allPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			//cout << "start to add segCandiNO: " << tmpSegCandiLoc << endl;
			//this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
			this->addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction(
				segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum, spliceJunctionDistanceMax);
		}	
	}

	void addNewSegGroupToCurrentPathInfoVec_cirRNA(Seg_Info* segInfo, int segGroupNO)
	{
		bool longSegBool = segInfo->checkSegLongOrNot(segGroupNO);
		//cout << "start to add segGroupNO: " << segGroupNO << endl;
		int currentPathNum = PathVec_seg.size();
		for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(segGroupNO); tmpSegCandiLoc++)
		{
			this->addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath_cirRNA(segInfo, segGroupNO, tmpSegCandiLoc, longSegBool, currentPathNum);
		}	
	}

	void copyOldPath_AddSegCandi_AddNewPath(int pathElementNO, int segGroupNO, int segCandiNO)
	{
		//copy path;
		vector< pair<int,int> > newPath;
		for(int tmpElementNO = 0; tmpElementNO < PathVec_seg[pathElementNO].size(); tmpElementNO++)
		{
			newPath.push_back((PathVec_seg[pathElementNO])[tmpElementNO]);
		}
		newPath.push_back(pair<int, int> (segGroupNO, segCandiNO));
		PathVec_seg.push_back(newPath);
		PathValidBoolVec.push_back(true);
		vector< pair<int,int> >().swap(newPath);
	}

	
	void copyOldPath_AddSegCandi_AddNewPath_fixGapAlso(int pathElementNO, int segGroupNO, int segCandiNO, int gapIndexNO)
	{
		//copy path;
		vector< pair<int,int> > newPath;
		for(int tmpElementNO = 0; tmpElementNO < PathVec_seg[pathElementNO].size(); tmpElementNO++)
		{
			newPath.push_back((PathVec_seg[pathElementNO])[tmpElementNO]);
		}
		newPath.push_back(pair<int, int> (segGroupNO, segCandiNO));
		PathVec_seg.push_back(newPath);

		vector<int> newPath_gapIndex;
		for(int tmpGapNO = 0; tmpGapNO < PathVec_seg_gapIndex[pathElementNO].size(); tmpGapNO++)
		{
			newPath_gapIndex.push_back((PathVec_seg_gapIndex[pathElementNO])[tmpGapNO]);
		}
		newPath_gapIndex.push_back(gapIndexNO);
		PathVec_seg_gapIndex.push_back(newPath_gapIndex);

		PathValidBoolVec.push_back(true);
		vector< pair<int,int> >().swap(newPath);
		vector<int>().swap(newPath_gapIndex);
	}

	int chooseBestCurrentPathForNewSegCandi(Seg_Info* segInfo, int segGroupNO, int segCandiNO, int currentPathNum, 
		int spliceJunctionDistanceMax)
	{
		//int bestCurrentPath = 1000000;
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO, spliceJunctionDistanceMax);
			if(tmpSegDistance <= minDistance_value)
			{
				minDistance_value = tmpSegDistance;
				minDistance_path = tmpPathElementNO;
			}	
		}
		return minDistance_path;
	}
	
	void generateBestCurrentPathVecForNewSegCandi(Seg_Info* segInfo, int segGroupNO, 
		int segCandiNO, int currentPathNum, vector<int>& targetPathVec, int spliceJunctionDistanceMax)
	{
		vector<int> bestPathVec_match;
		vector<int> bestPathVec_indel;
		vector<int> bestPathVec_shortSpliceJunction;
		vector<int> bestPathVec_longSpliceJunction;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO, spliceJunctionDistanceMax);
			if(tmpSegDistance == 0)
			{
				bestPathVec_match.push_back(tmpPathElementNO);
			}	
			else if (tmpSegDistance <= MAX_DELETION_LENGTH)
			{
				bestPathVec_indel.push_back(tmpPathElementNO);
			}
			else if (tmpSegDistance <= MAX_SHORT_SPLICE_DISTANCE)
			{
				bestPathVec_shortSpliceJunction.push_back(tmpPathElementNO);
			}
			else if(tmpSegDistance <= spliceJunctionDistanceMax)
			{
				bestPathVec_longSpliceJunction.push_back(tmpPathElementNO);
			}
			else
			{}
		}

		if(bestPathVec_match.size() > 0)
		{
			for(int tmp = 0; tmp < bestPathVec_match.size(); tmp++)
			{
				targetPathVec.push_back(bestPathVec_match[tmp]);
			}
		}
		else if(bestPathVec_indel.size() > 0)
		{
			for(int tmp = 0; tmp < bestPathVec_indel.size(); tmp++)
			{
				targetPathVec.push_back(bestPathVec_indel[tmp]);
			}
		}
		else if(bestPathVec_shortSpliceJunction.size() > 0)
		{
			for(int tmp = 0; tmp < bestPathVec_shortSpliceJunction.size(); tmp++)
			{
				targetPathVec.push_back(bestPathVec_shortSpliceJunction[tmp]);
			}
		}
		else if(bestPathVec_longSpliceJunction.size() > 0)
		{
			for(int tmp = 0; tmp < bestPathVec_longSpliceJunction.size(); tmp++)
			{
				targetPathVec.push_back(bestPathVec_longSpliceJunction[tmp]);
			}			
		}
		else
		{}

		//return minDistance_path;
		return;
	}

	void generateBestCurrentPathVecForNewSegCandi_fixGapAlso(Seg_Info* segInfo, int segGroupNO, 
		int segCandiNO, int currentPathNum, vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
		const string& readSeq_inProcess, Index_Info* indexInfo,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		//cout << " start to generateBestCurrentPathVecForNewSegCandi_fixGapAlso " << endl;
		vector< pair<int,int> > possiblePathVec_match;
		vector< pair<int,int> > possiblePathVec_insertion;
		vector< pair<int,int> > possiblePathVec_deletion;
		vector< pair<int,int> > possiblePathVec_spliceJunction;
		vector< pair<int,int> > possiblePathVec_cirRNA;

		vector<int> possiblePathVec_match_pathIndex;
		vector<int> possiblePathVec_insertion_pathIndex;
		vector<int> possiblePathVec_deletion_pathIndex;		
		vector<int> possiblePathVec_spliceJunction_pathIndex;
		vector<int> possiblePathVec_cirRNA_pathIndex;

		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			//cout << "tmpPathLastSegGroupNO: " << tmpPathLastSegGroupNO << endl;
			//cout << "tmpPathLastSegCandiNO: " << tmpPathLastSegCandiNO << endl;
			int tmpSegDistance = segInfo->distanceBetweenSegment_signed(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO, spliceJunctionDistanceMax);
			//cout << "tmpSegDistance: " << tmpSegDistance << endl;
			if(tmpSegDistance == 0)
			{
				possiblePathVec_match.push_back( pair<int, int>(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO) );
				possiblePathVec_match_pathIndex.push_back(tmpPathElementNO);
			}
			else if(tmpSegDistance < 0)
			{
				if((0 - tmpSegDistance) <= MAX_INSERTION_LENGTH)
				{
					possiblePathVec_insertion.push_back( pair<int, int>(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO) );
					possiblePathVec_insertion_pathIndex.push_back(tmpPathElementNO);
				}
				else
				{
					#ifdef DETECT_CIRCULAR_RNA
					if((0 - tmpSegDistance) <= spliceJunctionDistanceMax)
					{
						//cout << "pushBack2candidateCirRNAgap" << endl;
						possiblePathVec_cirRNA.push_back( pair<int, int>(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO) );
						possiblePathVec_cirRNA_pathIndex.push_back(tmpPathElementNO);
					}
					#endif
				}
			}	
			else if (tmpSegDistance <= MAX_DELETION_LENGTH)
			{
				possiblePathVec_deletion.push_back( pair<int, int>(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO) );
				possiblePathVec_deletion_pathIndex.push_back(tmpPathElementNO);
			}
			else if(tmpSegDistance <= spliceJunctionDistanceMax)
			{
				possiblePathVec_spliceJunction.push_back( pair<int, int>(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO) );
				possiblePathVec_spliceJunction_pathIndex.push_back(tmpPathElementNO);
			}
			else
			{}
		}
		
		//cout << "start to check match_gap_fixed_bool or not ..." << endl;	
		bool match_gap_fixed_bool = this->generateBestTargetPathVec_fixGapAlso_match(
			segInfo, segGroupNO, segCandiNO, targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_match, possiblePathVec_match_pathIndex, indexInfo, readSeq_inProcess,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

		//cout << "end of checking match_gap_fixed_bool: " << match_gap_fixed_bool << endl;
		if(match_gap_fixed_bool)
		{
			return;
		}
			
		bool indel_gap_fixed_bool = this->generateBestTargetPathVec_fixGapAlso_indel(
			segInfo, segGroupNO, segCandiNO, targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_insertion, possiblePathVec_insertion_pathIndex,
			possiblePathVec_deletion, possiblePathVec_deletion_pathIndex, indexInfo, readSeq_inProcess,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
				
		if(indel_gap_fixed_bool)
		{
			return;
		}
			
		bool spliceJunction_gap_fixed_bool = this->generateBestTargetPathVec_fixGapAlso_spliceJunction(
			segInfo, segGroupNO, segCandiNO, targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_spliceJunction, possiblePathVec_spliceJunction_pathIndex, indexInfo, readSeq_inProcess,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
		if(spliceJunction_gap_fixed_bool)
		{
			return;
		}

		#ifdef DETECT_CIRCULAR_RNA
		bool cirRNA_gap_fixed_bool = this->generateBestTargetPathVec_fixGapAlso_cirRNA(
			segInfo, segGroupNO, segCandiNO, targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_cirRNA, possiblePathVec_cirRNA_pathIndex, indexInfo, readSeq_inProcess,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
		if(cirRNA_gap_fixed_bool)
			return;
		#endif

		return;
	}

	bool generateBestTargetPathVec_fixGapAlso_match(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_match, 
			vector<int>& possiblePathVec_match_pathIndex, Index_Info* indexInfo, const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		bool match_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		vector<int> tmpFixGapVec_match_index;
		for(int tmp = 0; tmp < possiblePathVec_match.size(); tmp++)
		{
			//cout << "try to fix: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			//cout << "try to fix: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;

			int tmpPathLastSegGroupNO = possiblePathVec_match[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_match[tmp].second;

			//cout << "try to fix match: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			//cout << "gapIndex: " << currentFixGapVec_index << endl;

			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				//cout << "error in generateBestTargetPathVec_fixGapAlso fixMatch " 
 				//	<< tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
 				//return false;
 				//cout << "found in map ! fixed or unfixed: " << fixGapVec[tmpFoundIndexInFixGapMap].first
 				//	<< " found index: " << tmpFoundIndexInFixGapMap << endl;
 				tmpFixGapVec_match_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//int currentFixGapVec_index = fixGapVec.size();

 			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess,
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);

 			tmpFixGapVec_match_index.push_back(currentFixGapVec_index);
 			//cout << "not found in map ! after trying to fix, fixed or not: " << fixGapVec[currentFixGapVec_index].first << endl;
			currentFixGapVec_index ++;
		}

		for(int tmp = 0; tmp < tmpFixGapVec_match_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_match_index[tmp];
			int tmpPathIndex = possiblePathVec_match_pathIndex[tmp];
			bool tmpMatchGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpMatchGapFixed_bool)
			{
				match_gap_fixed_bool = true;
				targetPathVec.push_back(tmpPathIndex);
				targetPathVec_gapIndex.push_back(tmpGapIndex);
			}
		}
		return match_gap_fixed_bool;
	}

	bool generateBestTargetPathVec_fixGapAlso_indel(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_insertion, 
			vector<int>& possiblePathVec_insertion_pathIndex,
			vector< pair<int,int> >& possiblePathVec_deletion, 
			vector<int>& possiblePathVec_deletion_pathIndex, Index_Info* indexInfo, 
			const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		bool indel_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		
		vector<int> tmpFixGapVec_insertion_index;
		for(int tmp = 0; tmp < possiblePathVec_insertion.size(); tmp++)
		{
			//cout << "try to fix: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;

			int tmpPathLastSegGroupNO = possiblePathVec_insertion[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_insertion[tmp].second;

			//cout << "try to fix insertion: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(
 					tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
 					segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				//cout << "error in generateBestTargetPathVec_fixGapAlso fix insertion " 
 				//	<< tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
 				tmpFixGapVec_insertion_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//int currentFixGapVec_index = fixGapVec.size();
 			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess,
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
 				spliceJunctionDistanceMax);
 			tmpFixGapVec_insertion_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}
		vector<int> tmpFixGapVec_deletion_index;
		for(int tmp = 0; tmp < possiblePathVec_deletion.size(); tmp++)
		{
			//cout << "try to fix: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;

			int tmpPathLastSegGroupNO = possiblePathVec_deletion[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_deletion[tmp].second;

			//cout << "try to fix deletion: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(
 					tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
 					segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				//cout << "error in generateBestTargetPathVec_fixGapAlso fix deletion " 
 				//	<< tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
 				//return false;
 				tmpFixGapVec_deletion_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//int currentFixGapVec_index = fixGapVec.size();
 			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess,
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
 				spliceJunctionDistanceMax);
 			tmpFixGapVec_deletion_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}

		for(int tmp = 0; tmp < tmpFixGapVec_insertion_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_insertion_index[tmp];
			int tmpPathIndex = possiblePathVec_insertion_pathIndex[tmp];
			bool tmpInsertionGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpInsertionGapFixed_bool)
			{
				indel_gap_fixed_bool = true;
				targetPathVec.push_back(tmpPathIndex);
				targetPathVec_gapIndex.push_back(tmpGapIndex);
			}
		}
		for(int tmp = 0; tmp < tmpFixGapVec_deletion_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_deletion_index[tmp];
			int tmpPathIndex = possiblePathVec_deletion_pathIndex[tmp];
			bool tmpDeletionGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpDeletionGapFixed_bool)
			{
				indel_gap_fixed_bool = true;
				targetPathVec.push_back(tmpPathIndex);
				targetPathVec_gapIndex.push_back(tmpGapIndex);
			}
		}
		return indel_gap_fixed_bool;
	}

	bool generateBestTargetPathVec_fixGapAlso_spliceJunction(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_spliceJunction, 
			vector<int>& possiblePathVec_spliceJunction_pathIndex, 
			Index_Info* indexInfo, const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		bool spliceJunction_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		vector<int> tmpFixGapVec_spliceJunction_index;
		for(int tmp = 0; tmp < possiblePathVec_spliceJunction.size(); tmp++)
		{
			int tmpPathLastSegGroupNO = possiblePathVec_spliceJunction[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_spliceJunction[tmp].second;
 
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				tmpFixGapVec_spliceJunction_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				

			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess, 
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
 			tmpFixGapVec_spliceJunction_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}

		spliceJunction_gap_fixed_bool = this->selectBestSpliceJunctionCase(
			segInfo,
			targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_spliceJunction, possiblePathVec_spliceJunction_pathIndex, 
			indexInfo, tmpFixGapVec_spliceJunction_index
			);

		return spliceJunction_gap_fixed_bool;
	}

	bool generateBestTargetPathVec_fixGapAlso_cirRNA(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_cirRNA, 
			vector<int>& possiblePathVec_cirRNA_pathIndex, 
			Index_Info* indexInfo, const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		//cout << "start to generateBestTargetPathVec_fixGapAlso_cirRNA( ..." << endl;
		//cout << "possiblePathVec_cirRNA.size(): " << possiblePathVec_cirRNA.size() << endl;
		bool cirRNA_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		vector<int> tmpFixGapVec_cirRNA_index;
		for(int tmp = 0; tmp < possiblePathVec_cirRNA.size(); tmp++)
		{
			int tmpPathLastSegGroupNO = possiblePathVec_cirRNA[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_cirRNA[tmp].second;
 			//cout << "tmpPathLastSegGroupNO: " << tmpPathLastSegGroupNO << endl;
 			//cout << "tmpPathLastSegCandiNO: " << tmpPathLastSegCandiNO << endl;
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				tmpFixGapVec_cirRNA_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//cout << "start to pushBackAndFixGap ...." << endl;
			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess, 
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
 			tmpFixGapVec_cirRNA_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}
		//cout << "tmpFixGapVec_cirRNA_index.size(): " << tmpFixGapVec_cirRNA_index.size() << endl;
		cirRNA_gap_fixed_bool = this->selectBestCirRnaCase(
			segInfo,
			targetPathVec, targetPathVec_gapIndex,
			possiblePathVec_cirRNA, possiblePathVec_cirRNA_pathIndex, 
			indexInfo, tmpFixGapVec_cirRNA_index
			);

		return cirRNA_gap_fixed_bool;
	}

	bool selectBestSpliceJunctionCase(
			Seg_Info* segInfo,
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_spliceJunction, 
			vector<int>& possiblePathVec_spliceJunction_pathIndex, 
			Index_Info* indexInfo, vector<int>& tmpFixGapVec_spliceJunction_index
			)
	{
		bool selectBestSpliceJunctionCase_bool = false;

		vector<int> validPathVec;
		vector<int> validPathVec_gapIndex;
		//vector< pair<int, int> > validPathVec_lastSeg;
		vector<int> validPathVec_flankStringCase;
		vector<int> validPathVec_spliceJunctionDistance;

		for(int tmp = 0; tmp < tmpFixGapVec_spliceJunction_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_spliceJunction_index[tmp];
			int tmpPathIndex = possiblePathVec_spliceJunction_pathIndex[tmp];
			bool tmpSpliceJunctionGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpSpliceJunctionGapFixed_bool)
			{
				selectBestSpliceJunctionCase_bool = true;
				validPathVec.push_back(tmpPathIndex);
				validPathVec_gapIndex.push_back(tmpGapIndex);
				//validPathVec_lastSeg.push_back(possiblePathVec_spliceJunction[tmp]);
				int tmpSpliceJunctionDistance;
				int tmpFlankStringCase = this->getTmpFlankStringCaseFromJumpCodeVec(
					segInfo, possiblePathVec_spliceJunction[tmp].first, possiblePathVec_spliceJunction[tmp].second,
					//tmpSpliceJunctionJumpCodeVec[tmp],
					(fixGapVec[tmpGapIndex].second).second,
					indexInfo, tmpSpliceJunctionDistance
					);
				validPathVec_flankStringCase.push_back(tmpFlankStringCase);
				validPathVec_spliceJunctionDistance.push_back(tmpSpliceJunctionDistance);
			}
		}

		if(selectBestSpliceJunctionCase_bool)
		{
			/*
			this->getFinalTargetPathAndGap_unique_nearestCanonicalSJ(
				validPathVec, validPathVec_gapIndex,
				validPathVec_flankStringCase, validPathVec_spliceJunctionDistance,
				targetPathVec, targetPathVec_gapIndex);*/
			this->getFinalTargetPathAndGap_multi_shortSJ(
				validPathVec, validPathVec_gapIndex,
				validPathVec_flankStringCase, validPathVec_spliceJunctionDistance,
				targetPathVec, targetPathVec_gapIndex);
		}

		return selectBestSpliceJunctionCase_bool;
	}

	bool selectBestCirRnaCase(
			Seg_Info* segInfo,
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_cirRNA, 
			vector<int>& possiblePathVec_cirRNA_pathIndex, 
			Index_Info* indexInfo, vector<int>& tmpFixGapVec_cirRNA_index
			)
	{
		bool selectBestCirRnaCase_bool = false;

		vector<int> validPathVec;
		vector<int> validPathVec_gapIndex;
		//vector< pair<int, int> > validPathVec_lastSeg;
		vector<int> validPathVec_flankStringCase;
		vector<int> validPathVec_cirRnaDistance;

		for(int tmp = 0; tmp < tmpFixGapVec_cirRNA_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_cirRNA_index[tmp];
			int tmpPathIndex = possiblePathVec_cirRNA_pathIndex[tmp];
			bool tmpCirRnaGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpCirRnaGapFixed_bool)
			{
				selectBestCirRnaCase_bool = true;
				validPathVec.push_back(tmpPathIndex);
				validPathVec_gapIndex.push_back(tmpGapIndex);
				//validPathVec_lastSeg.push_back(possiblePathVec_spliceJunction[tmp]);
				int tmpCirRnaDistance;
				int tmpFlankStringCase = this->getTmpFlankStringCaseFromJumpCodeVec(
					segInfo, possiblePathVec_cirRNA[tmp].first, possiblePathVec_cirRNA[tmp].second,
					//tmpSpliceJunctionJumpCodeVec[tmp],
					(fixGapVec[tmpGapIndex].second).second,
					indexInfo, tmpCirRnaDistance
					);
				validPathVec_flankStringCase.push_back(tmpFlankStringCase);
				validPathVec_cirRnaDistance.push_back(tmpCirRnaDistance);
			}
		}

		if(selectBestCirRnaCase_bool)
		{
			/*
			this->getFinalTargetPathAndGap_unique_nearestCanonicalSJ(
				validPathVec, validPathVec_gapIndex,
				validPathVec_flankStringCase, validPathVec_spliceJunctionDistance,
				targetPathVec, targetPathVec_gapIndex);*/
			this->getFinalTargetPathAndGap_multi_shortCirRNA(
				validPathVec, validPathVec_gapIndex,
				validPathVec_flankStringCase, validPathVec_cirRnaDistance,
				targetPathVec, targetPathVec_gapIndex);
		}

		return selectBestCirRnaCase_bool;
	}

	void getFinalTargetPathAndGap_unique_nearestCanonicalSJ(
		vector<int>& validPathVec, vector<int>& validPathVec_gapIndex,
		vector<int>& validPathVec_flankStringCase, 
		vector<int>& validPathVec_spliceJunctionDistance,
		vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex)
	{
		int tmp_SJ_distance_min = 600000;
		int tmp_SJ_distance_min_canonical = 600000;

		int tmp_SJ_distance_min_indice = -1;
		int tmp_SJ_distance_min_canonical_indice = -1;	
		
			for(int tmp = 0; tmp < validPathVec_flankStringCase.size(); tmp++)
			{
				int tmpSJdistance = validPathVec_spliceJunctionDistance[tmp];
				if(tmpSJdistance < tmp_SJ_distance_min)
				{
					tmp_SJ_distance_min_indice = tmp;
					tmp_SJ_distance_min = tmpSJdistance;
				}

				if(validPathVec_flankStringCase[tmp] >= 5)
				{
					if(tmpSJdistance < tmp_SJ_distance_min_canonical)
					{
						tmp_SJ_distance_min_canonical_indice = tmp;
						tmp_SJ_distance_min_canonical = tmpSJdistance;
					}
				}
			}

			if(tmp_SJ_distance_min_canonical_indice >= 0)
			{
				targetPathVec.push_back(validPathVec[tmp_SJ_distance_min_canonical_indice]);
				targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp_SJ_distance_min_canonical_indice]);
			}
			else if(tmp_SJ_distance_min_indice >= 0)
			{
				targetPathVec.push_back(validPathVec[tmp_SJ_distance_min_indice]);
				targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp_SJ_distance_min_indice]);
			}
			else
			{
				//selectBestSpliceJunctionCase_bool = false;
				cout << "error in selectBestSpliceJunctionCase path_info.h " << endl;
			}		
	}

	void getFinalTargetPathAndGap_multi_shortSJ(
		vector<int>& validPathVec, vector<int>& validPathVec_gapIndex,
		vector<int>& validPathVec_flankStringCase, 
		vector<int>& validPathVec_spliceJunctionDistance,
		vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex)
	{
		for(int tmp = 0; tmp < validPathVec_flankStringCase.size(); tmp++)
		{
			int tmpSJdistance = validPathVec_spliceJunctionDistance[tmp];
			if(tmpSJdistance <= MAX_SHORT_SPLICE_DISTANCE)
			{
				targetPathVec.push_back(validPathVec[tmp]);
				targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp]);
			}
		}
	
		if(targetPathVec.size() > 0)
			return;
	
		for(int tmp = 0; tmp < validPathVec_flankStringCase.size(); tmp++)
		{
			targetPathVec.push_back(validPathVec[tmp]);
			targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp]);		
		}		
		return;
	}

	void getFinalTargetPathAndGap_multi_shortCirRNA(
		vector<int>& validPathVec, vector<int>& validPathVec_gapIndex,
		vector<int>& validPathVec_flankStringCase, 
		vector<int>& validPathVec_CirRnaDistance,
		vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex)
	{
		for(int tmp = 0; tmp < validPathVec_flankStringCase.size(); tmp++)
		{
			int tmpSJdistance = 0 - validPathVec_CirRnaDistance[tmp];
			if(tmpSJdistance <= MAX_SHORT_SPLICE_DISTANCE)
			{
				targetPathVec.push_back(validPathVec[tmp]);
				targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp]);
			}
		}
	
		if(targetPathVec.size() > 0)
			return;
	
		for(int tmp = 0; tmp < validPathVec_flankStringCase.size(); tmp++)
		{
			targetPathVec.push_back(validPathVec[tmp]);
			targetPathVec_gapIndex.push_back(validPathVec_gapIndex[tmp]);		
		}		
		return;
	}

	int getTmpFlankStringCaseFromJumpCodeVec(
		Seg_Info* segInfo, int lastSegGroupNO, int lastSegCandiNO,
		vector<Jump_Code>& tmpSpliceJunctionJumpCodeVec, Index_Info* indexInfo, int& tmpSpliceJunctionDistance
		)
	{
		int firstMatchLength = tmpSpliceJunctionJumpCodeVec[0].len;
		int spliceJunctionDistance = tmpSpliceJunctionJumpCodeVec[1].len;
		
		tmpSpliceJunctionDistance = spliceJunctionDistance;

		unsigned int donerSitePos = segInfo->returnSegmentMapPos_end(lastSegGroupNO, lastSegCandiNO)
			+ firstMatchLength;
		unsigned int acceptorSitePos = donerSitePos + spliceJunctionDistance + 1;

		//string flankString = (indexInfo->chromString).substr(donerSitePos, 2) 
		//	+ (indexInfo->chromString).substr(acceptorSitePos - 3, 2);
		string flankString = indexInfo->returnChromStringSubstr(donerSitePos + 1, 2)
			+ indexInfo->returnChromStringSubstr(acceptorSitePos - 2, 2);
		return this->returnFlankStringCaseFromFlankString(flankString);
	}

	int returnFlankStringCaseFromFlankString(const string& flankString)
	{
		if(flankString == "GTAG")
		{
			return 5;
		}
		else if(flankString == "CTAC")
		{
			return 6;
		}
		else if(flankString == "ATAC")
		{
			return 1;
		}
		else if(flankString == "GTAT")
		{
			return 2;
		}
		else if(flankString == "CTGC")
		{
			return 3;
		}
		else if(flankString == "GCAG")
		{
			return 4;
		}
		else
		{
			return 0;
		}
	}



	bool generateBestTargetPathVec_fixGapAlso_shortSpliceJunction(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_shortSpliceJunction, 
			vector<int>& possiblePathVec_shortSpliceJunction_pathIndex, Index_Info* indexInfo, const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		bool shortSpliceJunction_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		vector<int> tmpFixGapVec_shortSpliceJunction_index;
		for(int tmp = 0; tmp < possiblePathVec_shortSpliceJunction.size(); tmp++)
		{
			//cout << "try to fix: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;

			int tmpPathLastSegGroupNO = possiblePathVec_shortSpliceJunction[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_shortSpliceJunction[tmp].second;
 
			//cout << "try to fix shortSpliceJunction: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				//cout << "error in generateBestTargetPathVec_fixGapAlso fixshortSpliceJunction " 
 				//	<< tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
 				//return false;
 				tmpFixGapVec_shortSpliceJunction_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//int currentFixGapVec_index = fixGapVec.size();
 			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess,
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
 			tmpFixGapVec_shortSpliceJunction_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}

		// select the best shortSpliceJunction ...
		for(int tmp = 0; tmp < tmpFixGapVec_shortSpliceJunction_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_shortSpliceJunction_index[tmp];
			int tmpPathIndex = possiblePathVec_shortSpliceJunction_pathIndex[tmp];
			bool tmpShortSpliceJunctionGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpShortSpliceJunctionGapFixed_bool)
			{
				shortSpliceJunction_gap_fixed_bool = true;
				targetPathVec.push_back(tmpPathIndex);
				targetPathVec_gapIndex.push_back(tmpGapIndex);
			}
		}
		return shortSpliceJunction_gap_fixed_bool;
	}

	bool generateBestTargetPathVec_fixGapAlso_longSpliceJunction(
			Seg_Info* segInfo, int segGroupNO, int segCandiNO, //int currentPathNum, 
			vector<int>& targetPathVec, vector<int>& targetPathVec_gapIndex,
			vector< pair<int,int> >& possiblePathVec_longSpliceJunction, 
			vector<int>& possiblePathVec_longSpliceJunction_pathIndex, Index_Info* indexInfo, const string& readSeq_inProcess,
			bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
			int spliceJunctionDistanceMax
			)
	{
		bool longSpliceJunction_gap_fixed_bool = false;
		
		int currentFixGapVec_index = fixGapVec.size();
		vector<int> tmpFixGapVec_longSpliceJunction_index;
		for(int tmp = 0; tmp < possiblePathVec_longSpliceJunction.size(); tmp++)
		{
			int tmpPathLastSegGroupNO = possiblePathVec_longSpliceJunction[tmp].first;
			int tmpPathLastSegCandiNO = possiblePathVec_longSpliceJunction[tmp].second;

			//cout << "try to fix longSpliceJunction: " << tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
			int tmpFoundIndexInFixGapMap;
 			if(this->foundInFixGapMap_index(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO, tmpFoundIndexInFixGapMap)) // to test: should never happen
 			{
 				//cout << "error in generateBestTargetPathVec_fixGapAlso fixlongSpliceJunction: " 
 				//	<< tmpPathLastSegGroupNO << " " << tmpPathLastSegCandiNO << " " << segGroupNO << " " << segCandiNO << endl;
 				//return false;
 				tmpFixGapVec_longSpliceJunction_index.push_back(tmpFoundIndexInFixGapMap);
 				continue;
 			}				
 			//int currentFixGapVec_index = fixGapVec.size();
 			this->pushBackAndFixGap(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, segGroupNO, segCandiNO,
 				currentFixGapVec_index, segInfo, indexInfo, readSeq_inProcess,
 				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax);
 			tmpFixGapVec_longSpliceJunction_index.push_back(currentFixGapVec_index);
			currentFixGapVec_index ++;
		}

		// select the best longSpliceJunction ...
		for(int tmp = 0; tmp < tmpFixGapVec_longSpliceJunction_index.size(); tmp++)
		{
			int tmpGapIndex = tmpFixGapVec_longSpliceJunction_index[tmp];
			int tmpPathIndex = possiblePathVec_longSpliceJunction_pathIndex[tmp];
			bool tmpLongSpliceJunctionGapFixed_bool = fixGapVec[tmpGapIndex].first;
			if(tmpLongSpliceJunctionGapFixed_bool)
			{
				longSpliceJunction_gap_fixed_bool = true;
				targetPathVec.push_back(tmpPathIndex);
				targetPathVec_gapIndex.push_back(tmpGapIndex);
			}
		}
		return longSpliceJunction_gap_fixed_bool;
	}

	int minDistanceWithCurrentPath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, int currentPathNum, int* minSegNumGap,
		int spliceJunctionDistanceMax)
	{
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		int tmpMinSegNumGap = 100;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO, spliceJunctionDistanceMax);
			int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
			
			if(tmpSegNumGap <= 0)
			{
				continue;
			}

			if(tmpSegDistance < spliceJunctionDistanceMax)
			{
				if(tmpSegNumGap < tmpMinSegNumGap)
				{
					tmpMinSegNumGap = tmpSegNumGap;
					minDistance_value = tmpSegDistance;
					minDistance_path = tmpPathElementNO;
				}
				else if(tmpSegNumGap == tmpMinSegNumGap)
				{
					if(tmpSegDistance < minDistance_value)
					{
						minDistance_value = tmpSegDistance;
						minDistance_path = tmpPathElementNO;
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
		(*minSegNumGap) = tmpMinSegNumGap;
		return minDistance_value;		
	}

	int minDistanceWithCurrentPath_cirRNA(Seg_Info* segInfo, int segGroupNO, int segCandiNO, int currentPathNum, int* minSegNumGap)
	{
		int minDistance_value = 900000;
		int minDistance_path = 1000;
		int tmpMinSegNumGap = 100;
		for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
		{
			int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
			int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
			int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
			int tmpSegDistance = segInfo->distanceBetweenSegment_cirRNA(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
				segGroupNO, segCandiNO);
			int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
			
			if(tmpSegNumGap <= 0)
			{
				continue;
			}

			if(tmpSegDistance < MAX_SPLICE_DISTANCE_PHASE1)
			{
				if(tmpSegNumGap < tmpMinSegNumGap)
				{
					tmpMinSegNumGap = tmpSegNumGap;
					minDistance_value = tmpSegDistance;
					minDistance_path = tmpPathElementNO;
				}
				else if(tmpSegNumGap == tmpMinSegNumGap)
				{
					if(tmpSegDistance < minDistance_value)
					{
						minDistance_value = tmpSegDistance;
						minDistance_path = tmpPathElementNO;
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
		(*minSegNumGap) = tmpMinSegNumGap;
		return minDistance_value;		
	}

	void addNewSegCandiToCurrentPathInfoVec_matchIndelUniqueSpliceMultiPath_cirRNA(Seg_Info* segInfo, 
		int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	{
		bool relatedToSomePath = false;
		int minSegNumGap = 0;
		int minDistance = this->minDistanceWithCurrentPath_cirRNA(segInfo, segGroupNO, segCandiNO, currentPathNum, &minSegNumGap);
		
		//int currentPathNum = PathVec_seg.size();
		if(minDistance <= MAX_DELETION_LENGTH)
		{
			//cout << "PathVec_seg.size(): " << PathVec_seg.size() << endl;

			for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO++)
			{
				//cout << tmpPathElementNO << endl;
				//cout << " ... PathVec_seg.size(): " << PathVec_seg.size() << endl;
				int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
				int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
				int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
				int tmpSegDistance = segInfo->distanceBetweenSegment_cirRNA(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
					segGroupNO, segCandiNO);
				int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
				if((tmpSegDistance <= MAX_DELETION_LENGTH)&&(tmpSegNumGap <= minSegNumGap))
				{
					//cout << "tmpPathElementNO: " << tmpPathElementNO << endl;
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);
					relatedToSomePath = true;
				}
			}
		}
		else
		{
			for(int tmpPathElementNO = 0; tmpPathElementNO < currentPathNum; tmpPathElementNO ++)
			{

				int tmpPathElementVecSize = (PathVec_seg[tmpPathElementNO]).size();
				int tmpPathLastSegGroupNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).first;
				int tmpPathLastSegCandiNO = ((PathVec_seg[tmpPathElementNO])[tmpPathElementVecSize-1]).second;
				int tmpSegDistance = segInfo->distanceBetweenSegment_cirRNA(tmpPathLastSegGroupNO, tmpPathLastSegCandiNO, 
					segGroupNO, segCandiNO);
				int tmpSegNumGap = segGroupNO - tmpPathLastSegGroupNO;
				if((tmpSegDistance < MAX_SPLICE_DISTANCE_PHASE1)&&(tmpSegNumGap <= minSegNumGap))
				{
					relatedToSomePath = true;
					
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);					
				}
				/*if(relatedOrNot)
				{
					relatedToSomePath = true;
					
					PathValidBoolVec[tmpPathElementNO] = false;
					this->copyOldPath_AddSegCandi_AddNewPath(tmpPathElementNO, segGroupNO, segCandiNO);
					relatedToSomePath = true;
					//(PathVec_seg[tmpPathElementNO]).push_back(pair<int,int> (segGroupNO, segCandiNO));
				}*/
			}
		}

		if((!relatedToSomePath)&&
			(longSegBool))
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (segGroupNO, segCandiNO) );	
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);		
			vector< pair<int,int> >().swap(tmpPathElementVec);
		}
	}
	
	void addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction(
		Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum,
		int spliceJunctionDistanceMax)
	{
		vector<int> targetPathVec;

		this->generateBestCurrentPathVecForNewSegCandi(
			segInfo, segGroupNO, segCandiNO, currentPathNum, targetPathVec, spliceJunctionDistanceMax);
		
		if(targetPathVec.size() > 0)
		{
			for(int tmp = 0; tmp < targetPathVec.size(); tmp++)
			{
				int tmpTargetPath = targetPathVec[tmp];
				PathValidBoolVec[tmpTargetPath] = false;
				this->copyOldPath_AddSegCandi_AddNewPath(tmpTargetPath, segGroupNO, segCandiNO);			
			}
		}
		else if(longSegBool)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (segGroupNO, segCandiNO) );	
			PathVec_seg.push_back(tmpPathElementVec);
			PathValidBoolVec.push_back(true);		
			vector< pair<int,int> >().swap(tmpPathElementVec);
		}
		else
		{

		}
	}

	void addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction_fixGapAlso(
		Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum, 
		Index_Info* indexInfo, const string& readSeq_inProcess,
		bool annotation_provided_bool, bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax)
	{
		//cout << endl << "addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction" << endl;
		//cout << "segGroupNO: " << segGroupNO << endl;
		//cout << "segCandiNO: " << segCandiNO << endl;

		vector<int> targetPathVec;
		vector<int> targetPathVec_gapIndex;
		
		this->generateBestCurrentPathVecForNewSegCandi_fixGapAlso(
			segInfo, segGroupNO, segCandiNO, currentPathNum, targetPathVec, targetPathVec_gapIndex, 
			readSeq_inProcess, indexInfo, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
			spliceJunctionDistanceMax);

		if(targetPathVec.size() > 0)
		{
			for(int tmp = 0; tmp < targetPathVec.size(); tmp++)
			{
				int tmpTargetPath = targetPathVec[tmp];
				int tmpTargetPath_gapIndex = targetPathVec_gapIndex[tmp];
				PathValidBoolVec[tmpTargetPath] = false;
				this->copyOldPath_AddSegCandi_AddNewPath_fixGapAlso(
					tmpTargetPath, segGroupNO, segCandiNO, tmpTargetPath_gapIndex);			
			}
		}
		else if(longSegBool)
		{
			vector< pair<int,int> > tmpPathElementVec;
			tmpPathElementVec.push_back( pair<int,int> (segGroupNO, segCandiNO) );	
			PathVec_seg.push_back(tmpPathElementVec);

			vector<int> tmpPathElementVec_gapIndex;
			PathVec_seg_gapIndex.push_back(tmpPathElementVec_gapIndex);

			PathValidBoolVec.push_back(true);		
			vector< pair<int,int> >().swap(tmpPathElementVec);
			vector<int>().swap(tmpPathElementVec_gapIndex);
		}
		else
		{

		}
		//cout << "end of addNewSegCandiToCurrentPathInfoVec_prefer_Match_Indel_SmallSpliceJunction_fixGapAlso ..." << endl << endl;
	}
	/*
	void addNewSegCandiToCurrentPathInfoVec_uniquePath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	void addNewSegCandiToCurrentPathInfoVec_allPath(Seg_Info* segInfo, int segGroupNO, int segCandiNO, bool longSegBool, int currentPathNum)
	*/
	
	bool getPossiPathFromSeg(Seg_Info* segInfo, int spliceJunctionDistanceMax)
	{
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec(segInfo, tmpSegNO, spliceJunctionDistanceMax);
		}

		possiPathExists = true;
		return possiPathExists;
	}


	bool getPossiPathFromSeg_incompleteHead(Seg_Info* segInfo, int spliceJunctionDistanceMax)
	{
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec(segInfo, tmpSegNO, spliceJunctionDistanceMax);
		}

		possiPathExists = true;
		return possiPathExists;
	}

	/*
	bool getPossiPathFromSeg_fixOneEndUnmapped(Seg_Info* segInfo)
	{
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec_localGreedyMapping(segInfo, tmpSegNO);
		}

		possiPathExists = true;
		return possiPathExists;
	}*/

	bool getPossiPathFromSeg_fixOneEndUnmapped(Seg_Info* segInfo, int spliceJunctionDistanceMax)
	{
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec_localGreedyMapping(segInfo, tmpSegNO, spliceJunctionDistanceMax);
		}

		possiPathExists = true;
		return possiPathExists;
	}


	bool getPossiPathFromSeg_incompleteTail(Seg_Info* segInfo, int spliceJunctionDistanceMax)
	{

		bool possiPathExists = false;
		int firstLongSegNO = 0;//segInfo->getFirstLongSegNO();
			
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec_localGreedyMapping(segInfo, tmpSegNO, spliceJunctionDistanceMax);
		}

		possiPathExists = true;
		return possiPathExists;
	}

	bool getPossiPathFromSeg_cirRNA(Seg_Info* segInfo)
	{
		if(segInfo->returnRepeatRegion_index() >= 0)
		{
			repeatRegion_index = segInfo->returnRepeatRegion_index();
			return true;
		}

		bool possiPathExists = false;
		int firstLongSegNO = segInfo->getFirstLongSegNO();
		
		if(firstLongSegNO < 0)
		{
			//cout << "error in getPossiPathFromSeg " << endl;
			return false;
		}
		
			for(int tmpSegCandiLoc = 0; tmpSegCandiLoc < segInfo->returnSegmentAlignNum(firstLongSegNO); tmpSegCandiLoc++)
			{
				vector< pair<int,int> > tmpPathElementVec;
				tmpPathElementVec.push_back( pair<int,int> (firstLongSegNO, tmpSegCandiLoc) );
				PathVec_seg.push_back(tmpPathElementVec);
				PathValidBoolVec.push_back(true);
				vector< pair<int,int> > ().swap(tmpPathElementVec);
			}

		for(int tmpSegNO = firstLongSegNO + 1; tmpSegNO < segInfo->returnSegmentNum(); tmpSegNO ++)
		{
			if(//(!segInfo->checkSegLongOrNot(tmpSegNO)) && 
				(segInfo->returnSegmentAlignNum(tmpSegNO)) > CANDALILOC)
				continue;
			this->addNewSegGroupToCurrentPathInfoVec_cirRNA(segInfo, tmpSegNO);
		}

		possiPathExists = true;
		return possiPathExists;
	}

	string possiPathStr()
	{
		string possiPathStr = "Path Info: \n";
		for(int tmpPath = 0; tmpPath < PathVec_seg.size(); tmpPath++)
		{
			possiPathStr = possiPathStr + "... Path " + int_to_str(tmpPath+1) + ": ";
			for(int tmpPathElement = 0; tmpPathElement < PathVec_seg[tmpPath].size(); tmpPathElement++)
			{
				int tmpSegGroupNO = ((PathVec_seg[tmpPath])[tmpPathElement]).first;
				int tmpSegCandiNO = ((PathVec_seg[tmpPath])[tmpPathElement]).second;
				possiPathStr = possiPathStr + int_to_str(tmpSegGroupNO+1) + "," + int_to_str(tmpSegCandiNO+1) + "--";
			}
			if(PathValidBoolVec[tmpPath])
			{
				possiPathStr += " -- true";
			}
			else
			{
				possiPathStr += " -- false";
			}
			possiPathStr += "\n";
		}
		return possiPathStr;
	}

};

#endif
