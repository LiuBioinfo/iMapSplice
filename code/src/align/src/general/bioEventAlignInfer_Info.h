// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef BIOEVENTALIGNINFER_INFO_H
#define BIOEVENTALIGNINFER_INFO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include <sstream>

#include "../constantDefinitions.h"
#include "splice_info.h"
#include "fixDoubleAnchorMatch_info.h"

using namespace std;

class BioEventAlignInfer_Info
{
private:
	int chrNameInt;
	int bioEventDonerEndPos;
	int bioEventLen,
	string bioEventType,

	int pathNum_backward;
	int pathNum_forward;

	vector< vector<Jump_Code> > pathVec_backward;
	vector< vector<Jump_Code> > pathVec_forward;

	// Note: support numbers for paths are not that accurate as some alignment with 
	// few number bases beside SJ may comply with muliple paths, but only the first
	// complied one got extended.
	vector< int > pathSupportNumVec_backward; 
	vector< int > pathSupportNumVec_forward;

	int supportNum;

public:
	BioEventAlignInfer_Info()
	{
		supportNum = 0;
		pathNum_backward = 0;
		pathNum_forward = 0;		
	}

	void initiateBioEventAlignInferInfo(
		int mapChrNameInt, int tmpBioEventDonerEndPos, 
		int tmpBioEventLen, const string& tmpBioEventType, 
		vector<Jump_Code>& tmpBioEventJumpCodeVec_backward, 
		vector<Jump_Code>& tmpBioEventJumpCodeVec_forward)
	{
		chrNameInt = mapChrNameInt;
		bioEventDonerEndPos = tmpBioEventDonerEndPos;
		bioEventLen = tmpBioEventLen;
		bioEventType = tmpBioEventType;

		pathVec_backward.push_back(tmpBioEventJumpCodeVec_backward);
		pathSupportNumVec_backward.push_back(1);
		pathNum_backward = 1;	

		pathVec_forward.push_back(tmpBioEventJumpCodeVec_forward);
		pathSupportNumVec_forward.push_back(1);
		pathNum_forward = 1;		
	}

	void updateOldBioEventAlignInferInfo(
		vector<Jump_Code>& tmpBioEventJumpCodeVec_backward,
		vector<Jump_Code>& tmpBioEventJumpCodeVec_forward, 
		int maxReadBaseNum)
	{
		pair<int,int> checkWithOldPath_backward_relation_pair 
			= this->checkWithOldPathVec_backward(tmpBioEventJumpCodeVec_backward);
		int compliedOldPathIndex_backward = checkWithOldPath_backward_relation_pair.second;
		int checkWithOldPath_backward_relation = checkWithOldPath_backward_relation_pair.first;			
		if(checkWithOldPath_backward_relation == -1) // no need to update path
		{
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
		}
		else if(checkWithOldPath_backward_relation == -2) // new path
		{
			pathVec_backward.push_back(tmpBioEventJumpCodeVec_backward);
			pathSupportNumVec_backward.push_back(1);
			pathNum_backward ++;
		}
		else if(checkWithOldPath_backward_relation == 0) // update
		{
			pathVec_backward[compliedOldPathIndex_backward] = tmpBioEventJumpCodeVec_backward;
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
		}
		else
		{
			cout << "checkWithOldPath_backward_relation: " << checkWithOldPath_backward_relation << endl;
			cout << "error in checkWithOldPath_backward_relation ..." << endl;
			exit(1);
		}

		pair<int,int> checkWithOldPath_forward_relation_pair
			= this->checkWithOldPathVec_forward(tmpBioEventJumpCodeVec_forward);
		int compliedOldPathIndex_forward = checkWithOldPath_forward_relation_pair.second;
		int checkWithOldPath_forward_relation = checkWithOldPath_forward_relation_pair.first; 
		if(checkWithOldPath_forward_relation == -1) // no need to update path
		{
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
		}
		else if(checkWithOldPath_forward_relation == -2) // new path
		{
			pathVec_forward.push_back(tmpBioEventJumpCodeVec_forward);
			pathSupportNumVec_forward.push_back(1);
			pathNum_forward ++;
		}
		else if(checkWithOldPath_forward_relation == 0) // update
		{
			pathVec_forward[compliedOldPathIndex_forward] = tmpBioEventJumpCodeVec_forward;
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
		}
		else
		{
			cout << "checkWithOldPath_forward_relation: " << checkWithOldPath_forward_relation << endl;
			cout << "error in checkWithOldPath_forward_relation ..." << endl;
			exit(1);
		}

		this->updateSupportNum();
	}

	pair<int,int> checkWithOldPathVec_backward(
		vector<Jump_Code>& tmpBioEventJumpCodeVec_backward)
	{
		// check with the old paths and return the relation checking results:
		// Output:
		//	<0, indexInOldPath(>=0)>, -- comply with an old path and need to extend
		//  or
		//	<-1, indexInOldPath(>=0)> -- comply with an old path but no need to extend
		//	or
		//	<-2, -2> -- new path	
		// Methods:				 
		// 1. comply with one of the old paths
		//		1.1 longer than the old path
		//			1.1.1 bases num in old path < ALIGNINFER_BASE_NUM_MAX -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX), 
		//				 return <0,indexInOldPath>
		//			1.1.2 bases num in old path == ALIGNINFER_BASE_NUM_MAX -- skip it 
		//				 return <-1,indexInOldPath>;
		//		1.2 shorter than the old path or the same with it -- skip it
		//				return <-1,indexInOldPath>;
		// 2. got a new path -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX)
		//			return <-2,-2>
		for(int tmp = 0; tmp < pathVec_backward.size(); tmp++)
		{
			int checkOldPath_relation 
				= this->checkWithOldPath_backward(newSJjumpCodeVec_backward, tmp);
			if(checkOldPath_relation == 0)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(-1,tmp);
			}
			else if(checkOldPath_relation == 1)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(0,tmp);
			}
			else
			{}
		}
		return pair<int,int>(-2,-2);
	}

	int checkWithOldPath_backward(
		vector<Jump_Code>& newSJjumpCodeVec_backward, int index_oldPathVec_backward)
		// return: 1 -- comply with old path, need to update
		// 		   0 -- comply with old path, no need to update
		//        -1 -- new path
	{
		int newSJjumpCodeVecSize = newSJjumpCodeVec_backward.size();
		int oldSJjumpCodeVecSize = (pathVec_backward[index_oldPathVec_backward]).size();
		int minSJjumpCodeVecSize = newSJjumpCodeVecSize;
		if(oldSJjumpCodeVecSize < minSJjumpCodeVecSize)
			minSJjumpCodeVecSize = oldSJjumpCodeVecSize;

		for(int tmp = 0; tmp < minSJjumpCodeVecSize-1; tmp++) // check first N-1 JumpCode
		{
			int newSJjumpCodeLength = newSJjumpCodeVec_backward[tmp].len;
			string newSJjumpCodeType = newSJjumpCodeVec_backward[tmp].type;
			int oldSJjumpCodeLength = (pathVec_backward[index_oldPathVec_backward])[tmp].len;
			string oldSJjumpCodeType = (pathVec_backward[index_oldPathVec_backward])[tmp].type;			

			if( (newSJjumpCodeLength != oldSJjumpCodeLength) || (newSJjumpCodeType != oldSJjumpCodeType))
				return -1;
		}

		// check last shared JumpCode
		int newSJjumpCodeLength_lastShared = newSJjumpCodeVec_backward[minSJjumpCodeVecSize-1].len;
		string newSJjumpCodeType_lastShared = newSJjumpCodeVec_backward[minSJjumpCodeVecSize-1].type;
		int oldSJjumpCodeLength_lastShared = (pathVec_backward[index_oldPathVec_backward])[minSJjumpCodeVecSize-1].len;
		string oldSJjumpCodeType_lastShared = (pathVec_backward[index_oldPathVec_backward])[minSJjumpCodeVecSize-1].type;	
		if(newSJjumpCodeType_lastShared == oldSJjumpCodeType_lastShared)
		{
			if(newSJjumpCodeVecSize < oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return -1;
			}
			else if(newSJjumpCodeVecSize == oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return 1;				
			}
			else // newSJjumpCodeVecSize > oldSJjumpCodeVecSize
			{
				if(newSJjumpCodeLength_lastShared < oldSJjumpCodeLength_lastShared)
					return -1;
				else
					return 1;				
			}
		}
		else // last shared JumpCode inconsistent
		{
			return -1;
		}
	}

	pair<int,int> checkWithOldPathVec_forward(
		vector<Jump_Code>& newSJjumpCodeVec_forward)
	{
		// check with the old paths and return the relation checking results:
		// Output:
		//	<0, indexInOldPath(>=0)>, -- comply with an old path and need to extend
		//  or
		//	<-1, indexInOldPath(>=0)> -- comply with an old path but no need to extend
		//	or
		//	<-2, -2> -- new path	
		// Methods:				 
		// 1. comply with one of the old paths
		//		1.1 longer than the old path
		//			1.1.1 bases num in old path < ALIGNINFER_BASE_NUM_MAX -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX), 
		//				 return <0,indexInOldPath>
		//			1.1.2 bases num in old path == ALIGNINFER_BASE_NUM_MAX -- skip it 
		//				 return <-1,indexInOldPath>;
		//		1.2 shorter than the old path or the same with it -- skip it
		//				return <-1,indexInOldPath>;
		// 2. got a new path -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX)
		//			return <-2,-1>
		for(int tmp = 0; tmp < pathVec_forward.size(); tmp++)
		{
			int checkOldPath_relation 
				= this->checkWithOldPath_forward(newSJjumpCodeVec_forward, tmp);
			if(checkOldPath_relation == 0)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(-1,tmp);
			}
			else if(checkOldPath_relation == 1)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(0,tmp);
			}
			else
			{}
		}
		return pair<int,int>(-2,-2);
	}

	int checkWithOldPath_forward(
		vector<Jump_Code>& newSJjumpCodeVec_forward, int index_oldPathVec_forward)
		// return: 1 -- comply with old path, need to update
		// 		   0 -- comply with old path, no need to update
		//        -1 -- new path
	{
		int newSJjumpCodeVecSize = newSJjumpCodeVec_forward.size();
		int oldSJjumpCodeVecSize = (pathVec_forward[index_oldPathVec_forward]).size();
		int minSJjumpCodeVecSize = newSJjumpCodeVecSize;
		if(oldSJjumpCodeVecSize < minSJjumpCodeVecSize)
			minSJjumpCodeVecSize = oldSJjumpCodeVecSize;

		for(int tmp = 0; tmp < minSJjumpCodeVecSize-1; tmp++) // check first N-1 JumpCode
		{
			int newSJjumpCodeLength = newSJjumpCodeVec_forward[tmp].len;
			string newSJjumpCodeType = newSJjumpCodeVec_forward[tmp].type;
			int oldSJjumpCodeLength = (pathVec_forward[index_oldPathVec_forward])[tmp].len;
			string oldSJjumpCodeType = (pathVec_forward[index_oldPathVec_forward])[tmp].type;			

			if( (newSJjumpCodeLength != oldSJjumpCodeLength) || (newSJjumpCodeType != oldSJjumpCodeType))
				return -1;
		}

		// check last shared JumpCode
		int newSJjumpCodeLength_lastShared = newSJjumpCodeVec_forward[minSJjumpCodeVecSize-1].len;
		string newSJjumpCodeType_lastShared = newSJjumpCodeVec_forward[minSJjumpCodeVecSize-1].type;
		int oldSJjumpCodeLength_lastShared = (pathVec_forward[index_oldPathVec_forward])[minSJjumpCodeVecSize-1].len;
		string oldSJjumpCodeType_lastShared = (pathVec_forward[index_oldPathVec_forward])[minSJjumpCodeVecSize-1].type;	
		if(newSJjumpCodeType_lastShared == oldSJjumpCodeType_lastShared)
		{
			if(newSJjumpCodeVecSize < oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return -1;
			}
			else if(newSJjumpCodeVecSize == oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return 1;				
			}
			else // newSJjumpCodeVecSize > oldSJjumpCodeVecSize
			{
				if(newSJjumpCodeLength_lastShared < oldSJjumpCodeLength_lastShared)
					return -1;
				else
					return 1;				
			}
		}
		else // last shared JumpCode inconsistent
		{
			return -1;
		}
	}


	void updateSupportNum()
	{
		supportNum ++;
	}
};