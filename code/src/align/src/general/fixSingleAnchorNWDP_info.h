// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXSINGLEANCHORNWDP_INFO_H
#define FIXSINGLEANCHORNWDP_INFO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "nw_DP.h"

using namespace std;

class FixSingleAnchor_NWDP_Info
{
private:
	//******************** Final Results *************//
	vector<Jump_Code> gapJumpCodeVec;
	vector<int> gapMismatchPosVec;
	vector<char> gapMismatchCharVec;

	int penalty_mismatch;
	int penalty_insertion;
	int penalty_deletion;
public:
	FixSingleAnchor_NWDP_Info()
	{
		penalty_mismatch = 1;
		penalty_insertion = 1;
		penalty_deletion = 1;
	}

	int returnBestGapJumpCodeVecSize()
	{
		return gapJumpCodeVec.size();
	}

	int returnBestGapMismatchPosVecSize()
	{
		return gapMismatchPosVec.size();
	}
	int returnBestGapMismatchCharVecSize()
	{
		return gapMismatchCharVec.size();
	}

	int returnBestGapMismatchPos(int index)
	{
		return gapMismatchPosVec[index];
	}
	char returnBestGapMismatchChar(int index)
	{
		return gapMismatchCharVec[index];
	}

	void copyMismatchPosVec2TargetVec(vector<int>& targetVec)
	{
		for(int tmp = 0; tmp < gapMismatchPosVec.size(); tmp++)
		{
			targetVec.push_back(gapMismatchPosVec[tmp]);
		}
	}

	void copyMismatchCharVec2TargetVec(vector<char>& targetVec)
	{
		for(int tmp = 0; tmp < gapMismatchCharVec.size(); tmp++)
		{
			targetVec.push_back(gapMismatchCharVec[tmp]);
		}
	}


	void doNWDP(string& tmpReadSubSeq,
		string& tmpChrSubSeq) // 1-end open NWDP, 
	{
		nw_oneEndOpen(tmpReadSubSeq, tmpChrSubSeq, 
			gapJumpCodeVec, gapMismatchPosVec, gapMismatchCharVec);
	}	

	void doNWDP_withMismatchJumpCode(string& tmpReadSubSeq,
		string& tmpChrSubSeq) // 1-end open NWDP, 
	{
		nw_oneEndOpen_withMismatchJumpCode(tmpReadSubSeq, tmpChrSubSeq, 
			gapJumpCodeVec, gapMismatchPosVec, gapMismatchCharVec);
	}

	void doNWDP_reverse(string& tmpReadSubSeq,
		string& tmpChrSubSeq) // 1-end open NWDP, 
	{
		nw_oneEndOpen_reverse(tmpReadSubSeq, tmpChrSubSeq, 
			gapJumpCodeVec, gapMismatchPosVec, gapMismatchCharVec);
	}	

	void doNWDP_withMismatchJumpCode_reverse(string& tmpReadSubSeq,
		string& tmpChrSubSeq) // 1-end open NWDP, 
	{
		nw_oneEndOpen_withMismatchJumpCode_reverse(tmpReadSubSeq, tmpChrSubSeq, 
			gapJumpCodeVec, gapMismatchPosVec, gapMismatchCharVec);
	}

	void copyJumpCodeVec2TargetVec(vector<Jump_Code>& targetVec)
	{
		int vecSize = this->returnBestGapJumpCodeVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			targetVec.push_back(gapJumpCodeVec[tmp]);
		}
	}

	int getPenalty()
	{
		int tmpPenalty = 0;
		for(int tmp = 0; tmp < gapJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = gapJumpCodeVec[tmp].type;
			int tmpJumpCodeLen = gapJumpCodeVec[tmp].len;
			if(tmpJumpCodeType == "I")
			{
				tmpPenalty += tmpJumpCodeLen*penalty_insertion;
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpPenalty += tmpJumpCodeLen*penalty_deletion;
			}
		}
		tmpPenalty += gapMismatchPosVec.size()*penalty_mismatch;
		return tmpPenalty;
	}

	int getLength()
	{
		int tmpLength = 0;
		for(int tmp = 0; tmp < gapJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = gapJumpCodeVec[tmp].type;
			int tmpJumpCodeLen = gapJumpCodeVec[tmp].len;
			if((tmpJumpCodeType == "X")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "I"))
			{
				tmpLength += tmpJumpCodeLen;
			}

		}
		return tmpLength;
	}

	bool fixedOrNot()
	{
		int penalty = this->getPenalty();
		int length = this->getLength();
		if(penalty > (1 + length/5))
			return false;
		else
			return true;
	}
	bool fixedOrNot(int length)
	{
		int penalty = this->getPenalty();
		if(penalty > (1 + length/5))
			return false;
		else
			return true;
	}	
};

#endif