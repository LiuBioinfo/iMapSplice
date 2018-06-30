// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORNWDP_INFO_H
#define FIXDOUBLEANCHORNWDP_INFO_H

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

class FixDoubleAnchor_NWDP_Info
{
private:
	//******************** Final Results *************//
	vector<Jump_Code> gapJumpCodeVec;
	vector<int> gapMismatchPosVec;
	vector<char> gapMismatchCharVec;

	//******************** Inter Results ************//
	int mismatchPos_interval_min;
public:
	FixDoubleAnchor_NWDP_Info()
	{
		mismatchPos_interval_min = 3;
	}

	int getMinMismatchInterval(vector<int>& gapMismatchPosVec_sorted)
	{
		for(int tmp = 0; tmp < gapMismatchPosVec.size(); tmp++)
		{
			gapMismatchPosVec_sorted.push_back(gapMismatchPosVec[tmp]);
		}				
		sort(gapMismatchPosVec_sorted.begin(), gapMismatchPosVec_sorted.end());

		int tmpMinMismatchPosInterval = 100;
		for(int tmp = 0; tmp < gapMismatchPosVec_sorted.size()-1; tmp++)
		{
			int tmpIndex_1 = tmp;
			int tmpIndex_2 = tmp + 1;
			int tmpMismatchPos_1 = gapMismatchPosVec_sorted[tmpIndex_1];
			int tmpMismatchPos_2 = gapMismatchPosVec_sorted[tmpIndex_2];
			int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
			if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
				tmpMinMismatchPosInterval = tmpMismatchPosInterval;
		}
		return tmpMinMismatchPosInterval;
	}

		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> gapMismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(gapMismatchPosVec_sorted);
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < gapMismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = gapMismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = gapMismatchPosVec_sorted[tmpIndex_2];
				int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
				if(tmpMismatchPosInterval < minimumMismatchGap)
				{
					cout << "some thing wrong in getMinMismatchInterval_2" << endl;
					exit(1);
				}
				else if(tmpMismatchPosInterval == minimumMismatchGap)
				{	
					if(!checkedMinimumGap_bool)
						checkedMinimumGap_bool = true;
					else
						return tmpMismatchPosInterval;
				}
				else
				{
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval_2)
						tmpMinMismatchPosInterval_2 = tmpMismatchPosInterval;
				}
			}
			return tmpMinMismatchPosInterval_2;
		}

	bool filterInterResults()
	{
		int mismatchPosVecSize = gapMismatchPosVec.size();
		if(mismatchPosVecSize >= 3)
		{
			int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
			if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
				return false;
		}
		return true;
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

	void copyMismatchPos2TargetVec_startLoc(vector<int>& targetVec, int startLoc)
	{
		int vecSize = this->returnBestGapMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = this->returnBestGapMismatchPos(tmp) + startLoc - 1;
			targetVec.push_back(tmpMismatchPos);	
		}			
	}

	void copyMismatchPos2TargetVec(vector<int>& targetVec)
	{
		int vecSize = this->returnBestGapMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = this->returnBestGapMismatchPos(tmp);
			targetVec.push_back(tmpMismatchPos);	
		}	
	}
	void copyMismathChar2TargetVec(vector<char>& targetVec)
	{
		int vecSize = this->returnBestGapMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = this->returnBestGapMismatchChar(tmp);
			targetVec.push_back(tmpMismatchChar);
		}
	}

	void copyJumpCodeVec2TargetVec(vector<Jump_Code>& targetVec)
	{
		int vecSize = this->returnBestGapJumpCodeVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			targetVec.push_back(gapJumpCodeVec[tmp]);
		}
	}

	bool fix_NWDP(
		int toFix_read_start, int toFix_read_end,
		int toFix_chrom_start, int toFix_chrom_end,
		const string& readSeq_inProcess,
		Index_Info* indexInfo, int chrom_name_int,
		int max_allowed_mismatchNum)
	{
		//return false;
		//cout << "start in fix_NWDP" << endl;
		string readSeqToProcess 
			= readSeq_inProcess.substr(toFix_read_start-1, toFix_read_end - toFix_read_start + 1);
		string chromSeqToProcess 
			= indexInfo->returnChromStrSubstr(chrom_name_int, toFix_chrom_start, toFix_chrom_end - toFix_chrom_start + 1);
		//cout << "readSeqToProcess: " << readSeqToProcess << endl;
		//cout << "chromSeqToProcess: " << chromSeqToProcess << endl;

		this->doNWDP(readSeqToProcess, chromSeqToProcess);
		int fixedGap_mismatchNum = this->returnBestGapMismatchPosVecSize();

		//cout << "fixedGap_mismatchNum: " << fixedGap_mismatchNum << endl;
		if(fixedGap_mismatchNum > max_allowed_mismatchNum)
			return false;
		else
		{
			bool filter_bool = this->filterInterResults();
			if(filter_bool)	
				return true;
			else
				return false;
		}
	}

	void doNWDP(string& seq_1, string& seq_2) // Needleman-Wunsch algorithm for global alignment 
	{
		nw(seq_1, seq_2, gapJumpCodeVec, gapMismatchPosVec, gapMismatchCharVec);
	}
};

#endif