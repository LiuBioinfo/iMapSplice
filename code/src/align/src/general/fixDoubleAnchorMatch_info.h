// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORMATCH_INFO_H
#define FIXDOUBLEANCHORMATCH_INFO_H

//#include "nw_DP.h"

class FixDoubleAnchor_Match_Info
{
private:
	vector<int> mismatchPosInReadVec;
	vector<char> mismatchCharVec;
	int mismatchNum;
	int mismatchPos_interval_min;
public:
	FixDoubleAnchor_Match_Info()
	{
		mismatchNum = 0;
		mismatchPos_interval_min = 3;
	}

		int getMinMismatchInterval(vector<int>& mismatchPosInReadVec_sorted)
		{
				for(int tmp = 0; tmp < mismatchPosInReadVec.size(); tmp++)
				{
					mismatchPosInReadVec_sorted.push_back(mismatchPosInReadVec[tmp]);
				}				
				sort(mismatchPosInReadVec_sorted.begin(), mismatchPosInReadVec_sorted.end());

				int tmpMinMismatchPosInterval = 100;
				for(int tmp = 0; tmp < mismatchPosInReadVec_sorted.size()-1; tmp++)
				{
					int tmpIndex_1 = tmp;
					int tmpIndex_2 = tmp + 1;
					int tmpMismatchPos_1 = mismatchPosInReadVec_sorted[tmpIndex_1];
					int tmpMismatchPos_2 = mismatchPosInReadVec_sorted[tmpIndex_2];
					int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
						tmpMinMismatchPosInterval = tmpMismatchPosInterval;
				}
				return tmpMinMismatchPosInterval;
		}
		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> mismatchPosInReadVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(mismatchPosInReadVec_sorted);
			//cout << "minimumMismatchGap: " << minimumMismatchGap << endl;
			//cout << "start to check sorted mismatchPosVec" << endl;
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < mismatchPosInReadVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = mismatchPosInReadVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = mismatchPosInReadVec_sorted[tmpIndex_2];
				//cout << "mismatchInterval-" << tmp+1 << ": " << tmpMismatchPos_1 << "~" << tmpMismatchPos_2 << endl;
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
			//cout << "start to do filterInterResults ...." << endl;
			int mismatchPosVecSize = mismatchPosInReadVec.size();
			//cout << "mismatchPosVecSize: " << mismatchPosVecSize << endl;
			if(mismatchPosVecSize >= 3)
			{
				int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
				//cout << "tmpMinMismatchPosInterval_2: " << tmpMinMismatchPosInterval << endl; 
				if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
					return false;
			}
			return true;
		}

	int returnMismatchNum()
	{
		return mismatchNum;
	}

	void copyMismatchPos2TargetVec(vector<int>& targetVec)
	{
		for(int tmp = 0; tmp < mismatchPosInReadVec.size(); tmp++)
		{
			targetVec.push_back(mismatchPosInReadVec[tmp]);
		}
	}
	void copyMismatchChar2TargetVec(vector<char>& targetVec)
	{
		for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
		{
			targetVec.push_back(mismatchCharVec[tmp]);
		}
	}
	int returnMismatchPosInRead(int mismatchPosVec_index)
	{
		return mismatchPosInReadVec[mismatchPosVec_index];
	}
	int returnMismatchPos(int mismatchPosVec_index)
	{
		return mismatchPosInReadVec[mismatchPosVec_index];
	}
	int returnMismatchPosInReadVecSize()
	{
		return mismatchPosInReadVec.size();
	}
	int returnMismatchPosVecSize()
	{
		return mismatchPosInReadVec.size();
	}
	char returnMismatchChar(int mismatchCharVec_index)
	{
		return mismatchCharVec[mismatchCharVec_index];
	}
	int returnMismatchCharVecSize()
	{
		return mismatchCharVec.size();
	}

	bool fixMatch(const string& s1, const string& s2, int max_mismatch, int startPosInRead)
	{
		if(STORE_MISMATCH_POS)
		{
			if(STORE_MISMATCH_CHA)
			{
				return this->fixMatch_MismatchNum_MismatchPos_MismatchChar(s1, s2, max_mismatch, startPosInRead);
			}
			else
			{
				return this->fixMatch_MismatchNum_MismatchPos(s1, s2, max_mismatch, startPosInRead);
			}
		}
		else
		{
			return this->fixMatch_MismatchNum(s1, s2, max_mismatch);
		}
	}

	bool fixMatch_MismatchNum(const string& s1, const string& s2, int max_mismatch)
	{
		bool match_bool = true;
		if(s1.length() != s2.length())
		{
			cout << "s1.length != s2.length in score_string DoubleAnchorScore.h";
			match_bool = false;
		}
		//int mismatchNum = 0;
		for(int tmp = 0; tmp < s1.length(); tmp++)
		{
			if(s1[tmp] != s2[tmp])
			{	
				mismatchNum++;
				if(mismatchNum > max_mismatch)
				{
					match_bool = false; // not match
					break;
				}
			}
		}
		//num_mismatch = mismatchNum;
		return match_bool;
	}

	bool fixMatch_MismatchNum_MismatchPos(const string& s1, const string& s2, int max_mismatch, int startPosInRead)
	{
		bool match_bool = true;
		if(s1.length() != s2.length())
		{
			cout << "s1.length != s2.length in score_string DoubleAnchorScore.h";
			match_bool = false;
		}
		//int mismatchNum = 0;
		for(int tmp = 0; tmp < s1.length(); tmp++)
		{
			if(s1[tmp] != s2[tmp])
			{	
				int tmpMismatchPos = startPosInRead + tmp;
				mismatchPosInReadVec.push_back(tmpMismatchPos);
				mismatchNum++;
				if(mismatchNum > max_mismatch)
				{
					match_bool = false; // not match
					break;
				}
			}
		}
		//num_mismatch = mismatchNum;
		return match_bool;
	}	

	bool fixMatch_MismatchNum_MismatchPos_MismatchChar(const string& s1/*read string*/, const string& s2/*chrom string*/, int max_mismatch, int startPosInRead)
	{
		bool match_bool = true;
		if(s1.length() != s2.length())
		{
			cout << "s1.length != s2.length in score_string DoubleAnchorScore.h";
			match_bool = false;
		}
		//int mismatchNum = 0;
		for(int tmp = 0; tmp < s1.length(); tmp++)
		{
			if(s1[tmp] != s2[tmp])
			{	
				int tmpMismatchPos = startPosInRead + tmp;
				//cout << "tmpMismatchPos: " << tmpMismatchPos << endl;
				mismatchPosInReadVec.push_back(tmpMismatchPos);
				mismatchCharVec.push_back(s2[tmp]);
				mismatchNum++;
				if(mismatchNum > max_mismatch)
				{
					match_bool = false; // not match
					break;
				}
			}
		}
		if(match_bool)
		{
			bool filter_bool = this->filterInterResults();
			if(filter_bool)
				return true;
			else 
				return false;
		}
		else
		{
			return false;
		}
		//num_mismatch = mismatchNum;
		//return match_bool;
	}	
};

#endif