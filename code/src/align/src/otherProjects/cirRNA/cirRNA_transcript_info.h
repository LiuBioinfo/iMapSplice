// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef CIRRNA_TRANSCRIPT_INFO_H
#define CIRRNA_TRANSCRIPT_INFO_H

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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

class CirRNA_Transcript_Info
{
private:
	int chrNameInt;
	int backSpliceAcceptorStart_normalTranscriptStart_pos;
	int backSpliceDonerEnd_normalTranscriptEnd_pos;
	vector< pair<int,int> > innerSJsitePosPairVec;
	vector< pair<int,int> > exonPosPairVec;
	int readCount;
public:
	CirRNA_Transcript_Info()
	{
		readCount = 0;
	}

	int return_backSpliceAcceptorStart_normalTranscriptStart_pos()
	{
		return backSpliceAcceptorStart_normalTranscriptStart_pos;
	}

	int return_backSpliceDonerEnd_normalTranscriptEnd_pos()
	{
		return backSpliceDonerEnd_normalTranscriptEnd_pos;
	}

	void readCountIncrement()
	{
		readCount ++;
	}

	void initiateCirRNAtranscriptInfo(
		int tmpChrNameInt, int tmpBackSplice_left,
		int tmpBackSplice_right)
	{
		chrNameInt = tmpChrNameInt;
		backSpliceAcceptorStart_normalTranscriptStart_pos
			= tmpBackSplice_left;
		backSpliceDonerEnd_normalTranscriptEnd_pos
			= tmpBackSplice_right;
	}

	bool generateInnerSJvec(vector< pair<int,int> >& toAddInnerSJsitePosPairVec)
	{
		vector<int> tmpToAddSJstartPosVec;
		vector<int> tmpToAddSJendPosVec;
		for(int tmp = 0; tmp < toAddInnerSJsitePosPairVec.size(); tmp++)
		{
			tmpToAddSJstartPosVec.push_back(toAddInnerSJsitePosPairVec[tmp].first);
			tmpToAddSJendPosVec.push_back(toAddInnerSJsitePosPairVec[tmp].second);
		}
		sort(tmpToAddSJstartPosVec.begin(), tmpToAddSJstartPosVec.end());
		sort(tmpToAddSJendPosVec.begin(), tmpToAddSJendPosVec.end());
		for(int tmp = 0; tmp < tmpToAddSJstartPosVec.size(); tmp++)
		{
			int tmpToAddSJstartPos = tmpToAddSJstartPosVec[tmp];
			int tmpToAddSJendPos = tmpToAddSJendPosVec[tmp];
			bool tmpToAddSJsitePair_exists = false;
			for(int tmp2 = 0; tmp2 < toAddInnerSJsitePosPairVec.size(); tmp2++)
			{
				int toAddInnerSJsitePosPair_start = toAddInnerSJsitePosPairVec[tmp2].first;
				int toAddInnerSJsitePosPair_end = toAddInnerSJsitePosPairVec[tmp2].second;
				if((tmpToAddSJstartPos == toAddInnerSJsitePosPair_start)&&
					(tmpToAddSJendPos == toAddInnerSJsitePosPair_end))
				{
					tmpToAddSJsitePair_exists = true;
					break;
				}
			}
			if(!tmpToAddSJsitePair_exists)
				return false;
			innerSJsitePosPairVec.push_back(pair<int,int>(tmpToAddSJstartPos, tmpToAddSJendPos));
		}
	}

	void generateExonVec()
	{
		int exonNum = innerSJsitePosPairVec.size() + 1;
		if(exonNum == 1)
		{
			// the only exon			
			exonPosPairVec.push_back(pair<int,int>(
				backSpliceAcceptorStart_normalTranscriptStart_pos,
				backSpliceDonerEnd_normalTranscriptEnd_pos));
		}
		else if(exonNum == 2)
		{
			// 1st exon
			exonPosPairVec.push_back(pair<int,int>(
				backSpliceAcceptorStart_normalTranscriptStart_pos,
				innerSJsitePosPairVec[0].first));
			// 2nd exon
			exonPosPairVec.push_back(pair<int,int>(
				innerSJsitePosPairVec[0].second,
				backSpliceDonerEnd_normalTranscriptEnd_pos));
		}
		else
		{
			// 1st exon
			exonPosPairVec.push_back(pair<int,int>(
				backSpliceAcceptorStart_normalTranscriptStart_pos,
				innerSJsitePosPairVec[0].first));
			// other exons
			for(int tmpSJindex = 0; tmpSJindex <= innerSJsitePosPairVec.size()-2; tmpSJindex++)
			{
				int tmpExonStartPos = innerSJsitePosPairVec[tmpSJindex].second;
				int tmpExonEndPos = innerSJsitePosPairVec[tmpSJindex+1].first;
				exonPosPairVec.push_back(pair<int,int>(
					tmpExonStartPos, tmpExonEndPos));
			}
			// last exon
			exonPosPairVec.push_back(pair<int,int>(
				innerSJsitePosPairVec[innerSJsitePosPairVec.size()-1].second,
				backSpliceDonerEnd_normalTranscriptEnd_pos));
		}		
	}

	string returnCirRNAtranscriptReadCountStr(Index_Info* indexInfo)
	{
		string tmpStr;
		string chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string backSpliceAcceptorStart_normalTranscriptStart_pos_str
			= int_to_str(backSpliceAcceptorStart_normalTranscriptStart_pos);
		string backSpliceDonerEnd_normalTranscriptEnd_pos_str
			= int_to_str(backSpliceDonerEnd_normalTranscriptEnd_pos);		
		string readCountStr = int_to_str(readCount);
		tmpStr = chrNameStr + "_" + backSpliceAcceptorStart_normalTranscriptStart_pos_str
			+ "_" + backSpliceDonerEnd_normalTranscriptEnd_pos_str + "\t"
			+ readCountStr;
		return tmpStr;
	}

	string returnCirRNAtranscriptInfoStr(Index_Info* indexInfo)
	{
		string tmpStr;
		string chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		int exonNum = innerSJsitePosPairVec.size() + 1;
		string exonNumStr = int_to_str(exonNum);
		string backSpliceAcceptorStart_normalTranscriptStart_pos_str
			= int_to_str(backSpliceAcceptorStart_normalTranscriptStart_pos);
		string backSpliceDonerEnd_normalTranscriptEnd_pos_str
			= int_to_str(backSpliceDonerEnd_normalTranscriptEnd_pos);
		tmpStr = chrNameStr + "\t" + backSpliceAcceptorStart_normalTranscriptStart_pos_str + "\t"
			+ backSpliceDonerEnd_normalTranscriptEnd_pos_str + "\t" 
			+ int_to_str(readCount) + "\t"
			+ exonNumStr + "\t";
		if(exonNum == 1)
		{
			// the only exon			
			tmpStr = tmpStr + int_to_str(backSpliceAcceptorStart_normalTranscriptStart_pos)
				+ "," + int_to_str(backSpliceDonerEnd_normalTranscriptEnd_pos);
		}
		else if(exonNum == 2)
		{
			// 1st exon
			tmpStr = tmpStr + int_to_str(backSpliceAcceptorStart_normalTranscriptStart_pos)
				+ "," + int_to_str(innerSJsitePosPairVec[0].first) + "\t";
			// 2nd exon
			tmpStr = tmpStr	+ int_to_str(innerSJsitePosPairVec[0].second) + ","
				+ int_to_str(backSpliceDonerEnd_normalTranscriptEnd_pos);
		}
		else
		{
			// 1st exon
			tmpStr = tmpStr + int_to_str(backSpliceAcceptorStart_normalTranscriptStart_pos)
				+ "," + int_to_str(innerSJsitePosPairVec[0].first) + "\t";
			// other exons
			for(int tmpSJindex = 0; tmpSJindex <= innerSJsitePosPairVec.size()-2; tmpSJindex++)
			{
				int tmpExonStartPos = innerSJsitePosPairVec[tmpSJindex].second;
				int tmpExonEndPos = innerSJsitePosPairVec[tmpSJindex+1].first;
				tmpStr = tmpStr + int_to_str(tmpExonStartPos) + ","
					+ int_to_str(tmpExonEndPos) + "\t";
			}
			// last exon
			tmpStr = tmpStr	+ int_to_str(innerSJsitePosPairVec[innerSJsitePosPairVec.size()-1].second) + ","
				+ int_to_str(backSpliceDonerEnd_normalTranscriptEnd_pos);
		}
		return tmpStr;
	}
};
#endif