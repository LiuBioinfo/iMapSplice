// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNEVENTSET_INFO_H
#define ALIGNEVENTSET_INFO_H

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
#include "alignEvent_info.h"

using namespace std;

typedef map< pair<int,int>, int> SJmap; // <donerEnd, acceptorStart>, index_SJmapIndexVec

class EventEntry_Info
{
private:
	string chrName;
	//int chrPos;
	int startChrPos;
	int endChrPos;
	int supportNum;

	vector< vector<Jump_Code> > followingJumpCodeVec;
	vector< vector< pair<int, pair<char,char> > > > followingJumpCodeMismatchVec;

	//vector<>
public:
	EventEntry_Info()
	{

	}

	int returnFollowingJumpCodeMismatchPos(int index_level_1, int index_level_2)
	{
		return (followingJumpCodeMismatchVec[index_level_1])[index_level_2].first;
	}

	char returnFollowingJumpCodeMismatchCharInRead(int index_level_1, int index_level_2)
	{
		return ((followingJumpCodeMismatchVec[index_level_1])[index_level_2].second).first;
	}	

	char returnFollowingJumpCodeMismatchQualityInRead(int index_level_1, int index_level_2)
	{
		return ((followingJumpCodeMismatchVec[index_level_1])[index_level_2].second).second;
	}	

	void initiateEventEntryInfo(const string& chrNameStr, int startChrPosInt, int endChrPosInt, 
		vector<Jump_Code>& tmpFollowingJumpCodeVec, vector< pair<int, pair<char,char> > >& tmpSJfollowingSegsMismatchPosVec,
		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool)
	{
		chrName = chrNameStr;
		startChrPos = startChrPosInt;
		endChrPos = endChrPosInt;
		supportNum = 1;
		if(recordFollowingJumpCode_bool)
			followingJumpCodeVec.push_back(tmpFollowingJumpCodeVec);
		if(recordFollowingSegsMismatch_bool)
			followingJumpCodeMismatchVec.push_back(tmpSJfollowingSegsMismatchPosVec);
	}

	void duplicateSJfound(vector<Jump_Code>& tmpFollowingJumpCodeVec, 
		vector< pair<int, pair<char,char> > >& tmpSJfollowingSegsMismatchPosVec,
		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool)
	{
		supportNum++;
		if(recordFollowingJumpCode_bool)
			followingJumpCodeVec.push_back(tmpFollowingJumpCodeVec);
		if(recordFollowingSegsMismatch_bool)
			followingJumpCodeMismatchVec.push_back(tmpSJfollowingSegsMismatchPosVec);
	}

	string returnEventEntryStr(const string& junc_name, bool output_followingJumpCode_bool, bool output_followingJumpCodeMismatch_bool)
	{
		string str;
		str = chrName + "\t" + int_to_str(startChrPos) + "\t" + int_to_str(endChrPos)
			+ "\t" + "JUNC_" + junc_name + "\t" + int_to_str(supportNum) 
			+ "\t" + "*" + "\t";
		if(output_followingJumpCode_bool)
		{
			for(int tmp = 0; tmp < followingJumpCodeVec.size(); tmp++)
			{
				for(int tmp2 = 0; tmp2 < followingJumpCodeVec[tmp].size(); tmp2 ++)
				{
					string tmpJumpCodeType = (followingJumpCodeVec[tmp])[tmp2].type;
					int tmpJumpCodeLen = (followingJumpCodeVec[tmp])[tmp2].len;
					str = str + int_to_str(tmpJumpCodeLen) + tmpJumpCodeType; 
				}
				/*if(output_followingJumpCodeMismatch_bool)
				{
					str += "-";
					for(int tmpMismatchPosIndex = 0; tmpMismatchPosIndex < followingJumpCodeMismatchVec.size(); tmpMismatchPosIndex++)
					{
						str = str + int_to_str(this->returnFollowingJumpCodeMismatchPos(tmp, tmpMismatchPosIndex)) + "/";
					}
					str += "-";
					for(int tmpMismatchPosIndex = 0; tmpMismatchPosIndex < followingJumpCodeMismatchVec.size(); tmpMismatchPosIndex++)
					{
						stringstream tmpCharStrStream;
						tmpCharStrStream << this->returnFollowingJumpCodeMismatchCharInRead(tmp, tmpMismatchPosIndex);
						string tmpCharStr;
						tmpCharStrStream >> tmpCharStr;
						str = str + tmpCharStr + "/";
					}
					str += "-";
					for(int tmpMismatchPosIndex = 0; tmpMismatchPosIndex < followingJumpCodeMismatchVec.size(); tmpMismatchPosIndex++)
					{
						stringstream tmpCharStrStream;
						tmpCharStrStream << this->returnFollowingJumpCodeMismatchQualityInRead(tmp, tmpMismatchPosIndex);
						string tmpCharStr;
						tmpCharStrStream >> tmpCharStr;
						str = str + tmpCharStr + "/";
					}
				}*/
				str += ",";
			}
		}
		return str;
	}
};

class AlignEventSet_Info
{
private:
	vector<SJmap> SJmapVec;
	vector<EventEntry_Info> eventInfoVec;
public:
	AlignEventSet_Info(int chrNum)
	{
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			SJmap newSJmap;
			SJmapVec.push_back(newSJmap);
		}
	}

	void insertEventInfo(AlignEvent_Info* alignEventInfo, Index_Info* indexInfo, 
		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool, 
		bool fasta_fastq_bool)
	{
		this->insertEventInfo_SJ(alignEventInfo, indexInfo, 
			recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool,
			fasta_fastq_bool);
	}

	void insertEventInfo_SJ(AlignEvent_Info* alignEventInfo, Index_Info* indexInfo, 
		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool,
		bool fasta_fastq_bool)
	{
		string tmpSJchr = alignEventInfo->returnMapChrStr();
		vector< pair<int,int> > tmpSJposVec;
		vector< vector<Jump_Code> > tmpSJfollowingJumpCodeVec;
		vector< vector < pair<int, pair<char,char> > > > tmpSJfollowingSegsMismatchPosVec; 

 		alignEventInfo->generateSJposVec(tmpSJposVec, tmpSJfollowingJumpCodeVec, tmpSJfollowingSegsMismatchPosVec,
 			recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool, fasta_fastq_bool);
 		for(int tmp = 0; tmp < tmpSJposVec.size(); tmp++)
 		{
 			int tmpDonerEndPos = tmpSJposVec[tmp].first;
 			int tmpAcceptorStartPos = tmpSJposVec[tmp].second;
 			this->insertSJ(tmpSJchr, tmpDonerEndPos, tmpAcceptorStartPos, indexInfo, tmpSJfollowingJumpCodeVec[tmp],
 				tmpSJfollowingSegsMismatchPosVec[tmp],
 				recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool, fasta_fastq_bool);
 		}
 	}

 	void insertSJ(const string& chrName, int donerEndPos, int acceptorStartPos, 
 		Index_Info* indexInfo, vector<Jump_Code>& tmpFollowingJumpCodeVec,
 		vector< pair<int, pair<char,char> > >& tmpSJfollowingSegsMismatchPosVec,
 		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool, bool fasta_fastq_bool)
 	{
 		int chrNameInt = indexInfo->convertStringToInt(chrName);
 		if(chrNameInt >= 0)
 		{
 			SJmap::iterator tmpMapIter = SJmapVec[chrNameInt].find(pair<int,int>(donerEndPos, acceptorStartPos));
 			if(tmpMapIter == SJmapVec[chrNameInt].end()) // new SJ
 			{
 				int tmpEventEntryVecSize = eventInfoVec.size();
 				SJmapVec[chrNameInt].insert(pair< pair<int,int>, int >
 					(pair<int,int>(donerEndPos, acceptorStartPos) , tmpEventEntryVecSize));
 				EventEntry_Info newEventEntryInfo;
 				newEventEntryInfo.initiateEventEntryInfo(chrName, donerEndPos, acceptorStartPos, tmpFollowingJumpCodeVec,
 					tmpSJfollowingSegsMismatchPosVec, recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool);
 				eventInfoVec.push_back(newEventEntryInfo);
 			}
 			else // found SJ
 			{
 				int index_eventInfoVec = tmpMapIter->second;
 				eventInfoVec[index_eventInfoVec].duplicateSJfound(tmpFollowingJumpCodeVec, tmpSJfollowingSegsMismatchPosVec,
 					recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool);
 			}
 		}
 	}

 	void outputEventSet(ofstream& output_ofs, bool output_followingJumpCode_bool, bool output_followingJumpCodeMismatch_bool)
 	{
 		for(int tmp = 0; tmp < eventInfoVec.size(); tmp++)
 		{	
 			string tmpEventInfoStr = eventInfoVec[tmp].returnEventEntryStr(int_to_str(tmp+1), 
 				output_followingJumpCode_bool, output_followingJumpCodeMismatch_bool);
 			output_ofs << tmpEventInfoStr << endl;
 		}
 	}
};

#endif