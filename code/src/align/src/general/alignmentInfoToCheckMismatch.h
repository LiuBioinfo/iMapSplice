// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// required: alignment file should be in SAM format

#ifndef ALIGNMENTINFOTOCHECKMISMATCH_H
#define ALIGNMENTINFOTOCHECKMISMATCH_H

#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <map>
#include "general/index_info.h"

using namespace std;

class AlignmentInfoToCheckMismatch_Info
{
private:
	string mapChrName;
	int mapChrPos;
	vector<Jump_Code> cigarStringJumpCodeVec;

public:
	AlignmentInfoToCheckMismatch_Info
	{}

	void getSJposVecFromSAM_InsertIntoSJmap(const string& samStr, Index_Info* indexInfo)
	{
		// get SJposVec from SAM string
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 10; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}

		mapChrName = samFieldVec[2];

		string mapChrPosStr = samFieldVec[3];
		mapChrPos = atoi(mapChrPosStr.c_str());
		
		string cigarString = samFieldVec[5];		
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		
		readSeq = samFieldVec[9];
	}		

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNID";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	int getMismatchBesideSJ()
	{
		vector<int> SJposInRead

		string mappedSubstrInRead = getMappedSubstrInRead()
		string mappedSubstrInChr = getMappedSubstrInChr()
	
		set<int> mismatchPosInReadSet;
		int mappedBasesNum = mappedSubstrInRead.length();
		for(int tmp = 0; tmp < mappedBasesNum; tmp++)
		{
			if(mappedSubstrInRead.at(tmp) != mappedSubstrInChr.at(tmp))
				mismatchPosInReadSet.insert(tmp);
		}

		

	}
};

#endif