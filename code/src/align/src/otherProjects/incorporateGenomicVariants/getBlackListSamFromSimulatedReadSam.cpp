// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
#include "./general/SNPhash_info.h"

using namespace std;

int getEndLocInReadOfSpecificJumpCode(
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	if(jumpCodeIndex < 0)
		return 0;
	int tmpEndLocInRead = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "M")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "I")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "D")
		{
		}
		else if(tmpJumpCodeType == "N")
		{
		}
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return tmpEndLocInRead;
}	

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int tmpEndPos = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{}
		else if(tmpJumpCodeType == "M")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
		{}
		else if(tmpJumpCodeType == "D")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "N")
			tmpEndPos += tmpJumpCodeLength;
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return (tmpEndPos + startPos-1);
}

void generateExonLocInReadPosInChr(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	vector<int>& endLocVecInRead, vector<int>& endPosVecInChr, vector<int>& lenVec)
{
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp ++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
		int tmpJumpCodeIndex = tmp;
		if(tmpJumpCodeType == "M")
		{
			int tmpEndLocInRead = getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmpJumpCodeIndex);
			int tmpEndPosInChr = getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec, tmpJumpCodeIndex);
			endLocVecInRead.push_back(tmpEndLocInRead);
			endPosVecInChr.push_back(tmpEndPosInChr);
			lenVec.push_back(tmpJumpCodeLength);
		}
	}
}

void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpJumpCodeLength;
	string tmpJumpCodeType;

	int jumpCodeStartPosInCigarStr = 0;
	int jumpCodeEndPosInCigarStr;
		
	string candidateJumpCodeType = "SMNIDX";
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

bool returnPosInChrWithLocInRead(int startPos, string& cigarString, int tmpLocInRead,
	int& correspondingPosInChr)
{
	vector<Jump_Code> cigarStringJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);

	vector<int> endLocVecInRead;
	vector<int> endPosVecInChr;
	vector<int> lenVec;
	generateExonLocInReadPosInChr(startPos, cigarStringJumpCodeVec, endLocVecInRead, endPosVecInChr, lenVec);
	int exonNum = lenVec.size();
	for(int tmp = 0; tmp < exonNum; tmp++)
	{
		int startLocInRead_exon = endLocVecInRead[tmp] - lenVec[tmp] + 1;
		int endLocInRead_exon = endLocVecInRead[tmp];
		if((tmpLocInRead >= startLocInRead_exon)&&(tmpLocInRead <= endLocInRead_exon))
		{
			int tmpDistanceToEndInExon = endLocInRead_exon - tmpLocInRead;
			correspondingPosInChr = endPosVecInChr[tmp] - tmpDistanceToEndInExon;
			return true;
		}
	}
	return false;
}

bool getForOrRcmFromFlag(int flagInt)
{
	if(flagInt & 0x10)
		return false;
	else
		return true;
}

bool extractGTandBlackListSAMfromStr_different(string& tmpSamStr, string& chrName_gt, 
	int& chrPos_gt, string& cigarString_gt, int& SNPloc_gt, string& chrName_sam, 
	int& chrPos_sam, string& cigarString_sam, bool& forOrRcm_sam, int& correspondingPosInChr)
{
	int tabLoc_1 = tmpSamStr.find("\t");
	int tabLoc_2 = tmpSamStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSamStr.find("\t", tabLoc_2 + 1);		
	int tabLoc_4 = tmpSamStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpSamStr.find("\t", tabLoc_4 + 1);
	int tabLoc_6 = tmpSamStr.find("\t", tabLoc_5 + 1);
	string gtStr = tmpSamStr.substr(0, tabLoc_1);
	string flagStr = tmpSamStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	int flagInt = atoi(flagStr.c_str());
	if(flagInt & 0x4)
		return false;
	forOrRcm_sam = getForOrRcmFromFlag(flagInt);
	chrName_sam = tmpSamStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	string chrPosStr = tmpSamStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	chrPos_sam = atoi(chrPosStr.c_str());
	cigarString_sam = tmpSamStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);

	int lineLoc_1 = gtStr.find("_");
	int lineLoc_2 = gtStr.find("_", lineLoc_1 + 1);
	int lineLoc_3 = gtStr.find("_", lineLoc_2 + 1);
	chrName_gt = gtStr.substr(0, lineLoc_1);
	string chrPosStr_gt = gtStr.substr(lineLoc_1 + 1, lineLoc_2 - lineLoc_1 - 1);
	chrPos_gt = atoi(chrPosStr_gt.c_str());
	cigarString_gt = gtStr.substr(lineLoc_2 + 1, lineLoc_3 - lineLoc_2 - 1);
	string SNPlocStr_gt = gtStr.substr(lineLoc_3 + 1);
	SNPloc_gt = atoi(SNPlocStr_gt.c_str());

	bool returnPosInChrWithLocInRead_bool = returnPosInChrWithLocInRead(
		chrPos_sam, cigarString_sam, SNPloc_gt, correspondingPosInChr);
	if(!returnPosInChrWithLocInRead_bool)
		return false;

	if((chrName_gt == chrName_sam)&&(chrPos_gt == chrPos_sam)&&(cigarString_gt == cigarString_sam))
		return false;
	else
		return true;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSAMforSimulatedRead outputCandiBlackListSAMpath" << endl;
		exit(1);
	}
	string inputSAMforSimulatedRead = argv[1];
	string outputCandiBlackListSAMpath = argv[2];
	ifstream SAM_ifs(inputSAMforSimulatedRead.c_str());
	ofstream candiBlackListSAM_ofs(outputCandiBlackListSAMpath.c_str());
	while(!SAM_ifs.eof())
	{
		string tmpSamStr;
		getline(SAM_ifs, tmpSamStr);
		if(tmpSamStr == "")
			break;
		if(tmpSamStr.substr(0,1) == "@")
			continue;
		string chrName_gt;
		int chrPos_gt;
		string cigarString_gt;
		int SNPloc_gt;
		string chrName_sam;
		int chrPos_sam;
		string cigarString_sam;
		bool forOrRcm_sam;
		int correspondingPosInChr;
		bool candiBlackListSAMfoundBool = extractGTandBlackListSAMfromStr_different(
			tmpSamStr, chrName_gt, chrPos_gt, cigarString_gt, SNPloc_gt,
			chrName_sam, chrPos_sam, cigarString_sam, forOrRcm_sam, correspondingPosInChr);
		if(candiBlackListSAMfoundBool)
		{
			string forOrRcmStr_sam = "+";
			if(!forOrRcm_sam)
				forOrRcmStr_sam = "-";
			candiBlackListSAM_ofs << chrName_gt << "\t" << chrPos_gt << "\t" << cigarString_gt << "\t"
				<< chrName_sam << "\t" << chrPos_sam << "\t" << cigarString_sam << "\t"
				<< forOrRcmStr_sam << "\t" << SNPloc_gt << "\t" << correspondingPosInChr << endl;
		}
	}
	candiBlackListSAM_ofs.close();
	SAM_ifs.close();
	return 0;
}