// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to count reads number within and overlap each region in a list.
//input region list format: chrName buildRefName(like hg19) regionType(like exon) startPos endPos ....... 
//output_detailed readCount format: < original string in input region list > 
//	+ "\t" + unique_within + "\t" + unique_overlap + "\t" + unique_total
//  + "\t" + multi_within + "\t" + multi_overlap + "\t" + multi_total
//  + "\t" + within_total + "\t" overlap_total + "\t" + total
//output_unique readCount format: < original string in input region list > 
// + "\t" + unique_within + "\t" + unique_overlap + "\t" + unique_total
//output_unique_multi readCount format: < original string in input region list > 
// + "\t" + unique_within + "\t" + unique_overlap + "\t" + unique_total
// + "\t" + multi_within + "\t" + multi_overlap + "\t" + multi_total + "\t" + total
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
//#include <omp.h>
#include "../general/read_block_test.h"
#include "../general/bwtmap_info.h"
#include "../general/DoubleAnchorScore.h"
#include "../general/sbndm.h"
#include "../general/splice_info.h"

using namespace std;

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool uniqueOrNot(string& IHfieldStr)
{
	int IHfieldStrLen = IHfieldStr.length();
	int IHnumInt = atoi((IHfieldStr.substr(5)).c_str());
	if(IHnumInt == 1)
		return true;
	else
		return false;
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

void getChrNamePosFromRegionStr(
	string& tmpRegion_chrName, 
	int& tmpRegion_startPos, int& tmpRegion_endPos,
	string& tmpRegionStr)
{
	vector<string> regionFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 5; tmp++)
	{
		int tabLoc = tmpRegionStr.find("\t", startLoc);
		string tmpRegionField = tmpRegionStr.substr(startLoc, tabLoc-startLoc);
		regionFieldVec.push_back(tmpRegionField);
		startLoc = tabLoc + 1;
	}	
	tmpRegion_chrName = regionFieldVec[0];
	tmpRegion_startPos = atoi((regionFieldVec[3]).c_str());
	tmpRegion_endPos = atoi((regionFieldVec[4]).c_str());
}

void getMapInfoFromSamStr(
	bool& tmpSAM_mapOrNot_bool,
	bool& tmpSAM_uniqueOrNot_bool, string& tmpSAM_chrName,
	int& tmpSAM_startMapPos, int& tmpSAM_endMapPos, 
	string& samStr)
{			
	vector<string> samFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 13; tmp++)
	{
		int tabLoc = samStr.find("\t", startLoc);
		string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
		samFieldVec.push_back(tmpSamField);
		startLoc = tabLoc + 1;
	}
	string flagStr = samFieldVec[1];
	int flagInt = atoi(flagStr.c_str());
	string mapChrPosStr = samFieldVec[3];
	string cigarString = samFieldVec[5];
	vector<Jump_Code> cigarStringJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
	string IHfieldStr = samFieldVec[12];

	tmpSAM_mapOrNot_bool = mappedOrNot(flagInt);
	tmpSAM_uniqueOrNot_bool = uniqueOrNot(IHfieldStr);
	tmpSAM_chrName = samFieldVec[2];
	tmpSAM_startMapPos = atoi(mapChrPosStr.c_str());
	tmpSAM_endMapPos = getEndPosOfSpecificJumpCode(tmpSAM_startMapPos,
		cigarStringJumpCodeVec, cigarStringJumpCodeVec.size()-1);
}

void checkSamWithOverlapCertainRegionOrNot(
	bool& withinOrNot_bool, bool& overlapOrNot_bool,
	string& tmpRegion_chrName, int tmpRegion_startPos, int tmpRegion_endPos,
	string& tmpSAM_chrName, int tmpSAM_startMapPos, int tmpSAM_endMapPos)
{
	if(tmpRegion_chrName != tmpSAM_chrName)
	{
		withinOrNot_bool = false;
		overlapOrNot_bool = false;
	}
	else
	{
		if((tmpSAM_startMapPos >= tmpRegion_endPos)||(tmpSAM_endMapPos <= tmpRegion_startPos))
		{
			withinOrNot_bool = false;
			overlapOrNot_bool = false;
		}
		else if((tmpSAM_startMapPos >= tmpRegion_startPos)&&(tmpSAM_endMapPos <= tmpRegion_endPos))
		{
			withinOrNot_bool = true;
			overlapOrNot_bool = true;
		}
		else if(((tmpSAM_startMapPos >= tmpRegion_startPos)&&(tmpSAM_startMapPos <= tmpRegion_endPos))
			||
			((tmpSAM_endMapPos >= tmpRegion_startPos)&&(tmpSAM_endMapPos <= tmpRegion_endPos)))
		{
			withinOrNot_bool = false;
			overlapOrNot_bool = true;
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputRegionListFilePath inputSAMpath outputFolder" << endl;
		exit(1);
	}
	string inputRegionListFilePath = argv[1];
	string inputSAMpath = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_log = outputFolderStr + "log.txt";
	string output_uniqueAlignmentCount = outputFolderStr + "uniqueAlignmentCount.txt";
	string output_totalAlignmentCount = outputFolderStr + "totalAlignmentCount.txt";
	string output_detailedCountInfo = outputFolderStr + "detailedCountInfo.txt";
	ofstream log_ofs(output_log.c_str());
	ofstream uniqueAlignmentCount_ofs(output_uniqueAlignmentCount.c_str());
	ofstream totalAlignmentCount_ofs(output_totalAlignmentCount.c_str());
	ofstream detailedCountInfo_ofs(output_detailedCountInfo.c_str());
	log_ofs << "inputRegionListFilePath: " << inputRegionListFilePath << endl;
	log_ofs << "inputSAMpath: " << inputSAMpath << endl;

	ifstream regionList_ifs(inputRegionListFilePath.c_str());
	while(!regionList_ifs.eof())
	{
		string tmpRegionStr;
		getline(regionList_ifs, tmpRegionStr);
		string tmpRegion_chrName;// = getChrNameFromRegionStr(tmpRegionStr);
		int tmpRegion_startPos;// = getStartPosFromRegionStr(tmpRegionStr);
		int tmpRegion_endPos;// = getEndPosFromRegionStr(tmpRegionStr);
		getChrNamePosFromRegionStr(tmpRegion_chrName, 
			tmpRegion_startPos, tmpRegion_endPos,
			tmpRegionStr);

		int tmp_unique_within_readNum = 0;
		int tmp_unique_overlap_readNum = 0;
		int tmp_multi_within_readNum = 0;
		int tmp_multi_overlap_readNum = 0;

		ifstream sam_ifs(inputSAMpath.c_str());
		while(!sam_ifs.eof())
		{
			string tmpSamStr;
			getline(sam_ifs, tmpSamStr);
			if(sam_ifs.eof())
				break;
			if(tmpSamStr.at(0) == '@')
				continue;
			bool tmpSAM_mapOrNot_bool;
			bool tmpSAM_uniqueOrNot_bool;
			string tmpSAM_chrName;
			int tmpSAM_startMapPos;
			int tmpSAM_endMapPos;
			getMapInfoFromSamStr(tmpSAM_mapOrNot_bool,
				tmpSAM_uniqueOrNot_bool, tmpSAM_chrName,
				tmpSAM_startMapPos, tmpSAM_endMapPos, tmpSamStr);
			bool withinOrNot_bool;
			bool overlapOrNot_bool;
			if(tmpSAM_mapOrNot_bool)
			{
				checkSamWithOverlapCertainRegionOrNot(
					withinOrNot_bool, overlapOrNot_bool,
					tmpRegion_chrName, tmpRegion_startPos, tmpRegion_endPos,
					tmpSAM_chrName, tmpSAM_startMapPos, tmpSAM_endMapPos);
				if(withinOrNot_bool)
				{
					if(tmpSAM_uniqueOrNot_bool)
						tmp_unique_within_readNum ++;
					else
						tmp_multi_within_readNum ++;
				}
				else if(overlapOrNot_bool)
				{
					if(tmpSAM_uniqueOrNot_bool)
						tmp_unique_overlap_readNum ++;
					else
						tmp_multi_overlap_readNum ++;
				}
				else
				{}
			}
			else
			{}
		}
		sam_ifs.close();

		int total_unique_readNum = tmp_unique_within_readNum
			+ tmp_unique_overlap_readNum;
		int total_multi_readNum = tmp_multi_within_readNum 
			+ tmp_multi_overlap_readNum;
		int total_within_readNum = tmp_unique_within_readNum
			+ tmp_multi_within_readNum;
		int total_overlap_readNum = tmp_unique_overlap_readNum
			+ tmp_multi_overlap_readNum;
		int total_readNum = tmp_unique_within_readNum 
			+ tmp_unique_overlap_readNum 
			+ tmp_multi_within_readNum 
			+ tmp_multi_overlap_readNum;
		uniqueAlignmentCount_ofs << tmpRegionStr 
			<< "\t" << tmp_unique_within_readNum
			<< "\t" << tmp_unique_overlap_readNum
			<< "\t" << total_unique_readNum << endl;
		totalAlignmentCount_ofs << tmpRegionStr
			<< "\t" << tmp_unique_within_readNum
			<< "\t" << tmp_unique_overlap_readNum
			<< "\t" << total_unique_readNum 
			<< "\t" << tmp_multi_within_readNum
			<< "\t" << tmp_multi_overlap_readNum
			<< "\t" << total_multi_readNum
			<< "\t" << total_readNum << endl;
		detailedCountInfo_ofs << tmpRegionStr
			<< "\t" << tmp_unique_within_readNum
			<< "\t" << tmp_unique_overlap_readNum
			<< "\t" << total_unique_readNum 
			<< "\t" << tmp_multi_within_readNum
			<< "\t" << tmp_multi_overlap_readNum
			<< "\t" << total_multi_readNum
			<< "\t" << total_within_readNum
			<< "\t" << total_overlap_readNum
			<< "\t" << total_readNum << endl;
	}
	regionList_ifs.eof();
	log_ofs.close();
	uniqueAlignmentCount_ofs.close();
	totalAlignmentCount_ofs.close();
	detailedCountInfo_ofs.close();
	return 0;
}

