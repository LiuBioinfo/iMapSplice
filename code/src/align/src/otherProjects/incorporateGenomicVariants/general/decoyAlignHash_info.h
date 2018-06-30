// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef DECOYALIGNHASH_INFO_H
#define DECOYALIGNHASH_INFO_H

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

#include "decoyAlign_info.h"

using namespace std;

typedef map<int, int> DecoyAlignInfoIndexMap;
typedef map<int, set<int> > DecoyAlignAreaMap;  //( areaNO = pos/100 ) intermediate hash to directly get all possible SNPs near a certain position

class DecoyAlignHash_Info
{
private: 
	vector< DecoyAlignInfoIndexMap > decoyAlignInfoIndexMapVec;
	vector< DecoyAlignAreaMap > decoyAlignAreaMapVec;

	int areaSize;
public:
	vector<DecoyAlign_Info> decoyAlignInfoVec;

	DecoyAlignHash_Info()
	{
		areaSize = 1000;
	}

	int searchAndReturnDecoyAlignInfoVecIndex(int tmpChrNameInt, int tmpDecoyAlignPos)
	{
		DecoyAlignInfoIndexMap::iterator tmpDecoyAlignInfoIndexMapIter 
			= decoyAlignInfoIndexMapVec[tmpChrNameInt].find(tmpDecoyAlignPos);
		if(tmpDecoyAlignInfoIndexMapIter != decoyAlignInfoIndexMapVec[tmpChrNameInt].end()) // found
			return (tmpDecoyAlignInfoIndexMapIter->second);
		else
			return -1;
	}

	void initiate_decoyAlignInfoIndexMapVec_decoyAlignAreaMapVec(int chromNum)
	{
		this->initiateDecoyAlignInfoIndexMapVec(chromNum);
		this->initiateDecoyAlignAreaMapVec(chromNum);
	}

	void initiateDecoyAlignInfoIndexMapVec(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			DecoyAlignInfoIndexMap tmpDecoyAlignInfoIndexMap;
			decoyAlignInfoIndexMapVec.push_back(tmpDecoyAlignInfoIndexMap);
		}
	}

	void initiateDecoyAlignAreaMapVec(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			DecoyAlignAreaMap tmpDecoyAlignAreaMap;
			decoyAlignAreaMapVec.push_back(tmpDecoyAlignAreaMap);
		}
	}	

	void insertDecoyAlignPos2areaHash(int tmpChrNameInt, int tmpPos)
	{
		int tmpDecoyAlignAreaNO = (int)(tmpPos/areaSize);
		DecoyAlignAreaMap::iterator tmpDecoyAlignAreaMapIter;
		tmpDecoyAlignAreaMapIter = decoyAlignAreaMapVec[tmpChrNameInt].find(tmpDecoyAlignAreaNO);
		if(tmpDecoyAlignAreaMapIter == decoyAlignAreaMapVec[tmpChrNameInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(tmpPos);
			decoyAlignAreaMapVec[tmpChrNameInt].insert(pair<int, set<int> > (tmpDecoyAlignAreaNO, newPosSet));
		}
		else
		{
			if((tmpDecoyAlignAreaMapIter->second).find(tmpPos) == (tmpDecoyAlignAreaMapIter->second).end())
				(tmpDecoyAlignAreaMapIter->second).insert(tmpPos);
			else
			{}
		}
	}

	bool mappedOrNot(int tmpFlag)
	{
		if(tmpFlag & 0x4)
			return false;
		else
			return true;
	}

	bool forOrRcm(int tmpFlag)
	{
		if(tmpFlag & 0x10)
			return false;
		else
			return true;		
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
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	bool returnBaseMapPosInChr(int& toReturnBaseMapPosInChr, int tmpSNPlocInRead, int tmpSimulatedReadLength, int tmpSimulatedReadAlignment_flag, 
		int tmpSimulatedReadAlignment_chrPos, vector<Jump_Code>& tmpSimulatedReadAlignment_jumpCodeVec)
	{
		int tmpSNPlocInRead_final = tmpSNPlocInRead;
		bool forOrRcmMapBool = this->forOrRcm(tmpSimulatedReadAlignment_flag);
		if(!forOrRcmMapBool)
			tmpSNPlocInRead_final = tmpSimulatedReadLength - tmpSNPlocInRead + 1;

		int tmpSimulatedReadAlignment_jumpCodeVecSize = tmpSimulatedReadAlignment_jumpCodeVec.size();
		if(tmpSimulatedReadAlignment_jumpCodeVecSize == 1)
		{
			if((tmpSNPlocInRead_final >= 1)&&(tmpSNPlocInRead_final <= tmpSimulatedReadLength))
			{
				toReturnBaseMapPosInChr = (tmpSimulatedReadAlignment_chrPos + tmpSNPlocInRead_final - 1);
				return true;
			}
			else
			{
				cout << "error in returnBaseMapPosInChr, tmpSNPlocInRead_final: " << tmpSNPlocInRead_final << endl;
				exit(1);
				return false;
			}
		}
		else
		{
			vector<int> endLocVecInRead;
			vector<int> endMapPosVecInChr;
			for(int tmp = 0; tmp < tmpSimulatedReadAlignment_jumpCodeVecSize; tmp++)
			{
				int tmpEndLocInRead = this->getEndLocInReadOfSpecificJumpCode(tmpSimulatedReadAlignment_jumpCodeVec, tmp);
				int tmpEndPosInChr = this->getEndPosOfSpecificJumpCode(tmpSimulatedReadAlignment_chrPos,
					tmpSimulatedReadAlignment_jumpCodeVec, tmp);
				endLocVecInRead.push_back(tmpEndLocInRead);
				endMapPosVecInChr.push_back(tmpEndPosInChr);
			}

			vector< pair<int,int> > locPairVecInRead;
			locPairVecInRead.push_back(pair<int,int>(1, endLocVecInRead[0]));
			for(int tmp = 1; tmp < tmpSimulatedReadAlignment_jumpCodeVec.size(); tmp++)
			{
				int tmpSeg_startLocInRead = endLocVecInRead[tmp-1] + 1;
				int tmpSeg_endLocInRead = endLocVecInRead[tmp];
				locPairVecInRead.push_back(pair<int,int>(tmpSeg_startLocInRead, tmpSeg_endLocInRead));
			}

			for(int tmp = 0; tmp < tmpSimulatedReadAlignment_jumpCodeVecSize; tmp++)
			{
				string tmpJumpCodeType = tmpSimulatedReadAlignment_jumpCodeVec[tmp].type;
				if(tmpJumpCodeType == "M")
				{
					int tmpSeg_startLocInRead = locPairVecInRead[tmp].first;
					int tmpSeg_endLocInRead = locPairVecInRead[tmp].second;
					if((tmpSNPlocInRead_final >= tmpSeg_startLocInRead)&&(tmpSNPlocInRead_final <= tmpSeg_endLocInRead))
					{
						int tmpSeg_endPosInChr = endMapPosVecInChr[tmp];
						toReturnBaseMapPosInChr = tmpSeg_endPosInChr - (tmpSeg_endLocInRead - tmpSNPlocInRead_final);
						return true;
					}
				}
			}
			return false;
		}	
	}

	bool initiateAndAddDecoyAlignInfo_fromSimulatedTransReadAlignResultStr(string& tmpSamStr, 
		string& tmpChrName_gt, int tmpChrPos_gt, string& tmpCigarString_gt, 
		int tmpSNPlocInRead, int tmpSimulatedReadLength, int tmpSNPposInChr, 
		string& tmpTranscriptID, int tmpSNPlocInTranscript, int tmpSimulatedReadAlignment_flag,
		string& tmpSimulatedReadAlignment_chrName, int tmpSimulatedReadAlignment_chrPos, 
		string& tmpSimulatedReadAlignment_cigarString, Index_Info* indexInfo)
	{
		int tmpChrNameInt_gt = indexInfo->convertStringToInt(tmpChrName_gt);
		int tmpSimulatedReadAlignment_chrNameInt = indexInfo->convertStringToInt(tmpSimulatedReadAlignment_chrName);
		if((tmpChrNameInt_gt < 0)||(tmpSimulatedReadAlignment_chrNameInt < 0))
			return false;
		bool mappedOrNot = this->mappedOrNot(tmpSimulatedReadAlignment_flag);
		if(!mappedOrNot)
			return false;
		//vector<Jump_Code> jumpCodeVec_gt;
		vector<Jump_Code> tmpSimulatedReadAlignment_jumpCodeVec;
		//this->cigarString2jumpCodeVec(tmpCigarString_gt, jumpCodeVec_gt);
		this->cigarString2jumpCodeVec(tmpSimulatedReadAlignment_cigarString, tmpSimulatedReadAlignment_jumpCodeVec);
		
		int tmpSNPcorrespondingChrNameIntInSimulatedReadAlignment = tmpSimulatedReadAlignment_chrNameInt;
		int tmpSNPcorrespondingPosInSimulatedReadAlignment;
		bool tmpSNPcorrespondingPosInSimulatedReadAlignmentBool = this->returnBaseMapPosInChr(
			tmpSNPcorrespondingPosInSimulatedReadAlignment,
			tmpSNPlocInRead, tmpSimulatedReadLength, tmpSimulatedReadAlignment_flag, 
			tmpSimulatedReadAlignment_chrPos, tmpSimulatedReadAlignment_jumpCodeVec);
		if(!tmpSNPcorrespondingPosInSimulatedReadAlignmentBool)
			return false;

		DecoyAlignInfoIndexMap::iterator decoyAlignInfoMapIter = decoyAlignInfoIndexMapVec[tmpSimulatedReadAlignment_chrNameInt].find(
				tmpSNPcorrespondingPosInSimulatedReadAlignment);
		if(decoyAlignInfoMapIter == decoyAlignInfoIndexMapVec[tmpSimulatedReadAlignment_chrNameInt].end()) // not found
		{
			// add pos to area map
			this->insertDecoyAlignPos2areaHash(tmpSimulatedReadAlignment_chrNameInt, tmpSNPcorrespondingPosInSimulatedReadAlignment);
			// add decoyAlignInfo to DecoyAlignInfoIndexMap
			DecoyAlign_Info tmpNewDecoyAlignInfo;
			tmpNewDecoyAlignInfo.initiateAndAddNewDecoyAlignElement(tmpSamStr,
				tmpSNPcorrespondingChrNameIntInSimulatedReadAlignment, tmpSNPcorrespondingPosInSimulatedReadAlignment,
				tmpChrNameInt_gt, tmpChrPos_gt, tmpCigarString_gt, tmpSNPlocInRead, 
				tmpSimulatedReadLength, tmpSNPposInChr, tmpTranscriptID, tmpSNPlocInTranscript, tmpSimulatedReadAlignment_flag,
				tmpSimulatedReadAlignment_chrNameInt, tmpSimulatedReadAlignment_chrPos, tmpSimulatedReadAlignment_cigarString);
			decoyAlignInfoVec.push_back(tmpNewDecoyAlignInfo);
			int tmpDecoyAlignInfoIndex = decoyAlignInfoVec.size() - 1;
			decoyAlignInfoIndexMapVec[tmpSimulatedReadAlignment_chrNameInt].insert(pair<int,int>(
				tmpSNPcorrespondingPosInSimulatedReadAlignment, tmpDecoyAlignInfoIndex));
		}
		else // found, update corresponding decoyAlignInfo
		{
			// update decoyAlignInfo to DecoyAlignInfoIndexMap
			int decoyAlignInfoIndex = decoyAlignInfoMapIter->second;
			decoyAlignInfoVec[decoyAlignInfoIndex].addNewDecoyAlignElement(tmpSamStr,
				tmpChrNameInt_gt, tmpChrPos_gt, tmpCigarString_gt, tmpSNPlocInRead, 
				tmpSimulatedReadLength, tmpSNPposInChr, tmpTranscriptID, tmpSNPlocInTranscript, tmpSimulatedReadAlignment_flag,
				tmpSimulatedReadAlignment_chrNameInt, tmpSimulatedReadAlignment_chrPos, tmpSimulatedReadAlignment_cigarString);
		}
		return true;
	}

	void extractGtAndSimulatedReadAlignInfo(string& tmpSamStr, 
		string& tmpChrName_gt, int& tmpChrPos_gt, string& tmpCigarString_gt, 
		int& tmpSNPlocInRead, int& tmpSimulatedReadLength, int& tmpSNPposInChr, 
		string& tmpTranscriptID, int& tmpSNPlocInTranscript, int& tmpSimulatedReadAlignment_flag,
		string& tmpSimulatedReadAlignment_chrName, int& tmpSimulatedReadAlignment_chrPos, 
		string& tmpSimulatedReadAlignment_cigarString)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = tmpSamStr.find("\t", startLoc);
			string tmpSamField = tmpSamStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}

		// extract info of simulatedReadAlignment
		string tmpSimulatedReadAlignment_flagStr = samFieldVec[1];
		tmpSimulatedReadAlignment_flag 
			= atoi(tmpSimulatedReadAlignment_flagStr.c_str());
		tmpSimulatedReadAlignment_chrName = samFieldVec[2];
		string tmpSimulatedReadAlignment_chrPosStr = samFieldVec[3];
		tmpSimulatedReadAlignment_chrPos 
			= atoi(tmpSimulatedReadAlignment_chrPosStr.c_str());
		tmpSimulatedReadAlignment_cigarString = samFieldVec[5];

		// extract info of ground truth (correct and original alignments of reads)
		string tmpGTinfoStr = samFieldVec[0];
		vector<string> gtInfoFieldVec;
		startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpGTinfoStr.find(":", startLoc);
			string tmpGTinfoField = tmpGTinfoStr.substr(startLoc, tabLoc-startLoc);
			gtInfoFieldVec.push_back(tmpGTinfoField);
			startLoc = tabLoc + 1;
		}
		gtInfoFieldVec.push_back(tmpGTinfoStr.substr(startLoc));
		tmpChrName_gt = gtInfoFieldVec[0];
		string tmpChrPosStr_gt = gtInfoFieldVec[1];
		tmpChrPos_gt = atoi(tmpChrPosStr_gt.c_str());
		tmpCigarString_gt = gtInfoFieldVec[2];
		string tmpSNPlocInReadStr = gtInfoFieldVec[3];
		tmpSNPlocInRead = atoi(tmpSNPlocInReadStr.c_str());
		string tmpSimulatedReadLengthStr = gtInfoFieldVec[4];
		tmpSimulatedReadLength = atoi(tmpSimulatedReadLengthStr.c_str());
		string tmpSNPposInChrStr = gtInfoFieldVec[5];
		tmpSNPposInChr = atoi(tmpSNPposInChrStr.c_str());
		tmpTranscriptID = gtInfoFieldVec[6];
		string tmpSNPlocInTranscriptStr = gtInfoFieldVec[7];
		tmpSNPlocInTranscript = atoi(tmpSNPlocInTranscriptStr.c_str());
	}

	void generateDecoyAlignHash_fromSimulatedTransReadAlignResults(string inputSimulatedTransReadAlignResults, Index_Info* indexInfo)
	{
		ifstream sam_ifs(inputSimulatedTransReadAlignResults.c_str());
		while(!sam_ifs.eof())
		{
			string tmpSamStr;
			getline(sam_ifs, tmpSamStr);
			if(tmpSamStr == "")
				break;
			if(tmpSamStr.substr(0,1) == "@")
				continue;
			string tmpChrName_gt;
			int tmpChrPos_gt;
			string tmpCigarString_gt;
			int tmpSNPlocInRead;
			int tmpSimulatedReadLength;
			int tmpSNPposInChr;
			string tmpTranscriptID;
			int tmpSNPlocInTranscript;
			int tmpSimulatedReadAlignment_flag;
			string tmpSimulatedReadAlignment_chrName;
			int tmpSimulatedReadAlignment_chrPos;
			string tmpSimulatedReadAlignment_cigarString;
			this->extractGtAndSimulatedReadAlignInfo(tmpSamStr, tmpChrName_gt, 
				tmpChrPos_gt, tmpCigarString_gt, tmpSNPlocInRead, tmpSimulatedReadLength, 
				tmpSNPposInChr, tmpTranscriptID, tmpSNPlocInTranscript, 
				tmpSimulatedReadAlignment_flag, tmpSimulatedReadAlignment_chrName, 
				tmpSimulatedReadAlignment_chrPos, tmpSimulatedReadAlignment_cigarString);
			if((tmpChrName_gt == tmpSimulatedReadAlignment_chrName)
				&&(tmpChrPos_gt == tmpSimulatedReadAlignment_chrPos)
				&&(tmpCigarString_gt == tmpSimulatedReadAlignment_cigarString))
				continue;
			else
			{
				this->initiateAndAddDecoyAlignInfo_fromSimulatedTransReadAlignResultStr(tmpSamStr, 
					tmpChrName_gt, tmpChrPos_gt, tmpCigarString_gt, tmpSNPlocInRead, tmpSimulatedReadLength, 
					tmpSNPposInChr, tmpTranscriptID, tmpSNPlocInTranscript, tmpSimulatedReadAlignment_flag,
					tmpSimulatedReadAlignment_chrName, tmpSimulatedReadAlignment_chrPos, 
					tmpSimulatedReadAlignment_cigarString, indexInfo);
			}
		}
		sam_ifs.close();
	}

	void output_decoySNPcorrespondingMapPos(string& outputFilePath, Index_Info* indexInfo)
	{
		ofstream decoySNPmapPos_ofs(outputFilePath.c_str());
		int decoyAlignInfoVecSize = decoyAlignInfoVec.size();
		for(int tmp = 0; tmp < decoyAlignInfoVecSize; tmp++)
		{
			int tmpDecoyAlignInfo_SNPmapChrNameInt = decoyAlignInfoVec[tmp].return_decoyAlign_chrNameInt();
			string tmpDecoyAlignInfo_SNPmapChrName = indexInfo->returnChrNameStr(tmpDecoyAlignInfo_SNPmapChrNameInt);
			int tmpDecoyAlignInfo_SNPmapChrPos =  decoyAlignInfoVec[tmp].return_decoyAlign_chrPos();
			decoySNPmapPos_ofs << tmpDecoyAlignInfo_SNPmapChrName << "\t" << tmpDecoyAlignInfo_SNPmapChrPos;// << endl;
			string tmpCorrespondingSNPinfoStr = decoyAlignInfoVec[tmp].returnSNPinfoVecStr(indexInfo);
			decoySNPmapPos_ofs << "\t" << tmpCorrespondingSNPinfoStr << endl;
		}
		decoySNPmapPos_ofs.close();
	}

	void output_decoyAlignment(string& outputFilePath, Index_Info* indexInfo)
	{
		ofstream decoyAlignment_ofs(outputFilePath.c_str());
		int decoyAlignInfoVecSize = decoyAlignInfoVec.size();
		for(int tmp = 0; tmp < decoyAlignInfoVecSize; tmp++)
		{
			string tmpDecoyAlignmentStr = decoyAlignInfoVec[tmp].return_decoyAlignmentStr(indexInfo);
			decoyAlignment_ofs << tmpDecoyAlignmentStr << endl;
		}
		decoyAlignment_ofs.close();
	}
};

#endif