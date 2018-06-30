// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SMALLEXONHASH_INFO_H
#define SMALLEXONHASH_INFO_H

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

#include "splice_info.h"

#define SMALL_EXON_LEN_MAX 15

using namespace std;

typedef map<int, int> SmallExonLenMap; // <len, index>
typedef map<int, SmallExonLenMap> SmallExonChrMapPosMap;

class SmallExon_Info
{
private:
	int chrNameInt;
	int chrMapPos;
	int len;
	int intronSizeLast;
	int intronSizeNext;
	int supportNum;
public:
	SmallExon_Info()
	{}

	void initiate_smallExonInfo(int tmpChrNameInt,
		int tmpChrMapPos, int tmpLen, 
		int tmpIntronSizeLast, int tmpIntronSizeNext)
	{
		chrNameInt = tmpChrNameInt;
		chrMapPos = tmpChrMapPos;
		len = tmpLen;
		intronSizeLast = tmpIntronSizeLast;
		intronSizeNext = tmpIntronSizeNext;
		supportNum = 1;
	}

	string returnSmallExonInfoStr(Index_Info* indexInfo, int tmpSmallExonIndex)
	{
		string str;
		str = indexInfo->returnChrNameStr(chrNameInt) + "\t" + int_to_str(chrMapPos)
			+ "\t" + int_to_str(len) 
			+ "\tSmallExon_" + int_to_str(tmpSmallExonIndex)
			+ "\t" + int_to_str(supportNum) 
			+ "\t" + int_to_str(intronSizeLast)
			+ "\t" + int_to_str(intronSizeNext);
		return str;
	}

	void updateWithFoundNewAlignSupportSameSmallExon()
	{
		supportNum ++;
	}
};

class SmallExonHash_Info
{
private:
	vector<SmallExonChrMapPosMap> smallExonChrMapPosMapVec;
	vector< SmallExon_Info > smallExonInfoVec; // vec< pair< start_pos, exon_len > >

	int currentSmallExonInfoVecIndex;
public:
	SmallExonHash_Info()
	{
		currentSmallExonInfoVecIndex = 0;
	}

	void initiateSmallExonHashInfo(int chrNum)
	{
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			SmallExonChrMapPosMap newSmallExonChrMapPosMap;
			smallExonChrMapPosMapVec.push_back(newSmallExonChrMapPosMap);
		}
		currentSmallExonInfoVecIndex = 0;
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;	
	}

	void insertSmallExonFromAlignmentFileVec(vector<string>& alignmentFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			//cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertSmallExonFromAlignmentFile(alignmentFileVec[tmp], indexInfo);
		}
	}

	void insertSmallExonFromAlignmentFile(string& alignmentFile, Index_Info* indexInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			if(sam_ifs.eof())
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getSmallExonInfoFromSAM_InsertIntoSmallExonInfoHash(tmpAlignStr, indexInfo);				
		}
		sam_ifs.close();
	}

	void getSmallExonInfoFromSAM_InsertIntoSmallExonInfoHash(const string& samStr, Index_Info* indexInfo)
	{
		// get SJposVec from SAM string
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];

		//cout << "read: " << samFieldVec[0] << endl;
		  // cout << endl <<  "mapChrNameStr: " << mapChrNameStr << endl;
		  // cout << "mapChrPos: " << mapChrPos << endl;
		  // cout << "cigarString: " << cigarString << endl;

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		  //for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		  //{
		  //	cout << "jumpCode: " << cigarStringJumpCodeVec[tmp].toString() << endl;
		  //}
		 //cout << "smallExon ... " << endl;
		vector< pair<int,int> > tmpSmallExonInfoPairVec; // <smallExonChrMapPos, smallExonLen>
		vector< pair<int,int> > tmpBesideSmallExonIntronSizePairVec; // <lastIntronSize, nextIntronSize>   

		this->generateSmallExonInfoPairVec(mapChrPos, 
			tmpSmallExonInfoPairVec,
			tmpBesideSmallExonIntronSizePairVec,
			cigarStringJumpCodeVec);
		//cout << "start to insert small exon into Hash" << endl;
		this->insert2SmallExonHashWithFoundSmallExonVec(mapChrNameInt, tmpSmallExonInfoPairVec,
			tmpBesideSmallExonIntronSizePairVec, indexInfo);
	}

	void insert2SmallExonHashWithFoundSmallExonVec(int mapChrNameInt, 
		vector< pair<int,int> >& tmpSmallExonInfoPairVec,
		vector< pair<int,int> >& tmpBesideSmallExonIntronSizePairVec,
		Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpSmallExonInfoPairVec.size(); tmp++)
		{
			//cout << "tmpSj: " << tmp+ 1 << " start: " << tmpSJposPairVec[tmp].first << " end: " << tmpSJposPairVec[tmp].second << endl;
			this->insert2SmallExonHashWithFoundSmallExon(mapChrNameInt, 
				tmpSmallExonInfoPairVec[tmp].first, 
				tmpSmallExonInfoPairVec[tmp].second, 
				tmpBesideSmallExonIntronSizePairVec[tmp].first, 
				tmpBesideSmallExonIntronSizePairVec[tmp].second, 
				indexInfo);
		}
	}

	void insert2SmallExonHashWithFoundSmallExon(int mapChrNameInt, int tmpSmallExonChrMapPos, 
		int tmpSmallExonLen, int tmpLastIntronSize, int tmpNextIntronSize, Index_Info* indexInfo)
	{
		SmallExonChrMapPosMap::iterator smallExonChrMapPosMapIter 
			= smallExonChrMapPosMapVec[mapChrNameInt].find(tmpSmallExonChrMapPos);
		if( smallExonChrMapPosMapIter == (smallExonChrMapPosMapVec[mapChrNameInt]).end() ) // new chrMapPos to insert
		{
			//cout << "not found" << endl;
			SmallExon_Info newSmallExonInfo;
			newSmallExonInfo.initiate_smallExonInfo(mapChrNameInt, tmpSmallExonChrMapPos, 
				tmpSmallExonLen, tmpLastIntronSize, tmpNextIntronSize);
			SmallExonLenMap newSmallExonLenMap;
			newSmallExonLenMap.insert(pair<int, int>(
				tmpSmallExonLen, currentSmallExonInfoVecIndex));

			smallExonChrMapPosMapVec[mapChrNameInt].insert(
				pair< int, SmallExonLenMap> (
					tmpSmallExonChrMapPos, newSmallExonLenMap));
			smallExonInfoVec.push_back(newSmallExonInfo);
			currentSmallExonInfoVecIndex ++;
		}
		else // old chrMapPos found,
		{
			SmallExonLenMap::iterator smallExonLenMapIter 
				= (smallExonChrMapPosMapIter->second).find(tmpSmallExonLen);
			//cout << "found" << endl;
			if(smallExonLenMapIter == (smallExonChrMapPosMapIter->second).end()) // new exon length found, insert it
			{
				SmallExon_Info newSmallExonInfo;
				newSmallExonInfo.initiate_smallExonInfo(mapChrNameInt, tmpSmallExonChrMapPos, 
					tmpSmallExonLen, tmpLastIntronSize, tmpNextIntronSize);				
				(smallExonChrMapPosMapIter->second).insert(pair<int, int>(
					tmpSmallExonLen, currentSmallExonInfoVecIndex));
				smallExonInfoVec.push_back(newSmallExonInfo);
				currentSmallExonInfoVecIndex ++;
			}
			else
			{	
				// fix me, for now only update support number
				smallExonInfoVec[smallExonLenMapIter->second].updateWithFoundNewAlignSupportSameSmallExon();
			}
		}
	}

	void generateSmallExonInfoPairVec(int mapChrPos, 
		vector< pair<int,int> >& tmpSmallExonInfoPairVec, 
		vector< pair<int,int> >& tmpBesideSmallExonIntronSizePairVec,
		vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		//cout << "generateSmallExonInfoPairVec starts ..." << endl;
		//cout << "cigarStringJumpCodeVec.size() - 1: " << cigarStringJumpCodeVec.size() - 1 << endl;
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		for(int tmp = 1; tmp < cigarStringJumpCodeVecSize - 1; tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;

			if((tmpJumpCodeType == "M")&&(tmpJumpCodeLength <= SMALL_EXON_LEN_MAX))
			{
				string tmpLastJumpCodeType = cigarStringJumpCodeVec[tmp-1].type;
				string tmpNextJumpCodeType = cigarStringJumpCodeVec[tmp+1].type;
				if((tmpLastJumpCodeType == "N")&&(tmpNextJumpCodeType == "N"))
				{
					int tmpLastJumpCodeEndChrMapPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, tmp-1);
					int tmpSmallExonChrMapPos = tmpLastJumpCodeEndChrMapPos + 1;
					tmpSmallExonInfoPairVec.push_back(pair<int,int>(tmpSmallExonChrMapPos, tmpJumpCodeLength));
					int tmpLastJumpCodeLen = cigarStringJumpCodeVec[tmp-1].len;
					int tmpNextJumpCodeLen = cigarStringJumpCodeVec[tmp+1].len;
					tmpBesideSmallExonIntronSizePairVec.push_back(pair<int,int>(tmpLastJumpCodeLen, tmpNextJumpCodeLen));					
				}
			}
		}		
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

	string returnSmallExonHashInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to gather Hash info:" << endl;
		string tmpStr;
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < smallExonInfoVec.size(); tmp++)
		{
			string tmpSmallExonInfoStr = smallExonInfoVec[tmp].returnSmallExonInfoStr(indexInfo, tmp+1);
			tmpStr = tmpStr + tmpSmallExonInfoStr + "\n";
		}
		//cout << "end of gathering Hash info:" << endl;
		return tmpStr;	
	}

	void outputSmallExonHashInfo(Index_Info* indexInfo, const string& outputSmallExonHashFile)
	{
		ofstream outputSmallExonHash_ofs(outputSmallExonHashFile.c_str());
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < smallExonInfoVec.size(); tmp++)
		{
			string tmpSmallExonInfoStr = smallExonInfoVec[tmp].returnSmallExonInfoStr(indexInfo, tmp+1);
			outputSmallExonHash_ofs << tmpSmallExonInfoStr << endl;
		}
		outputSmallExonHash_ofs.close();
	}
};





#endif