// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef BIOEVENTALIGNINFERHASH_INFO_H
#define BIOEVENTALIGNINFERHASH_INFO_H

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

#include "alignInfer_info.h"

using namespace std;

typedef map<int, int> BioEventLenMap; // <bioEventLen, index>, bioEventLen < 0 if bioEvent == "Insertion"
typedef map<int, BioEventLenMap> BioEventAlignInferMap; //<donerEndPos, >

class BioEventAlignInferHash_Info
{
private: 
	vector<BioEventAlignInferMap> bioEventAlignInferMapVec;

	vector<BioEventAlignInfer_Info> bioEventAlignInferInfoVec;

	int currentBioEventAlignInferInfoVecIndex;
public:
	BioEventAlignInferHash_Info()
	{
		currentBioEventAlignInferInfoVecIndex = 0;
	} 

	void initiateBioEventAlignInferHashInfo(int chrNum)
	{
		// initiate bioEventAlignInferMap vec
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			BioEventAlignInferMap newBioEventAlignInferMap;
			bioEventAlignInferMapVec.push_back(newBioEventAlignInferMap);
		}
		currentAlignInferInfoVecIndex = 0;
	}

	void insertBioEventAlignInferInfoFromAlignmentFileVec(vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int maxReadBaseNumInPathStructure)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertBioEventAlignInferInfoFromAlignmentFile(
				alignmentFileVec[tmp], indexInfo, maxReadBaseNumInPathStructure);
		}
	}

	void insertBioEventAlignInferInfoFromAlignmentFile(
		string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure)
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
			this->getBioEventAlignInferInfoFromSAM_InsertIntoBioEventAlignInferHash(
				tmpAlignStr, indexInfo, maxReadBaseNumInPathStructure);				
		}
		sam_ifs.close();
	}	

	void getBioEventAlignInferInfoFromSAM_InsertIntoBioEventAlignInferHash(
		const string& samStr, Index_Info* indexInfo, int maxReadBaseNum)
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

		// cout << endl <<  "mapChrNameStr: " << mapChrNameStr << endl;
		// cout << "mapChrPos: " << mapChrPos << endl;
		// cout << "cigarString: " << cigarString << endl;

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		//if(cigarStringJumpCodeVec[0].type == "I")
		//	cout << "Alignment: " << endl << samStr << endl;
		// for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		// {
		// 	cout << "jumpCode: " << cigarStringJumpCodeVec[tmp].toString() << endl;
		// }

		vector< pair<int,int> > tmpBioEventPosLenPairVec;
		vector< string > tmpBioEventTypeVec;
		vector< int > tmpBioEventIndexVec_cigarStringJumpCodeVec;
		this->generateBioEventPosLenPairVecFromJumpCodeVec(
			mapChrPos, cigarStringJumpCodeVec, 
			tmpBioEventPosLenPairVec, 
			tmpBioEventTypeVec,
			tmpBioEventIndexVec_cigarStringJumpCodeVec);
		// for(int tmp = 0; tmp < tmpBioEventPosLenPairVec.size(); tmp++)
		// {
		// 	cout << "bioEvent: " << tmpBioEventPosLenPairVec[tmp].first 
		//		<< " - " << tmpBioEventPosLenPairVec[tmp].second 
		//		<< " " << tmpBioEventTypeVec[tmp] << endl;
		// } 

		vector< vector<Jump_Code> > tmpBioEventJumpCodeVecVec_backward;
		vector< vector<Jump_Code> > tmpBioEventJumpCodeVecVec_forward;
		this->generateBioEventJumpCodeVec_maxReadBaseNumInPathStruct(
			cigarStringJumpCodeVec, tmpBioEventIndexVec_cigarStringJumpCodeVec,
			tmpBioEventJumpCodeVecVec_backward, tmpBioEventJumpCodeVecVec_forward,
			maxReadBaseNum);

		this->insertBioEventJumpCodeVecVec_maxReadBaseNumInPathStruct(
			mapChrNameInt, 
			tmpBioEventPosLenPairVec,
			tmpBioEventTypeVec,
			tmpBioEventJumpCodeVecVec_backward, 
			tmpBioEventJumpCodeVecVec_forward, 
			indexInfo, maxReadBaseNum);
	}

	void insertBioEventJumpCodeVecVec_maxReadBaseNumInPathStruct(
			int mapChrNameInt, 
			vector< pair<int,int> >& tmpBioEventPosLenPairVec,
			vector< string >& tmpBioEventTypeVec,
			vector< vector<Jump_Code> >& tmpBioEventJumpCodeVecVec_backward, 
			vector< vector<Jump_Code> >& tmpBioEventJumpCodeVecVec_forward, 
			Index_Info* indexInfo, int maxReadBaseNum)
	{
		for(int tmp = 0; tmp < tmpBioEventPosLenPairVec.size(); tmp++)
		{
			//cout << "tmpSj: " << tmp+ 1 << " start: " << tmpSJposPairVec[tmp].first << " end: " << tmpSJposPairVec[tmp].second << endl;
			this->insertBioEventJumpCodeVec_maxReadBaseNumInPathStruct(
				mapChrNameInt, 
				tmpBioEventPosLenPairVec[tmp].first,
				tmpBioEventPosLenPairVec[tmp].second, 
				tmpBioEventTypeVec[tmp],
				tmpBioEventJumpCodeVecVec_backward[tmp], 
				tmpBioEventJumpCodeVecVec_forward[tmp], 
				indexInfo, maxReadBaseNum);
		}
	}

	void insertBioEventJumpCodeVec_maxReadBaseNumInPathStruct(
		int mapChrNameInt, 
		int tmpBioEventDonerEndPos, 
		int tmpBioEventLen,
		const string& tmpBioEventType, 
		vector<Jump_Code>& tmpBioEventJumpCodeVec_backward, 
		vector<Jump_Code>& tmpBioEventJumpCodeVec_forward, 
		Index_Info* indexInfo, int maxReadBaseNum)
	{
		int tmpFinalBioEventLen = tmpBioEventLen;
		if(tmpBioEventType == "I")
			tmpFinalBioEventLen = 0 - tmpBioEventLen;

		BioEventAlignInferMap::iterator bioEventAlignInferMapIter 
			= bioEventAlignInferMapVec[mapChrNameInt].find(tmpBioEventDonerEndPos);
		if( bioEventAlignInferMapIter  == (bioEventAlignInferMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			BioEventAlignInfer_Info tmpBioEventAlignInferInfo;
			tmpBioEventAlignInferInfo.initiateBioEventAlignInferInfo(
				mapChrNameInt, 
				tmpBioEventDonerEndPos, tmpBioEventLen, tmpBioEventType, 
				tmpBioEventJumpCodeVec_backward, tmpBioEventJumpCodeVec_forward);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			BioEventLenMap tmpBioEventLenMap;
			tmpBioEventLenMap.insert(pair<int, int>(
				tmpFinalBioEventLen, currentBioEventAlignInferInfoVecIndex));

			bioEventAlignInferMapVec[mapChrNameInt].insert(
				pair< int, BioEventLenMap> (
					tmpBioEventDonerEndPos, tmpBioEventLenMap));
			bioEventAlignInferInfoVec.push_back(tmpBioEventAlignInferInfo);
			currentBioEventAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			BioEventLenMap::iterator bioEventLenMapIter 
				= (bioEventAlignInferMapIter->second).find(tmpFinalBioEventLen);
			//cout << "found" << endl;
			if(bioEventLenMapIter == (bioEventAlignInferMapIter->second).end()) // new acceptor pos found, insert it
			{
				BioEventAlignInfer_Info tmpBioEventAlignInferInfo;
				tmpBioEventAlignInferInfo.initiateBioEventAlignInferInfo(
					mapChrNameInt, 
					tmpBioEventDonerEndPos, tmpBioEventLen, tmpBioEventType,
					tmpBioEventJumpCodeVec_backward, tmpBioEventJumpCodeVec_forward);				
				(bioEventAlignInferMapIter->second).insert(pair<int, int>(
					tmpFinalBioEventLen, currentBioEventAlignInferInfoVecIndex));
				bioEventAlignInferInfoVec.push_back(tmpAlignInferInfo);
				currentBioEventAlignInferInfoVecIndex ++;
			}
			else
			{	
				bioEventAlignInferInfoVec[bioEventLenMapIter->second].updateOldAlignInfo(
					tmpBioEventJumpCodeVec_backward, tmpBioEventJumpCodeVec_forward, 
					maxReadBaseNum);
				//alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
			}
		}
	}



	void generateBioEventJumpCodeVec_maxReadBaseNumInPathStruct(
			vector<Jump_Code>& cigarStringJumpCodeVec,
			vector< int >& tmpBioEventIndexVec_cigarStringJumpCodeVec,
			vector< vector<Jump_Code> >& tmpBioEventJumpCodeVecVec_backward,
			vector< vector<Jump_Code> >& tmpBioEventJumpCodeVecVec_forward, 
			int maxReadBaseNum)
	{
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		//cout << "cigarStringJumpCodeVecSize: " << cigarStringJumpCodeVecSize << endl;
		for(int tmp = 0; tmp < tmpBioEventIndexVec_cigarStringJumpCodeVec.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpBioEventIndex_cigarStringJumpCodeVec = tmpBioEventIndexVec_cigarStringJumpCodeVec[tmp];
			//cout << "tmpSJindex_cigarStringJumpCodeVec " << tmpSJindex_cigarStringJumpCodeVec << endl;
			vector<Jump_Code> tmpBioEventJumpCodeVec_backward;
			vector<Jump_Code> tmpBioEventJumpCodeVec_forward;

			this->generateBioEventJumpCodeVec_backward_maxReadBaseNumInPathStruct(
				cigarStringJumpCodeVec, tmpBioEventIndex_cigarStringJumpCodeVec,
				tmpBioEventJumpCodeVec_backward, maxReadBaseNum);
			this->generateBioEventJumpCodeVec_forward_maxReadBaseNumInPathStruct(
				cigarStringJumpCodeVec, tmpBioEventIndex_cigarStringJumpCodeVec,
				tmpBioEventJumpCodeVec_forward, maxReadBaseNum);

			tmpBioEventJumpCodeVecVec_backward.push_back(tmpBioEventJumpCodeVec_backward);
			tmpBioEventJumpCodeVecVec_forward.push_back(tmpBioEventJumpCodeVec_forward);
		}
	}	

	void generateBioEventJumpCodeVec_backward_maxReadBaseNumInPathStruct(
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		int tmpBioEventIndex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpBioEventJumpCodeVec_backward, 
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		for(int tmpJumpCodeIndex = tmpBioEventIndex_cigarStringJumpCodeVec-1; 
			tmpJumpCodeIndex >= 0; tmpJumpCodeIndex--)
		{
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpBioEventJumpCodeVec_backward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpBioEventJumpCodeVec_backward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpBioEventJumpCodeVec_backward.push_back(lastJumpCode);
				return;
			}
		}
	}

	void generateBioEventJumpCodeVec_forward_maxReadBaseNumInPathStruct(
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		int tmpBioEventIndex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpBioEventJumpCodeVec_forward,
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		for(int tmpJumpCodeIndex = tmpBioEventIndex_cigarStringJumpCodeVec+1; 
			tmpJumpCodeIndex < cigarStringJumpCodeVecSize; tmpJumpCodeIndex++)
		{
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpBioEventJumpCodeVec_forward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpBioEventJumpCodeVec_forward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpBioEventJumpCodeVec_forward.push_back(lastJumpCode);
				return;
			}
		}
	}



	void generateBioEventPosLenPairVecFromJumpCodeVec(int tmpMapChrPos, 
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector< pair<int,int> >& tmpBioEventPosLenPairVec,
		vector< string >& tmpBioEventTypeVec,
		vector<int>& tmpBioEventIndexVec_cigarStringJumpCodeVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
			if((tmpJumpCodeType == "N")||(tmpJumpCodeType == "I")||(tmpJumpCodeType == "D"))
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(
					tmpMapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
				//int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
				tmpBioEventPosLenPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpJumpCodeLength));
				tmpBioEventTypeVec.push_back(tmpJumpCodeType);
				tmpBioEventIndexVec_cigarStringJumpCodeVec.push_back(tmp);
			}
		}	
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
};



