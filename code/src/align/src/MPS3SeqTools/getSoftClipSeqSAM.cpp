// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to check mappedLength's distribution of primary alignments
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

bool primaryOrNot(int tmpFlag)
{
	if(tmpFlag & 0x100)
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

int getAlignmentEndPos(int startMapPos, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpEndPos = startMapPos - 1;
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
		if((tmpJumpCodeType == "M")||(tmpJumpCodeType == "N")||(tmpJumpCodeType == "D"))
		{
			tmpEndPos += tmpJumpCodeLength;
		}
	}
	return tmpEndPos;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSAMfile outputSoftClipSeqAlignmentFile " << endl;
		exit(1);
	}

	string inputSAMfileStr = argv[1];
	string outputSoftClipSeqAlignmentFileStr = argv[2];
	ifstream inputSAM_ifs(inputSAMfileStr.c_str());
	ofstream outputSoftClipSAM_ofs(outputSoftClipSeqAlignmentFileStr.c_str());

	while(!(inputSAM_ifs.eof()))
	{
		string samStr;
		getline(inputSAM_ifs, samStr);
		if(inputSAM_ifs.eof())
		 	break;
		if(samStr.at(0) == '@')
			continue;

		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 10; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string readName = samFieldVec[0];
		string readName_softClip_head = readName + ".softClip.head";
		string readName_softClip_tail = readName + ".softClip.tail";	
		string flagStr = samFieldVec[1];
		int tmpFlag = atoi(flagStr.c_str());
		string mapChrNameStr = samFieldVec[2];
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];			
		string readSeq = samFieldVec[9];
		int readLength = readSeq.length();
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(!mappedOrNot_bool)
			continue;

		int headSoftClipLength = 0;
		int tailSoftClipLength = 0;
		vector<Jump_Code> tmpJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);

		int mapChrPos_end = getAlignmentEndPos(mapChrPos, tmpJumpCodeVec);

		int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
		string firstJumpCodeType = tmpJumpCodeVec[0].type;
		int firstJumpCodeLength = tmpJumpCodeVec[0].len;
		string lastJumpCodeType = tmpJumpCodeVec[tmpJumpCodeVecSize-1].type;
		int lastJumpCodeLength = tmpJumpCodeVec[tmpJumpCodeVecSize-1].len;
		
		if(firstJumpCodeType == "S")
		{
			headSoftClipLength = firstJumpCodeLength;
			outputSoftClipSAM_ofs << readName_softClip_head << "\t0\t" << mapChrNameStr 
				<< "\t" << mapChrPos-headSoftClipLength << "\t255\t" << headSoftClipLength 
				<< "M\t*\t0\t0\t" << readSeq.substr(0,headSoftClipLength) 
				<< "\t*\tIH:i:1\tIH:i:1" << endl;  
		}
		if(lastJumpCodeType == "S")
		{
			tailSoftClipLength = lastJumpCodeLength;
			outputSoftClipSAM_ofs << readName_softClip_tail << "\t0\t" << mapChrNameStr 
				<< "\t" << mapChrPos_end+1 << "\t255\t" << tailSoftClipLength 
				<< "M\t*\t0\t0\t" << readSeq.substr(readLength-tailSoftClipLength, tailSoftClipLength) 
				<< "\t*\tIH:i:1\tIH:i:1" << endl;  		
		}

	}
	inputSAM_ifs.close();
	outputSoftClipSAM_ofs.close();
	return 0;
}