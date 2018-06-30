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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
#include "../../general/alignInferJunctionHash_info.h"
using namespace std;

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

int getChrNameIntFromSamStr(const string& samStr, Index_Info* indexInfo)
{
	vector<string> samFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = samStr.find("\t", startLoc);
		string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
		samFieldVec.push_back(tmpSamField);
		startLoc = tabLoc + 1;
	}
	string mapChrNameStr = samFieldVec[2];
	int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
	//string mapChrPosStr = samFieldVec[3];
	//int mapChrPos = atoi(mapChrPosStr.c_str());
	//string cigarString = samFieldVec[5];	
	return mapChrNameInt;
}
int getEndMapPosFromSamStr(const string& samStr, Index_Info* indexInfo)
{
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
	//cout << "mapChrPos: " << mapChrPos << endl;
	string cigarString = samFieldVec[5];	
	//cout << "cigarString: " << cigarString << endl;
	vector<Jump_Code> cigarJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, cigarJumpCodeVec);
	//cout << "cigarJumpCodeVec.size(): " << cigarJumpCodeVec.size() << endl;
	int mapChrPos_end = getEndPosOfSpecificJumpCode(mapChrPos, cigarJumpCodeVec, cigarJumpCodeVec.size()-1);
	return mapChrPos_end;
}
int getStartMapPosFromSamStr(const string& samStr, Index_Info* indexInfo)
{
	vector<string> samFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 4; tmp++)
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
	//string cigarString = samFieldVec[5];	
	return mapChrPos;
}


int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndex inputBackSpliceJuncFile inputSamFile outputBackSpliceFileWithEncompassingSupNum" << endl;
		exit(1);
	}

	int bufferSizeToSearchForBackSpliceEncompassed = 10000;

	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	cout << "finish loading chromosomes" << endl;

	// generating alignInferJuncHash
	string inputJuncFile = argv[2];
	cout << "start to insert SJs from JuncFile into alignInferJuncHash" << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePos_supportNum_backSpliceOnly(
		inputJuncFile, indexInfo);

	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);

	//input alignment to generate encompassing read numbers for each back splice junction
	cout << "start to generate encompassing read num" << endl;
	string inputSamFile = argv[3];
	ifstream sam_ifs(inputSamFile.c_str());
	while(!sam_ifs.eof())
	{
		string samStr_1, samStr_2;
		getline(sam_ifs, samStr_1);
		if(samStr_1 == "")
			break;
		if(samStr_1.substr(0,1) == "@")
			continue;
		getline(sam_ifs, samStr_2);
		//cout << "samStr_1: " << samStr_1 << endl;
		//cout << "samStr_2: " << samStr_2 << endl;		
		int tmpChrNameInt = getChrNameIntFromSamStr(samStr_1, indexInfo);
		//cout << "tmpChrNameInt: " << tmpChrNameInt << endl;
		int startMapPos_1 = getStartMapPosFromSamStr(samStr_1, indexInfo);
		//cout << "startMapPos_1: " << startMapPos_1 << endl;
		int endMapPos_1 = getEndMapPosFromSamStr(samStr_1, indexInfo);
		//cout << "endMapPos_1: " << endMapPos_1 << endl;
		int startMapPos_2 = getStartMapPosFromSamStr(samStr_2, indexInfo);
		//cout << "startMapPos_2: " << startMapPos_2 << endl;
		int endMapPos_2 = getEndMapPosFromSamStr(samStr_2, indexInfo);
		//cout << "endMapPos_2: " << endMapPos_2 << endl;
		if((startMapPos_1 <= startMapPos_2)&&(endMapPos_1 <= endMapPos_2)) // coordinates meet colinear alignment standard
		{}
		else // possible back splice alignment
		{
			int tmpCandiBackSpliceStartPos_leftMost = endMapPos_1 + 1;
			int tmpCandiBackSpliceStartPos_rightMost = endMapPos_1 + bufferSizeToSearchForBackSpliceEncompassed;
			int tmpCandiBackSpliceEndPos_leftMost = startMapPos_2 - bufferSizeToSearchForBackSpliceEncompassed;
			int tmpCandiBackSpliceEndPos_rightMost = startMapPos_2 - 1;
			if(tmpCandiBackSpliceStartPos_rightMost <= tmpCandiBackSpliceEndPos_leftMost) // no space for back splice
			{}
			else
			{
				vector< pair<int,int> > tmpCandiSJposPairVec;
				SJhashInfo->searchForCandiSJwithinAreas(tmpChrNameInt,
					tmpCandiBackSpliceStartPos_leftMost, tmpCandiBackSpliceStartPos_rightMost,
					tmpCandiBackSpliceEndPos_leftMost, tmpCandiBackSpliceEndPos_rightMost, 
					tmpCandiSJposPairVec);
				int tmpEncompassingDistance_min = 1000000;
				int tmpEncompassingDistance_min_index = -1;
				for(int tmp = 0; tmp < tmpCandiSJposPairVec.size(); tmp++)
				{
					int tmpBackSpliceStartPos = tmpCandiSJposPairVec[tmp].first;
					int tmpBackSpliceEndPos = tmpCandiSJposPairVec[tmp].second;
					int tmpEncompassingDistance = tmpBackSpliceStartPos - endMapPos_1
						+ startMapPos_2 - tmpBackSpliceEndPos;
					if(tmpEncompassingDistance < tmpEncompassingDistance_min)
					{
						tmpEncompassingDistance_min = tmpEncompassingDistance;
						tmpEncompassingDistance_min_index = tmp;
					}
				}
				if(tmpEncompassingDistance_min_index >= 0)
				{
					int tmpFinalBackSpliceStartPos = tmpCandiSJposPairVec[tmpEncompassingDistance_min_index].first;
					int tmpFinalBackSpliceEndPos = tmpCandiSJposPairVec[tmpEncompassingDistance_min_index].second;
					int alignInferJuncIndex = alignInferJunctionHashInfo->searchAndReturnAlignInferInfoVecIndex(
						tmpChrNameInt, tmpFinalBackSpliceStartPos, tmpFinalBackSpliceEndPos);
					if(alignInferJuncIndex >= 0)
						alignInferJunctionHashInfo->updateEncompassingNum(alignInferJuncIndex);
					else
					{
						cout << "tmp back splice not found in alignInferJuncHash ...." << endl;
						exit(1);
					}
				}
			}
		}
	}
	sam_ifs.close();

	// output all back splice junctions with encompassing read numbers
	string outputFile = argv[4];
	alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_chrNamePos_supportNum_encompassingNum(
		indexInfo, outputFile);

	delete indexInfo;
	parameter_ifs.close();
	return 0;
}