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

void cigarString2jumpCodeVec(string& jumpCodeStr, 
	vector<Jump_Code>& cigarStringJumpCodeVec)
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

int getEndPosOfSpecificJumpCode(int startPos, 
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
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

void generateMatchPosPairVecFromJumpCodeVec(int startMapPos, 
	vector<Jump_Code>& cigarStringJumpCodeVec,
	vector< pair<int,int> >& matchPosPairVec)
{
	//vector<int> matchIndexVec_cigarStringJumpCodeVec;
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		if(tmpJumpCodeType == "M")
		{
			int lastJumpCodeIndex = tmp - 1;
			int currentJumpCodeIndex = tmp;
			int tmpSegStartPos, tmpSegEndPos;
			if(lastJumpCodeIndex < 0)
				tmpSegStartPos = startMapPos;
			else
				tmpSegStartPos = getEndPosOfSpecificJumpCode(startMapPos,
					cigarStringJumpCodeVec, lastJumpCodeIndex) + 1;
			tmpSegEndPos = getEndPosOfSpecificJumpCode(startMapPos,
				cigarStringJumpCodeVec, currentJumpCodeIndex);
			matchPosPairVec.push_back(pair<int,int>(tmpSegStartPos, tmpSegEndPos));
		}
	}		
}

void extractChrNamePosFromAsmStr(
	string& tmpAsmStr, string& tmpChrName, 
	int& tmpStartPos, int& tmpEndPos, string& tmpOtherStr)
{
	int tabLoc_1 = tmpAsmStr.find("\t");
	int tabLoc_2 = tmpAsmStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpAsmStr.find("\t", tabLoc_2 + 1);
	tmpChrName = tmpAsmStr.substr(0, tabLoc_1);
	string tmpStartPosStr = tmpAsmStr.substr(
		tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEndPosStr = tmpAsmStr.substr(
		tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	tmpStartPos = atoi(tmpStartPosStr.c_str());
	tmpEndPos = atoi(tmpEndPosStr.c_str());
	tmpOtherStr = tmpAsmStr.substr(tabLoc_3 + 1);
}

void getChrNameMapPosPairVecFromSamStr(string& samStr, 
	string& tmpSam_chrName, vector< pair<int,int> >& tmpSam_mapPosPairVec)
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
	tmpSam_chrName = samFieldVec[2];
	string tmpSam_startPosStr = samFieldVec[3];
	int tmpSam_startPos = atoi(tmpSam_startPosStr.c_str());
	string tmpSam_cigarString = samFieldVec[5];
	vector<Jump_Code> tmpCigarStringJumpCodeVec;
	cigarString2jumpCodeVec(tmpSam_cigarString, tmpCigarStringJumpCodeVec);
	generateMatchPosPairVecFromJumpCodeVec(tmpSam_startPos, 
		tmpCigarStringJumpCodeVec, tmpSam_mapPosPairVec);
}

void addReadCountToCoveredAsm(string& tmpSam_chrName, 
	vector< pair<int,int> >& tmpSam_mapPosPairVec,
	vector<string>& asmVec_chrName, vector<int>& asmVec_startPos, 
	vector<int>& asmVec_endPos, vector<int>& asmVec_readCount)
{
	for(int tmpAsm = 0; tmpAsm < asmVec_chrName.size(); tmpAsm ++)
	{
		string tmpAsmChrName = asmVec_chrName[tmpAsm];
		if(tmpSam_chrName == tmpAsmChrName)
		{
			int tmpAsmStartPos = asmVec_startPos[tmpAsm];
			int tmpAsmEndPos = asmVec_endPos[tmpAsm];
			for(int tmpSeg = 0; tmpSeg < tmpSam_mapPosPairVec.size(); tmpSeg++)
			{
				int tmpSegStartPos = tmpSam_mapPosPairVec[tmpSeg].first;
				int tmpSegEndPos = tmpSam_mapPosPairVec[tmpSeg].second;
				if(!((tmpSegStartPos > tmpAsmEndPos)||(tmpSegEndPos < tmpAsmStartPos)))
				{
					// overlapped with asm
					asmVec_readCount[tmpAsm] ++;
					break;
				}
			}	
		}
	}
}


int main(int argc, char** argv)
{
	if(argc <= 3)
	{
		cout << "Executable inputASMfile outputAsmFileWithReadCount inputSAM_1 (inputSAM_2 ...) " << endl;
		exit(1);
	}
	string inputASMfile = argv[1];
	string outputAsmFileWithReadCount = argv[2];
	int inputSamFileNum = argc - 3;
	vector<string> inputSamFileVec;
	for(int tmp = 3; tmp < argc; tmp++)
		inputSamFileVec.push_back(argv[tmp]);

	vector<string> asmVec_chrName;
	vector<int> asmVec_startPos;
	vector<int> asmVec_endPos;
	vector<int> asmVec_readCount;
	vector<string> asmVec_otherStr;
	cout << "start to extract asm info from asm file" << endl;
	ifstream asm_ifs(inputASMfile.c_str());
	string headerStr;
	getline(asm_ifs, headerStr);
	while(!asm_ifs.eof())
	{
		string tmpAsmStr;
		getline(asm_ifs, tmpAsmStr);
		string tmpChrName;
		int tmpStartPos, tmpEndPos;
		string tmpOtherStr;
		extractChrNamePosFromAsmStr(tmpAsmStr, tmpChrName, 
			tmpStartPos, tmpEndPos, tmpOtherStr);
		asmVec_chrName.push_back(tmpChrName);
		asmVec_startPos.push_back(tmpStartPos);
		asmVec_endPos.push_back(tmpEndPos);
		asmVec_readCount.push_back(0);
		asmVec_otherStr.push_back(tmpOtherStr);
	}
	asm_ifs.eof();
	cout << "start to scan files to do read counts for each asm" << endl;
	for(int tmpSamFileIndex = 0; tmpSamFileIndex < inputSamFileNum; 
		tmpSamFileIndex ++)
	{
		string tmpSamFile = inputSamFileVec[tmpSamFileIndex];
		cout << endl << "tmpSamFile: " << tmpSamFile << endl;
		ifstream tmpSam_ifs(tmpSamFile.c_str());
		int tmpAlignmentNum = 0;
		while(!tmpSam_ifs.eof())
		{
			string tmpSamStr;
			getline(tmpSam_ifs, tmpSamStr);
			if(tmpSamStr == "")
				break;
			if(tmpSamStr.substr(0,1) == "@")
				continue;
			tmpAlignmentNum ++;
			int tmpTenThousandIndex = tmpAlignmentNum / 10000;
			if(tmpTenThousandIndex * 10000 == tmpAlignmentNum)
				cout << "readNum: " << tmpAlignmentNum << endl;
			string tmpSam_chrName;
			vector< pair<int,int> > tmpSam_mapPosPairVec;
			getChrNameMapPosPairVecFromSamStr(tmpSamStr, 
				tmpSam_chrName, tmpSam_mapPosPairVec);
			addReadCountToCoveredAsm(tmpSam_chrName, tmpSam_mapPosPairVec,
				asmVec_chrName, asmVec_startPos, asmVec_endPos, asmVec_readCount);
		}
		tmpSam_ifs.eof();
	}
	cout << "start to output asm file with read count " << endl;
	cout << "outputAsmFileWithReadCount: " << outputAsmFileWithReadCount << endl;
	ofstream asmWithReadCount_ofs(outputAsmFileWithReadCount.c_str());
	for(int tmpAsm = 0; tmpAsm < asmVec_chrName.size(); tmpAsm ++)
	{
		//cout << "tmpAsm: " << tmpAsm << endl;
		int tmpAsmID = tmpAsm + 1;
		asmWithReadCount_ofs << asmVec_chrName[tmpAsm] << "\t"
			<< asmVec_startPos[tmpAsm] << "\t" << asmVec_endPos[tmpAsm] << "\t"
			<< asmVec_readCount[tmpAsm] << "\tASM_" << tmpAsmID << "\t"
			<< asmVec_otherStr[tmpAsm] << endl; 
	}
	asmWithReadCount_ofs.close();
	cout << "end of counting reads over ASMs" << endl;
	return 0;
}