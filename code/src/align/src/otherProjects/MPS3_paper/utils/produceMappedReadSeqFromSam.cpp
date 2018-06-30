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

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"

using namespace std;

bool bothEndsMapped(int tmpFlag)
{
	if(tmpFlag & 0x2)
		return true;
	else
		return false;
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
int getSoftClippedHeadLengthFromCigarString(string& tmpCigarStringStr)
{
	vector<Jump_Code> tmpJumpCodeVec;
	//cout << "cigarString " << cigarString << endl;
	cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);	
	if(tmpJumpCodeVec[0].type == "S")
		return tmpJumpCodeVec[0].len;
	else
		return 0; 
}	

int getSoftClippedTailLengthFromCigarString(string& tmpCigarStringStr)
{
	vector<Jump_Code> tmpJumpCodeVec;
	//cout << "cigarString " << cigarString << endl;
	cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);
	int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
	if(tmpJumpCodeVec[tmpJumpCodeVecSize-1].type == "S")
		return tmpJumpCodeVec[tmpJumpCodeVecSize-1].len;
	else
		return 0; 	
}


int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSAM outputFasta" << endl;
		exit(1);
	}
	string inputSAM = argv[1];
	string outputFasta = argv[2];
	ifstream sam_ifs(inputSAM.c_str());
	ofstream fasta_ofs(outputFasta.c_str());
	cout << "start to read sam and generate fa" << endl;
	while(!sam_ifs.eof())
	{
		string samStr;
		getline(sam_ifs, samStr);
		if(samStr == "")
			break;
		if(samStr.at(0) == '@')
			continue;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string readName = samFieldVec[0];
		samFieldVec.push_back(samStr.substr(startLoc));
		string flagStr = samFieldVec[1];
		int flagInt = atoi(flagStr.c_str());
		string cigarString = samFieldVec[5];
		string readSeq = samFieldVec[9];
		string qualSeq = samFieldVec[10];
		bool mappedOrNot_bool = mappedOrNot(flagInt);
		bool primaryOrNot_bool = primaryOrNot(flagInt);
		bool bothEndsMapped_bool = bothEndsMapped(flagInt);
		if(mappedOrNot_bool)
		{
			if(bothEndsMapped_bool)
			{
				if(primaryOrNot_bool)
				{
					int readLength = readSeq.length();
					int softClippedHeadLength = getSoftClippedHeadLengthFromCigarString(cigarString);
					int softClippedTailLength = getSoftClippedTailLengthFromCigarString(cigarString);
					int startMappedBaseLocInRead = softClippedHeadLength + 1;
					int endMappedBaseLocInRead = readLength - softClippedTailLength;
					string tmpReadSubSeq = readSeq.substr(startMappedBaseLocInRead - 1, endMappedBaseLocInRead - startMappedBaseLocInRead + 1);
					string tmpQualSubSeq = qualSeq.substr(startMappedBaseLocInRead - 1, endMappedBaseLocInRead - startMappedBaseLocInRead + 1);
					fasta_ofs << "@" << readName << endl << tmpReadSubSeq << endl << "+" << endl << tmpQualSubSeq << endl;
				}
				else
				{}
			}
			else
			{}
		}
		else
		{}
	}
	sam_ifs.close();
	fasta_ofs.close();
	return 0;
}