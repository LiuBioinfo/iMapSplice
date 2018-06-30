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
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "../../general/read_block_test.h"
#include "../../general/splice_info.h"
using namespace std;

bool forOrRcm(int flag)
{
	if(flag & 0x10)
		return false;
	else
		return true;
}

void extractFieldVal_IH_XI(int& tmpIHint, int& tmpXIint,
	string& tmpUnpairedPEalignmentStr)
{
	int IHfieldLoc = tmpUnpairedPEalignmentStr.find("IH:i:");
	int XIfieldLoc = tmpUnpairedPEalignmentStr.find("XI:i:");
	int IHfieldLoc_nextTab = tmpUnpairedPEalignmentStr.find("\t", IHfieldLoc);
	int lineLength = tmpUnpairedPEalignmentStr.length();
	string IHintStr = tmpUnpairedPEalignmentStr.substr(IHfieldLoc + 5, IHfieldLoc_nextTab - IHfieldLoc - 5);
	string XIintStr = tmpUnpairedPEalignmentStr.substr(XIfieldLoc + 5, lineLength - XIfieldLoc - 5);
	tmpIHint = atoi(IHintStr.c_str());
	tmpXIint = atoi(XIintStr.c_str());
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

void getChrNameAndStartEndMapPosFromSamStr(
	string& tmpChrName, int& tmpStartMapPos, int& tmpEndMapPos,
	int& unfixedHeadLength, int& unfixedTailLength, bool& Nor_or_Rcm_bool,
	string& tmpSamStr)
{
	vector<string> tmpSamFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		tmpSamFieldVec.push_back(tmpSamStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	int flagInt = atoi(tmpSamFieldVec[1].c_str());
	Nor_or_Rcm_bool = forOrRcm(flagInt);

	tmpChrName = tmpSamFieldVec[2];
	string tmpStartMapPosStr = tmpSamFieldVec[3];
	tmpStartMapPos = atoi(tmpStartMapPosStr.c_str());
	string cigarString = tmpSamFieldVec[5];
	vector<Jump_Code> tmpJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);
	tmpEndMapPos = getEndPosOfSpecificJumpCode(tmpStartMapPos, tmpJumpCodeVec, tmpJumpCodeVec.size()-1);

	if(tmpJumpCodeVec[0].type == "S")
		unfixedHeadLength = tmpJumpCodeVec[0].len;
	else
		unfixedHeadLength = 0;	

	if(tmpJumpCodeVec[tmpJumpCodeVec.size()-1].type == "S")
		unfixedTailLength = tmpJumpCodeVec[tmpJumpCodeVec.size()-1].len;
	else
		unfixedTailLength = 0;	
}

void getChrNameAndStartMapPosFromSamStr(
	string& tmpChrName, int& tmpStartMapPos, int& unfixedHeadLength, string& tmpSamStr)
{
	vector<string> tmpSamFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		tmpSamFieldVec.push_back(tmpSamStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	string cigarString = tmpSamFieldVec[5];
	vector<Jump_Code> tmpJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);	
	tmpChrName = tmpSamFieldVec[2];
	string tmpStartMapPosStr = tmpSamFieldVec[3];
	tmpStartMapPos = atoi(tmpStartMapPosStr.c_str());
	if(tmpJumpCodeVec[0].type == "S")
		unfixedHeadLength = tmpJumpCodeVec[0].len;
	else
		unfixedHeadLength = 0;
}

void getChrNameAndEndMapPosFromSamStr(
	string& tmpChrName, int& tmpEndMapPos, int& unfixedTailLength, string& tmpSamStr)
{
	vector<string> tmpSamFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		tmpSamFieldVec.push_back(tmpSamStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	tmpChrName = tmpSamFieldVec[2];
	string tmpStartMapPosStr = tmpSamFieldVec[3];
	int tmpStartMapPos = atoi(tmpStartMapPosStr.c_str());
	string cigarString = tmpSamFieldVec[5];
	vector<Jump_Code> tmpJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);
	tmpEndMapPos = getEndPosOfSpecificJumpCode(tmpStartMapPos, tmpJumpCodeVec, tmpJumpCodeVec.size()-1);
	if(tmpJumpCodeVec[tmpJumpCodeVec.size()-1].type == "S")
		unfixedTailLength = tmpJumpCodeVec[tmpJumpCodeVec.size()-1].len;
	else
		unfixedTailLength = 0;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputUnpairedPEalignmentPath outputFile" << endl;
		exit(1);	
	}

	string outputFile = argv[2];
	ofstream candiFusionRegion_ofs(outputFile.c_str());
	string inputUnpairedPEalignmentPath = argv[1];
	ifstream unpairedPEalignment_ifs(inputUnpairedPEalignmentPath.c_str());
	while(!unpairedPEalignment_ifs.eof())
	{
		string tmpUnpairedPEalignmentStr;
		getline(unpairedPEalignment_ifs, tmpUnpairedPEalignmentStr);
		if((unpairedPEalignment_ifs.eof())||(tmpUnpairedPEalignmentStr == ""))
			break;
		int tmpIHint, tmpXIint;
		vector<string> tmpEndSamStrVec_1;
		vector<string> tmpEndSamStrVec_2;
		extractFieldVal_IH_XI(tmpIHint, tmpXIint, tmpUnpairedPEalignmentStr);
		if(tmpIHint == 0)
		{}
		else
		{
			tmpEndSamStrVec_1.push_back(tmpUnpairedPEalignmentStr);
			for(int tmp = 1; tmp < tmpIHint; tmp++)
			{
				string tmpStr;
				getline(unpairedPEalignment_ifs, tmpStr);
				tmpEndSamStrVec_1.push_back(tmpStr);
			}
		}
		if(tmpXIint == 0)
		{
			string tmpStr;
			getline(unpairedPEalignment_ifs, tmpStr);
		}
		else
		{
			for(int tmp = 0; tmp < tmpXIint; tmp++)
			{
				string tmpStr;
				getline(unpairedPEalignment_ifs, tmpStr);
				tmpEndSamStrVec_2.push_back(tmpStr);
			}
		}

		if((tmpEndSamStrVec_1.size() == 1)&&(tmpEndSamStrVec_2.size() == 1))
		{
			string tmpChrName_1, tmpChrName_2;
			int tmpChrStartMapPos_1;
			int tmpChrEndMapPos_1;
			int tmpChrStartMapPos_2;
			int tmpChrEndMapPos_2;
			int unfixedHeadLength_1;
			int unfixedTailLength_1;
			int unfixedHeadLength_2;
			int unfixedTailLength_2;
			bool forOrRevBool_1;
			bool forOrRevBool_2;

			getChrNameAndStartEndMapPosFromSamStr(tmpChrName_1, tmpChrStartMapPos_1, tmpChrEndMapPos_1,
				unfixedHeadLength_1, unfixedTailLength_1, forOrRevBool_1, tmpEndSamStrVec_1[0]);

			getChrNameAndStartEndMapPosFromSamStr(tmpChrName_2, tmpChrStartMapPos_2, tmpChrEndMapPos_2,
				unfixedHeadLength_2, unfixedTailLength_2, forOrRevBool_2, tmpEndSamStrVec_2[0]);

			candiFusionRegion_ofs << tmpChrName_1 << "\t" << tmpChrStartMapPos_1 << "\t" 
				<< tmpChrEndMapPos_1 << "\t" << unfixedHeadLength_1 << "\t" << unfixedTailLength_1 << "\t";
			if(forOrRevBool_1)
				candiFusionRegion_ofs << "For\t";
			else
				candiFusionRegion_ofs << "Rev\t";

			candiFusionRegion_ofs << tmpChrName_2 << "\t" << tmpChrStartMapPos_2 << "\t" 
				<< tmpChrEndMapPos_2 << "\t" << unfixedHeadLength_2 << "\t" << unfixedTailLength_2 << "\t";
			if(forOrRevBool_2)
				candiFusionRegion_ofs << "For" << endl;
			else
				candiFusionRegion_ofs << "Rev" << endl;
		}
	}
	unpairedPEalignment_ifs.close();
	candiFusionRegion_ofs.close();
	return 0;
}