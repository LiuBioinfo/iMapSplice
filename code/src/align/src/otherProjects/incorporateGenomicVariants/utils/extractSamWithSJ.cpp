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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	int jumpCodeIndex)
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

void generateSJposVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	vector< pair<int,int> >& SJposPairVec)
{
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		if(tmpJumpCodeType == "N")
		{
			int lastJumpCodeIndex = tmp-1;
			int currentJumpCodeIndex = tmp;
			int tmpDonerEndPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
			int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
			SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
		}
	}		
}

void samStr2SJposPairVec(string& samStr, int& tmpSJchrNameInt,
	vector< pair<int,int> >& tmpSJposPairVec, Index_Info* indexInfo)
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
	tmpSJchrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
	string mapChrPosStr = samFieldVec[3];
	int mapChrPos = atoi(mapChrPosStr.c_str());
	string cigarString = samFieldVec[5];
	vector<Jump_Code> cigarStringJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
	generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, tmpSJposPairVec);
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndex inputSJlist inputSAM outputFile" << endl;
		exit(1);
	}
	cout << "loading indexInfo parameters ......" << endl;
	//log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	//parameter_ifs.close();
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	parameter_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;

	string juncHash_file_str = argv[2];
	AlignInferJunctionHash_Info* juncHash = new AlignInferJunctionHash_Info();
	juncHash->initiateAlignInferJunctionInfo(chromNum);
	juncHash->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_str, indexInfo);

	string inputSAM = argv[3];
	string outputFile = argv[4];
	ifstream sam_ifs(inputSAM.c_str());
	ofstream samWithSJ_ofs(outputFile.c_str());
	while(!sam_ifs.eof())
	{
		string tmpSamStr;
		getline(sam_ifs, tmpSamStr);
		if(tmpSamStr == "")
			break;
		if(tmpSamStr.substr(0,1) == "@")
			continue;
		int tmpSJchrNameInt;
		vector< pair<int,int> > tmpSJposPairVec;
		samStr2SJposPairVec(tmpSamStr, tmpSJchrNameInt, tmpSJposPairVec, indexInfo);
		int tmpSJnum = tmpSJposPairVec.size();
		vector< pair<int,int> > tmpSJposPairVec_found;
		for(int tmp = 0; tmp < tmpSJnum; tmp ++)
		{
			int tmpSJstartPos = tmpSJposPairVec[tmp].first;
			int tmpSJendPos = tmpSJposPairVec[tmp].second;
			//cout << "tmpSJstartPos: " << tmpSJstartPos << endl;
			//cout << "tmpSJendPos: " << tmpSJendPos << endl; 
 			int tmpSJ_index = juncHash->searchAndReturnAlignInferInfoVecIndex(
				tmpSJchrNameInt, tmpSJstartPos, tmpSJendPos);
			//cout << "tmpSJ_index: " << tmpSJ_index << endl;
			if(tmpSJ_index >= 0)
			{
				tmpSJposPairVec_found.push_back(pair<int,int>(tmpSJstartPos, tmpSJendPos));
			}
		}
		int tmpSJposPairVec_found_num = tmpSJposPairVec_found.size();
		if(tmpSJposPairVec_found_num > 0)
		{
			string tmpSJchrName = indexInfo->returnChrNameStr(tmpSJchrNameInt);
			samWithSJ_ofs << tmpSamStr << "\t";
			for(int tmp = 0; tmp < tmpSJposPairVec_found_num; tmp++)
			{
				samWithSJ_ofs << tmpSJchrName << ":" 
					<< tmpSJposPairVec_found[tmp].first << "-"
					<< tmpSJposPairVec_found[tmp].second << ","; 
			}
			samWithSJ_ofs << endl;
		}
	}
	samWithSJ_ofs.close();
	sam_ifs.close();
	delete juncHash;
	delete indexInfo;
	parameter_ifs.close();	
	return 0;
}