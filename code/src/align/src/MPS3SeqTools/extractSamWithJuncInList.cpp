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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"

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

void generateSJposVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	vector< pair<int,int> >& SJposPairVec, vector<int>& SJindexVec_cigarStringJumpCodeVec)
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
			SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
		}
	}		
}	

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputJuncListFile inputSamFile outputFolder" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string inputJuncListFile = argv[2];
	AlignInferJunctionHash_Info* juncHashInfo = new AlignInferJunctionHash_Info();	
	juncHashInfo->initiateAlignInferJunctionInfo(chromNum);
	juncHashInfo->insertJuncFromJuncFile_chrNamePosOnly(inputJuncListFile, indexInfo);

	string outputFolderStr = argv[4];
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string outputPrefix = outputFolderStr + "/";
	string outputSamWithJuncInList_ori = outputPrefix + "withJuncInList_ori.sam";
	string outputSamWithJuncInList_modified = outputPrefix + "withJuncInList_modified.sam";
	string outputSamWithOutJuncInList_ori = outputPrefix + "withOutJuncInList_ori.sam";
	string outputSam_final = outputPrefix + "final.sam";
	ofstream withJuncInList_ori_ofs(outputSamWithJuncInList_ori.c_str());
	ofstream withOutJuncInList_ori_ofs(outputSamWithOutJuncInList_ori.c_str());
	ofstream withJuncInList_modified_ofs(outputSamWithJuncInList_modified.c_str());
	ofstream final_ofs(outputSam_final.c_str());

	string inputSamFile = argv[3];
	ifstream oriSam_ifs(inputSamFile.c_str());
	while(!oriSam_ifs.eof())
	{
		string samStr;
		getline(oriSam_ifs, samStr);
		if(samStr == "")
			break;
		if(samStr.substr(0,1) == "@")
		{
			final_ofs << samStr << endl;
			continue;
		}
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string readNameStr = samFieldVec[0];
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];
		string readSeq = samFieldVec[9];
		string readQualSeq = samFieldVec[10];
		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		bool someJuncExistInAlignerInferJuncHashBool
			= juncHashInfo->someJuncInVecExistInAlignInferJuncHash(
				mapChrNameInt, tmpSJposPairVec);
		if(someJuncExistInAlignerInferJuncHashBool)
		{
			withJuncInList_ori_ofs << samStr << endl;
			string modifiedSamStr = readNameStr + "\t4\t*\t0\t0\t*\t*\t0\t0\t" 
				+ readSeq + "\t" + readQualSeq + "\t" + "IH:i:0\tHI:i:0"; 
			withJuncInList_modified_ofs << modifiedSamStr << endl;
			final_ofs << modifiedSamStr << endl;
		}
		else
		{
			withOutJuncInList_ori_ofs << samStr << endl;
			final_ofs << samStr << endl;
		}
	}
	oriSam_ifs.close();
	withJuncInList_ori_ofs.close();
	withOutJuncInList_ori_ofs.close();
	withJuncInList_modified_ofs.close();
	final_ofs.close();
	delete juncHashInfo;
	delete indexInfo;
	return 0;
}