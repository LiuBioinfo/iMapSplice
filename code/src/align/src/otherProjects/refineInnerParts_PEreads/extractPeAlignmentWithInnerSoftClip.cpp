// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// used to separate forward alignments and rev alignments
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

using namespace std;

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool alignmentForwardOrRev(int tmpFlag)
{
	if(tmpFlag & 0x10)
		return false;
	else
		return true;
}	

void cigarString2jumpCodeVec(vector<string>& alignJumpCodeTypeVec, string& jumpCodeStr)
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
			alignJumpCodeTypeVec.push_back(tmpJumpCodeType);
			jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
		}
	}
}

bool headSoftClippedOrNot(string& jumpCodeStr)
{
	vector<string> tmpJumpCodeTypeVec;
	cigarString2jumpCodeVec(tmpJumpCodeTypeVec, jumpCodeStr);
	if(tmpJumpCodeTypeVec[0] == "S")
		return true;
	else
		return false;
}

bool tailSoftClippedOrNot(string& jumpCodeStr)
{
	vector<string> tmpJumpCodeTypeVec;
	cigarString2jumpCodeVec(tmpJumpCodeTypeVec, jumpCodeStr);
	if(tmpJumpCodeTypeVec[tmpJumpCodeTypeVec.size()-1] == "S")
		return true;
	else
		return false;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputPeSamFile outputInnerSoftClipSamFile" << endl;
		exit(1);
	}
	string inputSAMfileStr = argv[1];
	//string outputFolderStr = argv[2];
	//outputFolderStr += "/";
	//string mkdir = "mkdir -p " + outputFolderStr;
	//system(mkdir.c_str());
	//string output_unmap_file_str = outputFolderStr + "unmap.sam";
	//string output_for_file_str = outputFolderStr + "for.sam";
	//string output_rev_file_str = outputFolderStr + "rev.sam";
	//ofstream unmap_ofs(output_unmap_file_str.c_str());
	//ofstream for_ofs(output_for_file_str.c_str());
	//ofstream rev_ofs(output_rev_file_str.c_str());
	string outputInnerSoftClipPeAlignmentFileStr = argv[2];
	ofstream innerSoftClip_ofs(outputInnerSoftClipPeAlignmentFileStr.c_str());

	ifstream tmpInputSAM_ifs(inputSAMfileStr.c_str());
	while(!(tmpInputSAM_ifs.eof()))
	{
		string samStr_for;
		getline(tmpInputSAM_ifs, samStr_for);
		if(tmpInputSAM_ifs.eof()||(samStr_for == ""))
		 	break;
		if(samStr_for.at(0) == '@')
			continue;
		vector<string> samFieldVec_for;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr_for.find("\t", startLoc);
			string tmpSamField = samStr_for.substr(startLoc, tabLoc-startLoc);
			samFieldVec_for.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_for.push_back(samStr_for.substr(startLoc));
		string samStr_rcm;
		getline(tmpInputSAM_ifs, samStr_rcm);
		vector<string> samFieldVec_rcm;
		startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr_rcm.find("\t", startLoc);
			string tmpSamField = samStr_rcm.substr(startLoc, tabLoc-startLoc);
			samFieldVec_rcm.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_rcm.push_back(samStr_rcm.substr(startLoc));		

		string flagStr_for = samFieldVec_for[1];
		int tmpFlag_for = atoi(flagStr_for.c_str());
		bool mappedOrNot_for_bool = mappedOrNot(tmpFlag_for);
		string flagStr_rcm = samFieldVec_rcm[1];
		int tmpFlag_rcm = atoi(flagStr_rcm.c_str());
		bool mappedOrNot_rcm_bool = mappedOrNot(tmpFlag_rcm);

		if(mappedOrNot_for_bool && mappedOrNot_rcm_bool)
		{
			string cigarString_for = samFieldVec_for[5];
			string cigarString_rcm = samFieldVec_rcm[5];
			bool tailClip_for_bool = tailSoftClippedOrNot(cigarString_for);
			bool headClip_rcm_bool = headSoftClippedOrNot(cigarString_rcm);
			if(tailClip_for_bool || headClip_rcm_bool)
			{	
				innerSoftClip_ofs << samStr_for << endl;
				innerSoftClip_ofs << samStr_rcm << endl;
			}
			else
			{}
		}
		else
		{
		}
	}
	innerSoftClip_ofs.close();
	tmpInputSAM_ifs.close();
	return 0;
}