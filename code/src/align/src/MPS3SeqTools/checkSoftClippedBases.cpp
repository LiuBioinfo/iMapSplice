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

int main(int argc, char** argv)
{
	if(argc < 6)
	{
		cout << "Executable outputSoftClippedBasesDistributionFile readLength_max read_num(NOT pair) Aligner_1_name SAM_file_1 (Aligner_2_name SAM_file_2 ...)" << endl;
		exit(1);
	}
	string outputFileStr = argv[1];
	ofstream output_ofs(outputFileStr.c_str());
	string outputFileStr_2plotInR = outputFileStr + ".2plotInR";
	ofstream output_2plotInR_ofs(outputFileStr_2plotInR.c_str());

	string readLengthMaxStr = argv[2];
	int readLengthMax = atoi(readLengthMaxStr.c_str());
	string readNumStr = argv[3];
	unsigned int readNum = atoi(readNumStr.c_str());

	vector<string> alignerNameVec;
	vector<string> alignmentFileVec;
	for(int tmp = 4; tmp < argc; tmp += 2)
	{
		string tmpAlignerName = argv[tmp];
		alignerNameVec.push_back(tmpAlignerName);
		string tmpAlignmentFileStr = argv[tmp+1];
		alignmentFileVec.push_back(tmpAlignmentFileStr);
	}

	int alignmentFileVecSize = alignmentFileVec.size();
	unsigned int* clippedBase_num = (unsigned int*)malloc(
		(readLengthMax) * alignmentFileVecSize * sizeof(unsigned int));
	for(int tmp = 0; tmp < readLengthMax * alignmentFileVecSize; tmp++)
		clippedBase_num[tmp] = 0;
	int mappedRead_num[alignmentFileVecSize];
	for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner++)
	{	
		//mappedRead_num[tmpAligner] = 0;
		/// count primary alignment number for each mapped length
		string tmpInputSAMfile = alignmentFileVec[tmpAligner];
		ifstream tmpInputSAM_ifs(tmpInputSAMfile.c_str());
		unsigned int unmapAlignNum = 0;
		unsigned int totalAlignNum = 0;
		unsigned int primaryAlignNum = 0;
		unsigned int secondaryAlignNum = 0;
		while(!(tmpInputSAM_ifs.eof()))
		{
			// if(inputSAMfile_ifs.eof())
			// 	break;
			string samStr;
			getline(tmpInputSAM_ifs, samStr);
			// if(inputSAMfile_ifs.eof())
			// 	break;
			if(samStr.at(0) == '@')
				continue;
			totalAlignNum ++;

			vector<string> samFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 5; tmp++)
			{
				int tabLoc = samStr.find("\t", startLoc);
				string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
				samFieldVec.push_back(tmpSamField);
				startLoc = tabLoc + 1;
			}
			samFieldVec.push_back(samStr.substr(startLoc));
			string flagStr = samFieldVec[1];
			int tmpFlag = atoi(flagStr.c_str());
			string mapChrNameStr = samFieldVec[2];
			string mapChrPosStr = samFieldVec[3];
			int mapChrPos = atoi(mapChrPosStr.c_str());
			string cigarString = samFieldVec[4];			
	
			bool mappedOrNot_bool = mappedOrNot(tmpFlag);
			bool primaryOrNot_bool = primaryOrNot(tmpFlag);
			
			if(!mappedOrNot_bool)
			{
				unmappedAlignmentNum ++;
			}
			else
			{
				if(primaryOrNot_bool)
				{
					primaryAlignmentNum ++;
				}
				else
				{
					secondaryAlignNum ++;
				}
			}

			if(primaryOrNot_bool)
			{
				int headSoftClipLength = 0;
				int tailSoftClipLength = 0;
				int softClipLength = 0;
				vector<Jump_Code> tmpJumpCodeVec;
				cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);
				int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
				string firstJumpCodeType = tmpJumpCodeVec[0].type;
				int firstJumpCodeLength = tmpJumpCodeVec[0].len;
				string lastJumpCodeType = tmpJumpCodeVec[tmpJumpCodeVecSize-1].type;
				int lastJumpCodeLength = tmpJumpCodeVec[tmpJumpCodeVecSize-1].len;
				if(firstJumpCodeType == "S")
				{
					headSoftClipLength = firstJumpCodeLength;
				}
				if(lastJumpCodeType == "S")
				{
					tailSoftClipLength = lastJumpCodeLength;
				}
				softClipLength = headSoftClipLength + tailSoftClipLength;
				if(softClipLength > 0)
				{	
					for(int tmp = 1; tmp <= softClipLength; tmp++)
					{
						clippedBase_num[readLengthMax * tmpAligner + tmp -1] ++;
					}
				}
			}
		}
		mappedRead_num[tmpAligner] = readNum - primaryAlignNum;
	}

	// output stats ....
	output_ofs << "clipped base";
	for(int tmp = 0; tmp < alignmentFileVecSize; tmp++)
	{
		output_ofs << "\t" << alignerNameVec[tmp]; 
	}
	output_ofs << endl;
	return 0;
}