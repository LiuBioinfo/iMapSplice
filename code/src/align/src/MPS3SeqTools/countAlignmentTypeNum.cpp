// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to check aligment type number (unsplied, single-spliced, multi-spliced) for primary ones
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

int getSJnumFromCigarString(string& tmpCigarStringStr)
{
	vector<Jump_Code> tmpJumpCodeVec;
	//cout << "cigarString " << cigarString << endl;
	cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);	
	int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
	int tmpSJnum = 0;
	for(int tmpJumpCodeIndex = 0; tmpJumpCodeIndex < tmpJumpCodeVecSize; tmpJumpCodeIndex++)
	{
		string tmpJumpCodeType = tmpJumpCodeVec[tmpJumpCodeIndex].type;
		if(tmpJumpCodeType == "N")
			tmpSJnum ++;
	}
	return tmpSJnum;
}

int main(int argc, char** argv)
{
	if(argc <= 4)
	{
		cout << "Executable TotalAlignmentNumber OutputFolder Aligner_1_name SAM_file_1 (Aligner_2_name SAM_file_2 ...)" << endl;
		exit(1);
	}

	string totalAlignmentNumberStr = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string outputAlignmentTypeNumFile = outputFolderStr + "alignmentType_num.txt";
	string outputAlignmentTypePercFile = outputFolderStr + "alignmentType_perc.txt";
	ofstream alignmentTypeNum_ofs(outputAlignmentTypeNumFile.c_str());
	ofstream alignmentTypePerc_ofs(outputAlignmentTypePercFile.c_str());

	int totalAlignmentNumber = atoi(totalAlignmentNumberStr.c_str());
	cout << "total alignments number: " << totalAlignmentNumber << endl;
	vector<string> alignerNameVec;
	vector<string> alignmentFileVec;
	for(int tmp = 3; tmp < argc; tmp += 2)
	{
		string tmpAlignerName = argv[tmp];
		alignerNameVec.push_back(tmpAlignerName);
		string tmpAlignmentFileStr = argv[tmp+1];
		alignmentFileVec.push_back(tmpAlignmentFileStr);
	}
	int alignmentFileVecSize = alignmentFileVec.size();
	cout << "alignment files number: " << alignmentFileVecSize << endl;
	int SJnumMax = 2;
	vector< vector<int> > alignmentTypeNumVec;
	for(int tmp = 0; tmp < alignmentFileVecSize; tmp++)
	{
		vector<int> tmpAlignerAlignmentTypeNumVec;
		for(int tmp2 = 0; tmp2 <= SJnumMax; tmp2++)
		{
			tmpAlignerAlignmentTypeNumVec.push_back(0);
		}
		alignmentTypeNumVec.push_back(tmpAlignerAlignmentTypeNumVec);
	}
	for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner++)
	{	
		cout << "tmpAligner: " << tmpAligner << endl;
		//mappedRead_num[tmpAligner] = 0;
		/// count primary alignment number for each mapped length
		string tmpInputSAMfile = alignmentFileVec[tmpAligner];
		ifstream tmpInputSAM_ifs(tmpInputSAMfile.c_str());
		while(!(tmpInputSAM_ifs.eof()))
		{
			string samStr;
			getline(tmpInputSAM_ifs, samStr);
			 if(tmpInputSAM_ifs.eof()||(samStr == ""))
			 	break;
			if(samStr.at(0) == '@')
				continue;
			//totalAlignNum ++;
			//cout << "samStr: " << samStr << endl;
			vector<string> samFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 6; tmp++)
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
			string cigarString = samFieldVec[5];			
	
			//cout << "tmpFlag" << tmpFlag << endl;
			bool mappedOrNot_bool = mappedOrNot(tmpFlag);
			bool primaryOrNot_bool = primaryOrNot(tmpFlag);
			//cout << "primary: " << primaryOrNot_bool << endl;
			//cout << "mapped: " << mappedOrNot_bool << endl;
			if(!mappedOrNot_bool)
			{
				//unmappedAlignmentNum ++;
			}
			else
			{
				if(primaryOrNot_bool)
				{
					int tmpAlignmentSJnum = getSJnumFromCigarString(cigarString);
					if(tmpAlignmentSJnum >= SJnumMax)
					{
						(alignmentTypeNumVec[tmpAligner])[SJnumMax] ++;
					}
					else
					{
						(alignmentTypeNumVec[tmpAligner])[tmpAlignmentSJnum] ++;
					}
				}
				else
				{
					//secondaryAlignNum ++;
				}
			}
		}
	}	

	for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner++)
	{
		int unsplicedAlignmentNum = (alignmentTypeNumVec[tmpAligner])[0];
		int singleSplicedAlignmentNum = (alignmentTypeNumVec[tmpAligner])[1];
		int multiSplicedAlignmentNum = (alignmentTypeNumVec[tmpAligner])[2]; 
		double unsplicedAlignmentPerc = ((double)unsplicedAlignmentNum/(double)totalAlignmentNumber) * 100;
		double singleSplicedAlignmentPerc = ((double)singleSplicedAlignmentNum/(double)totalAlignmentNumber) * 100;
		double multiSplicedAlignmentPerc = ((double)multiSplicedAlignmentNum/(double)totalAlignmentNumber) * 100;
		alignmentTypeNum_ofs << unsplicedAlignmentNum << "\t" << alignerNameVec[tmpAligner] << "\tunspliced" << endl;
		alignmentTypePerc_ofs << unsplicedAlignmentPerc << "\t" << alignerNameVec[tmpAligner] << "\tunspliced" << endl;
		alignmentTypeNum_ofs << singleSplicedAlignmentNum << "\t" << alignerNameVec[tmpAligner] << "\tsingleJunction" << endl;
		alignmentTypePerc_ofs << singleSplicedAlignmentPerc << "\t" << alignerNameVec[tmpAligner] << "\tsingleJunction" << endl;
		alignmentTypeNum_ofs << multiSplicedAlignmentNum << "\t" << alignerNameVec[tmpAligner] << "\tmultiJunction" << endl;
		alignmentTypePerc_ofs << multiSplicedAlignmentPerc << "\t" << alignerNameVec[tmpAligner] << "\tmultiJunction" << endl;		
	}

	return 0;
}
