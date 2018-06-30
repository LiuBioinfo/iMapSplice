// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to count aligment type number 
//(total unsplied, single-spliced, multi-spliced) for primary ones
//used to count soft clipped bases for primary ones
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
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"

using namespace std;

bool bothEndsMapped(int tmpFlag)
{
	if(tmpFlag & 0x2)
		return true;
	else
		return false;
}

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
	if(argc != 6)
	{
		cout << "Executable TotalReadNum ReadLength SamFile outputFolder methodName" << endl;
		exit(1);
	}
	string methodName = argv[5];
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	log_ofs << "start to generate files" << endl;
	string outputAlignmentTypeNumPercFile = outputFolderStr + "alignmentTypeNumPerc.txt";
	ofstream alignmentTypeNumPerc_ofs(outputAlignmentTypeNumPercFile.c_str());
	string outputUnmappedBaseNumPercFile = outputFolderStr + "UnmappedBaseNumPerc.txt";
	ofstream unmappedBaseNumPerc_ofs(outputUnmappedBaseNumPercFile.c_str());

	string totalAlignmentNumberStr = argv[1];
	int totalAlignmentNumber = atoi(totalAlignmentNumberStr.c_str());
	cout << "total alignments number: " << totalAlignmentNumber << endl;	
	log_ofs << "total alignments number: " << totalAlignmentNumber << endl;	
	string fixedReadLengthStr = argv[2];
	int fixedReadLength = atoi(fixedReadLengthStr.c_str());
	cout << "fixedReadLength: " << fixedReadLength << endl;
	log_ofs << "fixedReadLength: " << fixedReadLength << endl;
	
	string samPath = argv[3];
	log_ofs << "samPath: " << samPath << endl;
	int primary_total_alignmentNum = 0;
	int primary_unspliced_alignmentNum = 0;
	int primary_singleSpliced_alignmentNum = 0;
	int primary_multiSpliced_alignmentNum = 0;
	vector<int> mappedBaseNumVec_primaryAlignment;
	for(int tmp = 0; tmp < fixedReadLength; tmp++)
		mappedBaseNumVec_primaryAlignment.push_back(0);
	ifstream sam_ifs(samPath.c_str());
	while(!(sam_ifs.eof()))
	{
		string samStr;
		getline(sam_ifs, samStr);
		 if(sam_ifs.eof()||(samStr == ""))
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
		bool bothEndsMapped_bool = bothEndsMapped(tmpFlag);
		#ifdef INCLUDE_UNPAIRED
		bothEndsMapped_bool = true;
		#endif
		if(!mappedOrNot_bool)
		{}
		else
		{
			if(bothEndsMapped_bool)
			{
				if(primaryOrNot_bool)
				{
					primary_total_alignmentNum ++;
					int tmpAlignmentSJnum = getSJnumFromCigarString(cigarString);
					if(tmpAlignmentSJnum == 0)
						primary_unspliced_alignmentNum ++;
					else if(tmpAlignmentSJnum == 1)
						primary_singleSpliced_alignmentNum ++;
					else
						primary_multiSpliced_alignmentNum ++;
					int softClippedHeadLength = getSoftClippedHeadLengthFromCigarString(cigarString);
					int softClippedTailLength = getSoftClippedTailLengthFromCigarString(cigarString);
					int startMappedBaseLocInRead = softClippedHeadLength + 1;
					int endMappedBaseLocInRead = fixedReadLength - softClippedTailLength;
					for(int tmpBaseLoc = startMappedBaseLocInRead; 
							tmpBaseLoc <= endMappedBaseLocInRead; tmpBaseLoc++)
						mappedBaseNumVec_primaryAlignment[tmpBaseLoc-1] ++;
				}
				else
				{}
			}
		}
	}

	int unmappedAlignmentNum = totalAlignmentNumber - primary_total_alignmentNum;
	double primary_total_alignmentPerc = ((double)primary_total_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unspliced_alignmentPerc = ((double)primary_unspliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_singleSpliced_alignmentPerc = ((double)primary_singleSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_multiSpliced_alignmentPerc = ((double)primary_multiSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	alignmentTypeNumPerc_ofs << "Total_readNum:\t" << totalAlignmentNumber << endl;
	alignmentTypeNumPerc_ofs << "primary_total_alignmentNumPerc:\t" 
		<< primary_total_alignmentNum << "\t" << primary_total_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_unspliced_alignmentNumPerc:\t" 
		<< primary_unspliced_alignmentNum << "\t" << primary_unspliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_singleSpliced_alignmentNumPerc:\t" 
		<< primary_singleSpliced_alignmentNum << "\t" << primary_singleSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_multiSpliced_alignmentNumPerc:\t" 
		<< primary_multiSpliced_alignmentNum << "\t" << primary_multiSpliced_alignmentPerc << endl;
	for(int tmp = 0; tmp < fixedReadLength; tmp++)
	{
		int tmpBaseLoc = tmp + 1;
		int tmpMappedBaseNum = mappedBaseNumVec_primaryAlignment[tmp];
		int tmpUnmappedBaseNum = totalAlignmentNumber - mappedBaseNumVec_primaryAlignment[tmp];
		double tmpUnmappedBasePerc = ((double)tmpUnmappedBaseNum/(double)totalAlignmentNumber) * 100;
		unmappedBaseNumPerc_ofs << "Base[" << tmpBaseLoc << "]:\t" 
			<< tmpUnmappedBaseNum << "\t" << tmpUnmappedBasePerc << "\t" << methodName << endl;
	}			
	sam_ifs.close();
	alignmentTypeNumPerc_ofs.close();
	unmappedBaseNumPerc_ofs.close();
	return 0;
}