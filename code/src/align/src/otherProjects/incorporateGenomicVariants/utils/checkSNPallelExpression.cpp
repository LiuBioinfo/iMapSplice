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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../general/SNPhash_info.h"

using namespace std;

// inline char getCharRevComp(char ch)
// {
// 	int chInt = ch - 'A';
// 	static const char alphatChar[26] = {'T', 'N', 'G', 'N', 'N', 'N', 'C',
// 		'N', 'N', 'N', 'N', 'N', 'N', 'N',
// 		'N', 'N', 'N', 'N', 'N', 'A',
// 		'N', 'N', 'N', 'N', 'N', 'N'};
// 	return alphatChar[chInt];
// }

// string getRcmSeq(const string& readSeq)
// {
// 	int readSeqLength = readSeq.length();

// 	char readRcmSeqChar[readSeqLength];

// 	readRcmSeqChar[0] = getCharRevComp(readSeq.at(readSeqLength-1));

// 	for(int tmp = 1; tmp < readSeqLength; tmp ++)
// 	{
// 		readRcmSeqChar[tmp] = getCharRevComp((readSeq.at(readSeqLength - tmp - 1)));
// 	}
// 	string rcmSeq = readRcmSeqChar;
// 	return rcmSeq.substr(0, readSeqLength);
// }

int getEndLocInReadOfSpecificJumpCode(
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	if(jumpCodeIndex < 0)
		return 0;
	int tmpEndLocInRead = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "M")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "I")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "D")
		{
		}
		else if(tmpJumpCodeType == "N")
		{
		}
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return tmpEndLocInRead;
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

void generateExonLocInReadPosInChr(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	vector<int>& endLocVecInRead, vector<int>& endPosVecInChr, vector<int>& lenVec)
{
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp ++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
		int tmpJumpCodeIndex = tmp;
		if(tmpJumpCodeType == "M")
		{
			int tmpEndLocInRead = getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmpJumpCodeIndex);
			int tmpEndPosInChr = getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec, tmpJumpCodeIndex);
			endLocVecInRead.push_back(tmpEndLocInRead);
			endPosVecInChr.push_back(tmpEndPosInChr);
			lenVec.push_back(tmpJumpCodeLength);
		}
	}
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


bool parseSam2chrNamePosCigarString(string& tmpSamStr, int& chrNameInt, int& chrMapPos, 
	string& cigarString, string& readSeq, Index_Info* indexInfo, bool BeersSamOrNot)
{
	int tabLoc_1 = tmpSamStr.find("\t");
	int tabLoc_2 = tmpSamStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSamStr.find("\t", tabLoc_2 + 1);		
	int tabLoc_4 = tmpSamStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpSamStr.find("\t", tabLoc_4 + 1);
	int tabLoc_6 = tmpSamStr.find("\t", tabLoc_5 + 1);
	int tabLoc_7 = tmpSamStr.find("\t", tabLoc_6 + 1);		
	int tabLoc_8 = tmpSamStr.find("\t", tabLoc_7 + 1);
	int tabLoc_9 = tmpSamStr.find("\t", tabLoc_8 + 1);
	int tabLoc_10 = tmpSamStr.find("\t", tabLoc_9 + 1);
	string tmpChrName = tmpSamStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	chrNameInt = indexInfo->convertStringToInt(tmpChrName);
	if(chrNameInt < 0)
		return false;
	string tmpChrPosStr = tmpSamStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	chrMapPos = atoi(tmpChrPosStr.c_str());
	if(chrMapPos < 0)
		return false;
	cigarString = tmpSamStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
	if(cigarString == "*")
		return false;
	readSeq = tmpSamStr.substr(tabLoc_9 + 1, tabLoc_10 - tabLoc_9 - 1);
	if(BeersSamOrNot)
	{
		string tmpFlagStr = tmpSamStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		if(tmpFlagStr == "0")
		{}
		else if(tmpFlagStr == "16")
		{
			readSeq = getRcmSeq(readSeq);
		}
		else
		{
			cout << "error ! flag in Beers Sam: " << tmpFlagStr << endl;
			cout << "exiting ......" << endl;
			exit(1);
		}
	}
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputSNP inputSAM outputFile BeersSamOrNot" << endl;
		exit(1);
	}
	bool BeersSamOrNot = false;
	string BeersSamOrNotStr = argv[5];
	if((BeersSamOrNotStr == "true")||(BeersSamOrNotStr == "TRUE")||(BeersSamOrNotStr == "True"))
		BeersSamOrNot = true;
	else if((BeersSamOrNotStr == "false")||(BeersSamOrNotStr == "FALSE")||(BeersSamOrNotStr == "False"))
		BeersSamOrNot = false;
	else
	{
		cout << "BeersSamOrNot should be set as true or false" << endl;
		cout << "exiting ......" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputSNPpath = argv[2];
	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_formattedSNPfile(inputSNPpath, indexInfo);

	int totalSNPnum = tmpSNPhashInfo.returnSNPnum();
	vector< int > snpRefBaseReadCountVec;
	vector< int > snpAlterBaseReadCountVec;
	vector< int > snpOtherBaseReadCountVec; 
	vector< int > snpReadCountVec;
	for(int tmp = 0; tmp < totalSNPnum; tmp++)
	{
		snpRefBaseReadCountVec.push_back(0);
		snpAlterBaseReadCountVec.push_back(0);
		snpOtherBaseReadCountVec.push_back(0);
		snpReadCountVec.push_back(0);
	}

	string inputSAMpath = argv[3];
	ifstream SAM_ifs(inputSAMpath.c_str());

	int tmpReadNum = 0;
	while(!SAM_ifs.eof())
	{
		string tmpSamStr;
		getline(SAM_ifs, tmpSamStr);
		if(tmpSamStr == "")
			break;
		if(tmpSamStr.substr(0,1) == "@")
			continue;
		tmpReadNum ++;
		int tmpThousandIndex = tmpReadNum / 100000;
		if(tmpReadNum == tmpThousandIndex * 100000)
			cout << "Processed Read #: " << tmpReadNum << endl;
		int chrNameInt, chrMapPos;
		string cigarString, readSeq;
		//cout << "start to parse " << endl;
		bool parseBool = parseSam2chrNamePosCigarString(tmpSamStr, chrNameInt, chrMapPos, cigarString, readSeq, indexInfo, BeersSamOrNot);
		//cout << "parseBool: " << parseBool << endl;
		if(!parseBool)
			continue;
		//cout << "start to do cigarString2jumpCodeVec" << endl;
		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		int chrMapPos_end = getEndPosOfSpecificJumpCode(chrMapPos, cigarStringJumpCodeVec, cigarStringJumpCodeVecSize - 1);
		//cout << "start to do snp search " << endl;
		bool SNPexistsBool = tmpSNPhashInfo.searchSNPvecWithinRegion(chrNameInt, chrMapPos, chrMapPos_end);
		//cout << "SNPexistsBool: " << SNPexistsBool << endl;
		if(!SNPexistsBool)
			continue;
		vector<int> endLocVecInRead;
		vector<int> endPosVecInChr;
		vector<int> lenVec;
		//cout << "start to do generateExonLocInReadPosInChr " << endl;
		generateExonLocInReadPosInChr(chrMapPos, cigarStringJumpCodeVec, endLocVecInRead, endPosVecInChr, lenVec);		
		//cout << "start to check each exon ......" << endl;
		//cout << "tmpSamStr: " << tmpSamStr << endl;
		for(int tmp = 0; tmp < endLocVecInRead.size(); tmp++)
		{
			int tmpEndLocInRead = endLocVecInRead[tmp];
			int tmpEndPosInChr = endPosVecInChr[tmp];
			int tmpLen = lenVec[tmp];
			int tmpStartPosInChr = tmpEndPosInChr - tmpLen + 1;
			// cout << "tmpEndLocInRead: " << tmpEndLocInRead << endl;
			// cout << "tmpEndPosInChr: " << tmpEndPosInChr << endl;
			// cout << "tmpLen: " << tmpLen << endl;
			// cout << "start to do SNPexistsBool_thisExon" << endl;
			bool SNPexistsBool_thisExon = tmpSNPhashInfo.searchSNPvecWithinRegion(chrNameInt, tmpStartPosInChr, tmpEndPosInChr);
			//cout << "SNPexistsBool_thisExon: " << SNPexistsBool_thisExon << endl;
			if(!SNPexistsBool_thisExon)
				continue;
			for(int tmpBase = 0; tmpBase < tmpLen; tmpBase ++)
			{
				int tmpPosInChr = tmpStartPosInChr + tmpBase;
				int tmpLocInRead = tmpEndLocInRead - tmpLen + 1 + tmpBase;
				int tmpSNPindex = tmpSNPhashInfo.searchAndReturnSNPinfoVecIndex(chrNameInt, tmpPosInChr);
				if(tmpSNPindex < 0)
					continue;
				else
				{
					string tmpSNPrefBase = tmpSNPhashInfo.returnSNP_referBase(tmpSNPindex);
					string tmpSNPalterBase = tmpSNPhashInfo.returnSNP_alterBase(tmpSNPindex);					
					string tmpBase = readSeq.substr(tmpLocInRead - 1, 1);
					snpReadCountVec[tmpSNPindex] ++;
					if(tmpBase == tmpSNPrefBase)
						snpRefBaseReadCountVec[tmpSNPindex] ++;
					else if(tmpBase == tmpSNPalterBase)
						snpAlterBaseReadCountVec[tmpSNPindex] ++;
					else
						snpOtherBaseReadCountVec[tmpSNPindex] ++;
				}
			}
		}
	}	
	SAM_ifs.close();

	string outputFilePath = argv[4];
	ofstream ASE_ofs(outputFilePath.c_str());
	for(int tmp = 0; tmp < totalSNPnum; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo.returnSNP_chrPos(tmp);
		string tmpSNP_refBase = tmpSNPhashInfo.returnSNP_referBase(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo.returnSNP_alterBase(tmp);
		int tmpTotalReadCount = snpReadCountVec[tmp];
		int tmpRefReadCount = snpRefBaseReadCountVec[tmp];
		int tmpAlterReadCount = snpAlterBaseReadCountVec[tmp];
		int tmpOtherReadCount = snpOtherBaseReadCountVec[tmp];
		ASE_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_refBase << "\t" 
			<< tmpSNP_alterBase << "\t" << tmpTotalReadCount << "\t" 
			<< tmpRefReadCount << "\t" << tmpAlterReadCount << "\t" 
			<< tmpOtherReadCount << endl;
	}
	ASE_ofs.close();
	return 0;
}