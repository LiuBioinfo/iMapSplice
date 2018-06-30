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
	cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);	
	if(tmpJumpCodeVec[0].type == "S")
		return tmpJumpCodeVec[0].len;
	else
		return 0; 
}	

int getSoftClippedTailLengthFromCigarString(string& tmpCigarStringStr)
{
	vector<Jump_Code> tmpJumpCodeVec;
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
		cout << "Executable inputSAM outputFile thread_num toolName phred33or64" << endl;
		exit(1);
	}
	string toolName = argv[4];
	string thread_num_str = argv[3];
	int thread_num = atoi(thread_num_str.c_str());
	omp_set_num_threads(thread_num);

	string phred33or64_str = argv[5];
	bool phred33or64_bool = true;
	if((phred33or64_str == "phred33")||(phred33or64_str == "Phred33"))
		phred33or64_bool = true;
	else if((phred33or64_str == "phred64")||(phred33or64_str == "Phred64"))
		phred33or64_bool = false;
	else
	{
		cout << "invalid phred33or64 format ! " << endl; 
		exit(1);
	}
	cout << "start to initiate mappedBaseNumWithEachQualVecVec......" << endl;
	vector< vector<long long> > mappedBaseNumWithEachQualVecVec;
	int qualNum_max = 43;
	for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
	{	
		vector<long long> tmpMappedBaseNumWithEachQualVec;
		for(int tmp = 0; tmp < qualNum_max; tmp++)
			tmpMappedBaseNumWithEachQualVec.push_back(0);
		mappedBaseNumWithEachQualVecVec.push_back(tmpMappedBaseNumWithEachQualVec);
	}

	cout << "start to process SAM file ......" << endl;
	int primary_total_alignmentNum = 0;

	string inputSAM = argv[1];
	ifstream sam_ifs(inputSAM.c_str());

	bool EndOfRecord = false;
	int realRecordNum;
	int readPairNum = 0;
	int readTotalNum = 0;

	int normalRecordNum = 500000;
	vector<string> samStrVec(normalRecordNum);
	for(int tmpTurn = 0; ; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		//log_ofs << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
			if(sam_ifs.eof())
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    				
			}
			string samStr;
			getline(sam_ifs, samStr);
			if(samStr == "")
			{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    				
			}
		 	samStrVec[recordNumTmp] = samStr;
		}		
		readTotalNum += realRecordNum;

		omp_set_num_threads(thread_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			string samStr;
			samStr = samStrVec[tmpOpenMP];
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
						primary_total_alignmentNum ++;
						int tmpReadSeqLength = readSeq.length();
						int softClippedHeadLength = getSoftClippedHeadLengthFromCigarString(cigarString);
						int softClippedTailLength = getSoftClippedTailLengthFromCigarString(cigarString);
						int startMappedBaseLocInRead = softClippedHeadLength + 1;
						int endMappedBaseLocInRead = tmpReadSeqLength - softClippedTailLength;
						for(int tmpBaseLoc = startMappedBaseLocInRead; 
							tmpBaseLoc <= endMappedBaseLocInRead; tmpBaseLoc++)
						{
							char tmpBaseQualInQualSeq = qualSeq.at(tmpBaseLoc - 1);
							int tmpBaseQualQint;
							if(phred33or64_bool)
								tmpBaseQualQint = tmpBaseQualInQualSeq - '!';
							else
								tmpBaseQualQint = tmpBaseQualInQualSeq - '@';						
							(mappedBaseNumWithEachQualVecVec[threadNO])[tmpBaseQualQint] ++;
						}
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
	}
	string outputFileStr = argv[2];
	ofstream baseQual_ofs(outputFileStr.c_str());
	for(int tmp = 0; tmp < qualNum_max; tmp++)
	{
		long long int tmpBaseNumWithThisQual = 0;
		for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
			tmpBaseNumWithThisQual += (mappedBaseNumWithEachQualVecVec[tmpThread])[tmp];
		baseQual_ofs << tmp << "\t" << tmpBaseNumWithThisQual << "\t" << toolName << endl;
	}	
	baseQual_ofs.close();
	sam_ifs.close();
	return 0;
}