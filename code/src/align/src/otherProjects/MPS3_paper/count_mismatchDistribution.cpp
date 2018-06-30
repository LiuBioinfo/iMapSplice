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


//#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"

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
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "M")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "D")
		{}
		else if(tmpJumpCodeType == "N")
		{}
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

void generateMismatchPosVec(vector<int>& tmpMismatchPosVec, vector<int>& tmpMismatchBaseQualVec, int mapChrNameInt, 
	int mapChrPos, string& cigarString, string& readSeq, string& qualSeq, int readLength_max, 
	Index_Info* indexInfo, bool phred33or64_bool)
{
	vector<Jump_Code> cigarStringJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
	int jumpCodeVecSize = cigarStringJumpCodeVec.size();

	vector<int> mismatchPosVec_raw;
	for(int tmp = 0; tmp < jumpCodeVecSize; tmp++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		if(tmpJumpCodeType != "M")
			continue;
		int tmpEndLocInRead = getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmp);
		int tmpEndPosInChr = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, tmp);
		int tmpSegLen = cigarStringJumpCodeVec[tmp].len;
		int tmpStartLocInRead = tmpEndLocInRead - tmpSegLen + 1;
		int tmpStartPosInChr = tmpEndPosInChr - tmpSegLen + 1;

		string tmpSeqInChr = indexInfo->returnChromStrSubstr(mapChrNameInt, tmpStartPosInChr, tmpSegLen);
		string tmpSeqInRead = readSeq.substr(tmpStartLocInRead - 1, tmpSegLen);
		for(int tmp2 = 0; tmp2 < tmpSegLen; tmp2++)
		{
			if(tmpSeqInChr.at(tmp2) != tmpSeqInRead.at(tmp2))
			{
				int tmpMismatchBasePosInRead = tmp2 + tmpStartLocInRead;
				mismatchPosVec_raw.push_back(tmpMismatchBasePosInRead);
				char tmpMismatchBaseQualInQualSeq = qualSeq.at(tmpMismatchBasePosInRead - 1);
				int tmpMismatchBaseQualQint;
				if(phred33or64_bool)
					tmpMismatchBaseQualQint = tmpMismatchBaseQualInQualSeq - '!';
				else
					tmpMismatchBaseQualQint = tmpMismatchBaseQualInQualSeq - '@';
				tmpMismatchBaseQualVec.push_back(tmpMismatchBaseQualQint);
			}
		}
	}

	int tmpReadSeqLength = readSeq.length();
	if(tmpReadSeqLength == readLength_max)
	{
		for(int tmp = 0; tmp < mismatchPosVec_raw.size(); tmp++)
			tmpMismatchPosVec.push_back(mismatchPosVec_raw[tmp]);
	}
	else
	{
		for(int tmp = 0; tmp < mismatchPosVec_raw.size(); tmp++)
		{
			int tmpMismatchPos_raw = mismatchPosVec_raw[tmp];
			int tmpDistance_toReadStart = tmpMismatchPos_raw - 1;
			int tmpDistance_toReadEnd = readLength_max - tmpMismatchPos_raw;
			if(tmpDistance_toReadStart <= tmpDistance_toReadEnd)
				tmpMismatchPosVec.push_back(tmpMismatchPos_raw);
			else
				tmpMismatchPosVec.push_back(tmpMismatchPos_raw + readLength_max - tmpReadSeqLength);
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputIndexFolderPath readLengthMax inputSAM outputFolder thread_num toolName phred33or64" << endl;
		exit(1);
	}
	string toolName = argv[6];
	string thread_num_str = argv[5];
	int thread_num = atoi(thread_num_str.c_str());
	omp_set_num_threads(thread_num);

	string phred33or64_str = argv[7];
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
	string readLengthStr = argv[2];
	int readLength_max = atoi(readLengthStr.c_str());
	cout << "start to initiate mismatchNumVecVecs......" << endl;
	vector< vector<int> > mismatchNumInEachReadVecVec;
	vector< vector<int> > mismatchNumInEachBasePosVecVec;
	vector< vector<int> > mismatchNumWithEachQualVecVec;
	int qualNum_max = 43;
	for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
	{	
		vector<int> tmpMismatchNumInEachBasePosVec;
		for(int tmp = 0; tmp < readLength_max; tmp++)
			tmpMismatchNumInEachBasePosVec.push_back(0);
		mismatchNumInEachBasePosVecVec.push_back(tmpMismatchNumInEachBasePosVec);
		vector<int> tmpMismatchNumInEachReadVec;
		for(int tmp = 0; tmp < readLength_max; tmp++)
			tmpMismatchNumInEachReadVec.push_back(0);
		mismatchNumInEachReadVecVec.push_back(tmpMismatchNumInEachReadVec);
		vector<int> tmpMismatchNumWithEachQualVec;
		for(int tmp = 0; tmp < qualNum_max; tmp++)
			tmpMismatchNumWithEachQualVec.push_back(0);
		mismatchNumWithEachQualVecVec.push_back(tmpMismatchNumWithEachQualVec);
	}

	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	log_ofs << "start to initiate indexInfo" << endl;
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;

	cout << "start to process SAM file ......" << endl;
	log_ofs << "start to process SAM file ......" << endl;

	int primary_total_alignmentNum = 0;

	string inputSAM = argv[3];
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
		log_ofs << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
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
		log_ofs << "start to counting mismatches in SAM, turn: " << tmpTurn + 1 << endl;
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
			string mapChrNameStr = samFieldVec[2];
			int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
			string mapChrPosStr = samFieldVec[3];
			int mapChrPos = atoi(mapChrPosStr.c_str());
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
						vector<int> tmpMismatchPosVec;
						vector<int> tmpMismatchBaseQualVec;
						generateMismatchPosVec(tmpMismatchPosVec, tmpMismatchBaseQualVec, mapChrNameInt, 
							mapChrPos, cigarString, readSeq, qualSeq, readLength_max, indexInfo, phred33or64_bool);
						int tmpMismatchPosVecSize = tmpMismatchPosVec.size();
						for(int tmp2 = 0; tmp2 < tmpMismatchPosVecSize; tmp2++)
						{
							int tmpMismatchPos = tmpMismatchPosVec[tmp2];
							int tmpMismatchNumInPos_ori = (mismatchNumInEachBasePosVecVec[threadNO])[tmpMismatchPos-1];
							int tmpMismatchNumInPos_updated = tmpMismatchNumInPos_ori + 1;
							(mismatchNumInEachBasePosVecVec[threadNO])[tmpMismatchPos-1] = tmpMismatchNumInPos_updated;
							int tmpMismatchBaseQual = tmpMismatchBaseQualVec[tmp2];
							int tmpMismatchNumWithQual_ori = (mismatchNumWithEachQualVecVec[threadNO])[tmpMismatchBaseQual];
							int tmpMismatchNumWithQual_updated = tmpMismatchNumWithQual_ori + 1;
							(mismatchNumWithEachQualVecVec[threadNO])[tmpMismatchBaseQual] = tmpMismatchNumWithQual_updated;
						}
						(mismatchNumInEachReadVecVec[threadNO])[tmpMismatchPosVecSize] ++;
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
		log_ofs << "end of counting mismatches in SAM ... turn: " << tmpTurn+1 << endl;
		log_ofs << "# of reads have been fixed: " << readTotalNum << endl;
	}

	string mismatchNumInEachPos_file = outputFolderStr + "mismatchNumInEachBase.txt";
	ofstream mismatchNum_ofs(mismatchNumInEachPos_file.c_str());
	for(int tmp = 0; tmp < readLength_max; tmp++)
	{
		int tmpBasePos = tmp + 1;
		int tmpBaseMismatchNum = 0;
		for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
			tmpBaseMismatchNum += (mismatchNumInEachBasePosVecVec[tmpThread])[tmp];
		mismatchNum_ofs << tmpBasePos << "\t" << tmpBaseMismatchNum << "\t" << toolName << endl;
	}
	mismatchNum_ofs.close();
	
	string mismatchNumWithEachQual_file = outputFolderStr + "mismatchNumWithEachQual.txt";
	ofstream mismatchNumWithEachQual_ofs(mismatchNumWithEachQual_file.c_str());
	for(int tmp = 0; tmp < qualNum_max; tmp++)
	{
		int tmpMismatchNumWithThisQual = 0;
		for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
			tmpMismatchNumWithThisQual += (mismatchNumWithEachQualVecVec[tmpThread])[tmp];
		mismatchNumWithEachQual_ofs << tmp << "\t" << tmpMismatchNumWithThisQual << "\t" << toolName << endl;
	}
	mismatchNumWithEachQual_ofs.close();

	string mismatchNumInEachRead_file = outputFolderStr + "mismatchNumInEachRead.txt";
	ofstream mismatchNum_inEachRead_ofs(mismatchNumInEachRead_file.c_str());
	for(int tmp = 0; tmp < readLength_max; tmp++)
	{
		int tmpReadNum = 0;
		for(int tmpThread = 0; tmpThread < thread_num; tmpThread++)
			tmpReadNum += (mismatchNumInEachReadVecVec[tmpThread])[tmp];
		mismatchNum_inEachRead_ofs << tmp << "\t" << tmpReadNum << "\t" << toolName << endl;
	}
	mismatchNum_inEachRead_ofs.close();
	sam_ifs.close();
	free(chrom);
	delete indexInfo;
	log_ofs.close();
	return 0;
}