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

#include "../constantDefinitions.h"
#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"

using namespace std;

void generateSAMfieldInfo(
	const string& tmpSamStr,
	int tmpSAM_chrNameInt, int tmpSAM_startPos, 
	vector<Jump_Code>& tmpSAM_cigarStringJumpCodeVec,
	const string& tmpSAM_readSeq, vector<int>& tmpSAM_mismatchPosVec, 
	vector<char>& tmpSAM_mismatchCharVec, Index_Info* indexInfo)
{

}

void compareReadSeq2ref_generateMismatchVec(
	int tmpSAM_chrNameInt, int tmpSAM_startPos,
	vector<Jump_Code>& tmpSAM_cigarStringJumpCodeVec,
	const string& tmpSAM_readSeq, 
	vector<int>& tmpMismatchPosVec_compareReadSeq2ref,
	vector<char>& tmpMismatchCharVec_compareReadSeq2ref, 
	Index_Info* indexInfo)
{}

bool compareGeneratedMismatchVec(
	vector<int>& tmpSAM_mismatchPosVec,
	vector<char>& tmpSAM_mismatchCharVec,
	vector<int>& tmpMismatchPosVec_compareReadSeq2ref,
	vector<char>& tmpMismatchCharVec_compareReadSeq2ref)
{
	if(tmpSAM_mismatchPosVec.size() != tmpMismatchCharVec_compareReadSeq2ref.size())
		return false;
	else if((tmpSAM_mismatchPosVec.size() != tmpMismatchPosVec_compareReadSeq2ref.size())
		||(tmpSAM_mismatchCharVec.size() != tmpMismatchCharVec_compareReadSeq2ref.size()))
		return false;
	else
	{}
	int mismatchNum = tmpSAM_mismatchPosVec.size();
	set<int> tmpSAM_mismatchPosSet;
	for(int tmp = 0; tmp < mismatchNum; tmp++)
	{
		int tmpMismatchPos = tmpSAM_mismatchPosVec[tmp];
		tmpSAM_mismatchPosSet.insert(tmpMismatchPos);
	}
	int tmpSAM_mismatchPosSet_size = tmpSAM_mismatchPosSet.size();
	if(mismatchNum != tmpSAM_mismatchPosSet_size)
		return false;
	set<int> tmpMismatchPosSet_compareReadSeq2ref;	
	for(int tmp = 0; tmp < mismatchNum; tmp++)
	{
		int tmpMismatchPos = tmpMismatchPosVec_compareReadSeq2ref[tmp];

	}



	for(int tmp = 0; tmp < mismatchNum; tmp++)
	{
		int tmpMismatchPos_compareReadSeq2ref 
			= tmpMismatchPosVec_compareReadSeq2ref[tmp];
	}
}

int main(int argc, char** argv)
{
	if(argc != )
	{
		cout << "Executable inputIndexInfoPath inputSAMpath ";
		cout << " outputInconsistentAlignmentFilePath" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputSAMpath = argv[2];
	string outputInconsistentAlignmentFilePath = argv[3];

	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	cout << "finish loading chromosomes" << endl;

	ifstream sam_ifs(inputSAMpath.c_str());
	ofstream toCheckSAM_ofs(outputInconsistentAlignmentFilePath.c_str());
	while(!sam_ifs.eof())
	{
		string tmpSamStr;
		getline(sam_ifs, tmpSamStr);
		if((sam_ifs.eof())||(tmpSamStr == ""))
			break;
		if(tmpSamStr.at(0) == '@')
			continue;
		bool mapOrNot_bool 
			= mapOrNotWithSamStr(tmpSamStr);
		if(!mapOrNot_bool)
			continue;
		int tmpSAM_chrNameInt;
		int tmpSAM_startPos;
		vector<Jump_Code> tmpSAM_cigarStringJumpCodeVec;
		string tmpSAM_readSeq;
		vector<int> tmpSAM_mismatchPosVec;
		vector<char> tmpSAM_mismatchCharVec;
		generateSAMfieldInfo(tmpSamStr,
			tmpSAM_chrNameInt, tmpSAM_startPos, 
			tmpSAM_cigarStringJumpCodeVec,
			tmpSAM_readSeq, tmpSAM_mismatchPosVec, 
			tmpSAM_mismatchCharVec, indexInfo);
		vector<int> tmpMismatchPosVec_compareReadSeq2ref;
		vector<char> tmpMismatchCharVec_compareReadSeq2ref;
		compareReadSeq2ref_generateMismatchVec(
			tmpSAM_chrNameInt, tmpSAM_startPos,
			tmpSAM_cigarStringJumpCodeVec,
			tmpSAM_readSeq, tmpMismatchPosVec_compareReadSeq2ref,
			tmpMismatchCharVec_compareReadSeq2ref, indexInfo);
		bool mismatchVecConsistentOrNot 
			= compareGeneratedMismatchVec(
				tmpSAM_mismatchPosVec, tmpSAM_mismatchCharVec
				tmpMismatchPosVec_compareReadSeq2ref,
				tmpMismatchCharVec_compareReadSeq2ref);
		if(mismatchVecConsistentOrNot)
		{}
		else
		{
			toCheckSAM_ofs << tmpSamStr << "\t";
			string RDfieldStr = "RD:Z:";
			string RPfieldStr = "RP:Z:";
			for(int tmp = 0; tmp < tmpMismatchPosVec_compareReadSeq2ref.size();
				tmp ++)
			{
				RDfieldStr = RDfieldStr 
					+ int_to_str(tmpMismatchPosVec_compareReadSeq2ref[tmp]) + ",";
				string tmpRPchar = tmpMismatchCharVec_compareReadSeq2ref[tmp];
				RPfieldStr = RPfieldStr
					+ tmpRPchar + ",";
			}
			toCheckSAM_ofs << RDfieldStr << "\t" << RPfieldStr << endl;
		}
	}
	sam_ifs.close();
	return 0;
}