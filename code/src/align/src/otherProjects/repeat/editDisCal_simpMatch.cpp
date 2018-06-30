// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "executable inputIndexFolderPath inputQuerySeq outputResutlsFile" << endl;
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
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	string querySeq = argv[2];
	string outputPath = argv[3];
	ofstream seqSim_ofs(outputPath.c_str());

	int querySeqLen = querySeq.length();
	int toCompareChrNameInt = 0;
	int startToCompareChrPos = 1000;
	int toCompareChrLength = indexInfo->returnChromLength(toCompareChrNameInt);
	int endToCompareChrPos = toCompareChrLength - 1000;

	//int tmpCheckNum = 0;
	int baseNum_differentThanThePreviousOne = 0;
	int baseNum_theSameAsThePreviousOne = 0;
	for(int tmpPos = startToCompareChrPos; tmpPos <= endToCompareChrPos; tmpPos++)
	{
		// tmpCheckNum ++;
		// string tmpChrSubStr = indexInfo->returnChromStrSubstr(toCompareChrNameInt,
		// 	tmpPos, querySeqLen);
		// int tmpMisMatchNum = 0;
		// for(int tmp = 0; tmp < querySeqLen; tmp++)
		// {
		// 	string tmpBaseInQuery = querySeq.substr(tmp, 1);
		// 	string tmpBaseInChr = tmpChrSubStr.substr(tmp, 1);
		// 	if(tmpBaseInQuery != tmpBaseInChr)
		// 		tmpMisMatchNum ++;
		// }
		// seqSim_ofs << tmpPos << "\t" << tmpMisMatchNum << endl;
		// int tmpHundredThousandNumIndex = tmpCheckNum / 100000;
		// if(tmpHundredThousandNumIndex * 100000 == tmpCheckNum)
		// 	cout << tmpCheckNum << " seqs processed" << endl;
		string lastBaseInChr = indexInfo->returnChromStrSubstr(toCompareChrNameInt,
			tmpPos-1, 1);
		string thisBaseInChr = indexInfo->returnChromStrSubstr(toCompareChrNameInt,
			tmpPos, 1);
		if(lastBaseInChr == thisBaseInChr)
			baseNum_theSameAsThePreviousOne ++;
		else
		{	
			baseNum_differentThanThePreviousOne ++;
			seqSim_ofs << thisBaseInChr;
		}
	}
	seqSim_ofs << endl;
	cout << "baseNum_theSameAsThePreviousOne" << baseNum_theSameAsThePreviousOne << endl;
	cout << "baseNum_differentThanThePreviousOne" << baseNum_differentThanThePreviousOne << endl;
	cout << "checkedBaseNumSum: " << baseNum_theSameAsThePreviousOne 
		+ baseNum_differentThanThePreviousOne << endl;
	cout << "endToCompareChrPos - startToCompareChrPos: " << endToCompareChrPos 
		- startToCompareChrPos << endl;
	return 0;
}