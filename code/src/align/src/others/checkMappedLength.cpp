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
#include "general/read_block_test.h"
#include "general/bwtmap_info.h"
#include "general/DoubleAnchorScore.h"
#include "general/sbndm.h"
#include "general/splice_info.h"
#include "general/checkMappedLengthInfo.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputSAMfile readLength_max outputMappedLengthDistributionInfo" << endl;
		exit(1);
	}

	string inputSAMfileStr = argv[1];
	ifstream inputSAMfile_ifs(inputSAMfileStr.c_str());
	string readLengthMaxStr = argv[2];
	int readLengthMax = atoi(readLengthMaxStr.c_str());
	string outputMappedLengthDistributionInfoStr = argv[3];
	ofstream mappedLengthDistributionInfo_ofs(outputMappedLengthDistributionInfoStr.c_str());

	/// count primary alignment number for each mapped length
	unsigned int* mappedLength_num = (unsigned int*)malloc(readLengthMax* sizeof(unsigned int));
	for(int tmp = 0; tmp < readLengthMax; tmp++)
		mappedLength_num[tmp] = 0;

	unsigned int unmapAlignNum = 0;
	unsigned int totalAlignNum = 0;
	unsigned int primaryAlignNum = 0;
	unsigned int secondaryAlignNum = 0;
	while(1)
	{
		if(inputSAMfile_ifs.eof())
			break;
		string tmpAlignStr;
		getline(inputSAMfile_ifs, tmpAlignStr);
		if(inputSAMfile_ifs.eof())
			break;
		if(tmpAlignStr.at(0) == '@')
			continue;
		totalAlignNum ++;
		CheckMappedLengthInfo* newMapLenInfo = new CheckMappedLengthInfo();
		int alignType = newMapLenInfo->checkMappedLength(tmpAlignStr);
		if(alignType == 0)
		{
			int tmpMappedLength = newMapLenInfo->returnMappedLength();
			mappedLength_num[tmpMappedLength] = mappedLength_num[tmpMappedLength] + 1;

			primaryAlignNum ++;
			totalAlignNum ++;
		}
		else if(alignType == 1)
		{
			int tmpMappedLength = 0;
			mappedLength_num[tmpMappedLength] = mappedLength_num[tmpMappedLength] + 1;

			unmapAlignNum ++;
			primaryAlignNum ++;
			totalAlignNum ++;
		}
		else if(alignType == 2)
		{
			totalAlignNum ++;
			secondaryAlignNum ++;
		}
		else
		{
			cout << "error in alignType" << endl;
			exit(1);
		}
		delete newMapLenInfo;
	}


	/// output mappedLength_num ///
	mappedLengthDistributionInfo_ofs << "readLength\tmapped_primaryAlignment_num" << endl;
	unsigned tmpTotalPrimaryAlignNum = 0;
	for(int tmp = 0; tmp < readLengthMax; tmp++)
	{
		tmpTotalPrimaryAlignNum += mappedLength_num[tmp];
		mappedLengthDistributionInfo_ofs << tmp+1 << "\t" << mappedLength_num[tmp] << endl;
	}
	mappedLengthDistributionInfo_ofs << "total primaryAlignNum: " << tmpTotalPrimaryAlignNum << endl;
	mappedLengthDistributionInfo_ofs << endl << endl << "************************************" << endl << endl;
	mappedLengthDistributionInfo_ofs << "readLength\tmapped_primaryAlignment_perc" << endl;
	for(int tmp = 0; tmp < readLengthMax; tmp++)
	{
		double tmpPerc = (mappedLength_num[tmp]*100)/tmpTotalPrimaryAlignNum;
		mappedLengthDistributionInfo_ofs << tmp+1 << "\t" << tmpPerc << endl;
	}
	//free(mappedLength_perc);
	free(mappedLength_num);
	mappedLengthDistributionInfo_ofs.close();
	inputSAMfile_ifs.close();
	return 0;
}