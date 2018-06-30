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

#include "../../../general/read_block_test.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputUniqueBreakPointFusionFastaPath outputReadsPathPrefix" << endl;
		exit(1);
	}
	int skippingStepNum = 5;
	int insertSize = 100;
	int readLength = 100;

	string inputUniqueBreakPointFusionFastaPath = argv[1];
	string outputReadsPathPrefix = argv[2];
	string outputReads_end1 = outputReadsPathPrefix + ".1.fa";
	string outputReads_end2 = outputReadsPathPrefix + ".2.fa";
	ifstream fusionFa_ifs(inputUniqueBreakPointFusionFastaPath.c_str());
	ofstream reads_end1_ofs(outputReads_end1.c_str());
	ofstream reads_end2_ofs(outputReads_end2.c_str());
	while(!fusionFa_ifs.eof())
	{
		string tmpFusionName;
		getline(fusionFa_ifs, tmpFusionName);
		cout << "tmpFusionName: " << tmpFusionName << endl;
		string tmpFusionSeq;
		getline(fusionFa_ifs, tmpFusionSeq);
		int tmpFusionSeqLength = tmpFusionSeq.length();
		int tmpReadStartLoc_1 = 1;
		int tmpReadEndLoc_1, tmpReadStartLoc_2, tmpReadEndLoc_2;
		while(1)
		{
			tmpReadEndLoc_1 = tmpReadStartLoc_1 + readLength - 1;
			tmpReadStartLoc_2 = tmpReadStartLoc_1 + readLength + insertSize;
			tmpReadEndLoc_2 = tmpReadStartLoc_2 + readLength - 1;
			if(tmpReadEndLoc_2 > tmpFusionSeqLength)
				break;
			else
			{
				string tmpReadName_1 = tmpFusionName 
					+ "_" + int_to_str(tmpFusionSeqLength)
					+ "_" + int_to_str(tmpReadStartLoc_1)
					+ "_" + int_to_str(tmpReadEndLoc_1) 
					+ "_" + int_to_str(tmpReadStartLoc_2)
					+ "_" + int_to_str(tmpReadEndLoc_2)  
					+ "/1";
				string tmpReadName_2 = tmpFusionName 
					+ "_" + int_to_str(tmpFusionSeqLength)
					+ "_" + int_to_str(tmpReadStartLoc_1)
					+ "_" + int_to_str(tmpReadEndLoc_1) 
					+ "_" + int_to_str(tmpReadStartLoc_2)
					+ "_" + int_to_str(tmpReadEndLoc_2)  
					+ "/2";
				string tmpReadSeq_1 = tmpFusionSeq.substr(
					tmpReadStartLoc_1 - 1, readLength);
				string tmpReadSeq_2 = tmpFusionSeq.substr(
					tmpReadStartLoc_2 - 1, readLength);
				reads_end1_ofs << tmpReadName_1 << endl 
					<< tmpReadSeq_1 << endl;
				reads_end2_ofs << tmpReadName_2 << endl
					<< tmpReadSeq_2 << endl;
			}
			tmpReadStartLoc_1 += 5;
		}
	}
	fusionFa_ifs.close();
	reads_end1_ofs.close();
	reads_end2_ofs.close();
	return 0;
}