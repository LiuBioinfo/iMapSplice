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
#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputSupNumThresFile inputFusionInfoFile outputFilteredFusionInfoFilePrefix" << endl;
		exit(1);
	}
	string inputSupNumThresFile = argv[1];
	string inputFusionInfoFile = argv[2];
	string outputFilteredFusionInfoFilePrefix = argv[3];
	string kept_path = outputFilteredFusionInfoFilePrefix + "_kept.txt";
	string filteredOut_path = outputFilteredFusionInfoFilePrefix + "_filteredOut.txt";

	ifstream supNumThres_ifs(inputSupNumThresFile.c_str());
	ifstream ori_ifs(inputFusionInfoFile.c_str());
	ofstream kept_ofs(kept_path.c_str());
	ofstream filteredOut_ofs(filteredOut_path.c_str());
	vector<int> supNumThresVec; // 0~6 -- for canonical fusion; 7~13 -- for noncanonical fusion;
	while(!supNumThres_ifs.eof())
	{
		string tmpSupNumThresStr;
		getline(supNumThres_ifs, tmpSupNumThresStr);
		if(tmpSupNumThresStr == "")
			break;
		int tmpSupNumThres = atoi(tmpSupNumThresStr.c_str());
		supNumThresVec.push_back(tmpSupNumThres);
	}
	supNumThres_ifs.close();
	while(!ori_ifs.eof())
	{
		string tmpFusStr;
		getline(ori_ifs, tmpFusStr);
		if(tmpFusStr == "")
			break;
		vector<string> tmpFusFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 16; tmp++)
		{
			int tabLoc = tmpFusStr.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpFusField = tmpFusStr.substr(startLoc, tabLoc-startLoc);
			tmpFusFieldVec.push_back(tmpFusField);
			startLoc = tabLoc + 1;
		}
		string tmpStrand_1 = tmpFusFieldVec[4];
		string tmpStrand_2 = tmpFusFieldVec[5];
		int tmpGlobalSplit = atoi(tmpFusFieldVec[9].c_str());
		int tmpRemapSplit = atoi(tmpFusFieldVec[10].c_str());
		int tmpCompletePairEncomp = atoi(tmpFusFieldVec[11].c_str());
		int tmpIncompletePairEncomp = atoi(tmpFusFieldVec[12].c_str());
		int tmpSplit = tmpGlobalSplit + tmpRemapSplit;
		int tmpEncomp = tmpCompletePairEncomp + tmpIncompletePairEncomp;
		int tmpTotal = tmpSplit + tmpEncomp;
		if(tmpStrand_1 + tmpStrand_2 == "NN") // noncanonical
		{
			if((tmpGlobalSplit >= supNumThresVec[7])&&(tmpRemapSplit >= supNumThresVec[8])
				&&(tmpCompletePairEncomp >= supNumThresVec[9])&&(tmpIncompletePairEncomp >= supNumThresVec[10])
				&&(tmpSplit >= supNumThresVec[11])&&(tmpEncomp >= supNumThresVec[12])
				&&(tmpTotal >= supNumThresVec[13]))
				kept_ofs << tmpFusStr << endl;
			else
				filteredOut_ofs << tmpFusStr << endl;
		}
		else // canonical
		{
			if((tmpGlobalSplit >= supNumThresVec[0])&&(tmpRemapSplit >= supNumThresVec[1])
				&&(tmpCompletePairEncomp >= supNumThresVec[2])&&(tmpIncompletePairEncomp >= supNumThresVec[3])
				&&(tmpSplit >= supNumThresVec[4])&&(tmpEncomp >= supNumThresVec[5])
				&&(tmpTotal >= supNumThresVec[6]))
				kept_ofs << tmpFusStr << endl;
			else
				filteredOut_ofs << tmpFusStr << endl;
		}
	}
	ori_ifs.close();
	kept_ofs.close();
	filteredOut_ofs.close();
	return 0;
}	