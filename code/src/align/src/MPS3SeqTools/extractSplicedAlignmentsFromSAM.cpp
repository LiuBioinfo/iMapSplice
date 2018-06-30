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

#include "../general/read_block_test.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSAM outputSpliceAlignmentsFilePath" << endl;
		exit(1);
	}
	string inputSAMpath = argv[1];
	ifstream sam_ifs(inputSAMpath.c_str());
	string outputSpliceAlignmentsFilePath = argv[2];
	ofstream splicedSAM_ofs(outputSpliceAlignmentsFilePath.c_str());
	while(!sam_ifs.eof())
	{
		string tmpSAMstr;
		getline(sam_ifs, tmpSAMstr);
		if(tmpSAMstr.substr(0,1) == "@")
			continue;
		vector<string> tmpSAMfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = tmpSAMstr.find("\t", startLoc);
			string tmpSAMfieldStr = tmpSAMstr.substr(startLoc, tabLoc - startLoc);
			tmpSAMfieldVec.push_back(tmpSAMfieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpCigarString = tmpSAMfieldVec[5];
		int N_loc = tmpCigarString.find("N");
		if(N_loc != string::npos)
			splicedSAM_ofs << tmpSAMstr << endl;
	}
	splicedSAM_ofs.close();
	sam_ifs.close();
	return 0;
}