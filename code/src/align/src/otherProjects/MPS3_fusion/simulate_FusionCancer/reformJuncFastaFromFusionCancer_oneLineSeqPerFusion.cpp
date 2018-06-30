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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputOriFusionJuncFastaFileFromFusionCancer";
		cout << " outputReformedFusionJuncFasta" << endl;
		exit(1);
	}
	string inputFastaPath = argv[1];
	string outputFastaPath = argv[2];
	ifstream fa_ifs(inputFastaPath.c_str());
	ofstream fusionSeq_ofs(outputFastaPath.c_str());

	//vector<string> tmpJuncSeqVec;
	string fusionName;
	string fusionSeq;
	getline(fa_ifs, fusionName);
	//cout << "fusionName: " << endl << fusionName <<endl;
	while(!fa_ifs.eof())
	{
		string tmpFaStr;
		getline(fa_ifs, tmpFaStr);
		//cout << "tmpFaStr: " << endl << tmpFaStr << endl;
		if((fa_ifs.eof()))
			break;
		if(tmpFaStr == "")
			continue;
		if(tmpFaStr.at(0) == '>')
		{
			fusionSeq_ofs << fusionName << endl
				<< fusionSeq << endl;
			fusionName = tmpFaStr;
			fusionSeq = "";
		}
		else
		{
			fusionSeq += tmpFaStr;
		}	
	}

	fa_ifs.close();
	fusionSeq_ofs.close();
	return 0;
}