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

using namespace std;

void extractXMfromRegularJuncStr(int& tmpXMmin, int& tmpXMmax, string& tmpFusionJuncStr)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 9; tmp++)
	{
		int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
		string tmpFusionJuncField = tmpFusionJuncStr.substr(startLoc, tabLoc-startLoc);
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncField);
		startLoc = tabLoc + 1;
	}
	tmpXMmin = atoi(tmpFusionJuncFieldVec[7].c_str());
	tmpXMmax = atoi(tmpFusionJuncFieldVec[8].c_str());
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputRegularJunc outputPrefix maxReadLength" << endl;
		exit(1);
	}
	string inputFusionJuncPath = argv[1];
	string outputPrefix = argv[2];
	string maxReadLengthStr = argv[3];
	int maxReadLength = atoi(maxReadLengthStr.c_str());
	vector<int> XM_juncNumVec_min;
	vector<int> XM_juncNumVec_max;
	for(int tmp = 0; tmp < maxReadLength + 1; tmp++)
	{
		XM_juncNumVec_min.push_back(0);
		XM_juncNumVec_max.push_back(0);
	}

	string output_minXM_path = outputPrefix + "_distribution_min.txt";
	string output_maxXM_path = outputPrefix + "_distribution_max.txt";
	ofstream minXM_ofs(output_minXM_path.c_str());
	ofstream maxXM_ofs(output_maxXM_path.c_str());
	ifstream fusionJunc_ifs(inputFusionJuncPath.c_str());
	while(!fusionJunc_ifs.eof())
	{
		string tmpFusionJuncStr;
		getline(fusionJunc_ifs, tmpFusionJuncStr);
		//cout << "tmpFusionStr: " << tmpFusionStr << endl;
		if(tmpFusionJuncStr == "")
			break;		
		int tmpXMmin, tmpXMmax;
		extractXMfromRegularJuncStr(tmpXMmin, tmpXMmax, tmpFusionJuncStr);
		XM_juncNumVec_min[tmpXMmin] ++;
		XM_juncNumVec_max[tmpXMmax] ++;
	}
	for(int tmp = 0; tmp < maxReadLength + 1; tmp++)
	{
		minXM_ofs << tmp << "\t" << XM_juncNumVec_min[tmp] << endl;
		maxXM_ofs << tmp << "\t" << XM_juncNumVec_max[tmp] << endl;		
	}	
	fusionJunc_ifs.close();
	minXM_ofs.close();
	maxXM_ofs.close();
	return 0;
}