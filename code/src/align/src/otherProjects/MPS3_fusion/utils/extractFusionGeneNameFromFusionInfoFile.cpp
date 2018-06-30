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
	if(argc != 3)
	{
		cout << "Executable inputFusionInfoFile outputFusionGeneNameFile" << endl;
		exit(1);
	}
	string inputFusionInfoFile = argv[1];
	string outputFusionGeneNameFile = argv[2];
	ifstream fusion_ifs(inputFusionInfoFile.c_str());
	ofstream geneName_ofs(outputFusionGeneNameFile.c_str());
	while(!fusion_ifs.eof())
	{
		string tmpFusStr;
		getline(fusion_ifs, tmpFusStr);
		if(tmpFusStr == "")
			break;
		vector<string> tmpFusFieldVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpFusStr.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpFusField = tmpFusStr.substr(startLoc, tabLoc-startLoc);
			tmpFusFieldVec.push_back(tmpFusField);
			startLoc = tabLoc + 1;
		}
		tmpFusFieldVec.push_back(tmpFusStr.substr(startLoc));
		int tmpFusFieldVecSize = tmpFusFieldVec.size();
		string gene_1 = tmpFusFieldVec[tmpFusFieldVecSize - 2];
		string gene_2 = tmpFusFieldVec[tmpFusFieldVecSize - 1];
		string genePair = gene_1 + gene_2;
		geneName_ofs << genePair << endl;
	}
	fusion_ifs.close();
	geneName_ofs.close();
	return 0;
}