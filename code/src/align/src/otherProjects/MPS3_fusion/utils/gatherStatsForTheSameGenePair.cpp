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
		cout << "Executable inputIndexFolderPath inputFusionJuncPath outputGatheredFusionGenePairPath" << endl;
		exit(1);
	}
	//string inputIndexFolderPath = argv[1];
	string inputFusionJuncPath = argv[2];
	string outputGatheredFusionGenePairPath = argv[3];

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "end of initiating indexInfo" << endl;

	
	

	ifstream fj_ifs(inputFusionJuncPath.c_str());
	ofstream gatheredFusionInfo_ofs(outputGatheredFusionGenePairPath.c_str());


	fj_ifs.close();
	gatheredFusionInfo_ofs.close();
	return 0;
}