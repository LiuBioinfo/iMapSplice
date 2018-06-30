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

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputGroundTruthFusionListFile outputFusionGeneNamePair" << endl;
		exit(1);
	}
	string inputFusionList = argv[1];
	string outputFusionList = argv[2];
	ifstream fus_ifs(inputFusionList.c_str());
	ofstream fus_ofs(outputFusionList.c_str());
	while(!fus_ifs.eof())
	{
		string tmpFusStr;
		getline(fus_ifs, tmpFusStr);
		cout << "tmpFusStr" << tmpFusStr << endl;
		if(tmpFusStr == "")
			break;
		cout << "tmpFusStr" << tmpFusStr << endl;
		int tab_1 = tmpFusStr.find("\t");
		int tab_2 = tmpFusStr.find("\t", tab_1 + 1);
		int tab_3 = tmpFusStr.find("\t", tab_2 + 1);
		int tab_4 = tmpFusStr.find("\t", tab_3 + 1);
		string gene1 = tmpFusStr.substr(tab_2 + 1, tab_3 - tab_2 - 1);
		string gene2 = tmpFusStr.substr(tab_3 + 1, tab_4 - tab_3 - 1);
		fus_ofs << gene1 << "," << gene2 << "," << endl;
	}
	fus_ifs.close();
	fus_ofs.close();
	return 0;
}