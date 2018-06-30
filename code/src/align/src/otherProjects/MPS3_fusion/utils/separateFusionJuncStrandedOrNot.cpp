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
		cout << "Executable inputFusionJuncPath outputPrefix" << endl;
		exit(1);
	}
	string inputFusionJuncPath = argv[1];
	string outputPrefix = argv[2];
	string output_stranded = outputPrefix + ".stranded";
	string output_nonStranded = outputPrefix + ".nonStranded";

	ifstream fj_ifs(inputFusionJuncPath.c_str());
	ofstream stranded_ofs(output_stranded.c_str());
	ofstream nonStranded_ofs(output_nonStranded.c_str());
	while(!fj_ifs.eof())
	{
		string fjStr;
		getline(fj_ifs, fjStr);
		if(fjStr == "")
			break;
		vector<string> fjFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = fjStr.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpFjField = fjStr.substr(startLoc, tabLoc-startLoc);
			fjFieldVec.push_back(tmpFjField);
			startLoc = tabLoc + 1;
		}
		string tmpStrand_1 = fjFieldVec[4];
		string tmpStrand_2 = fjFieldVec[5];
		string tmpStrand = tmpStrand_1 + tmpStrand_2;
		if((tmpStrand == "++")||(tmpStrand == "--")||(tmpStrand == "+-")||(tmpStrand == "-+"))
			stranded_ofs << fjStr << endl;
		else
			nonStranded_ofs << fjStr << endl;
	}
	fj_ifs.close();
	stranded_ofs.close();
	nonStranded_ofs.close();
	return 0;
}