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
#include <hash_map>
#include <map>
#include <set>
//#include <omp.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable readNum_raw readNum_trimmed inputUnmapPerc outputUnmapPerc" << endl;
		exit(1);
	}
	string readNum_raw_str = argv[1];
	int readNum_raw = atoi(readNum_raw_str.c_str());
	string readNum_trimmed_str = argv[2];
	int readNum_trimmed = atoi(readNum_trimmed_str.c_str());
	string inputUnmapPerc_file = argv[3];
	string outputUnmapPerc_file = argv[4];
	ifstream unmapPerc_ifs(inputUnmapPerc_file.c_str());
	ofstream unmapPerc_ofs(outputUnmapPerc_file.c_str());
	while(!unmapPerc_ifs.eof())
	{
		string tmpUnmapPercStr;
		getline(unmapPerc_ifs, tmpUnmapPercStr);
		if(tmpUnmapPercStr == "")
			break;
		int tabLoc_1 = tmpUnmapPercStr.find("\t");
		int tabLoc_2 = tmpUnmapPercStr.find("\t", tabLoc_1 + 1);
		string basePosStr = tmpUnmapPercStr.substr(0, tabLoc_1);
		string toolNameStr = tmpUnmapPercStr.substr(tabLoc_2 + 1);
		string unmapPercStr_ori = tmpUnmapPercStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		double unmapPerc_ori = atof(unmapPercStr_ori.c_str());
		cout << "unmapPerc_ori: " << unmapPerc_ori << endl;
		double unmapPerc_updated = 100 - (100.00 - unmapPerc_ori) * ((double)readNum_trimmed/(double)readNum_raw); 
		cout << "unmapPerc_updated: " << unmapPerc_updated << endl; 
		unmapPerc_ofs << basePosStr << "\t" << unmapPerc_updated << "\t" << toolNameStr << endl;
	}
	unmapPerc_ifs.close();
	unmapPerc_ofs.close();
	return 0;
}