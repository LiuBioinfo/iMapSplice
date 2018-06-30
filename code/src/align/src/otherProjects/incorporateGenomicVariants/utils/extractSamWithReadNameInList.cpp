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

//#include "../../../general/read_block_test.h"
//#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputReadNameListFile inputSamFile outputSamWithReadNameInList " << endl;
		exit(1);
	}
	string inputReadNameListFile = argv[1];
	string inputSamFile = argv[2];
	string outputSamWithReadNameInList = argv[3];
	cout << "start to build readNameSet ......" << endl;
	set<string> readNameSet;
	ifstream readNameList_ifs(inputReadNameListFile.c_str());
	while(!readNameList_ifs.eof())
	{
		string tmpReadNameStr;
		getline(readNameList_ifs, tmpReadNameStr);
		if(tmpReadNameStr == "")
			break;
		readNameSet.insert(tmpReadNameStr);
	}
	readNameList_ifs.close();
	cout << "start to scan sam file and extract alignments with read name in list" << endl;
	ifstream sam_ifs(inputSamFile.c_str());
	ofstream samWithReadNameInList_ofs(outputSamWithReadNameInList.c_str());
	while(!sam_ifs.eof())
	{
		string tmpSamStr;
		getline(sam_ifs, tmpSamStr);
		if(tmpSamStr == "")
			break;
		int tabLoc_1st = tmpSamStr.find("\t");
		string tmpReadName = tmpSamStr.substr(0, tabLoc_1st);
		if(readNameSet.find(tmpReadName) != readNameSet.end())
			samWithReadNameInList_ofs << tmpSamStr << endl;
	}
	samWithReadNameInList_ofs.close();
	sam_ifs.close();
	return 0;
}