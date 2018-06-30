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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputReadNameListFile inputFaFile outputFaWithReadNameInList " << endl;
		exit(1);
	}
	string inputReadNameListFile = argv[1];
	string inputFaFile = argv[2];
	string outputFaWithReadNameInList = argv[3];
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
	ifstream fa_ifs(inputFaFile.c_str());
	ofstream faWithReadNameInList_ofs(outputFaWithReadNameInList.c_str());
	while(!fa_ifs.eof())
	{
		string tmpFaStr_1;
		getline(fa_ifs, tmpFaStr_1);
		if(tmpFaStr_1 == "")
			break;
		string tmpFaStr_2;
		getline(fa_ifs, tmpFaStr_2);
		int oriReadNameLength = tmpFaStr_1.length();
		string tmpReadName = tmpFaStr_1.substr(1, oriReadNameLength - 3);
		if(readNameSet.find(tmpReadName) != readNameSet.end())
			faWithReadNameInList_ofs << tmpFaStr_1 << endl << tmpFaStr_2 << endl;
	}
	faWithReadNameInList_ofs.close();
	fa_ifs.close();
	return 0;
}