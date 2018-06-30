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
#include <map>
#include <hash_map>
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputHISATjuncFile outputSTARjuncFile" << endl;
		exit(1);
	}

	string inputHISATjuncFileStr = argv[1];
	string outputSTARjuncFileStr = argv[2];

	FILE* fp_HISATjuncFile = fopen(inputHISATjuncFileStr.c_str(), "r");
	ofstream STARjuncFile_ofs(outputSTARjuncFileStr.c_str());

	string entryString;
	int tabLocation1;
	int tabLocation2;
	int tabLocation3;
	int tabLocation4;
	int tabLocation5;	
	char entry[500];

	int chrInt;
	int spliceStartPos_0_coordinate;
	int spliceEndPos_0_coordinate;
	int intronStartPos_1_coordinate;
	int intronEndPos_1_coordinate;

	string chrIntString;
	string spliceStartPosString;
	string spliceEndPosString;
	string strandString;

	while(!feof(fp_HISATjuncFile))
	{
		fgets(entry, sizeof(entry), fp_HISATjuncFile);
		entryString = entry;
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);

		chrIntString = entryString.substr(0, tabLocation1);
		spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		strandString = entryString.substr(tabLocation3+1, 1);

		//chrInt = covertStringToInt(chrIntString);
		spliceStartPos_0_coordinate = atoi(spliceStartPosString.c_str());
		spliceEndPos_0_coordinate = atoi(spliceEndPosString.c_str());
		intronStartPos_1_coordinate = spliceStartPos_0_coordinate + 1 + 1;
		intronEndPos_1_coordinate = spliceEndPos_0_coordinate;

		STARjuncFile_ofs << chrIntString << "\t" << intronStartPos_1_coordinate
			<< "\t" << intronEndPos_1_coordinate << "\t" << strandString << endl;
	}

	STARjuncFile_ofs.close();
	fclose(fp_HISATjuncFile);
	return 0;
}