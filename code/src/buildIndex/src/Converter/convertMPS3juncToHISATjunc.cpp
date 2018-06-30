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
	if(argc != 4)
	{
		cout << "Executable inputMPS3juncFile outputHISATjuncFile strandFieldNum" << endl;;
		exit(1);
	}

	string strandFieldNumStr = argv[3];
	int strandFieldNum = atoi(strandFieldNumStr.c_str());
	if((strandFieldNum != 4)&&(strandFieldNum != 5))
	{
		cout << "incorrect strandFieldNum: " << strandFieldNum << endl;
		cout << "exiting ......" << endl;
		exit(1);
	}

	string inputMPS3juncFileStr = argv[1];
	string outputHISATjuncFileStr = argv[2];

	//FILE* fp_MPS3juncFile = fopen(inputMPS3juncFileStr.c_str(), "r");
	ifstream MPS3juncFile_ifs(inputMPS3juncFileStr.c_str());
	ofstream HISATjuncFile_ofs(outputHISATjuncFileStr.c_str());

	//string entryString;
	int tabLocation1;
	int tabLocation2;
	int tabLocation3;
	int tabLocation4;
	int tabLocation5;	
	//char entry[500];

	int chrInt;
	int donerEndPos_1_coordinate;
	int acceptorStartPos_1_coordinate;
	int donerEndPos_0_coordinate;
	int acceptorStartPos_0_coordinate;

	string chrIntString;
	string donerEndPosString_1_coordinate;
	string acceptorStartPosString_1_coordinate;
	string strandString;

	while(!MPS3juncFile_ifs.eof())
	{
		//fgets(entry, sizeof(entry), fp_MPS3juncFile);
		//entryString = entry;
		string entryString;
		getline(MPS3juncFile_ifs, entryString);
		if(entryString == "")
			break;	
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);
		tabLocation4 = entryString.find('\t', tabLocation3+1);

		chrIntString = entryString.substr(0, tabLocation1);
		donerEndPosString_1_coordinate = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		acceptorStartPosString_1_coordinate = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		if(strandFieldNum == 5)
			strandString = entryString.substr(tabLocation4+1, 1);
		else
			strandString = entryString.substr(tabLocation3+1, 1);
		// cout << "entryString: " << entryString << endl;
		// cout << "chrIntString: " << chrIntString << endl;
		// cout << "donerEndPos_1_coordinate: " << donerEndPos_1_coordinate << endl;
		// cout << "acceptorStartPos_1_coordinate: " << acceptorStartPosString_1_coordinate << endl;
		// cout << "tabLocation3: " << tabLocation3 << endl;
		// cout << "tabLocation4: " << tabLocation4 << endl; 
		if(!((strandString == "+")||(strandString == "-")))
		{
			cout << "incorrect strand: " << strandString << endl;
			cout << "exiting ...... " << endl;
			exit(1);
		}
		//chrInt = covertStringToInt(chrIntString);
		donerEndPos_1_coordinate = atoi(donerEndPosString_1_coordinate.c_str());
		acceptorStartPos_1_coordinate = atoi(acceptorStartPosString_1_coordinate.c_str());
		donerEndPos_0_coordinate = donerEndPos_1_coordinate - 1;
		acceptorStartPos_0_coordinate = acceptorStartPos_1_coordinate - 1;

		HISATjuncFile_ofs << chrIntString << "\t" << donerEndPos_0_coordinate
			<< "\t" << acceptorStartPos_0_coordinate << "\t" << strandString << endl;
	}

	HISATjuncFile_ofs.close();
	MPS3juncFile_ifs.close();
	return 0;
}