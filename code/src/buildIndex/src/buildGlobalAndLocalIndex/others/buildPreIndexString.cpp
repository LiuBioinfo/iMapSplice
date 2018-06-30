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
#include <stdio.h>

using namespace std;

//typedef unsigned char BYTE;

int main(int argc, char**argv)
{
	if(argc < 2)
	{
		cout << "Executable <preIndexStringPrefix> <preIndexStringLength>" << endl;
		exit(0);
	}

	string preIndexStringPrefixStr = argv[1];
	string preIndexStringLengthStr = argv[2];
	int preIndexStringLength = atoi(preIndexStringLengthStr.c_str());
	
	string generateStringFilePrefix = preIndexStringPrefixStr + "_preIndexString";
	int stringLengthInHash = preIndexStringLength;

	//write stringfiles////////////////////////////////////////////////
	cout << "start to write string files" << endl;

	string baseStr[4] = {"A", "C", "G", "T"};
	string generateStringFile[stringLengthInHash];// = argv[1];
	//string generateStringFilePrefix = argv[1];
	generateStringFile[0] = generateStringFilePrefix + ".1";

	ofstream GenerateStringFile_ofs(generateStringFile[0].c_str());
	GenerateStringFile_ofs << "A" << endl << "C" << endl
		<< "G" << endl << "T";
	GenerateStringFile_ofs.close();
	string readString;
	char readChar[20];
	for(int tmp = 1; tmp < stringLengthInHash; tmp++)
	{
		char tmpChar[2];
		sprintf(tmpChar, "%d", tmp+1);
		string tmpString = tmpChar;
		generateStringFile[tmp] = generateStringFilePrefix + "." + tmpString;

		FILE *fp_in = fopen(generateStringFile[tmp-1].c_str(), "r");
		ofstream GenerateStringFile_ofs(generateStringFile[tmp].c_str());

    	while(!feof(fp_in))
    	{
    		fgets(readChar, sizeof(readChar), fp_in);
    		readString = readChar;
    		readString = readString.substr(0, tmp);
    		for(int tmpBase = 0; tmpBase < 4; tmpBase ++)
    		{
    			string resultString = readString + baseStr[tmpBase];
    			GenerateStringFile_ofs << resultString << endl;
    		}
    	}  fclose(fp_in);
    	GenerateStringFile_ofs.close();
    	remove(generateStringFile[tmp-1].c_str());
	}

	return 0;
}