// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
//#include <omp.h>

int read_num = 0; // to calculate the read num;

//#include "seg_info.h"

#define PreIndexSize 268435456

using namespace std;  

int baseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 100,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};

unsigned int getPreIndexNO(const string& readPreStr)
{
	int preIndexStrSize = readPreStr.length();

	unsigned int preIndexNO = 0;

	int baseForCount = 1;

	for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
	{
		//char tmpChar = readPreStr.at(tmp);
		//unsigned int tmpArrayIndex = tmpChar - 'A';
		//baseChar2int(tmpChar);
		//preIndexNO = preIndexNO + baseChar2int(tmpChar) * baseForCount;
		preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
		baseForCount = baseForCount * 4;
	}		
	return preIndexNO;
}

int main(int argc, char**argv)
{
	cout << "argc = " << argc << endl;
    if(argc < 1)
	{
		cout << "Executable <preIndexArrayPre>" << endl;
		exit(0);
	}

	string preIndexFileStr = argv[1];//"/data/homes/lxauky/adSA_result/chrAll/result_0222/preIndexRecord";
	preIndexFileStr.append("_preIndexRecord");
	FILE* fp_in_preIndex = fopen(preIndexFileStr.c_str(), "r"); 

	string preIndexArrayPreStr = argv[1];

	string preIndexMapLengthArrayStr = preIndexArrayPreStr;
	preIndexMapLengthArrayStr.append("_MapLength"); 
	ofstream preIndexMapLengthArray_ofs(preIndexMapLengthArrayStr.c_str(), ios::binary);

	string preIndexIntervalStartArrayStr = preIndexArrayPreStr;
	preIndexIntervalStartArrayStr.append("_IntervalStart");
	ofstream preIndexIntervalStartArray_ofs(preIndexIntervalStartArrayStr.c_str(), ios::binary);

	string preIndexIntervalEndArrayStr = preIndexArrayPreStr;
	preIndexIntervalEndArrayStr.append("_IntervalEnd");
	ofstream preIndexIntervalEndArray_ofs(preIndexIntervalEndArrayStr.c_str(), ios::binary);

	int* preIndexMapLengthArray;
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int));

	unsigned int *preIndexIntervalStartArray;
    preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));

	unsigned int *preIndexIntervalEndArray;
    preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));
	
	char preIndexRecordChar[100];
	char preIndexStrChar[20];
	char preIndexMapLengthChar[10];
	char preIndexIntervalStartChar[20];
	char preIndexIntervalEndChar[20];

	int preIndexMapLengthInt;
	unsigned int preIndexIntervalStartInt;
	unsigned int preIndexIntervalEndInt;

	unsigned int preIndexNO = 0;
	//cout << "start to read preIndex file " << endl;
	//if(!feof(fp_in_preIndex))
	//Seg_Info* preIndexSegInfo = new Seg_Info();

	while(!feof(fp_in_preIndex))
	{
		//cout << "preIndexNO: " << preIndexNO;// << endl;
		fgets(preIndexRecordChar, sizeof(preIndexRecordChar), fp_in_preIndex);
		//cout << "record: " << preIndexRecordChar << endl;

		sscanf(preIndexRecordChar, "%s\t%s\t%s\t%s", preIndexStrChar, preIndexMapLengthChar,
			preIndexIntervalStartChar, preIndexIntervalEndChar);
		preIndexMapLengthInt = atoi(preIndexMapLengthChar);
		preIndexIntervalStartInt = strtoul(preIndexIntervalStartChar, NULL, 10);
		preIndexIntervalEndInt = strtoul(preIndexIntervalEndChar, NULL, 10);
		
		string tmpPreIndexStr = preIndexStrChar;
		preIndexNO = getPreIndexNO(tmpPreIndexStr);

		preIndexMapLengthArray[preIndexNO] = preIndexMapLengthInt;
		preIndexIntervalStartArray[preIndexNO] = preIndexIntervalStartInt;
		preIndexIntervalEndArray[preIndexNO] = preIndexIntervalEndInt;

 	}
 	//delete(preIndexSegInfo);
 	fclose(fp_in_preIndex);

 	preIndexMapLengthArray_ofs.write((const char*) preIndexMapLengthArray, PreIndexSize * sizeof(int));
 	preIndexIntervalStartArray_ofs.write((const char*) preIndexIntervalStartArray, PreIndexSize * sizeof(unsigned int));
 	preIndexIntervalEndArray_ofs.write((const char*) preIndexIntervalEndArray, PreIndexSize * sizeof(unsigned int));

 	//cout << "mapLengthArray[0] = " << preIndexMapLengthArray[0] << endl;
 	//cout << "mapLengthArray[PreIndexSize-2] = " << preIndexMapLengthArray[PreIndexSize-2] << endl;
 	//cout << "mapLengthArray[PreIndexSize-1] = " << preIndexMapLengthArray[PreIndexSize-1] << endl;

 	//cout << "preIndexIntervalStartArray[0] = " << preIndexIntervalStartArray[0] << endl;
 	//cout << "preIndexIntervalStartArray[PreIndexSize-2] = " << preIndexIntervalStartArray[PreIndexSize-2] << endl;
 	//cout << "preIndexIntervalStartArray[PreIndexSize-1] = " << preIndexIntervalStartArray[PreIndexSize-1] << endl;

 	//cout << "preIndexIntervalEndArray[0] = " << preIndexIntervalEndArray[0] << endl;
 	//cout << "preIndexIntervalEndArray[PreIndexSize-2] = " << preIndexIntervalEndArray[PreIndexSize-2] << endl;
 	//cout << "preIndexIntervalEndArray[PreIndexSize-1] = " << preIndexIntervalEndArray[PreIndexSize-1] << endl;

 	free(preIndexMapLengthArray);
 	free(preIndexIntervalStartArray);
 	free(preIndexIntervalEndArray);

	return 0;
}