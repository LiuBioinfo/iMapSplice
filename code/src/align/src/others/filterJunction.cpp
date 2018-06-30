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

#define SPLICE_MIN_LENGTH 50
using namespace std;

int main(int argc, char**argv)
{
    if(argc < 3)
	{
		cout << "Executable <InputSpliceJunction> <OutputSpliceJunction>"<< endl;
		exit(0);
	}

	char* Input1 = argv[1];
	char* Output = argv[2];
	string OutputStr = Output;
	string filterOutJunctionFile = OutputStr + ".filterOut";
	string filteredJunctionFile = OutputStr + ".filtered";
	ofstream filterOut_ofs(filterOutJunctionFile.c_str());
	ofstream filtered_ofs(filteredJunctionFile.c_str());
	/*string outputFile = argv[3];

	string outputFoundSpliceJunction = outputFile + ".found";
	ofstream outputFoundSpliceJunction_ofs(outputFoundSpliceJunction.c_str());
	string outputNotFoundSpliceJunction = outputFile + ".notFound";
	ofstream outputNotFoundSpliceJunction_ofs(outputNotFoundSpliceJunction.c_str());*/

	FILE *fp_in1 = fopen(Input1, "r");
	string entryString;
	//int tabLocation1;
	//int tabLocation2;
	//int tabLocation3;
	char entry[500], chromChar[10], //spliceStartPosChar[20], spliceEndPosChar[20], 
	junctionNameChar[20],
	//coverageChar[20], 
	strandChar, 
	//blockEndChar[20], blockStartChar[20], 
	itemRgbChar[20], 
	//blockCountChar[10], 
	blockSizesChar[20], 
	blockStartsChar[20], entropyChar[20], 
	//flankStringCaseChar, 
	flankStringChar[4], othersChar[200];
	int spliceStartPosInt, spliceEndPosInt, coverageInt, flankStringCaseInt, 
	blockEndInt, blockStartInt, blockCountInt;
	int chrInt;
	//int spliceStartPos;
	//int spliceEndPos;
	//string chrIntString;
	//string spliceStartPosString;
	//string spliceEndPosString;
	int totalNum = 0;
	int deletionNum = 0;
	int spliceNum = 0;
	int leftAfterFilter = 0, throwAfterFilter = 0;
	//int groundTruthSpliceNum = 0;
	//int toFindSpliceNum = 0;

	fgets(entry, sizeof(entry), fp_in1);
	while(!feof(fp_in1))
	{
		totalNum++;
		fgets(entry, sizeof(entry), fp_in1);
		//cout << entry << endl;
		sscanf(entry, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s",
			chromChar, &spliceStartPosInt, &spliceEndPosInt, junctionNameChar, &coverageInt,
			&strandChar, &blockEndInt, &blockStartInt, itemRgbChar, &blockCountInt, blockSizesChar,
			blockStartsChar, entropyChar, &flankStringCaseInt, flankStringChar, othersChar);
		//cout << "spliceStartPosInt = " << spliceStartPosInt << endl;
		if((spliceEndPosInt - spliceStartPosInt) < SPLICE_MIN_LENGTH)
		{
			deletionNum++;
			continue;
		}
		spliceNum++;
		if((coverageInt <= 1)&&
			(flankStringCaseInt != 5) && (flankStringCaseInt != 6)
			) 
		{
			throwAfterFilter ++;
			filterOut_ofs << entry ;
		}
		else
		{
			leftAfterFilter++;
			filtered_ofs << entry ; 
		}
	}
	filterOut_ofs << endl;
	filtered_ofs << endl;
	cout << "totalNum: " << totalNum << endl;
	cout << "spliceNum: " << spliceNum << endl << "deletionNum: " << deletionNum << endl;
	cout << "filterOutNum: " << throwAfterFilter << endl << "LeftNum: " << leftAfterFilter << endl; 
	return 0;
}