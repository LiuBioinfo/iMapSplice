// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to check mappedLength's distribution of primary alignments
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
/*
int getFlagFromSamStr(const string& tmpSamStr)
{
	int tmpFlag;

	return tmpFlag;
}*/

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool primaryOrNot(int tmpFlag)
{
	if(tmpFlag & 0x100)
		return false;
	else
		return true;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Exetuable inputSamFile outputFile!" << endl;
		exit(1);
	}

	string inputFileStr = argv[1];
	string outputFileStr = argv[2];
	string outputFileStr_primary = outputFileStr + ".primary";
	string outputFileStr_nonPrimary = outputFileStr + ".secondary";
	string outputFileStr_unmapped = outputFileStr + ".unmapped";
	string outputFileStr_log = outputFileStr + ".log";

	FILE *fp;
	fp = fopen(inputFileStr.c_str(), "r");
	ofstream output_primary_ofs(outputFileStr_primary.c_str());
	ofstream output_nonPrimary_ofs(outputFileStr_nonPrimary.c_str());
	ofstream output_unmmapped_ofs(outputFileStr_unmapped.c_str());
	ofstream output_log_ofs(outputFileStr_log.c_str());

	int oriAlignmentNum = 0;
	int primaryAlignmentNum = 0;
	int nonPrimaryAlignmentNum = 0;
	int unmappedAlignmentNum = 0;

	while(!feof(fp))
	{
		char tmpSamStrChar[2000], tmpReadNameChar[100], tmpFlagChar[3],
			tmpOthersChar[1500];
		fgets(tmpSamStrChar, sizeof(tmpSamStrChar),fp);
		if(feof(fp))
			break;	
		if(tmpSamStrChar[0] == '@')
			continue;
		oriAlignmentNum ++;

		sscanf(tmpSamStrChar,"%s\t%s\t%s", 
			tmpReadNameChar, tmpFlagChar, tmpOthersChar);  

		int tmpFlag = atoi(tmpFlagChar);
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		bool primaryOrNot_bool = primaryOrNot(tmpFlag);
		if(!mappedOrNot_bool)
		{
			unmappedAlignmentNum ++;
			output_unmmapped_ofs << tmpSamStrChar;
		}
		else
		{
			if(primaryOrNot_bool)
			{
				primaryAlignmentNum ++;
				output_primary_ofs << tmpSamStrChar;
			}
			else
			{
				nonPrimaryAlignmentNum ++;
				output_nonPrimary_ofs << tmpSamStrChar;
			}
		}
	}
	output_nonPrimary_ofs << endl;	
	output_unmmapped_ofs << endl;	
	output_log_ofs << endl;	

	output_log_ofs << "oriAlignmentNum: " << oriAlignmentNum << endl;
	output_log_ofs << "mappedAlignmentNum: " << primaryAlignmentNum + nonPrimaryAlignmentNum 
		<< "\tperc: " << (((double)(primaryAlignmentNum + nonPrimaryAlignmentNum))/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "\t****** primaryAlignmentNum: " << primaryAlignmentNum
		<< "\tperc: " << ((double)primaryAlignmentNum/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "\t****** nonPrimaryAlignmentNum: " << nonPrimaryAlignmentNum
		<< "\tperc: " << ((double)nonPrimaryAlignmentNum/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "unmappedAlignmentNum: " << unmappedAlignmentNum 
		<< "\tperc: " << ((double)unmappedAlignmentNum/(double)oriAlignmentNum)*100 << endl;

	fclose(fp);
	output_primary_ofs.close();
	output_nonPrimary_ofs.close();
	output_unmmapped_ofs.close();
	output_log_ofs.close();
	return 0;
}