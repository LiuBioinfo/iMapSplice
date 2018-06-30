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
#include <sstream>
//#include "../../general/read_block_test.h"
//#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputBeersGroundTruthAlignmentPath outputFlagUpdatedSAMpath" << endl;
		exit(1);
	}
	string inputBeersGroundTruthAlignmentPath = argv[1];
	string outputFlagUpdatedSAMpath = argv[2];

	ifstream oriSAM_ifs(inputBeersGroundTruthAlignmentPath.c_str());
	ofstream flagUpdatedSAM_ofs(outputFlagUpdatedSAMpath.c_str());

	while(!oriSAM_ifs.eof())
	{
		string tmpSamStr_1, tmpSamStr_2;
		getline(oriSAM_ifs, tmpSamStr_1);
		if(tmpSamStr_1 == "")
			break;
		getline(oriSAM_ifs, tmpSamStr_2);
		int firstTabLoc = tmpSamStr_1.find("\t");
		int secondTabLoc = tmpSamStr_1.find("\t", firstTabLoc + 1);
		string readName_1 = tmpSamStr_1.substr(0, firstTabLoc);
		string flagStr_1 = tmpSamStr_1.substr(firstTabLoc + 1, secondTabLoc - firstTabLoc - 1);
		string otherFieldStr_1 = tmpSamStr_1.substr(secondTabLoc + 1);
		firstTabLoc = tmpSamStr_2.find("\t");
		secondTabLoc = tmpSamStr_2.find("\t", firstTabLoc + 1);
		string readName_2 = tmpSamStr_2.substr(0, firstTabLoc);
		string flagStr_2 = tmpSamStr_2.substr(firstTabLoc + 1, secondTabLoc - firstTabLoc - 1);
		string otherFieldStr_2 = tmpSamStr_2.substr(secondTabLoc + 1);

		int readName_1_length = readName_1.length();
		int readName_2_length = readName_2.length();
		if((flagStr_1 == "0")&&(flagStr_2 == "16"))
		{
			flagUpdatedSAM_ofs << readName_1.substr(0, readName_1_length-2) << "\t99\t" << otherFieldStr_1 << endl;
			flagUpdatedSAM_ofs << readName_2.substr(0, readName_2_length-2) << "\t147\t" << otherFieldStr_2 << endl;
		}
		else if((flagStr_1 == "16")&&(flagStr_2 == "0"))
		{
			flagUpdatedSAM_ofs << readName_1.substr(0, readName_1_length-2) << "\t163\t" << otherFieldStr_1 << endl;
			flagUpdatedSAM_ofs << readName_2.substr(0, readName_2_length-2) << "\t83\t" << otherFieldStr_2 << endl;			
		}
		else
		{
			cout << "error !" << endl << "tmpSamStr_1: " << tmpSamStr_1 << endl
				<< "tmpSamStr_2: " << tmpSamStr_2 <<endl;
		}
	}
	oriSAM_ifs.close();
	flagUpdatedSAM_ofs.close();
	return 0;
}