// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// used to get size for each exon
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
	if(argc != 3)
	{
		cout << "Executable inputGTF outputExonSizeFile" << endl;
		exit(1);
	}

	string gtf_input = argv[1];
	string exonSize_output = argv[2];

	ifstream gtf_ifs(gtf_input.c_str());
	ofstream exonSize_ofs(exonSize_output.c_str());
	while(!(gtf_ifs.eof()))
	{
		string lineStr;
		getline(gtf_ifs, lineStr);
		 if(gtf_ifs.eof()||(lineStr == ""))
		 	break;
		vector<string> fieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = lineStr.find("\t", startLoc);
			string tmpField = lineStr.substr(startLoc, tabLoc-startLoc);
			fieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		if(fieldVec[2] == "exon")
		{
			int tmpStartPos = atoi(fieldVec[3].c_str());
			int tmpEndPos = atoi(fieldVec[4].c_str());
			int exonSize = tmpEndPos - tmpStartPos + 1;
			exonSize_ofs << fieldVec[0] << "\t" << fieldVec[1] << "\t"
				<< fieldVec[2] << "\t" << fieldVec[3] << "\t" 
				<< fieldVec[4] << "\t" << exonSize << endl;

		}
	} 
	return 0;
}