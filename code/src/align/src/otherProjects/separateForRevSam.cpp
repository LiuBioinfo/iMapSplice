// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// used to separate forward alignments and rev alignments
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

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool alignmentForwardOrRev(int tmpFlag)
{
	if(tmpFlag & 0x10)
		return false;
	else
		return true;
}	

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputPEsamFile outputFolder" << endl;
		exit(1);
	}
	string inputSAMfileStr = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_unmap_file_str = outputFolderStr + "unmap.sam";
	string output_for_file_str = outputFolderStr + "for.sam";
	string output_rev_file_str = outputFolderStr + "rev.sam";
	ofstream unmap_ofs(output_unmap_file_str.c_str());
	ofstream for_ofs(output_for_file_str.c_str());
	ofstream rev_ofs(output_rev_file_str.c_str());

	ifstream tmpInputSAM_ifs(inputSAMfileStr.c_str());
	while(!(tmpInputSAM_ifs.eof()))
	{
		string samStr;
		getline(tmpInputSAM_ifs, samStr);
		 if(tmpInputSAM_ifs.eof()||(samStr == ""))
		 	break;
		if(samStr.at(0) == '@')
			continue;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec.push_back(samStr.substr(startLoc));
		string flagStr = samFieldVec[1];
		int tmpFlag = atoi(flagStr.c_str());
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		bool forOrRev_bool = alignmentForwardOrRev(tmpFlag);
		if(mappedOrNot_bool)
		{
			if(forOrRev_bool)
				for_ofs << samStr << endl;
			else
				rev_ofs << samStr << endl;
		}
		else
		{
			unmap_ofs << samStr << endl;
		}
	}
	for_ofs.close();
	rev_ofs.close();
	tmpInputSAM_ifs.close();
	return 0;
}