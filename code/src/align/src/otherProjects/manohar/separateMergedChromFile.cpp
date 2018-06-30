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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputMergedChromFile outputChrNameListFile outputSeparatedChromFolder" << endl;
		exit(1);
	}

	string inputMergedChromFile = argv[1];
	string outputChrNameListFile = argv[2];
	string outputSeparatedChromFolder = argv[3];

	vector<string> chrNameVec;
	vector< vector<string> > chrSeqVecVec;
	cout << "start to read merged chrom file" << endl;
	ifstream mergedChrom_ifs(inputMergedChromFile);
	while(!mergedChrom_ifs.eof())
	{
		string tmpStr;
		getline(mergedChrom_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == ">")
		{	
			chrNameVec.push_back(tmpStr.substr(1));
			vector<string> tmpSeqVec;
			chrSeqVecVec.push_back(tmpSeqVec);
		}
		else
		{
			chrSeqVecVec[chrSeqVecVec.size()-1].push_back(tmpStr);
		}
	}

	cout << "chrFile #: " << chrSeqVecVec.size() << endl;
	cout << "start to output chrFiles " << endl;
	string outputFolderStr = outputSeparatedChromFolder;
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	ofstream chrNameList_ofs(outputChrNameListFile.c_str());	
	for(int tmp = 0; tmp < chrNameVec.size()-1; tmp++)
	{
		string tmpChrName = chrNameVec[tmp];
		string tmpChrFileName = outputFolderStr + tmpChrName + ".fa";
		chrNameList_ofs << tmpChrFileName << " \\" << endl;
		ofstream tmp_ofs(tmpChrFileName.c_str());
		tmp_ofs << ">" << chrNameVec[tmp] << endl;
		for(int tmpSeqLine = 0; tmpSeqLine < chrSeqVecVec[tmp].size(); tmpSeqLine++)
			tmp_ofs << (chrSeqVecVec[tmp])[tmpSeqLine] << endl;
		tmp_ofs.close();
	}
	string lastChrSeqFile = outputFolderStr + chrNameVec[chrNameVec.size()-1] + ".fa";
	chrNameList_ofs << lastChrSeqFile << endl;
	ofstream tmp_ofs(lastChrSeqFile.c_str());
	tmp_ofs << ">" << chrNameVec[chrNameVec.size()-1] << endl;
	for(int tmpSeqLine = 0; tmpSeqLine < chrSeqVecVec[chrNameVec.size()-1].size(); tmpSeqLine++)
		tmp_ofs << (chrSeqVecVec[chrNameVec.size()-1])[tmpSeqLine] << endl;	
	tmp_ofs.close();

	chrNameList_ofs.close();
	mergedChrom_ifs.close();
	return 0;
}