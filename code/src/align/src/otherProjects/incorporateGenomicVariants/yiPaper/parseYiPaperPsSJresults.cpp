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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputYiPaperPsSJresults outputFormattedPsSJ" << endl;
		exit(1);
	}
	string inputYiPaperPsSJresults = argv[1];
	string outputFormattedPsSJ = argv[2];
	ifstream yiPaperPsSJ_ifs(inputYiPaperPsSJresults.c_str());
	ofstream formattedPsSJ_ofs(outputFormattedPsSJ.c_str());
	while(!yiPaperPsSJ_ifs.eof())
	{
		string tmpStr;
		getline(yiPaperPsSJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpGeneName = tmpStr.substr(0, tabLoc_1);
		string tmpSJinfoStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpOtherStr = tmpStr.substr(tabLoc_2 + 1);
		int tmpCommaLoc = tmpSJinfoStr.find(":");
		int tmpCrossLineLoc = tmpSJinfoStr.find("-");
		string tmpSJ_chrNameStr = tmpSJinfoStr.substr(0, tmpCommaLoc);
		string tmpSJ_intronStartPosStr = tmpSJinfoStr.substr(tmpCommaLoc + 1, 
			tmpCrossLineLoc - tmpCommaLoc - 1);
		int tmpSJ_intronStartPos = atoi(tmpSJ_intronStartPosStr.c_str());
		int tmpSJ_chrStartPos = tmpSJ_intronStartPos - 1;
		string tmpSJ_intronEndPosStr = tmpSJinfoStr.substr(tmpCrossLineLoc + 1);
		int tmpSJ_intronEndPos = atoi(tmpSJ_intronEndPosStr.c_str());
		int tmpSJ_chrEndPos = tmpSJ_intronEndPos + 1;
		formattedPsSJ_ofs << tmpSJ_chrNameStr << "\t"
			<< tmpSJ_chrStartPos << "\t" << tmpSJ_chrEndPos << "\t"
			<< tmpGeneName << "\t" << tmpOtherStr << endl;
	}
	yiPaperPsSJ_ifs.close();
	formattedPsSJ_ofs.close();
	return 0;
}