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

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		cout << "Executable inputPseudogeneNameList inputFeatureCounts outputPseudogeneFeatureCounts" << endl;
		exit(1);
	}
	string inputPseudogeneNameList = argv[1];
	set<string> pseudoGeneNameSet;
	ifstream pseudoGene_ifs(inputPseudogeneNameList.c_str());
	while(!pseudoGene_ifs.eof())
	{
		string tmpGeneNameStr;
		getline(pseudoGene_ifs, tmpGeneNameStr);
		if(tmpGeneNameStr == "")
			break;
		pseudoGeneNameSet.insert(tmpGeneNameStr);
	}
	pseudoGene_ifs.close();
	
	string inputFeatureCounts = argv[2];
	string outputPseudogeneFeatureCounts = argv[3];
	ifstream fc_ifs(inputFeatureCounts.c_str());
	ofstream pfc_ofs(outputPseudogeneFeatureCounts.c_str());
	while(!fc_ifs.eof())
	{
		string tmpStr;
		getline(fc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpGene = tmpStr.substr(0, tabLoc);
		if(pseudoGeneNameSet.find(tmpGene) != pseudoGeneNameSet.end())
			pfc_ofs << tmpStr << endl;
	}
	fc_ifs.close();
	pfc_ofs.close();
	return 0;
}