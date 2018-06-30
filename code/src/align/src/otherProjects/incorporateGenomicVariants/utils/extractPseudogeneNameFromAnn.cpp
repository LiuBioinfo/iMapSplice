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
	if(argc < 3)
	{
		cout << "Executable inputAnn outputPseudogeneNameList" << endl;
		exit(1);
	}
	string inputAnn = argv[1];
	string outputPseudogeneNameList = argv[2];
	set<string> pseudoGeneNameSet;
	ifstream ann_ifs(inputAnn.c_str());
	ofstream pseudoGene_ofs(outputPseudogeneNameList.c_str());
	while(!ann_ifs.eof())
	{
		string tmpAnnStr;
		getline(ann_ifs, tmpAnnStr);
		if(tmpAnnStr == "")
			break;
		int gene_type_loc = tmpAnnStr.find("pseudogene");
		int gene_name_loc = tmpAnnStr.find("gene_name");
		if((gene_type_loc == string::npos)||(gene_name_loc == string::npos))
			continue;
		int nextComma = tmpAnnStr.find(";",gene_name_loc + 1);
		string tmpGeneName = tmpAnnStr.substr(gene_name_loc+11, nextComma - 2 - gene_name_loc - 11 + 1);
		pseudoGeneNameSet.insert(tmpGeneName);
	}
	for(set<string>::iterator tmpIter = pseudoGeneNameSet.begin(); tmpIter != pseudoGeneNameSet.end();
		tmpIter ++)
	{
		pseudoGene_ofs << (*tmpIter) << endl;
	}

	ann_ifs.close();
	pseudoGene_ofs.close();
	return 0;
}