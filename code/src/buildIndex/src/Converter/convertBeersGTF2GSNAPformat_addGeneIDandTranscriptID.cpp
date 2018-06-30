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
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputBeersGTFfile outputGSNAPgtfFile" << endl;;
		exit(1);
	}

	string inputBeersGTFfileStr = argv[1];
	string outputGSNAPgtfFileStr = argv[2];

	FILE* fp_BeersGTFfile = fopen(inputBeersGTFfileStr.c_str(), "r");
	ofstream outputGSNAPgtf_ofs(outputGSNAPgtfFileStr.c_str());

	string entryString;
	char entry[500];

	while(!feof(fp_BeersGTFfile))
	{
		fgets(entry, sizeof(entry), fp_BeersGTFfile);
		if(feof(fp_BeersGTFfile))
			break;
		entryString = entry;
		vector<string> gtfFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 8; tmp++)
		{
			int tabLoc = entryString.find("\t", startLoc);
			string tmpGtfFieldStr = entryString.substr(startLoc, tabLoc-startLoc);
			gtfFieldStrVec.push_back(tmpGtfFieldStr);
			startLoc = tabLoc + 1;
		}
		gtfFieldStrVec.push_back(entryString.substr(startLoc));
		string tmpGeneNameStr = gtfFieldStrVec[8];
		string tmpGeneNameNum = tmpGeneNameStr.substr(5);
		tmpGeneNameNum = tmpGeneNameNum.substr(0, tmpGeneNameNum.length()-1);
		string GSNAP_gene_transcript_ID = "gene_id \"" 
			+ tmpGeneNameNum + "\"; transcript_id \""
			+ tmpGeneNameNum + "\";";
		outputGSNAPgtf_ofs << gtfFieldStrVec[0] << "\t" << gtfFieldStrVec[1] << "\t"
			<< gtfFieldStrVec[2] << "\t" << gtfFieldStrVec[3] << "\t"
			<< gtfFieldStrVec[4] << "\t" << gtfFieldStrVec[5] << "\t"
			<< gtfFieldStrVec[6] << "\t" << gtfFieldStrVec[7] << "\t"
			<< GSNAP_gene_transcript_ID << endl; 	
	}
	outputGSNAPgtf_ofs.close();
	fclose(fp_BeersGTFfile);
	return 0;
}