// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

void extractFusionGeneNameBreakPointFrom_fusionFa1stLine(
	string& tmpFusionGeneName_1, string& tmpFusionGeneName_2,
	string& tmpFusionChrName_1, string& tmpFusionChrName_2,
	string& tmpFusionBreakPoint_1, string& tmpFusionBreakPoint_2,
	const string& tmpFusionFaStr_name)
{
	int loc_shortCrossLine = tmpFusionFaStr_name.find(")", 0) + 1;
	int loc_comma = tmpFusionFaStr_name.find(",", 0);
	string tmpFusionGeneInfoStr_1
		= tmpFusionFaStr_name.substr(1, loc_shortCrossLine-1);
	string tmpFusionGeneInfoStr_2
		= tmpFusionFaStr_name.substr(loc_shortCrossLine+1,
			loc_comma - loc_shortCrossLine - 1);
	int loc_leftPar_gene1 = tmpFusionGeneInfoStr_1.find("(", 0);
	int loc_leftPar_gene2 = tmpFusionGeneInfoStr_2.find("(", 0);
	tmpFusionGeneName_1 = tmpFusionGeneInfoStr_1.substr(0, loc_leftPar_gene1);
	tmpFusionGeneName_2 = tmpFusionGeneInfoStr_2.substr(0, loc_leftPar_gene2);
	int loc_straightLine_gene1 = tmpFusionGeneInfoStr_1.find("|", 0);
	int loc_straightLine_gene2 = tmpFusionGeneInfoStr_2.find("|", 0);
	int loc_colon_gene1 = tmpFusionGeneInfoStr_1.find(":", 0);
	int loc_colon_gene2 = tmpFusionGeneInfoStr_2.find(":", 0);
	int loc_rightPar_gene1 = tmpFusionGeneInfoStr_1.find(")", 0);
	int loc_rightPar_gene2 = tmpFusionGeneInfoStr_2.find(")", 0);
	tmpFusionChrName_1 = tmpFusionGeneInfoStr_1.substr(
		loc_straightLine_gene1 + 1, loc_colon_gene1 - loc_straightLine_gene1 - 1);
	tmpFusionChrName_2 = tmpFusionGeneInfoStr_2.substr(
		loc_straightLine_gene2 + 1, loc_colon_gene2 - loc_straightLine_gene2 - 1);
	tmpFusionBreakPoint_1 = tmpFusionGeneInfoStr_1.substr(
		loc_colon_gene1 + 1, loc_rightPar_gene1 - loc_colon_gene1 - 1);
	tmpFusionBreakPoint_2 = tmpFusionGeneInfoStr_2.substr(
		loc_colon_gene2 + 1, loc_rightPar_gene2 - loc_colon_gene2 - 1);
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputReformedJuncFastaPath outputUniqueBreakPointFusionPath " << endl;
		exit(1);
	}
	string inputReformedFusionFastaPath = argv[1];
	string outputUniqueFusionFastaPath = argv[2];
	ifstream reformedFusionFa_ifs(inputReformedFusionFastaPath.c_str());
	ofstream uniqueFusionFa_ofs(outputUniqueFusionFastaPath.c_str());
	
	// string lastFusion_geneName_1;
	// string lastFusion_geneName_2;
	// string lastFusion_chrName_1;
	// string lastFusion_chrName_2;
	// string lastFusion_breakPointStr_1;
	// string lastFusion_breakPointStr_2;

	set<string> existingFusionIDset;
	while(!reformedFusionFa_ifs.eof())
	{
		string tmpFusionFaStr_name;
		getline(reformedFusionFa_ifs, tmpFusionFaStr_name);
		if((reformedFusionFa_ifs.eof())||(tmpFusionFaStr_name == ""))
			break;
		string tmpFusionFaStr_seq;
		getline(reformedFusionFa_ifs, tmpFusionFaStr_seq);

		string tmpFusionGeneName_1;
		string tmpFusionGeneName_2;
		string tmpFusionChrName_1; 
		string tmpFusionChrName_2;
		string tmpFusionBreakPoint_1;
		string tmpFusionBreakPoint_2;	
		extractFusionGeneNameBreakPointFrom_fusionFa1stLine(
			tmpFusionGeneName_1, tmpFusionGeneName_2,
			tmpFusionChrName_1, tmpFusionChrName_2,
			tmpFusionBreakPoint_1, tmpFusionBreakPoint_2,				
			tmpFusionFaStr_name);
		string tmpFusionIDstr = tmpFusionGeneName_1 
			+ "_" + tmpFusionChrName_1 + "_" + tmpFusionBreakPoint_1
			+ "_" + tmpFusionGeneName_2
			+ "_" + tmpFusionChrName_2 + "_" + tmpFusionBreakPoint_2;
		//cout << "tmpFusionIDstr: " << endl << tmpFusionIDstr << endl;
		// if((tmpFusionGeneName_1 != lastFusion_geneName_1)
		// 	||(tmpFusionGeneName_2 != lastFusion_geneName_2)
		// 	||(tmpFusionChrName_1 != lastFusion_chrName_1)
		// 	||(tmpFusionChrName_2 != lastFusion_chrName_2)
		// 	||(tmpFusionBreakPoint_1 != lastFusion_breakPointStr_1)
		// 	||(tmpFusionBreakPoint_2 != lastFusion_breakPointStr_2))
		if(existingFusionIDset.find(tmpFusionIDstr) 
			== existingFusionIDset.end())
		{	
			//cout << "notFound" << endl;
			existingFusionIDset.insert(tmpFusionIDstr);
			// lastFusion_geneName_1 = tmpFusionGeneName_1;
			// lastFusion_geneName_2 = tmpFusionGeneName_2;
			// lastFusion_chrName_1 = tmpFusionChrName_1;
			// lastFusion_chrName_2 = tmpFusionChrName_2;
			// lastFusion_breakPointStr_1 = tmpFusionBreakPoint_1;
			// lastFusion_breakPointStr_2 = tmpFusionBreakPoint_2;
			uniqueFusionFa_ofs << ">" << tmpFusionIDstr << endl
				<< tmpFusionFaStr_seq << endl;
		}
		else
		{
			//cout << "found" << endl;
			//exit(1);
		}
	}
	reformedFusionFa_ifs.close();
	uniqueFusionFa_ofs.close();
	return 0;
}