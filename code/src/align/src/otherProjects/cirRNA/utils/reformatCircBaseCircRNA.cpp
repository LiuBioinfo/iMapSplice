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

#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputCircBaseReport outputReformatedCircRNAfile" << endl;
		exit(1);
	}	

	string inputFilePath = argv[1];
	string outputFilePath = argv[2];
	ifstream circBase_ifs(inputFilePath.c_str());
	ofstream reformatedCircRNA_ofs(outputFilePath.c_str());
	int tmpJuncID = 0;
	while(!circBase_ifs.eof())
	{
		string tmpStr;
		getline(circBase_ifs, tmpStr);
		if((circBase_ifs.eof())||(tmpStr == ""))
			break;	
		tmpJuncID ++;	
		int loc_1 = tmpStr.find(":");
		int loc_2 = tmpStr.find("-");
		int loc_3 = tmpStr.find("\t");
		string tmpStartPosStr = tmpStr.substr(loc_1 + 1, loc_2 - loc_1 - 1);
		string tmpEndPosStr = tmpStr.substr(loc_2 + 1, loc_3 - loc_2 - 1);
		string tmpChrName = tmpStr.substr(0, loc_1);
		reformatedCircRNA_ofs << tmpChrName << "\t" 
			<< tmpEndPosStr << "\t" << tmpStartPosStr << "\t" 
			<< "circJunc_" << tmpJuncID << endl;
	}
	circBase_ifs.close();
	reformatedCircRNA_ofs.close();
	return 0;
}
