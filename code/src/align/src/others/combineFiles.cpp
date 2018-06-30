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
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>

#include "read_block_test.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 3)
	{
		cout << "execution <inputFile_prefix> <Num_of_files> <outputFilePath>" << endl;
		exit(1);
	}

	string inputPrefix = argv[1];
	string numStr = argv[2];
	string outputFilePath = argv[3];
	int num_of_files = atoi(numStr.c_str());

	ofstream output_ofs(outputFilePath.c_str());

	for(int tmp = 1; tmp <= num_of_files; tmp ++)
	{
		string tmpFile = inputPrefix + "." + int_to_str(tmp);
		ifstream tmp_ifs(tmpFile.c_str());
		while(!tmp_ifs.eof())
		{
			string tmpLineStr;
			getline(tmp_ifs, tmpLineStr);
			output_ofs << tmpLineStr << endl;
		}
	}

	return 0;
}