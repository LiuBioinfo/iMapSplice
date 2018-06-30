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
	if(argc != 4)
	{
		cout << "Executable readLength_max inputFqFile outputReadNumDistributionAndTotalNumResults" << endl;
		exit(1);
	}
	string readLength_max_str = argv[1];
	int readLength_max = atoi(readLength_max_str.c_str());
	vector<int> readNumVec;
	for(int tmp = 0; tmp < readLength_max; tmp++)
		readNumVec.push_back(0);
	string inputFqFile = argv[2];
	string outputReadNumDistributionAndTotalNumResults = argv[3];

	ifstream fq_ifs(inputFqFile.c_str());
	while(!fq_ifs.eof())
	{
		string tmpStr;
		getline(fq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpStr_2, tmpStr_3, tmpStr_4;
		getline(fq_ifs, tmpStr_2);
		getline(fq_ifs, tmpStr_3);
		getline(fq_ifs, tmpStr_4);
		int tmpReadLength = tmpStr_2.length();
		readNumVec[tmpReadLength - 1]++;		
	}
	fq_ifs.close();

	ofstream readNum_ofs(outputReadNumDistributionAndTotalNumResults.c_str());
	for(int tmp = 0; tmp < readLength_max; tmp++)
		readNum_ofs << tmp+1 << "\t" << readNumVec[tmp] << endl;
	readNum_ofs.close();
	return 0;
}