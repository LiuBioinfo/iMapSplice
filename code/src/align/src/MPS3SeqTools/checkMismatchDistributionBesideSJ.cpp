// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//input: SAM file
//mismatch distribution beside SJ:
//mismatch[0...max_distaceInRead]

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

#include "general/read_block_test.h"
#include "general/index_info.h"
#include "general/splice_info.h"
#include "general/alignmentInfoToCheckMismatch.h"

using namespace std;

void checkMismatchBesideSJ_aligner(int tmpAligner_index, const string& alignerNameStr, 
	const string& inputSAMfileName, int mismatchDistance_max)
{
	int* mismatchDistanceArray = (int*)malloc(mismatchDistance_max * sizeof(int)); 
	for(int tmp = 0; tmp < mismatchDistance_max; tmp++)
	{
		mismatchDistanceArray[tmp] = 0;
	}

	ifstream inputSAM_ifs(inputSAMfileName.c_str());
	while(!inputSAM_ifs.eof())
	{
		string tmpSAMstr;
		getline(inputSAM_ifs, tmpSAMstr);
		if(inputSAM_ifs.eof())
			break;
		AlignmentInfoToCheckMismatch_Info tmpAlignInfo2checkMismatch;
		tmpAlignInfo2checkMismatch
	}
	inputSAM_ifs.close();
}


int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable mismatchDistributionOutputFile max_distaceInRead AlignerName_1 SAM_input_1 (AlignerName_2 SAM_input_2 ...)" << endl;
		exit(1);
	}

	vector<string> alignerNameVec;
	vector<string> inputSAMfileVec;
	cout << "*****************************************************" << endl;
	for(int tmp = 3; tmp < argc; )
	{
		string tmpAlignName = argv[tmp];
		string tmpInputSJfile = argv[tmp+1];
		cout << "tmpAligner: " << endl << tmpAlignName << endl;
		cout << "tmpInputSAMfile: " << endl << tmpInputSJfile << endl;		
		alignerNameVec.push_back(tmpAlignName);
		inputSAMfileVec.push_back(tmpInputSJfile);
		tmp += 2;
	}

	string outputFileStr = argv[1];
	ofstream output_ofs(outputFileStr.c_str());
	string maxDistanceStr = argv[2];
	int maxDistanceInRead = atoi(maxDistanceStr.c_str());

	//int* mismatchDistanceArray = (int*)malloc(maxDistanceInRead * (alignerNameVec.size())); // mismatch[ i * maxDistanceInRead + j ] jth pos mismatch Num of ith aligner 
	/*for(int tmp = 0; tmp < maxDistanceInRead * (alignerNameVec.size()); tmp++)
	{
		mismatchDistanceArray[tmp] = 0;
	}*/

	for(int tmpAligner = 0; tmpAligner < alignerNameVec.size(); tmpAligner++)
	{
		checkMismatchBesideSJ_aligner(tmpAligner, alignerNameVec[tmpAligner], inputSAMfileVec[tmpAligner], mismatchDistanceArray);
	}


	output_ofs.close();
	return 0;
}