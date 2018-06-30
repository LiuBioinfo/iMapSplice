// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: backSplice read count from all samples
// format: chrName, chrPos_right, chrPos_left, rawReadCount_1, rawReadCount_2, .....
// output: chrName, chrPos_right, chrPos_left, rawReadCount_1, rawReadCount_2, ... rawReadCount_N, 
//			normalizedReadCount_1, normalizedReadCount_2, normalizedReadCount_N
//          foldChange_1:2, foldChange_1:3

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
	if(argc < 5)
	{
		cout << "Executable inputRawReadCount outputNormalizedReadCountFoldChangePath ";
		cout << " normalizeFactor_1 normalizeFactor_2 ...." << endl;
		exit(1);
	}
	int sampleNum = argc - 3;
	string inputRawReadCountPath = argv[1];
	string outputNormalizedReadCountFoldChangePath = argv[2];
	vector<double> normalizeFactorVec;
	for(int tmpSample = 0; tmpSample < sampleNum; tmpSample++)
	{
		string tmpNormalizeFactorStr = argv[3+tmpSample];
		double tmpNormalizeFactorDouble = atof(tmpNormalizeFactorStr.c_str());
		cout << "tmpSample: " << tmpSample << " tmpNormalizeFactor: " << tmpNormalizeFactorDouble << endl;
		normalizeFactorVec.push_back(tmpNormalizeFactorDouble);
	}

	ifstream rawReadCount_ifs(inputRawReadCountPath.c_str());
	ofstream normalizedReadCount_foldChange_ofs(outputNormalizedReadCountFoldChangePath.c_str());
	while(!rawReadCount_ifs.eof())
	{
		string tmpBackSpliceSiteReadCountStr;
		getline(rawReadCount_ifs, tmpBackSpliceSiteReadCountStr);
		if((rawReadCount_ifs.eof())&&(tmpBackSpliceSiteReadCountStr == ""))
			break;
		vector<string> tmpBackSpliceSiteReadCountFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < sampleNum + 2; tmp++)
		{
			int tabLoc = tmpBackSpliceSiteReadCountStr.find("\t", startLoc);
			string tmpBackSpliceSiteReadCountFieldStr = tmpBackSpliceSiteReadCountStr.substr(
				startLoc, tabLoc - startLoc);
			tmpBackSpliceSiteReadCountFieldVec.push_back(tmpBackSpliceSiteReadCountFieldStr);
			//cout << "tmpBackSpliceSiteReadCountFieldStr: " << tmpBackSpliceSiteReadCountFieldStr << endl;
			startLoc = tabLoc + 1;
		}
		tmpBackSpliceSiteReadCountFieldVec.push_back(tmpBackSpliceSiteReadCountStr.substr(startLoc));
		// string tmpBackSplice_chrNameStr = tmpBackSpliceSiteReadCountFieldVec[0];
		// //cout << "tmpBackSplice_chrNameStr: " << tmpBackSplice_chrNameStr << endl;
		// string tmpBackSplice_chrPosStr_right = tmpBackSpliceSiteReadCountFieldVec[1];
		// //cout << "tmpBackSplice_chrPosStr_right: " << tmpBackSplice_chrPosStr_right << endl;
		// string tmpBackSplice_chrPosStr_left = tmpBackSpliceSiteReadCountFieldVec[2];
		// //cout << "tmpBackSplice_chrPosStr_left: " << tmpBackSplice_chrPosStr_left << endl;
		normalizedReadCount_foldChange_ofs	<< tmpBackSpliceSiteReadCountStr;		
		vector< double > tmpNormalizedReadCountVec;
		for(int tmpSample = 0; tmpSample < sampleNum; tmpSample++)
		{
			string tmpSampleRawReadCountStr = tmpBackSpliceSiteReadCountFieldVec[3+tmpSample];
			int tmpSampleRawReadCountInt = atoi(tmpSampleRawReadCountStr.c_str());
			double tmpSampleNormalizedReadCountDouble 
				= (double)tmpSampleRawReadCountInt * normalizeFactorVec[tmpSample];
			tmpNormalizedReadCountVec.push_back(tmpSampleNormalizedReadCountDouble);
			normalizedReadCount_foldChange_ofs << "\t" << tmpSampleNormalizedReadCountDouble;
		}
		// vector< double > tmpFoldChangeVec_normalizedReadCount;
		// for(int tmpSample = 1; tmpSample < sampleNum; tmpSample++)
		// {
		// 	double tmpFoldChange 
		// 		= tmpNormalizedReadCountVec[0]/tmpNormalizedReadCountVec[tmpSample];
		// 	tmpFoldChangeVec_normalizedReadCount.push_back(tmpFoldChange);
		// 	normalizedReadCount_foldChange_ofs << "\t" << tmpFoldChange;
		// }
		normalizedReadCount_foldChange_ofs << endl;
	}
	rawReadCount_ifs.close();
	normalizedReadCount_foldChange_ofs.close();
	return 0;
}