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
#include <math.h>
#include "../general/read_block_test.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc <= 4)
	{
		cout << "Executable inputBedtoolsResultsFile OutputFolder " << endl;
		cout << "DataSetPair_1_control Normalizer_1_control DataSetPair_1_compare Normalizer_1_compare" << endl;
		cout << "(DataSetPair_2_control Normalizer_2_control DataSetPair_2_compare Normalizer_2_compare ...)" << endl;
		exit(1);
	}
	int fieldNumBeforeDataField = 12;

	string inputBedtoolsResultsFileStr = argv[1];
	ifstream bedtoolsResults_ifs(inputBedtoolsResultsFileStr.c_str());
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	string outputFolderStr_coverageCompare = outputFolderStr + "coverageCompare_allData.txt";
	ofstream coverageCompare_ofs(outputFolderStr_coverageCompare.c_str());	
	string outputFolderStr_coverageCompare_log10 = outputFolderStr + "coverageCompare_allData_log10.txt";
	ofstream coverageCompare_log10_ofs(outputFolderStr_coverageCompare_log10.c_str());
	string outputFolderStr_coverageCompare_group = outputFolderStr + "coverageCompare_allData_group.txt";
	ofstream coverageCompare_group_ofs(outputFolderStr_coverageCompare_group.c_str());
	string outputFolderStr_coverageCompare_group_log10 = outputFolderStr + "coverageCompare_allData_group_log10.txt";
	ofstream coverageCompare_group_log10_ofs(outputFolderStr_coverageCompare_group_log10.c_str());

	vector<string> DataSetPair_control_vec;
	vector<double> Normalizer_control_vec;
	vector<string> DataSetPair_compare_vec;
	vector<double> Normalizer_compare_vec;
	for(int tmp = 3; tmp < argc; tmp += 4)
	{
		string tmpDataSetControlName = argv[tmp];
		DataSetPair_control_vec.push_back(tmpDataSetControlName);
		string tmpNormalizerControlStr = argv[tmp+1];
		double tmpNormalizerControlDouble = atof(tmpNormalizerControlStr.c_str());
		Normalizer_control_vec.push_back(tmpNormalizerControlDouble);
		string tmpDataSetCompareName = argv[tmp+2];
		DataSetPair_compare_vec.push_back(tmpDataSetCompareName);
		string tmpNormalizerCompareStr = argv[tmp+3];
		double tmpNormalizerCompareDouble = atof(tmpNormalizerCompareStr.c_str());
		Normalizer_compare_vec.push_back(tmpNormalizerCompareDouble);		
	}
	int DataSetNum = DataSetPair_control_vec.size() * 2;
	int groupNum = DataSetNum/2;
	if(DataSetNum != groupNum * 2)
	{
		cout << "error in data number, no way to separate as control and compare groups..." << endl;
		exit(1);
	}

	while(!(bedtoolsResults_ifs.eof()))
	{
		string bedtoolsStr;
		getline(bedtoolsResults_ifs, bedtoolsStr);
		if(bedtoolsResults_ifs.eof()||(bedtoolsStr == ""))
			break; 
	
		vector<string> bedtoolsFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < (fieldNumBeforeDataField + DataSetNum - 1); tmp++)
		{
			int tabLoc = bedtoolsStr.find("\t", startLoc);
			string tmpBedtoolsField = bedtoolsStr.substr(startLoc, tabLoc-startLoc);
			bedtoolsFieldVec.push_back(tmpBedtoolsField);
			startLoc = tabLoc + 1;			
		}
		bedtoolsFieldVec.push_back(bedtoolsStr.substr(startLoc));
		
		vector<double> dataSetCoverageVec_control_normalized;
		vector<double> dataSetCoverageVec_compare_normalized;
		for(int tmp = fieldNumBeforeDataField; tmp < bedtoolsFieldVec.size(); tmp+=2)
		{
			int dataSetPairIndex = (tmp-fieldNumBeforeDataField)/2;
			int tmpDataSetCoverage_control = atoi(bedtoolsFieldVec[tmp].c_str());
			int tmpDataSetCoverage_compare = atoi(bedtoolsFieldVec[tmp+1].c_str());
			double tmpDataSetCoverage_control_normalized 
				= ((double)tmpDataSetCoverage_control)*((double)Normalizer_control_vec[dataSetPairIndex]);
			double tmpDataSetCoverage_compare_normalized 
				= ((double)tmpDataSetCoverage_compare)*((double)Normalizer_compare_vec[dataSetPairIndex]);
			dataSetCoverageVec_control_normalized.push_back(tmpDataSetCoverage_control_normalized);
			dataSetCoverageVec_compare_normalized.push_back(tmpDataSetCoverage_compare_normalized);

			double tmpDataSetCoverage_control_normalized_log10 = log10(tmpDataSetCoverage_control_normalized);
			double tmpDataSetCoverage_compare_normalized_log10 = log10(tmpDataSetCoverage_compare_normalized);

			coverageCompare_ofs << tmpDataSetCoverage_control_normalized << "\t"
				<< tmpDataSetCoverage_compare_normalized << "\t";
			coverageCompare_group_ofs << tmpDataSetCoverage_control_normalized << "\t"
				<< tmpDataSetCoverage_compare_normalized << "\t"
				<< DataSetPair_control_vec[dataSetPairIndex] << "_vs_" 
				<< DataSetPair_compare_vec[dataSetPairIndex] << endl;
			coverageCompare_log10_ofs << tmpDataSetCoverage_control_normalized_log10 << "\t"
				<< tmpDataSetCoverage_compare_normalized_log10 << "\t";
			coverageCompare_group_log10_ofs << tmpDataSetCoverage_control_normalized_log10 << "\t"
				<< tmpDataSetCoverage_compare_normalized_log10 << "\t"
				<< DataSetPair_control_vec[dataSetPairIndex] << "_vs_" 
				<< DataSetPair_compare_vec[dataSetPairIndex] << endl;
		}
		coverageCompare_log10_ofs << endl;
		coverageCompare_ofs << endl;
	}	
	coverageCompare_ofs.close();
	coverageCompare_group_ofs.close();
	return 0;
}	