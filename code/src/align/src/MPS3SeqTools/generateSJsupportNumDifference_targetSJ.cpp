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

#include "general/read_block_test.h"
#include "general/index_info.h"
#include "general/splice_info.h"
#include "general/alignmentToJunc_supportNum.h"

using namespace std;

int getSJdistanceFromSJstr_gt(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);	

	string donerEndPosStr = SJstr.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
	//cout << "donerEndPosStr: " << donerEndPosStr << endl;
	string acceptorStartPosStr = SJstr.substr(tabLocation2+1, tabLocation3- tabLocation2-1);
	//cout << "acceptorStartPosStr: " << acceptorStartPosStr << endl;
	int donerEndPos = atoi(donerEndPosStr.c_str());
	//cout << "donerEndPos: " << donerEndPos << endl;
	int acceptorStartPos = atoi(acceptorStartPosStr.c_str());
	//cout << "acceptorStartPos: " << acceptorStartPos << endl;
	int SJdistance = acceptorStartPos - donerEndPos - 1;	
	return SJdistance;
}

int getSJdistanceFromSJstr(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);	

	string donerEndPosStr = SJstr.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
	//cout << "donerEndPosStr: " << donerEndPosStr << endl;
	string acceptorStartPosStr = SJstr.substr(tabLocation2+1, tabLocation3- tabLocation2-1);
	//cout << "acceptorStartPosStr: " << acceptorStartPosStr << endl;
	int donerEndPos = atoi(donerEndPosStr.c_str());
	//cout << "donerEndPos: " << donerEndPos << endl;
	int acceptorStartPos = atoi(acceptorStartPosStr.c_str());
	//cout << "acceptorStartPos: " << acceptorStartPos << endl;
	int SJdistance = acceptorStartPos - donerEndPos - 1;	
	//cout << "SJdistance: " << SJdistance << endl;
	return SJdistance;
}

int getSJsupportNumFromSJstr(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);
	int tabLocation4 = SJstr.find('\t', tabLocation3+1);
	int tabLocation5 = SJstr.find('\t', tabLocation4+1);
	string supportNumStr = SJstr.substr(tabLocation4+1, tabLocation5- tabLocation4-1);
	int supportNum = atoi(supportNumStr.c_str());
	return supportNum;
}

int main(int argc, char** argv)
{
	if(argc < 7)
	{
		cout << "Executable inputIndexFile outputSJsupportNumFile inputSJfileName inputSJfileCompared2 aligner_1 inputSamSJ_1 ..." << endl;
		exit(1);
	}
	int min_intron_size = 50;
	int max_intron_size = 300000;

	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());	
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	parameter_ifs.close();

	string outputSJsupportNumFile = argv[2];
	ofstream output_ofs(outputSJsupportNumFile.c_str());
	string outputSupportNumCompareResultsFile = outputSJsupportNumFile + ".2plotInR";
	ofstream output_2plotInR_ofs(outputSupportNumCompareResultsFile.c_str());

	string inputGroundTruthSJName = argv[3];
	string inputGroundTruthSJfile = argv[4];
	//ifstream inputGroundTruthSJ_ifs(inputGroundTruthSJfile.c_str());

	vector<string> alignerNameVec;
	vector<string> inputSJfileVec;
	vector<AlignmentToJunc_supportNum_Info*> alignment2juncInfoVec;
	cout << "*****************************************************" << endl;
	output_ofs << "chr\tdonerEndPos\tacceptorStartPos\t" << inputGroundTruthSJName;
	for(int tmp = 5; tmp < argc; )
	{
		string tmpAlignName = argv[tmp];
		string tmpInputSJfile = argv[tmp+1];
		cout << "tmpAligner: " << endl << tmpAlignName << endl;
		output_ofs << "\t" << tmpAlignName;
		cout << "tmpInputSJfile: " << endl << tmpInputSJfile << endl;
		AlignmentToJunc_supportNum_Info* tmpAlignment2JuncInfo = new AlignmentToJunc_supportNum_Info();
		int tmpJunc_validNum, tmpJunc_invalidNum;
		tmpAlignment2JuncInfo->initiateAlignmentToJuncInfo(chromNum);
		tmpAlignment2JuncInfo->getAlignmentToJuncInfoFromSJfile(tmpInputSJfile, indexInfo);
		cout << "tmpJuncSize: " << tmpAlignment2JuncInfo->returnJuncInfoSize() << endl;
		alignerNameVec.push_back(tmpAlignName);
		inputSJfileVec.push_back(tmpInputSJfile);
		alignment2juncInfoVec.push_back(tmpAlignment2JuncInfo);
		tmp += 2;
	}
	output_ofs << endl;	

	int alignerNum = alignerNameVec.size();

	FILE *fp_groundTruth = fopen(inputGroundTruthSJfile.c_str(), "r");
	int juncNum = 0;

	double* supportNumDiffPercSquareSum = (double*)malloc((alignerNum) * sizeof(int));;

	while(!feof(fp_groundTruth))
	{
		char entry[500];
		fgets(entry, sizeof(entry), fp_groundTruth);
		if(feof(fp_groundTruth))
			break;
		string entryString = entry;
		int tabLocation1 = entryString.find('\t', 0);
		int tabLocation2 = entryString.find('\t', tabLocation1+1);
		int tabLocation3 = entryString.find('\t', tabLocation2+1);
		int tabLocation4 = entryString.find('\t', tabLocation3+1);
		int tabLocation5 = entryString.find('\t', tabLocation4+1);

		string chrIntString = entryString.substr(0, tabLocation1);
		string spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		string spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		int chrInt = indexInfo->convertStringToInt(chrIntString);
		int spliceStartPos = atoi(spliceStartPosString.c_str());
		int spliceEndPos = atoi(spliceEndPosString.c_str());
		string supportNumString = entryString.substr(tabLocation4+1, tabLocation5-tabLocation4-1);
		int supportNum_groundTruth = atoi(supportNumString.c_str());
		output_ofs << chrIntString << "\t" << spliceStartPosString << "\t" 
			<< spliceEndPosString << "\t" << supportNumString;
		juncNum ++;

		int tmpSupportNum_1st 
			= alignment2juncInfoVec[0]->returnSupportNumInAlignmentToJuncInfo(
				chrInt, spliceStartPos, spliceEndPos);
		output_ofs << "\t" << tmpSupportNum_1st;
		for(int tmp = 1; tmp < alignment2juncInfoVec.size(); tmp++)
		{
			int tmpSupportNum 
				= alignment2juncInfoVec[tmp]->returnSupportNumInAlignmentToJuncInfo(
					chrInt, spliceStartPos, spliceEndPos);
			int tmpSupportNum_diff = tmpSupportNum - tmpSupportNum_1st;
			double tmpSupportNum_diff_perc;
			if(tmpSupportNum_1st == 0)
			{
				if(tmpSupportNum_diff > 0)
					tmpSupportNum_diff_perc = 1.0;
				else
					tmpSupportNum_diff_perc = 0.0;
			}
			else
			{
				tmpSupportNum_diff_perc = (double)tmpSupportNum_diff/(double)tmpSupportNum_1st;
				if(tmpSupportNum_diff_perc > 1.0)
					tmpSupportNum_diff_perc = 1.0;
				else if(tmpSupportNum_diff_perc < -1.0)
					tmpSupportNum_diff_perc = -1.0;
				else
				{}
			}
			supportNumDiffPercSquareSum[tmp] 
				= supportNumDiffPercSquareSum[tmp] + tmpSupportNum_diff_perc * tmpSupportNum_diff_perc;

			output_2plotInR_ofs << juncNum << "\t" << tmpSupportNum_diff_perc << "\t" << alignerNameVec[tmp] << endl;
			output_ofs << "\t" << tmpSupportNum;	
		}
		output_ofs << endl;
	}

	fclose(fp_groundTruth);
	output_ofs << endl << "supportNum SSE: " << endl;
	for(int tmp = 0; tmp < alignment2juncInfoVec.size(); tmp++)
	{	
		output_ofs << alignerNameVec[tmp] << ": " << supportNumDiffPercSquareSum[tmp]/(double)juncNum << endl;
		delete alignment2juncInfoVec[tmp];
	}
	output_2plotInR_ofs.close();
	output_ofs.close();
	return 0;
}