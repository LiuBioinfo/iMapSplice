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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

void assignWithMaxScore_vec(vector<int>&tmpScoreVec, vector<int>& tmpAssignedAluFamilyIndexVec)
{
	int tmpScore_max = -1;
	for(int tmp = 0; tmp < tmpScoreVec.size(); tmp++)
	{
		int tmpScore = tmpScoreVec[tmp];
		if(tmpScore > tmpScore_max)
		{
			tmpScore_max = tmpScore;
		}
	}
	for(int tmp = 0; tmp < tmpScoreVec.size(); tmp++)
	{
		int tmpScore = tmpScoreVec[tmp];
		if(tmpScore == tmpScore_max)
		{
			tmpAssignedAluFamilyIndexVec.push_back(tmp);;
		}
	}
}

int assignWithMaxScore(vector<int>& tmpScoreVec)
{
	int tmpScore_max = -1;
	int tmpIndex_max = -1;
	for(int tmp = 0; tmp < tmpScoreVec.size(); tmp++)
	{
		int tmpScore = tmpScoreVec[tmp];
		if(tmpScore > tmpScore_max)
		{
			tmpScore_max = tmpScore;
			tmpIndex_max = tmp;
		}
	}
	return tmpIndex_max;
}	

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputAluFamilyListPath segMapResultsFolder " << endl;
		exit(1);
	}
	string inputAluFamilyListPath = argv[1];
	string outputFolder = argv[2];
   	cout << "start to read alu family list" << endl;
   	vector<string> aluFamilyVec;
   	ifstream aluFamilyList_ifs(inputAluFamilyListPath.c_str());
   	while(!aluFamilyList_ifs.eof())
   	{
   		string tmpAluFamilyStr;
   		getline(aluFamilyList_ifs, tmpAluFamilyStr);
   		if(tmpAluFamilyStr == "")
   			break;
   		aluFamilyVec.push_back(tmpAluFamilyStr);
   	}
   	aluFamilyList_ifs.close();

   	cout << "start to read segMap resutls ......." << endl;

   	string tmpMergedOutputFolder = outputFolder + "/merged";
   	string cmd_mk_mergedFolder = "mkdir " + tmpMergedOutputFolder;
   	system(cmd_mk_mergedFolder.c_str());
   	//vector<string> aluFamilySegMapResultFileVec;
   	string toMergeScoreFile = "/toMergeScore.txt";
   	int aluFamilyNum = aluFamilyVec.size();
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   	{
   		string tmpAluFamilyName = aluFamilyVec[tmp];
   		string tmpAluFamilySegMapResultFilePath = outputFolder + "/" + tmpAluFamilyName + "/score.txt";
   		string tmpToMergeScoreFile = outputFolder + "/" + tmpAluFamilyName + toMergeScoreFile;
   		string cmd_cutToMergedScore = "cut -f 2 " + tmpAluFamilySegMapResultFilePath + " > " + tmpToMergeScoreFile;
   		system(cmd_cutToMergedScore.c_str());
   	}	

   	string toMerge_readName_file = tmpMergedOutputFolder + "/toMerge_readName.txt";
   	string cmd_cutToMergedReadName = "cut -f 1 " + outputFolder + "/" + aluFamilyVec[0] + "/score.txt > " + toMerge_readName_file;
   	system(cmd_cutToMergedReadName.c_str());

   	string cmd_paste_toMakeMergedScoreFile = "paste " + toMerge_readName_file;
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   	{
   		string tmpAluFamilyName = aluFamilyVec[tmp];
   		string tmpToMergeScoreFile = outputFolder + "/" + tmpAluFamilyName + toMergeScoreFile;
   		cmd_paste_toMakeMergedScoreFile += " ";
   		cmd_paste_toMakeMergedScoreFile += tmpToMergeScoreFile;
   	}
   	string mergedScore_file = tmpMergedOutputFolder + "/merged_score.txt";
   	cmd_paste_toMakeMergedScoreFile = cmd_paste_toMakeMergedScoreFile + " > " + mergedScore_file;
   	system(cmd_paste_toMakeMergedScoreFile.c_str());

   	cout << "start to process merged file to select the best alu family for each read and generate the final results" << endl;
   	string final_assignment_file = tmpMergedOutputFolder + "/read_assignment_final.txt";
   	ofstream final_assignment_ofs(final_assignment_file.c_str());
   	final_assignment_ofs << "Read_Name";
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   		final_assignment_ofs << "\t" << aluFamilyVec[tmp];
   	final_assignment_ofs << "\tAssigned_Family\tScore_Max" << endl;


   	vector<double> aluFamilyCountVec;
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   		aluFamilyCountVec.push_back(0);

   	ifstream mergedScore_ifs(mergedScore_file.c_str());
   	while(!mergedScore_ifs.eof())
   	{
   		string tmpStr;
   		getline(mergedScore_ifs, tmpStr);
   		if(tmpStr == "")
   			break;
   		int tabLoc_1 = tmpStr.find("\t");
   		string tmpReadName = tmpStr.substr(0, tabLoc_1);
   		vector<int> tmpScoreVec;
   		int startLoc = tabLoc_1 + 1;
   		for(int tmp = 1; tmp < aluFamilyNum; tmp++)
   		{
   			int tabLoc = tmpStr.find("\t", startLoc);
   			string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
   			int tmpScore = atoi(tmpField.c_str());
   			tmpScoreVec.push_back(tmpScore);
   			startLoc = tabLoc + 1;
   		}
   		string tmpLastScoreStr = tmpStr.substr(startLoc);
   		int tmpLastScore = atoi(tmpLastScoreStr.c_str());
   		tmpScoreVec.push_back(tmpLastScore);
   		vector<int> tmpAssignedAluFamilyIndexVec;
   		assignWithMaxScore_vec(tmpScoreVec, tmpAssignedAluFamilyIndexVec);
   		if(tmpAssignedAluFamilyIndexVec.size() == 0)
   		{
   			cout << "error in assignedAluFamilyIndex ! " << endl;
   			exit(1);
   		}
   		final_assignment_ofs << tmpStr << "\t";
   		double tmpScoreToAdd = 1;
   		if(tmpAssignedAluFamilyIndexVec.size() > 1)
   			tmpScoreToAdd = double(1) / tmpAssignedAluFamilyIndexVec.size();
   		for(int tmp = 0; tmp < tmpAssignedAluFamilyIndexVec.size(); tmp++)
   		{
   			int tmpAssignedAluFamilyIndex = tmpAssignedAluFamilyIndexVec[tmp];
   			final_assignment_ofs << aluFamilyVec[tmpAssignedAluFamilyIndex] << ",";
   			double originalScore = aluFamilyCountVec[tmpAssignedAluFamilyIndex];
   			double updatedScore = originalScore + tmpScoreToAdd;
   			aluFamilyCountVec[tmpAssignedAluFamilyIndex] = updatedScore;
   		}
   		int tmpAssignedAluFamilyIndex_1st = tmpAssignedAluFamilyIndexVec[0];
   		final_assignment_ofs << "\t" << tmpScoreVec[tmpAssignedAluFamilyIndex_1st] << endl;
   	}
   	mergedScore_ifs.close();
   	final_assignment_ofs.close();
   	
   	string final_aluFamilyReadCount_file = tmpMergedOutputFolder + "/readCount.txt";
   	ofstream readCount_ofs(final_aluFamilyReadCount_file.c_str());
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   	{
   		readCount_ofs << aluFamilyVec[tmp] << "\t" << aluFamilyCountVec[tmp] << endl;
   	}
   	readCount_ofs.close();
	return 0;
}