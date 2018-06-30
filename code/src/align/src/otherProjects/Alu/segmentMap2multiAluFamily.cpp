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

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputAluFamilyListPath inputMultiIndexFolderPath inputRead_1 inputRead_2 outputFolder threadNum" << endl;
		exit(1);
	}

	string inputAluFamilyListPath = argv[1];
	string inputMultiIndexFolderPath = argv[2]; inputMultiIndexFolderPath += "/";
	string inputRead_1 = argv[3];
	string inputRead_2 = argv[4];
	string outputFolder = argv[5];
	string threadNumStr = argv[6];

   	string mkdirOutputCommand = "mkdir -p " + outputFolder;
   	system(mkdirOutputCommand.c_str());	
   	string logStr = outputFolder + "/log.txt";
   	ofstream log_ofs(logStr.c_str());

   	log_ofs << "start to read alu family list " << endl;
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
   	cout << "start to run segMap on each alu family " << endl;
   	log_ofs << "start to run segMap on each alu family " << endl;   	
   	int aluFamilyNum = aluFamilyVec.size();
   	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
   	{
   		string tmpAluFamilyName = aluFamilyVec[tmp];
   		cout << "start to run segMap on aluFamily: " << tmpAluFamilyName << endl;
   		log_ofs << "start to run segMap on aluFamily: " << tmpAluFamilyName << endl;
   		string tmpAluFamilyIndex = inputMultiIndexFolderPath + tmpAluFamilyName;
   		string tmpOutputFolder = outputFolder + "/" + tmpAluFamilyName; 
   		string cmd_segMap2singlAluFamily = "./segmentMap2singleAluFamily "
   			+ tmpAluFamilyIndex + " " + inputRead_1 + " " + inputRead_2 + " "
   			+ tmpOutputFolder + " " + threadNumStr;
   		system(cmd_segMap2singlAluFamily.c_str());
   		cout << "end of running segMap on aluFamily: " << tmpAluFamilyName << endl;
   		log_ofs << "end of running segMap on aluFamily: " << tmpAluFamilyName << endl;
   	}

   	log_ofs.close();
	return 0;
}