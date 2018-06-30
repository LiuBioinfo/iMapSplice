// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <math>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sys/types.h>    
#include <dirent.h>    
#include <stdio.h>    
#include <errno.h>
#include <set>
#include "buildIndexParameter.h"

using namespace std;  

#define Range 6
#define LONGLCP 255

typedef unsigned char BYTE;

//set<string> chrNameSet;

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		cout << "Executable <InputChromosomesFolder> <outputIndexFolder>" << endl;
		exit(0);
	}	

	buildIndex_info* buildIndexInfo = new buildIndex_info();

	string InputChrFolder = argv[1];
	string OutputIndexFolder = argv[2];

	///////////// initiate index files //////////////
	//string outputIndexFileStr = OutputIndexFolder + "/WholeGenomeIndex";
	buildIndexInfo->InputChromFolderStr = InputChrFolder;
	buildIndexInfo->OutputIndexFolderStr = OutputIndexFolder;
	
	string parameter_file = buildIndexInfo->OutputIndexFolderStr + "_parameter";
	ofstream parameter_file_ofs(parameter_file.c_str());
	string chrom_file = buildIndexInfo->OutputIndexFolderStr + "_chromMerge";
	ofstream chrom_file_ofs(chrom_file.c_str(),ios::binary);	
	//buildIndexInfo->initiateAllOutputIndexFile();

	////////////// scan all chromosomes, generate index_parameter, index_chrom //////
	//set<string> chrNameSet;
	DIR *dp;    
    struct dirent *dirp;    
    if((dp=opendir(argv[1])) == NULL)    
    {
        printf("can't open %s",argv[1]);   
    }
    while ((dirp=readdir(dp))!=NULL)    
    {     
    	if (dirp->d_type==8)
	    {  
	    	(buildIndexInfo->chrNameStrSet).insert((dirp->d_name));
    	}
    }    
    closedir(dp);    

    //cout << "chromFileSetNum: " << (buildIndexInfo->chrNameStrSet).size() << endl;
    for(set<string>::iterator tmpIter = (buildIndexInfo->chrNameStrSet).begin(); 
   		 	tmpIter != (buildIndexInfo->chrNameStrSet).end(); tmpIter++)
    {   
		//cout << "chromNameStr: " << (*tmpIter) << endl;
		(buildIndexInfo->chrNameStrVec).push_back((*tmpIter));
	}

	char chrFileLine[100];
	string chrFileLineStr;
	int lineNum = 0;
	int tmpChromEndPosInGenome = 0;
	for(int tmpChromNum = 0; tmpChromNum < (buildIndexInfo->chrNameStrVec).size();
		tmpChromNum ++)
	{
		lineNum = 0;
		string tmpChromFileName = buildIndexInfo->InputChromFolderStr 
			+ "/" + (buildIndexInfo->chrNameStrVec)[tmpChromNum];

		FILE *fp_tmpChr = fopen(tmpChromFileName.c_str(), "r");
		//cout << endl << "readFileName " << tmpChromNum + 1 << ": " << tmpChromFileName << endl;

		fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);
		chrFileLineStr = chrFileLine; 
		fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);
		chrFileLineStr = chrFileLine;
		lineNum ++;
		int chrBaseNum = 0;
		while(!feof(fp_tmpChr))
		{
			lineNum ++;
			fgets(chrFileLine, sizeof(chrFileLine), fp_tmpChr);
		
			chrFileLineStr = chrFileLine;
			
			if(feof(fp_tmpChr))
			{
				if(chrBaseNum == 0)
				{
					chrBaseNum = 50;
				}
				(buildIndexInfo->chromLengthVec).push_back((lineNum-2)*50 + chrBaseNum);
				break;		
			}
			if(chrFileLineStr.length() != 51)
			{
				chrBaseNum = chrFileLineStr.length()-1;
			}	

		}
		fclose(fp_tmpChr);

	}
	
	buildIndexInfo->getChrEndPosVec();
	buildIndexInfo->outputChromSeq(chrom_file_ofs);
	buildIndexInfo->outputIndexParameter(parameter_file_ofs);

	return 0;
}