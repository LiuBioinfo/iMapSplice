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
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputFusionBreakPointHashPath inputSEsamFilePath outputFolderPath" << endl;
		exit(1);
	}
	int buffer_breakPointSearch = 4;

	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_log = outputFolderStr + "log.txt";
	ofstream log_ofs(output_log.c_str());

	string globalIndexStr = argv[1];
	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	log_ofs << "finish loading chromosomes" << endl;

	string outputFusionBreakPointPath = outputFolderStr + "fusionBreakPoint.txt";
	string inputFusionBreakPointHashPath = argv[2];
	FusionBreakPointHash_Info* fusionBreakPointHashInfo = new FusionBreakPointHash_Info();
	fusionBreakPointHashInfo->initiateWithChromNum(chromNum);
	fusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromFuionBreakPointFile(
		inputFusionBreakPointHashPath, indexInfo)
	// output fusionBreakPointInfo
	fusionBreakPointHashInfo->outputFusionBreakPointHashInfoStr(outputFusionBreakPointPath, indexInfo);

	string keptSAM_path = outputFolderStr + "kept.sam";
	string unfixedFusionSAM_path = outputFolderStr + "unfixedFusion.sam";
	string fixedFusionSAM_path = outputFolderStr + "fixedFusion.sam";
	ofstream keptSAM_ofs(keptSAM_path.c_str());
	ofstream unfixedFusionSAM_ofs(unfixedFusionSAM_path.c_str());
	ofstream fixedFusionSAM_ofs(fixedFusionSAM_path.c_str());

	string inputSEsamPath = argv[3];
	ifstream seSAM_ifs(inputSEsamPath.c_str());
	while(!seSAM_ifs.eof())
	{
		string tmpSAMstr;
		getline(seSAM_ifs, tmpSAMstr);
		if((seSAM_ifs.eof())||(tmpSAMstr == ""))
			break;
		if(tmpSAMstr.at(0) == '@')
			continue;
		IncompleteSeAlignment2remapAgainstFusionBreakPoint_Info* incompleteSEalignInfo 
			= new IncompleteSeAlignment2remapAgainstFusionBreakPoint_Info();
		bool initiateSeAlignInfoBool 
			= incompleteSEalignInfo->initiateWithSamStr_remapUnfixedEnd2fusionBreakPoint(
				tmpSAMstr, indexInfo);
		if(initiateSeAlignInfoBool)
		{
			bool Nor_or_Rcm_bool = incompleteSEalignInfo->NorOrRcmBool();
			bool headUnfixed_bool = incompleteSEalignInfo->headUnfixed();
			bool tailUnfixed_bool = incompleteSEalignInfo->tailUnfixed();
			if(Nor_or_Rcm_bool)
			{
				if(headUnfixed_bool)
				{
					int unfixedHeadLen = incompleteSEalignInfo->returnUnfixedHeadLength();
					int oriStartPos = incompleteSEalignInfo->returnStartPos();
				}
				if(tailUnfixed_bool)
				{
					int unfixedTailLen = incompleteSEalignInfo->returnUnfixedTailLength();
					int oriEndPos = incompleteSEalignInfo->returnEndPos();
				}
			}
			else
			{
				if(headUnfixed_bool)
				{
					int unfixedHeadLen = incompleteSEalignInfo->returnUnfixedHeadLength();

				}
				if(tailUnfixed_bool)
				{
					int unfixedTailLen = incompleteSEalignInfo->returnUnfixedTailLength();

				}
			}
		}	
		else
		{
			keptSAM_ofs << tmpSAMstr << endl;
		}
		delete incompleteSEalignInfo;	
	}

	seSAM_ifs.close();
	fusionBreakPointHashInfo.memoryFree();
	delete fusionBreakPointHashInfo;
	delete indexInfo;
	parameter_file_ifs.close();
	chrom_bit_file_ifs.close();
	keptSAM_ofs.close();
	unfixedFusionSAM_ofs.close();
	fixedFusionSAM_ofs.close();
	return 0;
}