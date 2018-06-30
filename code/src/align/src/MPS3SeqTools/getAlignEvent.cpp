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

#include "general/read_block_test.h"
#include "general/index_info.h"
#include "general/alignEvent_info.h"
#include "general/alignEventSet_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <InputSam> <OutputFilePrefix>" << endl;
		exit(0);
	}
	bool recordFollowingJumpCode_bool = true;
	bool recordFollowingSegsMismatch_bool = true;
	bool fasta_or_fastq_bool = false;

	string inputIndexStr = argv[1]; inputIndexStr += "/";	
	string inputSamStr = argv[2];
	ifstream inputAlignment_ifs(inputSamStr.c_str());

	string OutputFilePrefix = argv[3];
	string OutputFile_SJ = OutputFilePrefix + ".SJ";
	ofstream SJ_ofs(OutputFile_SJ.c_str());
	string log_file_str = OutputFilePrefix + ".process.log";
	ofstream log_ofs(log_file_str.c_str());


	string parameter_file = inputIndexStr; 
	parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	int readTotalNum = 0;

	int normalRecordNum = 1000000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	AlignEventSet_Info* alignEventSetInfo = new AlignEventSet_Info(indexInfo->returnChromNum());

	vector<string> alignmentStrVec(normalRecordNum);
	vector< AlignEvent_Info* > alignEventInfoVec(normalRecordNum);
	log_ofs << endl << "start to extract event from alignemnt" << endl << endl;
	for(tmpTurn = 0; ; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		realRecordNum = normalRecordNum;
		log_ofs << endl <<  "start to read alignment from SAM file, turn: " << tmpTurn+1 << endl;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; )
		{
    		if(inputAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		string tmpAlignStrLine;
    		getline(inputAlignment_ifs, tmpAlignStrLine); // readName_1
    		if(inputAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		} 
    		if(tmpAlignStrLine.at(0) == '@')
    		{
    			continue;
    		}
    		else
    		{
    			alignmentStrVec[recordNumTmp] = tmpAlignStrLine;
    			recordNumTmp++;
    		}    		 
		}

		log_ofs << "finish reading alignment form SAM file, turn: " << tmpTurn+1 << endl;
		log_ofs << "realRecordNum: " << realRecordNum << endl;
		readTotalNum += realRecordNum;

		log_ofs << "start to extract event from alignment, turn: " << tmpTurn+1 << endl;
		
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			AlignEvent_Info* alignEventInfo = new AlignEvent_Info();
			alignEventInfo->readAlignInfoFromSAM(alignmentStrVec[tmp], 
				recordFollowingSegsMismatch_bool, fasta_or_fastq_bool);
			alignEventInfoVec[tmp] = alignEventInfo;
		}	
		
		log_ofs << "finish extracting event from alignment, turn: " << tmpTurn+1 << endl;
		log_ofs << "start to insert event info into eventSet, turn: " << tmpTurn+1 << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			alignEventSetInfo->insertEventInfo(alignEventInfoVec[tmp], indexInfo, 
				recordFollowingJumpCode_bool, recordFollowingSegsMismatch_bool, fasta_or_fastq_bool);
		}

		log_ofs << "finish inserting all event info into event set, turn: " << tmpTurn+1 << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			delete alignEventInfoVec[tmp];
		}		
	}

	bool output_followingJumpCode_bool = true;
	bool output_followingJumpCodeMismatch_bool = true;
	log_ofs << "start to output eventSetInfo" << endl;
	alignEventSetInfo->outputEventSet(SJ_ofs, output_followingJumpCode_bool, output_followingJumpCodeMismatch_bool);
	log_ofs << "finish output eventSetInfo, turn: " << tmpTurn + 1 << endl << endl;;


	inputAlignment_ifs.close();
	SJ_ofs.close();
	log_ofs.close();
	parameter_file_ifs.close();	

	delete alignEventSetInfo;
	delete indexInfo;
	return 0;
}
