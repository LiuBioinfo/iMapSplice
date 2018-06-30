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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputGAFfilePath outputSJfilePath" << endl;
		exit(1);
	}

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	////////////////// transcript info loading .... //////////
	string transcript_file_path = argv[2];
	cout << "start to load transcriptInfo " << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type_GAF = "GAF";
	transcriptInfo->extractTranscript(transcript_file_ifs, indexInfo, transcript_type_GAF);

	AlignInferJunctionHash_Info* juncHash = new AlignInferJunctionHash_Info();
	juncHash->initiateAlignInferJunctionInfo(chromNum);

	int transcriptNum = transcriptInfo->returnTranscriptNum();
	for(int tmp = 0; tmp < transcriptNum; tmp++)
	{
		string tmpChrName = transcriptInfo->returnTranscriptChromName(tmp);
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		string tmpStrand = transcriptInfo->returnTranscriptStrand(tmp);
		int tmpTranscriptExonNum = transcriptInfo->returnTranscriptExonNum(tmp);
		int tmpTranscriptSJnum = tmpTranscriptExonNum - 1;
		for(int tmp2 = 0; tmp2 < tmpTranscriptSJnum; tmp2++)
		{
			int tmpUpstreamExonIndex = tmp2;
			int tmpDownStreamExonIndex = tmp2+1;
			int tmpSJ_doner = transcriptInfo->returnTranscriptExonEndPosWithExonIndex(tmp, tmpUpstreamExonIndex);
			int tmpSJ_acceptor = transcriptInfo->returnTranscriptExonStartPosWithExonIndex(tmp, tmpDownStreamExonIndex);
			juncHash->insertSJ_chrNamePos_strand(tmpChrNameInt, tmpSJ_doner, tmpSJ_acceptor, tmpStrand, indexInfo);
		}
	}

	string outputSJfilePath = argv[3];
	juncHash->outputAlignInferInfoHashInfo_chrNamePos_strand(indexInfo, outputSJfilePath);

	cout << "all jobs done ..." << endl;
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}