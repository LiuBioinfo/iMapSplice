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
#include "../../../general/alignInferJunctionHash_info.h"
#include "../general/SNPhash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexFolderPath inputFormattedSNPfile inputSJfile outputFilePath SJminSize SJmaxSize" << endl;
		return 0;
	}
	int toCountSNPnumMaxInSJanchor = 10;
	string SJminSizeStr = argv[5];
	string SJmaxSizeStr = argv[6];
	int min_intron_size_int = atoi(SJminSizeStr.c_str());
	int max_intron_size_int = atoi(SJmaxSizeStr.c_str());

	string anchorSNPstats_path = argv[4];
	ofstream anchorSNPstats_ofs(anchorSNPstats_path.c_str());	

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

	string inputSNPpath = argv[2];
	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_formattedSNPfile(inputSNPpath, indexInfo);

	string inputSJfile = argv[3];
	AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
	tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
	tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(inputSJfile, indexInfo);

	int totalJuncNum = tmpJuncHash->returnAlignInferInfoVecSize();
	int validJuncNum = 0;
	int invalidJuncNum = 0;
	vector<int> juncNumVec_withDifferentSNPnum;
	for(int tmp = 0; tmp < toCountSNPnumMaxInSJanchor + 1; tmp++)
		juncNumVec_withDifferentSNPnum.push_back(0);
	for(int tmp = 0; tmp < totalJuncNum; tmp++)
	{
		int tmpJunc_chrNameInt = tmpJuncHash->returnAlignInferInfo_chrNameInt(tmp);
		int tmpJunc_startPos = tmpJuncHash->returnAlignInferInfo_donerEndPos(tmp);
		int tmpJunc_endPos = tmpJuncHash->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpJunc_donerAnchorSize = tmpJuncHash->returnAnchorSizeMax_doner(tmp);
		int tmpJunc_acceptorAnchorSize = tmpJuncHash->returnAnchorSizeMax_acceptor(tmp);
		//cout << "tmpJunc_donerAnchorSize: " << tmpJunc_donerAnchorSize << endl;
		//cout << "tmpJunc_acceptorAnchorSize: " << tmpJunc_acceptorAnchorSize << endl;
		int tmpJunc_size = tmpJunc_endPos - tmpJunc_startPos - 1;
		if((tmpJunc_size > max_intron_size_int)||(tmpJunc_size < min_intron_size_int))
		{
			invalidJuncNum ++;
			continue;
		}
		validJuncNum ++;
		int donerAnchor_chrPos_start = tmpJunc_startPos - tmpJunc_donerAnchorSize + 1;
		int donerAnchor_chrPos_end = tmpJunc_startPos;
		int acceptorAnchor_chrPos_start = tmpJunc_endPos;
		int acceptorAnchor_chrPos_end = tmpJunc_endPos + tmpJunc_acceptorAnchorSize - 1;
		vector<int> SNPposVec_donerAnchor;
		vector<int> SNPposVec_accetprAnchor;
		tmpSNPhashInfo.returnSNPposVecWithinRegion(tmpJunc_chrNameInt, donerAnchor_chrPos_start, 
			donerAnchor_chrPos_end, SNPposVec_donerAnchor);
		tmpSNPhashInfo.returnSNPposVecWithinRegion(tmpJunc_chrNameInt, acceptorAnchor_chrPos_start,
			acceptorAnchor_chrPos_end, SNPposVec_accetprAnchor);
		int SNPpos_donerAnchor_num = SNPposVec_donerAnchor.size();
		int SNPpos_acceptorAnchor_num = SNPposVec_accetprAnchor.size();
		int SNPpos_bothAnchor_num = SNPpos_donerAnchor_num + SNPpos_acceptorAnchor_num;
		if(SNPpos_bothAnchor_num > toCountSNPnumMaxInSJanchor)
			juncNumVec_withDifferentSNPnum[toCountSNPnumMaxInSJanchor] ++;
		else
			juncNumVec_withDifferentSNPnum[SNPpos_bothAnchor_num] ++;
	}
	anchorSNPstats_ofs << endl << "totalJuncNum: " << totalJuncNum << endl << endl
		<< "validJuncNum: " << validJuncNum << endl 
		<< "invalidJuncNum: " << invalidJuncNum << endl << endl;

	for(int tmp = 0; tmp < toCountSNPnumMaxInSJanchor + 1; tmp++)
	{
		double tmpSnpNumJuncPerc = ((double)juncNumVec_withDifferentSNPnum[tmp] / (double)validJuncNum) * 100; 
		anchorSNPstats_ofs << "SNP[" << tmp << "]:\t" << juncNumVec_withDifferentSNPnum[tmp] << "\t" 
			<< tmpSnpNumJuncPerc << "%" << endl; 
	}
	anchorSNPstats_ofs.close();
	free(chrom);
	delete indexInfo;
	//log_ofs.close();
	return 0;
}