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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"
#include "../../../phase2/spliceJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <IndexInput> <inferJuncFile> <toCheckSJfile> <outputFolder> " << endl;
		// Note: actually inferFile should be exon files like gtf instead of junction file
		exit(1);
	}
	int simulatedAlignmentReadSeqLength = 50;
	int anchorSeqAroundEachSpliceSite = simulatedAlignmentReadSeqLength / 2;

	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());

	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());
	string output_candidateLariat_backSplice_SJ
		= outputFolderStr + "candidateLariat_backSplice.junc";
	ofstream candidateLariat_backSplice_SJ_ofs(
		output_candidateLariat_backSplice_SJ.c_str());

	string output_candidateLariat_backSplice_simAlignment
		= outputFolderStr + "candidateLariat_backSplice_simAlignment.sam";
	ofstream candidateLariat_backSplice_simAlignment_ofs(
		output_candidateLariat_backSplice_simAlignment.c_str());

	int chromNum = indexInfo->returnChromNum();
	string juncHash_file_1_str = argv[2];
	cout << "start to initiate alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_1_str, indexInfo);

	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
	int junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
	cout << "junctionNum_in_alignInferJuncHashInfo: " << junctionNum_in_alignInferJuncHashInfo << endl;

	string inputJuncFile = argv[3];
	ifstream junc_ifs(inputJuncFile.c_str());
	while(!(junc_ifs.eof()))
	{	
		string juncStr;
		getline(junc_ifs, juncStr);
		 if(junc_ifs.eof()||(juncStr == ""))
		 	break;

		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int lastTabLoc = juncStr.find("\t", startLoc);
		if(lastTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, lastTabLoc-startLoc));
		string tmpSJchrName = juncFieldVec[0];
		string tmpSJstartPosStr = juncFieldVec[1];
		string tmpSJendPosStr = juncFieldVec[2];
		string tmpSJnameStr = juncFieldVec[3];
		string tmpSJsupportNumStr = juncFieldVec[4];
		int tmpSJchrNameInt = indexInfo->convertStringToInt(tmpSJchrName);
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());
		int tmpSJsupportNumInt = atoi(tmpSJsupportNumStr.c_str());
		//cout << "tmpSJstartPosInt: " << tmpSJstartPosInt << endl;
		//cout << "tmpSJendPosInt: " << tmpSJendPosInt << endl;
		if(tmpSJstartPosInt > tmpSJendPosInt)
		{
			//cout << "tmpSJendPosInt: " << tmpSJendPosInt << endl;
			vector<int> tmpSJacceptorStartPosVec;
		 	bool tmpSJdonerPosFoundBool_exactPos = SJ->searchSJdonerEndPos(
		 		tmpSJchrNameInt, tmpSJendPosInt-1, tmpSJacceptorStartPosVec);
			bool tmpSJdonerPosFoundBool_exactPosWithin1base = SJ->searchSJdonerEndPos(
				tmpSJchrNameInt, tmpSJendPosInt-2, tmpSJacceptorStartPosVec);
			for(int tmp = 0; tmp < tmpSJacceptorStartPosVec.size(); tmp++)
			{
				//cout << "tmp: " << tmp << endl;
				int tmpSJacceptorStartPos = tmpSJacceptorStartPosVec[tmp];
				if(tmpSJstartPosInt < tmpSJacceptorStartPos)
				{
					candidateLariat_backSplice_SJ_ofs << juncStr << endl;
					int startPos_1stSegment = tmpSJendPosInt - 1 - anchorSeqAroundEachSpliceSite + 1;
					int endPos_1stSegment = tmpSJendPosInt - 1;
					int startPos_2ndSegment = tmpSJstartPosInt + 1;
					int endPos_2ndSegment = tmpSJstartPosInt + 1 + anchorSeqAroundEachSpliceSite - 1;
					string seq_1stSegment = indexInfo->returnChromStrSubstr(
						tmpSJchrNameInt, startPos_1stSegment, anchorSeqAroundEachSpliceSite);
					string seq_2ndSegment = indexInfo->returnChromStrSubstr(
						tmpSJchrNameInt, startPos_2ndSegment, anchorSeqAroundEachSpliceSite);
					string simulatedAlignmentReadSeq = seq_1stSegment + seq_2ndSegment;
					for(int tmp = 0; tmp < tmpSJsupportNumInt; tmp++)
					{
						string tmpReadSam = "seq." + tmpSJnameStr + "_" + int_to_str(tmp+1)
							+ "\t0\t" + tmpSJchrName + "\t" + int_to_str(startPos_1stSegment) + "\t255\t"
							+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" 
							+ int_to_str(startPos_2ndSegment - endPos_1stSegment - 1) + "N"
							+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" + "\t"
							+ "*\t0\t0\t" + simulatedAlignmentReadSeq + "\t*\tNM:i:0\tIH:i:1\tHI:i:1";
						candidateLariat_backSplice_simAlignment_ofs << tmpReadSam << endl;
					}					
					break;
				}
			}
		}
		else
		{}
	}	

	delete alignInferJunctionHashInfo;
	delete SJ;
	return 0;
}