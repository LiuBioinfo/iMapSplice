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
	//int simulatedAlignmentReadSeqLength = 50;
	//int anchorSeqAroundEachSpliceSite = simulatedAlignmentReadSeqLength / 2;

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

	vector<int> offsetVec;
	int offset_max = 5;
	offsetVec.push_back(0);
	for(int tmp = 1; tmp <= offset_max; tmp++)
	{
		offsetVec.push_back(tmp);
		offsetVec.push_back(0-tmp);
	}

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
			vector< pair<int,int> > tmpFoundSJposPairVec; 
			bool tmpSJdonerEndPosFound_bool = false;
			for(int tmp = 0; tmp < offsetVec.size(); tmp++)
			{
				vector<int> tmpSJacceptorStartPosVec;
				int tmpOffset = offsetVec[tmp];
				int tmpToSearchDonerEndPos = tmpSJendPosInt - 1 + tmpOffset;
				bool tmpToSearchDonerEndPos_foundBool = SJ->searchSJdonerEndPos(
					tmpSJchrNameInt, tmpToSearchDonerEndPos, tmpSJacceptorStartPosVec);
				if(tmpToSearchDonerEndPos)
				{	
					tmpSJdonerEndPosFound_bool = true;
					for(int tmp2 = 0; tmp2 < tmpSJacceptorStartPosVec.size(); tmp2++)
						tmpFoundSJposPairVec.push_back(pair<int,int>(tmpToSearchDonerEndPos, 
							tmpSJacceptorStartPosVec[tmp2]));
				}
			}
			if(tmpSJdonerEndPosFound_bool)
			{	
				vector< pair<int,int> > tmpFinalFoundSJposPairVec;
				for(int tmp = 0; tmp < tmpFoundSJposPairVec.size(); tmp++)
				{
					//cout << "tmp: " << tmp << endl;
					//int tmpSJdonerEndPos = tmpFoundSJposPairVec[tmp].first;
					int tmpSJacceptorStartPos = tmpFoundSJposPairVec[tmp].second;
					if(tmpSJstartPosInt < tmpSJacceptorStartPos + offset_max)
						tmpFinalFoundSJposPairVec.push_back(tmpFoundSJposPairVec[tmp]);
				}
				if(tmpFinalFoundSJposPairVec.size() > 0)
				{
					candidateLariat_backSplice_SJ_ofs << tmpSJchrName << "\t" << tmpSJstartPosInt << "\t" 
						<< tmpSJendPosInt << "\t" << tmpSJnameStr << "\t" << tmpSJsupportNumInt << "\t";
					for(int tmpFinalFoundSJindex = 0; tmpFinalFoundSJindex < tmpFinalFoundSJposPairVec.size(); tmpFinalFoundSJindex++)
						candidateLariat_backSplice_SJ_ofs << tmpFinalFoundSJposPairVec[tmpFinalFoundSJindex].first
							<< ":" << tmpFinalFoundSJposPairVec[tmpFinalFoundSJindex].second 
							<< ":" << tmpSJendPosInt - tmpFinalFoundSJposPairVec[tmpFinalFoundSJindex].first - 1 << ","; 
					candidateLariat_backSplice_SJ_ofs << endl; 
				}
			}
		}
		else
		{}
	}	
	junc_ifs.close();
	candidateLariat_backSplice_SJ_ofs.close();
	delete alignInferJunctionHashInfo;
	return 0;
}