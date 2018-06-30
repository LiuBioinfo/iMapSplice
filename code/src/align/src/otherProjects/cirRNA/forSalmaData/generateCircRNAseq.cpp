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

using namespace std;

int returnLeftMostSJdonerSite(vector< int >& donerVec)
{
	int leftMost = donerVec[0];
	if(donerVec.size() > 1)
	{
		for(int tmp = 1; tmp < donerVec.size(); tmp++)
		{
			int tmpDonerSite = donerVec[tmp];
			if(tmpDonerSite < leftMost)
				leftMost = tmpDonerSite;
		}
	}
	return leftMost;
}

int returnRightMostSJacceptorSite(vector< int >& acceptorVec)
{
	int rightMost = acceptorVec[0];
	if(acceptorVec.size() > 1)
	{
		for(int tmp = 1; tmp < acceptorVec.size(); tmp++)
		{
			int tmpAcceptorSite = acceptorVec[tmp];
			if(tmpAcceptorSite > rightMost)
				rightMost = tmpAcceptorSite;
		}
	}
	return rightMost;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolder inputJuncWithTranscriptIdFile inputBackSpliceJuncFile outputFile" << endl;
		exit(1);
	}
	cout << "loading indexInfo parameters ......" << endl;
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
	cout << "start to initiate annotated juncVec" << endl;
	string inputJuncWithTranscriptIdFile = argv[2];
	vector< vector< pair<int,int> > > juncPosPairVecVec;
	vector< vector< vector<string> > > juncTranscriptIdVecVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		vector< pair<int,int> > tmpJuncPosPairVec;
		vector< vector<string> > tmpJuncTranscriptIdVecVec;
		juncPosPairVecVec.push_back(tmpJuncPosPairVec);
		juncTranscriptIdVecVecVec.push_back(tmpJuncTranscriptIdVecVec);
	}
	cout << "start to load annotated juncs" << endl;
	ifstream juncWithTranscriptId_ifs(inputJuncWithTranscriptIdFile.c_str());
	while(!juncWithTranscriptId_ifs.eof())
	{
		string tmpStr;
		getline(juncWithTranscriptId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		string tmpJunc_chrName = tmpStr.substr(0, tabLoc_1);
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		string tmpJunc_startPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		string tmpJunc_endPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		string tmpJunc_transcriptIdVecStr = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		vector<string> transcriptIdFieldVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpJunc_transcriptIdVecStr.find(";", startLoc);
			if(tabLoc != string::npos)
			{
				string tmpTranscriptIdField = tmpJunc_transcriptIdVecStr.substr(startLoc, tabLoc - startLoc);
				transcriptIdFieldVec.push_back(tmpTranscriptIdField);
			}
			else
				break;
			startLoc = tabLoc + 1;
		}
		juncPosPairVecVec[tmpJunc_chrNameInt].push_back(pair<int,int>(tmpJunc_startPos, tmpJunc_endPos));
		juncTranscriptIdVecVecVec[tmpJunc_chrNameInt].push_back(transcriptIdFieldVec);
	}
	cout << "start to search for juncs within back splice junctions" << endl;
	string outputJuncWithBackSplice = argv[4];
	ofstream juncWithinBackSplice_ofs(outputJuncWithBackSplice.c_str());
	string inputBackSpliceJuncFile = argv[3];
	ifstream backSplice_ifs(inputBackSpliceJuncFile.c_str());
	while(!backSplice_ifs.eof())
	{
		string tmpStr;
		getline(backSplice_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpBackSplice_chrName = tmpStr.substr(0, tabLoc_1);
		string tmpBackSplice_donerStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpBackSplice_acceptorStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpBackSplice_chrNameInt = indexInfo->convertStringToInt(tmpBackSplice_chrName);
		int tmpBackSplice_doner = atoi(tmpBackSplice_donerStr.c_str());
		int tmpBackSplice_acceptor = atoi(tmpBackSplice_acceptorStr.c_str());
		int tmpRegion_start = tmpBackSplice_acceptor;
		int tmpRegion_end = tmpBackSplice_doner;
		vector<int> tmpWithinJuncVec_doner;
		vector<int> tmpWithinJuncVec_acceptor;
		vector< vector<string> > tmpWithinJuncVec_transcriptIdVecVec;
		int tmpChrJuncNum = juncPosPairVecVec[tmpBackSplice_chrNameInt].size();
		for(int tmp = 0; tmp < tmpChrJuncNum; tmp++)
		{
			int tmpChrJunc_doner = (juncPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].first;
			int tmpChrJunc_acceptor = (juncPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].second;
			if((tmpChrJunc_doner >= tmpRegion_start)&&(tmpChrJunc_acceptor <= tmpRegion_end))
			{
				tmpWithinJuncVec_doner.push_back(tmpChrJunc_doner);
				tmpWithinJuncVec_acceptor.push_back(tmpChrJunc_acceptor);
				tmpWithinJuncVec_transcriptIdVecVec.push_back((juncTranscriptIdVecVecVec[tmpBackSplice_chrNameInt])[tmp]);
			}
		}
		juncWithinBackSplice_ofs << "--------------------------------------------------------" << endl;
		juncWithinBackSplice_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr << "--" << tmpBackSplice_acceptorStr << ":" << endl;
		// for(int tmp = 0; tmp < tmpWithinJuncVec_doner.size(); tmp++)
		// {
		// 	juncWithinBackSplice_ofs << tmpWithinJuncVec_doner[tmp] << "\t" << tmpWithinJuncVec_acceptor[tmp] << "\t";
		// 	for(int tmpTranscriptIndex = 0; tmpTranscriptIndex < (tmpWithinJuncVec_transcriptIdVecVec[tmp]).size(); tmpTranscriptIndex ++)
		// 		juncWithinBackSplice_ofs << (tmpWithinJuncVec_transcriptIdVecVec[tmp])[tmpTranscriptIndex] << ";";
		// 	juncWithinBackSplice_ofs << endl;
		// }
		int withinJuncNum = tmpWithinJuncVec_doner.size();
		if(withinJuncNum >= 1)
		{
			int tmpLeftMost = returnLeftMostSJdonerSite(tmpWithinJuncVec_doner);
			int tmpRightMost = returnRightMostSJacceptorSite(tmpWithinJuncVec_acceptor);
			int tmpLeftDistance = tmpLeftMost - tmpBackSplice_acceptor + 1;
			int tmpRightDistance = tmpBackSplice_doner - tmpRightMost + 1;
			juncWithinBackSplice_ofs << "leftExonicSeqLen:\t" << tmpLeftDistance << "\trightExonicSeqLen:\t" << tmpRightDistance << endl;
		}
	}
	backSplice_ifs.close();
	juncWithinBackSplice_ofs.close();
	juncWithTranscriptId_ifs.close();
	return 0;
}