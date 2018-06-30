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

int main(int argc, char** argv)
{
	if(argc <= 4)
	{
		cout << "Executable inputIndexFolder outputFolder inputJuncFileNum SJ_1 (SJ_2 ...)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "outputFolder: " << argv[2] << endl;
	log_ofs << "inputJuncFileNum: " << argv[3] << endl;

	string output_doner_alterAcceptor_supNumInSample = outputFolderStr + "doner_alterAcceptor.txt";
	string output_acceptor_alterDoner_supNumInSample = outputFolderStr + "acceptor_alterDoner.txt";
	ofstream doner_alterAcceptor_supNum_ofs(output_doner_alterAcceptor_supNumInSample.c_str());
	ofstream acceptor_alterDoner_supNum_ofs(output_acceptor_alterDoner_supNumInSample.c_str());

	string inputJuncFileNumStr = argv[3];
	int inputJuncFileNum = atoi(inputJuncFileNumStr.c_str());
	int SJ_file_num_in_command = argc - 4;
	if(inputJuncFileNum != SJ_file_num_in_command)
	{
		cout << "inputJuncFileNum != SJ_file_num_in_command" << endl;
		exit(1);
	}
	vector<string> juncFileVec;
	for(int tmp = 0; tmp < inputJuncFileNum; tmp++)
	{
		string tmpJuncFile = argv[4 + tmp];
		juncFileVec.push_back(tmpJuncFile);
	}

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
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

	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < inputJuncFileNum; tmp++)
	{
		string tmpJuncFile = juncFileVec[tmp];
		cout << "start to do initiate_juncHash: " << tmp + 1 << endl;
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum_includeBackSplice(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}	

	AlignInferJunctionHash_Info* juncHash_merged = new AlignInferJunctionHash_Info();
	juncHash_merged->initiateAlignInferJunctionInfo(chromNum);
	juncHash_merged->insertJuncFromJuncFileVec_chrNamePosOnly(juncFileVec, indexInfo);

	cout << "start to initiate SJhashInfo" << endl;
	log_ofs << "start to initiate SJhashInfo" << endl;	
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	cout << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	log_ofs << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	juncHash_merged->convert2SJhashInfo(SJ, indexInfo);

	cout << "start to initiate donerSiteSetVec and acceptorSiteSetVec ..." << endl;
	vector< set<int> > donerSiteSetVec;
	vector< set<int> > acceptorSiteSetVec;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerSiteSet;
		set<int> tmpAcceptorSiteSet;
		donerSiteSetVec.push_back(tmpDonerSiteSet);
		acceptorSiteSetVec.push_back(tmpAcceptorSiteSet);
	}

	cout << "start to insert splice sites 2 donerSiteSetVec and acceptorSiteSetVec ..." << endl;
	int alignInferInfoVecSize = juncHash_merged->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
	{	
		int tmpSJ_chrNameInt = juncHash_merged->returnAlignInferInfo_chrNameInt(tmp);
		//string tmpSJ_chrNameStr = indexInfo->returnChrNameStr(tmpSJ_chrNameInt);
		int tmpSJ_donerEndPos = juncHash_merged->returnAlignInferInfo_donerEndPos(tmp);
		int tmpSJ_acceptorStartPos = juncHash_merged->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpSJ_supportNum = juncHash_merged->returnAlignInferInfo_supportNum(tmp);
		int tmpSJ_juncSize = tmpSJ_acceptorStartPos - tmpSJ_donerEndPos - 1;
		if(tmpSJ_juncSize >= 0)
			continue;
		donerSiteSetVec[tmpSJ_chrNameInt].insert(tmpSJ_donerEndPos);
		acceptorSiteSetVec[tmpSJ_chrNameInt].insert(tmpSJ_acceptorStartPos);
	}	

	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		cout << "tmpChr: " << tmpChr << endl;
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(set<int>::iterator tmpIntSetIter = donerSiteSetVec[tmpChr].begin();
			tmpIntSetIter != donerSiteSetVec[tmpChr].end(); tmpIntSetIter ++)
		{
			int tmpDonerSite = *tmpIntSetIter;
			vector<int> tmpAlterAcceptorSiteVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChr, tmpDonerSite, tmpAlterAcceptorSiteVec);
			for(int tmpAlterAcceptorSiteIndex = 0; tmpAlterAcceptorSiteIndex < tmpAlterAcceptorSiteVec.size();
				tmpAlterAcceptorSiteIndex ++)
			{
				int tmpAlterAcceptorSite = tmpAlterAcceptorSiteVec[tmpAlterAcceptorSiteIndex];
				if(tmpDonerSite <= tmpAlterAcceptorSite) // normal splice
					continue;
				vector<int> thisBackSpliceSupNumInSampleVec;
				vector<int> otherAlterSJsupNumInSampleVec;
				for(int tmpSample = 0; tmpSample < inputJuncFileNum; tmpSample ++)
				{
					int thisBackSplice_supNum = juncHashVec[tmpSample]->searchAndReturnAlignInferJuncHashSupNum(
						tmpChr, tmpDonerSite, tmpAlterAcceptorSite);
					int thisDonerSite_supNumTotal = 0;
					for(int tmpAlterAcceptorSiteIndex_2 = 0; tmpAlterAcceptorSiteIndex_2 < tmpAlterAcceptorSiteVec.size();
						tmpAlterAcceptorSiteIndex_2 ++)
					{
						int tmpAlterSJ_supNum = juncHashVec[tmpSample]->searchAndReturnAlignInferJuncHashSupNum(
							tmpChr, tmpDonerSite, tmpAlterAcceptorSiteVec[tmpAlterAcceptorSiteIndex_2]);
						thisDonerSite_supNumTotal += tmpAlterSJ_supNum;
					}
					int thisDonerSite_supNumOther = thisDonerSite_supNumTotal - thisBackSplice_supNum;
					thisBackSpliceSupNumInSampleVec.push_back(thisBackSplice_supNum);
					otherAlterSJsupNumInSampleVec.push_back(thisDonerSite_supNumOther);
				}
				doner_alterAcceptor_supNum_ofs << tmpChrName << "\t" << tmpDonerSite << "\t" << tmpAlterAcceptorSite << "\tDoner2alterAcceptor\t";
				for(int TMP = 0; TMP < tmpAlterAcceptorSiteVec.size(); TMP++)
					doner_alterAcceptor_supNum_ofs << tmpAlterAcceptorSiteVec[TMP] << ";";
				for(int tmpSample = 0; tmpSample < inputJuncFileNum; tmpSample ++)
					doner_alterAcceptor_supNum_ofs << "\t" << thisBackSpliceSupNumInSampleVec[tmpSample]
						<< "\t" << otherAlterSJsupNumInSampleVec[tmpSample];
				doner_alterAcceptor_supNum_ofs << endl;
			}	
		}
		for(set<int>::iterator tmpIntSetIter = acceptorSiteSetVec[tmpChr].begin();
			tmpIntSetIter != acceptorSiteSetVec[tmpChr].end(); tmpIntSetIter ++)
		{
			int tmpAcceptorSite = *tmpIntSetIter;
			vector<int> tmpAlterDonerSiteVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChr, tmpAcceptorSite, tmpAlterDonerSiteVec);
			for(int tmpAlterDonerSiteIndex = 0; tmpAlterDonerSiteIndex < tmpAlterDonerSiteVec.size();
				tmpAlterDonerSiteIndex ++)
			{
				int tmpAlterDonerSite = tmpAlterDonerSiteVec[tmpAlterDonerSiteIndex];
				if(tmpAlterDonerSite <= tmpAcceptorSite)
					continue;
				vector<int> thisBackSpliceSupNumInSampleVec;
				vector<int> otherAlterSJsupNumInSampleVec;
				for(int tmpSample = 0; tmpSample < inputJuncFileNum; tmpSample ++)
				{
					int thisBackSplice_supNum = juncHashVec[tmpSample]->searchAndReturnAlignInferJuncHashSupNum(
						tmpChr, tmpAlterDonerSite, tmpAcceptorSite);					
					int thisAcceptorSite_supNumTotal = 0;
					for(int tmpAlterDonerSiteIndex_2 = 0; tmpAlterDonerSiteIndex_2 < tmpAlterDonerSiteVec.size();
						tmpAlterDonerSiteIndex_2 ++)
					{
						int tmpAlterSJ_supNum = juncHashVec[tmpSample]->searchAndReturnAlignInferJuncHashSupNum(
							tmpChr, tmpAlterDonerSiteVec[tmpAlterDonerSiteIndex_2], tmpAcceptorSite);
						thisAcceptorSite_supNumTotal += tmpAlterSJ_supNum;
					}
					int thisAcceptorSite_supNumOther = thisAcceptorSite_supNumTotal - thisBackSplice_supNum;
					thisBackSpliceSupNumInSampleVec.push_back(thisBackSplice_supNum);
					otherAlterSJsupNumInSampleVec.push_back(thisAcceptorSite_supNumOther);					
				}
				acceptor_alterDoner_supNum_ofs << tmpChrName << "\t" << tmpAcceptorSite << "\t" << tmpAlterDonerSite << "\tAcceptor2alterDoner\t";
				for(int TMP = 0; TMP < tmpAlterDonerSiteVec.size(); TMP++)
					acceptor_alterDoner_supNum_ofs << tmpAlterDonerSiteVec[TMP] << ";";
				for(int tmpSample = 0; tmpSample < inputJuncFileNum; tmpSample ++)
					acceptor_alterDoner_supNum_ofs << "\t" << thisBackSpliceSupNumInSampleVec[tmpSample]
						<< "\t" << otherAlterSJsupNumInSampleVec[tmpSample];
				acceptor_alterDoner_supNum_ofs << endl;
			}	
		}		
	}		

	delete juncHash_merged;
	for(int tmp = 0; tmp < inputJuncFileNum; tmp++)
		delete juncHashVec[tmp];
	return 0;
}