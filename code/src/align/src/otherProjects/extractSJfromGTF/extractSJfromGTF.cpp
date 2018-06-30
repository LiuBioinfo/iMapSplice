// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: Gencode
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
#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputGTFfile outputFolder" << endl;
		exit(1);
	}
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());	
	string invalidChrGTF = outputFolderStr + "/invalidChr.gtf";
	string validChrNonExonGTF = outputFolderStr + "/validChr_nonExon.gtf";
	string validChrExonInvalidFormatGTF = outputFolderStr + "/validChr_exon_invalidFormat.gtf";
	string header = outputFolderStr + "/header";
	string validChrExonGTF = outputFolderStr + "/validChr_exon.gtf";
	ofstream invalidChrGTF_ofs(invalidChrGTF.c_str());
	ofstream validChrNonExonGTF_ofs(validChrNonExonGTF.c_str());
	ofstream validChrExonInvalidFormatGTF_ofs(validChrExonInvalidFormatGTF.c_str());
	ofstream header_ofs(header.c_str());
	ofstream validChrExonGTF_ofs(validChrExonGTF.c_str());

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	parameter_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;


	cout << "start to get exon info ..." << endl;	
	string inputGTFfile = argv[2];
	vector<int> chrNameIntVec_validChrExon;
	vector<int> startPosVec_validChrExon;
	vector<int> endPosVec_validChrExon;
	vector<string> strandVec_validChrExon;
	vector<string> transcriptIDvec_validChrExon;
	vector<string> geneNameVec_validChrExon;
	ifstream gtf_ifs(inputGTFfile.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == "#")
		{
			header_ofs << tmpStr << endl;
			continue;
		}
		vector<string> gtfFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpGtfField = tmpStr.substr(startLoc, tabLoc-startLoc);
			gtfFieldVec.push_back(tmpGtfField);
			startLoc = tabLoc + 1;
		}
		string tmpChrName = gtfFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			invalidChrGTF_ofs << tmpStr << endl;
			continue;
		}
		string tmpFeature = gtfFieldVec[2];
		int tmpStartPos = atoi((gtfFieldVec[3]).c_str());
		int tmpEndPos = atoi((gtfFieldVec[4]).c_str());
		string tmpStrand = gtfFieldVec[6];
		if(tmpFeature != "exon")
		{
			validChrNonExonGTF_ofs << tmpStr << endl;
			continue;
		}
		int tmpTranscriptIDlocInStr = tmpStr.find("transcript_id");
		int tmpGeneNameLocInStr = tmpStr.find("gene_name");
		if((tmpTranscriptIDlocInStr == string::npos)||(tmpGeneNameLocInStr == string::npos))
		{
			validChrExonInvalidFormatGTF_ofs << tmpStr << endl;
			continue;
		}
		int firstQuota_transcriptID = tmpStr.find("\"", tmpTranscriptIDlocInStr + 1);
		int secondQuota_transcriptID = tmpStr.find("\"", firstQuota_transcriptID + 1);
		int firstQuota_geneName = tmpStr.find("\"", tmpGeneNameLocInStr + 1);
		int secondQuota_geneName = tmpStr.find("\"", firstQuota_geneName + 1);
		if((firstQuota_transcriptID == string::npos)||(secondQuota_transcriptID == string::npos)
			||(firstQuota_geneName == string::npos)||(secondQuota_geneName == string::npos))
		{
			validChrExonInvalidFormatGTF_ofs << tmpStr << endl;
			continue;
		}
		validChrExonGTF_ofs << tmpStr << endl;
		string tmpTranscriptID = tmpStr.substr(firstQuota_transcriptID + 1,
			secondQuota_transcriptID - firstQuota_transcriptID - 1);
		string tmpGeneName = tmpStr.substr(firstQuota_geneName + 1,
			secondQuota_geneName - firstQuota_geneName - 1);
		chrNameIntVec_validChrExon.push_back(tmpChrNameInt);
		startPosVec_validChrExon.push_back(tmpStartPos);
		endPosVec_validChrExon.push_back(tmpEndPos);
		strandVec_validChrExon.push_back(tmpStrand);
		transcriptIDvec_validChrExon.push_back(tmpTranscriptID);
		geneNameVec_validChrExon.push_back(tmpGeneName);
	}
	gtf_ifs.close();
	invalidChrGTF_ofs.close();
	validChrNonExonGTF_ofs.close();
	validChrExonInvalidFormatGTF_ofs.close();
	header_ofs.close();
	validChrExonGTF_ofs.close();	


	cout << "start to generate splice junctions ......" << endl;
	log_ofs << "start to gene splice junctions ......" << endl;
	int totalExonNum = chrNameIntVec_validChrExon.size();
	cout << "totalExonNum: " << totalExonNum << endl;
	vector<int> SJchrNameIntVec;
	vector<int> SJstartPosVec;
	vector<int> SJendPosVec;
	vector<string> SJtranscriptIdVec;
	vector<string> SJgeneNameVec;
	vector<string> SJstrandVec;
	for(int tmp = 0; tmp < totalExonNum-1; tmp++)
	{
		int index_thisExon = tmp;
		int index_nextExon = tmp+1;
		int chrNameInt_thisExon = chrNameIntVec_validChrExon[index_thisExon];
		int chrNameInt_nextExon = chrNameIntVec_validChrExon[index_nextExon];
		int startPos_thisExon = startPosVec_validChrExon[index_thisExon];
		int startPos_nextExon = startPosVec_validChrExon[index_nextExon];
		int endPos_thisExon = endPosVec_validChrExon[index_thisExon];
		int endPos_nextExon = endPosVec_validChrExon[index_nextExon];		
		string strand_thisExon = strandVec_validChrExon[index_thisExon];
		string strand_nextExon = strandVec_validChrExon[index_nextExon];
		string transcriptID_thisExon = transcriptIDvec_validChrExon[index_thisExon];
		string transcriptID_nextExon = transcriptIDvec_validChrExon[index_nextExon];
		string geneName_thisExon = geneNameVec_validChrExon[index_thisExon];
		string geneName_nextExon = geneNameVec_validChrExon[index_nextExon];

		if(transcriptID_thisExon != transcriptID_nextExon)
			continue;
		if((chrNameInt_thisExon != chrNameInt_nextExon)
			||(strand_thisExon != strand_nextExon)
			||(geneName_thisExon != geneName_nextExon))
		{
			cout << "incorrect transcript structure of exons: " << endl << "transcriptID: " << transcriptID_thisExon << " & "
				<< transcriptID_nextExon << endl << "startPos_thisExon: " << startPos_thisExon << endl << "startPos_nextExon: " 
				<< startPos_nextExon << endl;
		}
		else
		{
			int tmpSJ_startPos, tmpSJ_endPos;
			if(endPos_thisExon < startPos_nextExon)
			{
				tmpSJ_startPos = endPos_thisExon;
				tmpSJ_endPos = startPos_nextExon;				
			}
			else if(endPos_nextExon < startPos_thisExon)
			{
				tmpSJ_startPos = endPos_nextExon;
				tmpSJ_endPos = startPos_thisExon;				
			}
			else
			{
				cout << "incorrect transcript structure of exons: " << endl << "transcriptID: " << transcriptID_thisExon << " & "
					<< transcriptID_nextExon << endl << "startPos_thisExon: "<< startPos_thisExon << endl << "startPos_nextExon: " 
					<< startPos_nextExon << endl;
			}

			SJchrNameIntVec.push_back(chrNameInt_thisExon);				
			SJstartPosVec.push_back(tmpSJ_startPos);
			SJendPosVec.push_back(tmpSJ_endPos);
			SJtranscriptIdVec.push_back(transcriptID_thisExon);
			SJstrandVec.push_back(strand_thisExon);
			SJgeneNameVec.push_back(geneName_thisExon);
		}
	}

	cout << "start to separate SJs according to chrNameInt " << endl;
	log_ofs << "start to separate SJs according to chrNameInt " << endl;
	vector< vector<int> > SJstartPosVecVec_separateByChrNameInt;  
	vector< vector<int> > SJendPosVecVec_separateByChrNameInt;
	vector< vector<string> > SJtranscriptIdVecVec_separateByChrNameInt;
	vector< vector<string> > SJgeneNameVecVec_separateByChrNameInt;
	vector< vector<string> > SJstrandVecVec_separateByChrNameInt;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{	
		vector<int> tmpStartPosVec;
		vector<int> tmpEndPosVec;
		vector<string> tmpTranscriptIdVec;
		vector<string> tmpGeneNameVec;
		vector<string> tmpStrandVec;
		SJstartPosVecVec_separateByChrNameInt.push_back(tmpStartPosVec);
		SJendPosVecVec_separateByChrNameInt.push_back(tmpEndPosVec);
		SJtranscriptIdVecVec_separateByChrNameInt.push_back(tmpTranscriptIdVec);
		SJgeneNameVecVec_separateByChrNameInt.push_back(tmpGeneNameVec);
		SJstrandVecVec_separateByChrNameInt.push_back(tmpStrandVec);
	}
	int totalSJnum = SJchrNameIntVec.size();
	cout << "total original SJ #: "<< totalSJnum << endl;
	log_ofs << "total original SJ #: "<< totalSJnum << endl;
	for(int tmp = 0; tmp < totalSJnum; tmp++)
	{
		int tmpSJchrNameInt = SJchrNameIntVec[tmp];
		SJstartPosVecVec_separateByChrNameInt[tmpSJchrNameInt].push_back(SJstartPosVec[tmp]);
		SJendPosVecVec_separateByChrNameInt[tmpSJchrNameInt].push_back(SJendPosVec[tmp]);
		SJtranscriptIdVecVec_separateByChrNameInt[tmpSJchrNameInt].push_back(SJtranscriptIdVec[tmp]);
		SJgeneNameVecVec_separateByChrNameInt[tmpSJchrNameInt].push_back(SJgeneNameVec[tmp]);
		SJstrandVecVec_separateByChrNameInt[tmpSJchrNameInt].push_back(SJstrandVec[tmp]);
	}

	cout << "start to generate non-duplicate splice junction final vec ......" << endl;
	log_ofs << "start to generate non-duplicate splice junction final vec ......" << endl;
	vector< vector<int> > SJstartPosVecVec_separateByChrNameInt_final; 
	vector< vector<int> > SJendPosVecVec_separateByChrNameInt_final;
	vector< vector<string> > SJtranscriptIdVecVec_separateByChrNameInt_final;
	vector< vector<string> > SJgeneNameVecVec_separateByChrNameInt_final;
	vector< vector<string> > SJstrandVecVec_separateByChrNameInt_final;	
	for(int tmp = 0; tmp < chromNum; tmp++)
	{	
		vector<int> tmpStartPosVec;
		vector<int> tmpEndPosVec;
		vector<string> tmpTranscriptIdVec;
		vector<string> tmpGeneNameVec;
		vector<string> tmpStrandVec;
		SJstartPosVecVec_separateByChrNameInt_final.push_back(tmpStartPosVec);
		SJendPosVecVec_separateByChrNameInt_final.push_back(tmpEndPosVec);
		SJtranscriptIdVecVec_separateByChrNameInt_final.push_back(tmpTranscriptIdVec);
		SJgeneNameVecVec_separateByChrNameInt_final.push_back(tmpGeneNameVec);
		SJstrandVecVec_separateByChrNameInt_final.push_back(tmpStrandVec);
	}
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		cout << endl << "start to remove duplicate SJs in " << tmpChrName << endl;
		log_ofs << endl << "start to remove duplicate SJs in " << tmpChrName << endl;
		int SJnumInChr = (SJstartPosVecVec_separateByChrNameInt[tmpChr]).size();
		cout << "SJ # " << SJnumInChr << endl;
		log_ofs << "SJ # " << SJnumInChr << endl;
		for(int tmp = 0; tmp < SJnumInChr; tmp++)
		{
			int tmpSJstartPos = (SJstartPosVecVec_separateByChrNameInt[tmpChr])[tmp];
			int tmpSJendPos = (SJendPosVecVec_separateByChrNameInt[tmpChr])[tmp];
			int currentFinalSJnum = SJstartPosVecVec_separateByChrNameInt_final[tmpChr].size();
			bool theSameWithOneSJinFinalVecBool = false;
			for(int tmp_inFinalVec = 0; tmp_inFinalVec < currentFinalSJnum; tmp_inFinalVec ++)
			{
				int tmpSJstartPos_inFinalVec = (SJstartPosVecVec_separateByChrNameInt_final[tmpChr])[tmp_inFinalVec];
				int tmpSJendPos_inFinalVec = (SJendPosVecVec_separateByChrNameInt_final[tmpChr])[tmp_inFinalVec];
				if((tmpSJstartPos == tmpSJstartPos_inFinalVec)&&(tmpSJendPos == tmpSJendPos_inFinalVec))
				{
					string tmpSJtranscriptIDstr_ori = (SJtranscriptIdVecVec_separateByChrNameInt_final[tmpChr])[tmp_inFinalVec];
					string tmpSJtranscriptIDstr_updated 
						= tmpSJtranscriptIDstr_ori + (SJtranscriptIdVecVec_separateByChrNameInt[tmpChr])[tmp] + ";";
					(SJtranscriptIdVecVec_separateByChrNameInt_final[tmpChr])[tmp_inFinalVec] = tmpSJtranscriptIDstr_updated;
					theSameWithOneSJinFinalVecBool = true;
					break;
				}
			}
			if(!theSameWithOneSJinFinalVecBool)
			{
				(SJstartPosVecVec_separateByChrNameInt_final[tmpChr]).push_back(tmpSJstartPos);
				(SJendPosVecVec_separateByChrNameInt_final[tmpChr]).push_back(tmpSJendPos);
				string tmpSJtranscriptIDtoAdd = (SJtranscriptIdVecVec_separateByChrNameInt[tmpChr])[tmp] + ";";
				(SJtranscriptIdVecVec_separateByChrNameInt_final[tmpChr]).push_back(tmpSJtranscriptIDtoAdd);
				string tmpGeneNameToAdd = (SJgeneNameVecVec_separateByChrNameInt[tmpChr])[tmp] + ";";
				(SJgeneNameVecVec_separateByChrNameInt_final[tmpChr]).push_back(tmpGeneNameToAdd);
				(SJstrandVecVec_separateByChrNameInt_final[tmpChr]).push_back((SJstrandVecVec_separateByChrNameInt[tmpChr])[tmp]);
			}					
		}
	}

	cout << "start to output final non-duplicate splice junction ......" << endl;
	log_ofs << "start to output final non-duplicate splice junction ......" << endl;
	string outputSJfile = outputFolderStr + "/SJ.txt";
	ofstream finalSJlist_ofs(outputSJfile.c_str());
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{	
		int finalSJnumInThisChr = SJstartPosVecVec_separateByChrNameInt_final[tmpChr].size();
		cout << endl << "chrName: " << indexInfo->returnChrNameStr(tmpChr) << endl;
		cout << "finalSJnumInThisChr: " << finalSJnumInThisChr << endl;
		log_ofs << endl << "chrName: " << indexInfo->returnChrNameStr(tmpChr) << endl;
		log_ofs << "finalSJnumInThisChr: " << finalSJnumInThisChr << endl;
		string tmpSJ_chrName = indexInfo->returnChrNameStr(tmpChr);
		for(int tmpSJ = 0; tmpSJ < finalSJnumInThisChr; tmpSJ++)
		{
			int tmpSJ_startPos = (SJstartPosVecVec_separateByChrNameInt_final[tmpChr])[tmpSJ];
			int tmpSJ_endPos = (SJendPosVecVec_separateByChrNameInt_final[tmpChr])[tmpSJ];
			string tmpSJ_transcriptID = (SJtranscriptIdVecVec_separateByChrNameInt_final[tmpChr])[tmpSJ];
			string tmpSJ_geneName = (SJgeneNameVecVec_separateByChrNameInt_final[tmpChr])[tmpSJ];
			string tmpSJ_strand = (SJstrandVecVec_separateByChrNameInt_final[tmpChr])[tmpSJ];
			finalSJlist_ofs << tmpSJ_chrName << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t"
				<< tmpSJ_strand << "\t" << tmpSJ_transcriptID << "\t" << tmpSJ_geneName << endl;				
		}
	}

	cout << "All jobs done !" << endl; 
	log_ofs << "All jobs done !" << endl;
	finalSJlist_ofs.close();
	log_ofs.close();
	return 0;
}