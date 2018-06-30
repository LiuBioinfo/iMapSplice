// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// inputAnn: /scratch/xli262/yizhang/Emma/new_annotation/102/ref_EquCab2.0_top_level.gff3
// inputChrNameAccession: /scratch/xli262/yizhang/Emma/new_annotation/102/chr_accessions.txt 
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

using namespace std;

void returnGeneSymbolAndIDwithTranscriptID(string& tmpTranscriptID, 
	string& tmpGeneSymbol, string& tmpGeneId,
	vector<string>& transcriptIdVec, vector<string>& geneIdVec, vector<string>& geneNameVec)
{
	int transcriptIdVecSize = transcriptIdVec.size();
	for(int tmp = 0; tmp < transcriptIdVecSize; tmp++)
	{
		if(tmpTranscriptID == transcriptIdVec[tmp])
		{
			tmpGeneSymbol = geneNameVec[tmp];
			tmpGeneId = geneIdVec[tmp];
		}
	}
}

void getGeneSymbolID(string& tmpGeneSymbolInfoStr, string& tmpGeneIdInfoStr, 
	string& tmpTranscriptInfoStr, vector<string>& transcriptIdVec, 
	vector<string>& geneIdVec, vector<string>& geneNameVec)
{
	vector<string> tmpTranscriptIdVec;
	//vector<string> tmpParseStr;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int nextSemiColon = tmpTranscriptInfoStr.find(";", startLoc);
		if(nextSemiColon == string::npos)
			break;
		else
			tmpTranscriptIdVec.push_back(tmpTranscriptInfoStr.substr(startLoc, nextSemiColon - startLoc));
		startLoc = nextSemiColon + 1;
	}

	int tmpTranscriptIdVecSize = tmpTranscriptIdVec.size();
	if(tmpTranscriptIdVecSize == 0)
		cout << "tmpTranscriptIdVecSize == 0" << endl;
	set<string> tmpGeneIdSet;
	set<string> tmpGeneNameSet;
	for(int tmp = 0; tmp < tmpTranscriptIdVecSize; tmp++)
	{
		string tmpGeneSymbol_ori, tmpGeneId_ori;
		returnGeneSymbolAndIDwithTranscriptID(tmpTranscriptIdVec[tmp],
			tmpGeneSymbol_ori, tmpGeneId_ori, transcriptIdVec, geneIdVec, geneNameVec);
		if(tmpGeneId_ori != "")
			tmpGeneIdSet.insert(tmpGeneId_ori);
		if(tmpGeneSymbol_ori != "")
			tmpGeneNameSet.insert(tmpGeneSymbol_ori);
	}

	for(set<string>::iterator setStrIter = tmpGeneIdSet.begin(); setStrIter != tmpGeneIdSet.end(); setStrIter ++)
	{
		tmpGeneIdInfoStr += (*setStrIter);
		tmpGeneIdInfoStr += ";";
	}
	for(set<string>::iterator setStrIter = tmpGeneNameSet.begin(); setStrIter != tmpGeneNameSet.end(); setStrIter ++)
	{		
		tmpGeneSymbolInfoStr += (*setStrIter);
		tmpGeneSymbolInfoStr += ";";
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputGff3 inputSJ outputSJwithGeneSymbolAccession" << endl;
		exit(1);
	}
	cout << "start to read gff3 file and generate transcriptIdVec, geneIdVec, geneNameVec" << endl;
	vector<string> transcriptIdVec;
	vector<string> geneIdVec;
	vector<string> geneNameVec;
	string inputGff3 = argv[1];
	ifstream gff3_ifs(inputGff3.c_str());
	while(!gff3_ifs.eof())
	{
		string tmpStr;
		getline(gff3_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == "#")
			continue;
		vector<string> gff3fieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpGff3Field = tmpStr.substr(startLoc, tabLoc-startLoc);
			gff3fieldVec.push_back(tmpGff3Field);
			startLoc = tabLoc + 1;
		}
		string tmpChrAccession = gff3fieldVec[0];
		string tmpFeature = gff3fieldVec[2];
		int tmpStartPos = atoi((gff3fieldVec[3]).c_str());
		int tmpEndPos = atoi((gff3fieldVec[4]).c_str());
		string tmpStrand = gff3fieldVec[6];
		//if(tmpFeature != "mRNA")
		//	continue;
		if(tmpStr.find("transcript_id=") == string::npos)
			continue;
		if((tmpFeature == "gene")||(tmpFeature == "exon")||(tmpFeature == "CDS"))
			continue;
		int tmpTranscriptIDlocInStr = tmpStr.find("transcript_id=");
		string tmpTranscriptID;
		if(tmpTranscriptIDlocInStr == string::npos)
			cout << "transcript ID not found in this exon: " << endl << tmpStr << endl;
		else
		{
			int nextSemicolon = tmpStr.find(";", tmpTranscriptIDlocInStr);
			if(nextSemicolon != string::npos)
				cout << "some other info after transcript id: " << endl << tmpStr << endl;
			else
				tmpTranscriptID = tmpStr.substr(tmpTranscriptIDlocInStr + 14);
		}

		int tmpGeneIDlocInStr = tmpStr.find("GeneID:");
		string tmpGeneID;
		if(tmpGeneIDlocInStr == string::npos)
			cout << "gene ID not found in this mRNA: " << endl << tmpStr << endl;
		else
		{
			int nextComma = tmpStr.find(",", tmpGeneIDlocInStr);
			if(nextComma == string::npos)
				cout << "Comma not found after tmpGeneIDlocInStr: " << endl  << tmpStr << endl;		
			else
				tmpGeneID = tmpStr.substr(tmpGeneIDlocInStr, nextComma - tmpGeneIDlocInStr);
		}

		int tmpGeneNameLocInStr = tmpStr.find("gene=");
		string tmpGeneName;
		if(tmpGeneNameLocInStr == string::npos)
			cout << "gene name not found in this mRNA: " << endl << tmpStr << endl;
		else
		{
			int nextSemicolon = tmpStr.find(";", tmpGeneNameLocInStr);
			if(nextSemicolon == string::npos)
				cout << "Semicolon not found after tmpGeneNameLocInStr: " << endl << tmpStr << endl;
			else
				tmpGeneName = tmpStr.substr(tmpGeneNameLocInStr + 5, nextSemicolon - tmpGeneNameLocInStr - 5);
		}

		transcriptIdVec.push_back(tmpTranscriptID);
		geneIdVec.push_back(tmpGeneID);
		geneNameVec.push_back(tmpGeneName);		
	}
	gff3_ifs.close();

	cout << "start to add gene symbol and accession to SJ" << endl;
	string outputSJwithGeneSymbolAccession = argv[3];
	ofstream SJwithGene_ofs(outputSJwithGeneSymbolAccession.c_str());
	string inputSJ = argv[2];
	ifstream SJ_ifs(inputSJ.c_str());
	while(!SJ_ifs.eof())
	{
		string tmpSJstr;
		getline(SJ_ifs, tmpSJstr);
		if(tmpSJstr == "")
			break;
		vector<string> SJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = tmpSJstr.find("\t", startLoc);
			string tmpSJfield = tmpSJstr.substr(startLoc, tabLoc-startLoc);
			SJfieldVec.push_back(tmpSJfield);
			startLoc = tabLoc + 1;
		}
		SJfieldVec.push_back(tmpSJstr.substr(startLoc));
		string tmpChrName = SJfieldVec[0];
		string tmpChrAccession = SJfieldVec[1];
		string tmpStartPosStr = SJfieldVec[2];
		string tmpEndPosStr = SJfieldVec[3];
		string tmpStrand = SJfieldVec[4];
		string tmpTranscriptInfoStr = SJfieldVec[5];
		if(tmpStrand == "+")
			SJwithGene_ofs << tmpChrName << "\t" << tmpChrAccession << "\t"  
				<< tmpStartPosStr << "\t" << tmpEndPosStr << "\t" 
				<< tmpStrand << "\t" << tmpTranscriptInfoStr << "\t";
		else
			SJwithGene_ofs << tmpChrName << "\t" << tmpChrAccession << "\t"  
				<< tmpEndPosStr << "\t" << tmpStartPosStr << "\t"
				<< tmpStrand << "\t" << tmpTranscriptInfoStr << "\t";
		string tmpGeneSymbol;
		string tmpGeneID;
		getGeneSymbolID(tmpGeneSymbol, tmpGeneID, tmpTranscriptInfoStr, 
			transcriptIdVec, geneIdVec, geneNameVec);
		SJwithGene_ofs << tmpGeneSymbol << "\t" << tmpGeneID << endl;
		if((tmpGeneSymbol == "")||(tmpGeneID == ""))
			cout << "tmpGeneSymbol or tmpGeneID is empty:" << endl << tmpSJstr << endl;
	}
	SJ_ifs.close();
	SJwithGene_ofs.close();
	return 0;
}