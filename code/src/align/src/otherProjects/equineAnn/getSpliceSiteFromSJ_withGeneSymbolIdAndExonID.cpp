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

using namespace std;

string returnExonIdInfoStr_posSearch(int pos, vector<int>& exonPosVec, vector<string>& exonIdVec)
{
	string tmpStr = "";
	for(int tmp = 0; tmp < exonPosVec.size(); tmp++)
	{
		int tmpPos = exonPosVec[tmp];
		if(pos == tmpPos)
		{
			tmpStr += exonIdVec[tmp];
			tmpStr += ";";
		}
	}
	return tmpStr;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputGff3 inputSJ outputSpliceSiteWithGeneSymbolIDandExonID" << endl;
		exit(1);
	}
	cout << "start to read gff3 file and generate gene2exonVecVec" << endl;
	vector<string> geneSymbolVec;
	vector< vector<int> > exonStartPosVecVec;
	vector< vector<int> > exonEndPosVecVec;
	vector< vector<string> > exonIdVecVec;
	string inputGff3 = argv[1];
	ifstream gff3_ifs(inputGff3.c_str());
	string lastGeneName = "";
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
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpGff3Field = tmpStr.substr(startLoc, tabLoc-startLoc);
			gff3fieldVec.push_back(tmpGff3Field);
			startLoc = tabLoc + 1;
		}
		//string tmpChrAccession = gff3fieldVec[0];
		string tmpFeature = gff3fieldVec[2];
		int tmpStartPos = atoi((gff3fieldVec[3]).c_str());
		int tmpEndPos = atoi((gff3fieldVec[4]).c_str());
		if(tmpFeature != "exon")
			continue;	
		if(tmpStr.find("transcript_id=") == string::npos)
			continue;
		int tmpIdLocInStr = tmpStr.find("ID=");
		string tmpId;
		if(tmpIdLocInStr == string::npos)
			cout << "ID= not found in this exon: " << endl << tmpStr << endl;
		else
		{
			int nextSemicolon = tmpStr.find(";", tmpIdLocInStr);
			if(nextSemicolon == string::npos)
				cout << "Semicolon not found after tmpIdLocInStr: " << endl << tmpStr << endl;
			else
				tmpId = tmpStr.substr(tmpIdLocInStr, nextSemicolon - tmpIdLocInStr);
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

		if(tmpGeneName == lastGeneName)
		{
			int currentGeneNameVecSize = geneSymbolVec.size();
			exonStartPosVecVec[currentGeneNameVecSize-1].push_back(tmpStartPos);
			exonEndPosVecVec[currentGeneNameVecSize-1].push_back(tmpEndPos);
			exonIdVecVec[currentGeneNameVecSize-1].push_back(tmpId);
		}
		else
		{
			lastGeneName = tmpGeneName;
			geneSymbolVec.push_back(tmpGeneName);
			vector<int> tmpStartPosVec;
			tmpStartPosVec.push_back(tmpStartPos);
			exonStartPosVecVec.push_back(tmpStartPosVec);
			vector<int> tmpEndPosVec;
			tmpEndPosVec.push_back(tmpEndPos);
			exonEndPosVecVec.push_back(tmpEndPosVec);
			vector<string> tmpIdVec;
			tmpIdVec.push_back(tmpId);
			exonIdVecVec.push_back(tmpIdVec);			
		}
	}
	gff3_ifs.close();

	cout << "start to parse SJ from file and generate donerSiteVec_ori and acceptorSiteVec_ori" << endl;
	vector<string> geneSymbolVec_spliceSite;
	vector<string> chrNameVec_spliceSite;
	vector<string> chrAccessionVec_spliceSite;
	vector< set<int> > siteSetVec_doner_ori;
	vector< set<int> > siteSetVec_acceptor_ori;
	vector<string> strandVec_spliceSite;
	vector<string> geneIdVec_spliceSite;

	string inputSJ = argv[2];
	ifstream SJ_ifs(inputSJ.c_str());
	string lastGeneName_spliceSite = "";
	while(!SJ_ifs.eof())
	{
		string tmpSJstr;
		getline(SJ_ifs, tmpSJstr);
		if(tmpSJstr == "")
			break;
		vector<string> SJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpSJstr.find("\t", startLoc);
			string tmpSJfield = tmpSJstr.substr(startLoc, tabLoc-startLoc);
			SJfieldVec.push_back(tmpSJfield);
			startLoc = tabLoc + 1;
		}
		SJfieldVec.push_back(tmpSJstr.substr(startLoc));
		string tmpChrName = SJfieldVec[0];
		string tmpChrAccession = SJfieldVec[1];
		string tmpDonerSiteStr = SJfieldVec[2];
		string tmpAcceptorSiteStr = SJfieldVec[3];
		string tmpStrand = SJfieldVec[4];
		string tmpTranscriptInfoStr = SJfieldVec[5];
		string tmpGeneSymbol = SJfieldVec[6];
		string tmpGeneID = SJfieldVec[7];
		int tmpDonerSite = atoi(tmpDonerSiteStr.c_str());
		int tmpAcceptorSite = atoi(tmpAcceptorSiteStr.c_str());
		if(tmpGeneSymbol == lastGeneName_spliceSite)
		{
			int currentVecSize = geneSymbolVec_spliceSite.size();
			siteSetVec_doner_ori[currentVecSize-1].insert(tmpDonerSite);
			siteSetVec_acceptor_ori[currentVecSize-1].insert(tmpAcceptorSite);
		}
		else
		{
			lastGeneName_spliceSite = tmpGeneSymbol;
			geneSymbolVec_spliceSite.push_back(tmpGeneSymbol);
			chrNameVec_spliceSite.push_back(tmpChrName);
			chrAccessionVec_spliceSite.push_back(tmpChrAccession);
			strandVec_spliceSite.push_back(tmpStrand);
			geneIdVec_spliceSite.push_back(tmpGeneID);
			set<int> tmpDonerSiteSet;
			set<int> tmpAcceptorSiteSet;
			tmpDonerSiteSet.insert(tmpDonerSite);
			tmpAcceptorSiteSet.insert(tmpAcceptorSite);
			siteSetVec_doner_ori.push_back(tmpDonerSiteSet);
			siteSetVec_acceptor_ori.push_back(tmpAcceptorSiteSet);
		}
	}
	SJ_ifs.close();

	cout << "start to assign exon ids to splice sites and output" << endl;
	int geneNumWithSpliceSite = geneSymbolVec_spliceSite.size();
	string outputSpliceSiteWithGeneSymbolIDandExonID = argv[3];
	string doner_output = outputSpliceSiteWithGeneSymbolIDandExonID + "_doner.txt";
	string acceptor_output = outputSpliceSiteWithGeneSymbolIDandExonID + "_acceptor.txt";
	ofstream doner_ofs(doner_output.c_str());
	ofstream acceptor_ofs(acceptor_output.c_str());
	for(int tmp = 0; tmp < geneNumWithSpliceSite; tmp++)
	{
		string tmpGeneSymbolInfoStr = geneSymbolVec_spliceSite[tmp];
		vector<int> tmpGeneSymbolIndexInGeneSymbolVecFromGFF3;
		int startLoc = 0;
		for(int tmpSemiColonIndex = 0; ; tmpSemiColonIndex++)
		{
			int nextSemiColon = tmpGeneSymbolInfoStr.find(";", startLoc);
			if(nextSemiColon == string::npos)
				break;
			else
			{
				string tmpGeneSymbol = tmpGeneSymbolInfoStr.substr(startLoc, nextSemiColon - startLoc);
				for(int tmp2 = 0; tmp2 < geneSymbolVec.size(); tmp2++)
				{
					string tmpGeneSymbolInVec = geneSymbolVec[tmp2];
					if(tmpGeneSymbol == tmpGeneSymbolInVec)
					{
						tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.push_back(tmp2);
						//break;
					}
				}
			}
			startLoc = nextSemiColon + 1;
		}
		if(tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.size() == 0)
		{
			cout << "no such gene symbol is found in the geneSymbolVec derived from gff3: " << tmpGeneSymbolInfoStr << endl;
			continue;
		}

		string tmpGeneIdStr = geneIdVec_spliceSite[tmp]; 
		string tmpChrName = chrNameVec_spliceSite[tmp];
		string tmpChrAccession = chrAccessionVec_spliceSite[tmp];
		string tmpStrand = strandVec_spliceSite[tmp];
		// check doner site
		for(set<int>::iterator intSetIter = siteSetVec_doner_ori[tmp].begin();
			intSetIter != siteSetVec_doner_ori[tmp].end(); intSetIter ++)
		{
			int tmpDonerSite = (*intSetIter);
			string tmpExonIdStr = "";
			if(tmpStrand == "+") // check exon end pos;
			{
				for(int tmp3 = 0; tmp3 < tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.size(); tmp3 ++)
				{
					int tmpIndexInGeneSymbolVec = tmpGeneSymbolIndexInGeneSymbolVecFromGFF3[tmp3];
					string tmpStr = returnExonIdInfoStr_posSearch(tmpDonerSite, 
						exonEndPosVecVec[tmpIndexInGeneSymbolVec], exonIdVecVec[tmpIndexInGeneSymbolVec]);
					tmpExonIdStr += tmpStr;
				}
			}
			else // check exon start pos
			{
				for(int tmp3 = 0; tmp3 < tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.size(); tmp3 ++)
				{
					int tmpIndexInGeneSymbolVec = tmpGeneSymbolIndexInGeneSymbolVecFromGFF3[tmp3];
					string tmpStr = returnExonIdInfoStr_posSearch(tmpDonerSite, 
						exonStartPosVecVec[tmpIndexInGeneSymbolVec], exonIdVecVec[tmpIndexInGeneSymbolVec]);
					tmpExonIdStr += tmpStr;
				}
			}
			if(tmpExonIdStr == "")
				cout << "tmpExonIdStr: " << tmpGeneSymbolInfoStr << endl;
			doner_ofs << tmpChrName << "\t" << tmpChrAccession << "\t" << tmpDonerSite << "\t" 
				<< tmpStrand << "\t" << tmpExonIdStr << "\t" << tmpGeneSymbolInfoStr << "\t" << tmpGeneIdStr << endl;
		}
		// check acceptor site
		for(set<int>::iterator intSetIter = siteSetVec_acceptor_ori[tmp].begin();
			intSetIter != siteSetVec_acceptor_ori[tmp].end(); intSetIter ++)
		{
			int tmpAcceptorSite = (*intSetIter);
			string tmpExonIdStr = "";
			if(tmpStrand == "+") // check exon start pos;
			{
				for(int tmp3 = 0; tmp3 < tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.size(); tmp3 ++)
				{
					int tmpIndexInGeneSymbolVec = tmpGeneSymbolIndexInGeneSymbolVecFromGFF3[tmp3];
					string tmpStr = returnExonIdInfoStr_posSearch(tmpAcceptorSite, 
						exonStartPosVecVec[tmpIndexInGeneSymbolVec], exonIdVecVec[tmpIndexInGeneSymbolVec]);
					tmpExonIdStr += tmpStr;
				}
			}
			else // check exon end pos
			{
				for(int tmp3 = 0; tmp3 < tmpGeneSymbolIndexInGeneSymbolVecFromGFF3.size(); tmp3 ++)
				{
					int tmpIndexInGeneSymbolVec = tmpGeneSymbolIndexInGeneSymbolVecFromGFF3[tmp3];
					string tmpStr = returnExonIdInfoStr_posSearch(tmpAcceptorSite, 
						exonEndPosVecVec[tmpIndexInGeneSymbolVec], exonIdVecVec[tmpIndexInGeneSymbolVec]);
					tmpExonIdStr += tmpStr;
				}
			}
			if(tmpExonIdStr == "")
				cout << "tmpExonIdStr: " << tmpGeneSymbolInfoStr << endl;
			acceptor_ofs << tmpChrName << "\t" << tmpChrAccession << "\t" << tmpAcceptorSite << "\t" 
				<< tmpStrand << "\t" << tmpExonIdStr << "\t" << tmpGeneSymbolInfoStr << "\t" << tmpGeneIdStr << endl; 
		}		
	}
	doner_ofs.close();
	acceptor_ofs.close();
	return 0;
}