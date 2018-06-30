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

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputChrNameAccessionFile inputGff3file outputSJ" << endl;
		exit(1);
	}
	string inputChrNameAccessionFile = argv[1];
	string inputGff3file = argv[2];
	string outputSJfile = argv[3];

	cout << "start to get chrNameAccessionPairVec" << endl;
	vector< pair<string, string> > chrNameAccessionPairVec;
	ifstream chrNameAccession_ifs(inputChrNameAccessionFile.c_str());
	while(!chrNameAccession_ifs.eof())
	{
		string tmpStr;
		getline(chrNameAccession_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpChrName = tmpStr.substr(0, tabLoc);
		string tmpAccession = tmpStr.substr(tabLoc + 1);
		string tmpChrName_final = "chr" + tmpChrName;
		chrNameAccessionPairVec.push_back(pair<string,string>(tmpAccession, tmpChrName_final));
	}
	chrNameAccession_ifs.close();
	cout << "start to get exon info ..." << endl;	
	vector<string> chrAccessionVec;
	vector<int> startPosVec;
	vector<int> endPosVec;
	vector<string> strandVec;
	vector<string> transcriptIDvec;
	ifstream gff3_ifs(inputGff3file.c_str());
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
		if(tmpFeature != "exon")
			continue;
		int tmpTranscriptIDlocInStr = tmpStr.find("transcript_id=");
		string tmpTranscriptID;
		if(tmpTranscriptIDlocInStr == string::npos)
		{
			cout << "transcript ID not found in this exon: " << endl
				<< tmpStr << endl;
			//continue;
		}
		else
		{
			int nextSemicolon = tmpStr.find(";", tmpTranscriptIDlocInStr);
			if(nextSemicolon != string::npos)
			{
				cout << "some other info after transcript id: " << endl 
					<< tmpStr << endl;
			}
			else
			{
				tmpTranscriptID = tmpStr.substr(tmpTranscriptIDlocInStr + 14);
				chrAccessionVec.push_back(tmpChrAccession);
				startPosVec.push_back(tmpStartPos);
				endPosVec.push_back(tmpEndPos);
				strandVec.push_back(tmpStrand);
				transcriptIDvec.push_back(tmpTranscriptID);		
			}
		}
	}
	gff3_ifs.close();
	cout << "start to generate splice junctions ......" << endl;
	int totalExonNum = chrAccessionVec.size();
	cout << "totalExonNum: " << totalExonNum << endl;
	vector<string> SJchrAccessionVec;
	vector<int> SJstartPosVec;
	vector<int> SJendPosVec;
	vector<string> SJtranscriptIdVec;
	vector<string> SJstrandVec;
	for(int tmp = 0; tmp < totalExonNum-1; tmp++)
	{
		int index_thisExon = tmp;
		int index_nextExon = tmp+1;
		string chrAccession_thisExon = chrAccessionVec[index_thisExon];
		string chrAccession_nextExon = chrAccessionVec[index_nextExon];
		int startPos_thisExon = startPosVec[index_thisExon];
		int startPos_nextExon = startPosVec[index_nextExon];
		int endPos_thisExon = endPosVec[index_thisExon];
		int endPos_nextExon = endPosVec[index_nextExon];		
		string strand_thisExon = strandVec[index_thisExon];
		string strand_nextExon = strandVec[index_nextExon];
		string transcriptID_thisExon = transcriptIDvec[index_thisExon];
		string transcriptID_nextExon = transcriptIDvec[index_nextExon];
		if(transcriptID_thisExon != transcriptID_nextExon)
			continue;
		if((chrAccession_thisExon != chrAccession_nextExon)
			||(strand_thisExon != strand_nextExon))
		{
			cout << "incorrect transcript structure of exons: " << endl 
				<< "transcriptID: " << transcriptID_thisExon << " & "
				<< transcriptID_nextExon << endl << "startPos_thisExon: "
				<< startPos_thisExon << endl << "startPos_nextExon: " 
				<< startPos_nextExon << endl;
		}
		else
		{
			int tmpSJ_startPos, tmpSJ_endPos;
			if(strand_thisExon == "+")
			{
				tmpSJ_startPos = endPos_thisExon;
				tmpSJ_endPos = startPos_nextExon;
			}
			if(strand_thisExon == "-")
			{
				tmpSJ_startPos = endPos_nextExon;
				tmpSJ_endPos = startPos_thisExon;
			}
			if(tmpSJ_startPos >= tmpSJ_endPos)
				cout << "incorrect SJ sites: " << "tmpSJ_startPos: " << tmpSJ_startPos
					<< " tmpSJ_endPos: " << tmpSJ_endPos << endl;
			else
			{
				SJchrAccessionVec.push_back(chrAccession_thisExon);				
				SJstartPosVec.push_back(tmpSJ_startPos);
				SJendPosVec.push_back(tmpSJ_endPos);
				SJtranscriptIdVec.push_back(transcriptID_thisExon);
				SJstrandVec.push_back(strand_thisExon);
			}
		}
	}
	cout << "start to generate splice junction final vec ......" << endl;
	int totalSJnum = SJstartPosVec.size();
	cout << "totalSJnum: " << totalSJnum << endl;
	vector<string> SJchrAccessionVec_final;
	vector<int> SJstartPosVec_final;
	vector<int> SJendPosVec_final;
	vector<string> SJtranscriptIdVec_final;
	vector<string> SJstrandVec_final;
	for(int tmp = 0; tmp < totalSJnum; tmp++)
	{
		cout << "tmpSJ: " << tmp << endl;
		string tmpSJchrAccession = SJchrAccessionVec[tmp];
		int tmpSJstartPos = SJstartPosVec[tmp];
		int tmpSJendPos = SJendPosVec[tmp];
		string tmpSJstrand = SJstrandVec[tmp];
		int currentFinalSJnum = SJchrAccessionVec_final.size();
		bool theSameWithOneSJinFinalVecBool = false;
		for(int tmp_inFinalVec = 0; tmp_inFinalVec < currentFinalSJnum; tmp_inFinalVec ++)
		{
			string tmpSJchrAccession_inFinalVec = SJchrAccessionVec_final[tmp_inFinalVec];
			int tmpSJstartPos_inFinalVec = SJstartPosVec_final[tmp_inFinalVec];
			int tmpSJendPos_inFinalVec = SJendPosVec_final[tmp_inFinalVec];
			string tmpSJstrand_inFinalVec = SJstrandVec_final[tmp_inFinalVec];
			if((tmpSJchrAccession == tmpSJchrAccession_inFinalVec)
				&&(tmpSJstartPos == tmpSJstartPos_inFinalVec)
				&&(tmpSJendPos == tmpSJendPos_inFinalVec)
				&&(tmpSJstrand == tmpSJstrand_inFinalVec))
			{
				string tmpSJtranscriptIDstr_ori = SJtranscriptIdVec_final[tmp_inFinalVec];
				string tmpSJtranscriptIDstr_updated 
					= tmpSJtranscriptIDstr_ori + SJtranscriptIdVec[tmp] + ";";
				SJtranscriptIdVec_final[tmp_inFinalVec] = tmpSJtranscriptIDstr_updated;
				theSameWithOneSJinFinalVecBool = true;
				break;
			}
		}
		if(!theSameWithOneSJinFinalVecBool)
		{
			SJchrAccessionVec_final.push_back(tmpSJchrAccession);
			SJstartPosVec_final.push_back(tmpSJstartPos);
			SJendPosVec_final.push_back(tmpSJendPos);
			string tmpSJtranscriptIDtoAdd = SJtranscriptIdVec[tmp] + ";";
			SJtranscriptIdVec_final.push_back(tmpSJtranscriptIDtoAdd);
			SJstrandVec_final.push_back(tmpSJstrand);
		}
	}
	cout << "start to output final splice junction ......" << endl;
	int totalFinalSJnum = SJchrAccessionVec_final.size();
	cout << "totalSJnum_final: " << totalFinalSJnum << endl;
	ofstream finalSJlist_ofs(outputSJfile.c_str());
	for(int tmp = 0; tmp < totalFinalSJnum; tmp++)
	{
		string tmpSJchrAccession_inFinalVec = SJchrAccessionVec_final[tmp];
		int tmpSJstartPos_inFinalVec = SJstartPosVec_final[tmp];
		int tmpSJendPos_inFinalVec = SJendPosVec_final[tmp];
		string tmpSJtranscriptIDstr = SJtranscriptIdVec_final[tmp];
		string tmpSJstrand = SJstrandVec_final[tmp];
		// search for accession ....
		string tmpSJchrName = "NotFound";
		for(int tmpAccession = 0; tmpAccession < chrNameAccessionPairVec.size(); tmpAccession ++)
		{
			string tmpChrAccessionStr = chrNameAccessionPairVec[tmpAccession].first;
			string tmpChrNameStr = chrNameAccessionPairVec[tmpAccession].second;
			if(tmpChrAccessionStr == tmpSJchrAccession_inFinalVec)
			{
				tmpSJchrName = tmpChrNameStr;
				break;
			}
		}
		finalSJlist_ofs << tmpSJchrName << "\t" << tmpSJchrAccession_inFinalVec << "\t"
			<< tmpSJstartPos_inFinalVec << "\t" << tmpSJendPos_inFinalVec << "\t"
			<< tmpSJstrand << "\t" << tmpSJtranscriptIDstr << endl;
	}
	finalSJlist_ofs.close();
	return 0;
}