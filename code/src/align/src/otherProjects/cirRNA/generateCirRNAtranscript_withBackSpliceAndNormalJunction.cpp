// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: 1 -- backSplice junctions;
// 		  2 -- normalSplice junctions
// output: cirRNA transcript and their read counts

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

#include "cirRNA_transcript_info.h"
#include "../../general/splice_info.h"

using namespace std;

// bool mappedOrNot(int tmpFlag)
// {
// 	if(tmpFlag & 0x4)
// 		return false;
// 	else
// 		return true;
// }

void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpJumpCodeLength;
	string tmpJumpCodeType;

	int jumpCodeStartPosInCigarStr = 0;
	int jumpCodeEndPosInCigarStr;
		
	string candidateJumpCodeType = "SMNIDX";
	while(1)
	{
		jumpCodeEndPosInCigarStr = 
			jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
		if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
			{break;}
		else
		{
			tmpJumpCodeLength = 
				atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
			tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
			cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
			jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
		}
	}
}

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int tmpEndPos = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{}
		else if(tmpJumpCodeType == "M")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
		{}
		else if(tmpJumpCodeType == "D")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "N")
			tmpEndPos += tmpJumpCodeLength;
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}	
	return (tmpEndPos + startPos-1);
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexPath inputBackSpliceJunction inputNormalSpliceJunction inputSAM outputCirRNAtranscript" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;
	int chromNum = indexInfo->returnChromNum();

	string inputBackSpliceJunctionPath = argv[2];
	ifstream backSpliceJunction_ifs(inputBackSpliceJunctionPath.c_str());
	string inputNormalSpliceJunctionPath = argv[3];
	ifstream normalSpliceJunction_ifs(inputNormalSpliceJunctionPath.c_str());
	string inputSAMpath = argv[4];
	ifstream sam_ifs(inputSAMpath.c_str());
	string outputCirRNAtranscriptPath = argv[5];
	string outputCirRNAtranscriptPath_transcriptInfo = outputCirRNAtranscriptPath + ".transcriptInfo";
	ofstream cirRNAtranscript_ofs(outputCirRNAtranscriptPath_transcriptInfo.c_str());
	string outputCirRNAtranscriptPath_readCount = outputCirRNAtranscriptPath + ".readCount";
	ofstream cirRNAtranscript_readCount_ofs(outputCirRNAtranscriptPath_readCount.c_str());

	vector< vector< CirRNA_Transcript_Info* > > cirRNAtranscriptInfoVecVec;
	vector< vector< pair<int,int> > > backSpliceSitePairVecVec;
	vector< vector< pair<int,int> > > normalSpliceSitePairVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		vector<CirRNA_Transcript_Info*> tmpCirRNAtranscriptInfoVec;
		cirRNAtranscriptInfoVecVec.push_back(tmpCirRNAtranscriptInfoVec);
		vector< pair<int,int> > tmpBackSpliceSitePairVec;
		backSpliceSitePairVecVec.push_back(tmpBackSpliceSitePairVec);
		vector< pair<int,int> > tmpNormalSpliceSitePairVec;
		normalSpliceSitePairVecVec.push_back(tmpNormalSpliceSitePairVec);
	}

	cout << "start to extract backSpliceJunc" << endl;
	while(!backSpliceJunction_ifs.eof())
	{
		string juncStr;
		getline(backSpliceJunction_ifs, juncStr);
		 if(backSpliceJunction_ifs.eof()||(juncStr == ""))
		 	break;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));
		string tmpSJchrName = juncFieldVec[0];
		string tmpSJstartPosStr = juncFieldVec[1];
		string tmpSJendPosStr = juncFieldVec[2];
		int tmpSJchrNameInt = indexInfo->convertStringToInt(tmpSJchrName);
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());
		if(tmpSJstartPosInt > tmpSJendPosInt)
			backSpliceSitePairVecVec[tmpSJchrNameInt].push_back(pair<int,int>(tmpSJstartPosInt, tmpSJendPosInt));
		else
		{
			cout << "error in backSpliceSJ pos, tmpSJstartPosInt <= tmpSJendPosInt ! ";
			cout << "  tmpSJstartPosInt: " << tmpSJstartPosInt << " tmpSJendPosInt: " << tmpSJendPosInt;
			exit(1);
		}		
	}
	backSpliceJunction_ifs.close();

	cout << "start to extract normalSpliceJunc" << endl;
	while(!normalSpliceJunction_ifs.eof())
	{
		string juncStr;
		getline(normalSpliceJunction_ifs, juncStr);
		 if(normalSpliceJunction_ifs.eof()||(juncStr == ""))
		 	break;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));
		string tmpSJchrName = juncFieldVec[0];
		string tmpSJstartPosStr = juncFieldVec[1];
		string tmpSJendPosStr = juncFieldVec[2];
		int tmpSJchrNameInt = indexInfo->convertStringToInt(tmpSJchrName);
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());
		if(tmpSJstartPosInt < tmpSJendPosInt)
			normalSpliceSitePairVecVec[tmpSJchrNameInt].push_back(pair<int,int>(tmpSJstartPosInt, tmpSJendPosInt));
		else
		{
			cout << "error in normalSpliceSJ pos, tmpSJstartPosInt >= tmpSJendPosInt !";
			cout << "  tmpSJstartPosInt: " << tmpSJstartPosInt << " tmpSJendPosInt: " << tmpSJendPosInt;
			exit(1);
		}		
	}
	normalSpliceJunction_ifs.close();

	cout << "start to generate cirRNA transcript" << endl;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		int tmpBackSpliceSitePairVecSize = backSpliceSitePairVecVec[tmpChr].size();
		int tmpNormalSpliceSitePairVecSize = normalSpliceSitePairVecVec[tmpChr].size();
		for(int tmpIndexInBackSpliceSitePairVec = 0; tmpIndexInBackSpliceSitePairVec < tmpBackSpliceSitePairVecSize; 
			tmpIndexInBackSpliceSitePairVec ++)
		{
			int tmpBackSpliceSite_left = (backSpliceSitePairVecVec[tmpChr])[tmpIndexInBackSpliceSitePairVec].second;
			int tmpBackSpliceSite_right = (backSpliceSitePairVecVec[tmpChr])[tmpIndexInBackSpliceSitePairVec].first;
			vector< pair<int,int> > tmpCandiNormalSpliceSitePairVec;
			for(int tmpIndexInNormalSpliceSitePairVec = 0; tmpIndexInNormalSpliceSitePairVec < tmpNormalSpliceSitePairVecSize;
				tmpIndexInNormalSpliceSitePairVec ++)
			{
				int tmpNormalSpliceSite_left = (normalSpliceSitePairVecVec[tmpChr])[tmpIndexInNormalSpliceSitePairVec].first;
				int tmpNormalSpliceSite_right = (normalSpliceSitePairVecVec[tmpChr])[tmpIndexInNormalSpliceSitePairVec].second;
				if((tmpNormalSpliceSite_left > tmpBackSpliceSite_left)&&(tmpNormalSpliceSite_right < tmpBackSpliceSite_right))
					tmpCandiNormalSpliceSitePairVec.push_back(pair<int,int>(tmpNormalSpliceSite_left, tmpNormalSpliceSite_right));
			}
			CirRNA_Transcript_Info* tmpCirRNAtranscriptInfo = new CirRNA_Transcript_Info();
			tmpCirRNAtranscriptInfo->initiateCirRNAtranscriptInfo(tmpChr, tmpBackSpliceSite_left, tmpBackSpliceSite_right);
			bool addInnerSJsitePair_bool = tmpCirRNAtranscriptInfo->generateInnerSJvec(tmpCandiNormalSpliceSitePairVec);
			if(addInnerSJsitePair_bool)
			{	
				tmpCirRNAtranscriptInfo->generateExonVec();
				cirRNAtranscriptInfoVecVec[tmpChr].push_back(tmpCirRNAtranscriptInfo);
			}
			else
				delete tmpCirRNAtranscriptInfo;
		}
	}

	cout << "start to parse and process SAM file" << endl;
	while(1)
	{
		if(sam_ifs.eof())
			break;
		string samStr;
		getline(sam_ifs, samStr);
		if(samStr == "")
			break;
		if(samStr.at(0) == '@')
			continue;
		//cout << "samStr: " << samStr << endl;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string flagStr = samFieldVec[1];
		int flagInt = atoi(flagStr.c_str());
		bool mappedOrNotBool = mappedOrNot(flagInt);
		if(!mappedOrNotBool)
			continue;
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];
		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		int mapChrPos_end = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, cigarStringJumpCodeVec.size()-1);
		int tmpChr = mapChrNameInt;
		// cout << "mapChrPos: " << mapChrPos << endl;
		// cout << "mapChrPos_end: " << mapChrPos_end << endl;
		// cout << "mapChrNameInt: " << mapChrNameInt << endl;
		for(int tmpCirRNAtranscriptIndex = 0; tmpCirRNAtranscriptIndex < cirRNAtranscriptInfoVecVec[tmpChr].size(); tmpCirRNAtranscriptIndex++)
		{
			int tmpTranscriptStartPos 
				= (cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex]->return_backSpliceAcceptorStart_normalTranscriptStart_pos();
			int tmpTranscriptEndPos 
				= (cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex]->return_backSpliceDonerEnd_normalTranscriptEnd_pos();
			if((mapChrPos >= tmpTranscriptStartPos)&&(mapChrPos <= tmpTranscriptEndPos)
				&&(mapChrPos_end >= tmpTranscriptStartPos)&&(mapChrPos_end <= tmpTranscriptEndPos))
			(cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex]->readCountIncrement();	
		}
	}

	cout << "start to output transcript info" << endl;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		for(int tmpCirRNAtranscriptIndex = 0; tmpCirRNAtranscriptIndex < cirRNAtranscriptInfoVecVec[tmpChr].size(); tmpCirRNAtranscriptIndex++)
		{
			string tmpCirRNAtranscriptInfoStr = (cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex]->returnCirRNAtranscriptInfoStr(indexInfo);
			cirRNAtranscript_ofs << tmpCirRNAtranscriptInfoStr << endl;
			string tmpCirRNAtranscriptReadCountStr = (cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex]->returnCirRNAtranscriptReadCountStr(indexInfo);
			cirRNAtranscript_readCount_ofs << tmpCirRNAtranscriptReadCountStr << endl;
			delete (cirRNAtranscriptInfoVecVec[tmpChr])[tmpCirRNAtranscriptIndex];
		}
	}

	delete indexInfo;
	cirRNAtranscript_ofs.close();
	return 0;
}