// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input -- original SAM format alignment results generated from MPS3_cirRNA
// output -- reformed alignments results: split back-splice alignments into several alignments,
// each alignment represents each segment

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
#include "../../../general/splice_info.h"

using namespace std;

string jumpCodeVec2Str(vector<Jump_Code>& cigarStringJumpCode)
{
	string str;
	for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
	{	
		str += cigarStringJumpCode[tmp].toString();
	}
	return str;
}

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

int getEndLocInReadOfSpecificJumpCode(
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	if(jumpCodeIndex < 0)
		return 0;
	int tmpEndLocInRead = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "M")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "D")
		{}
		else if(tmpJumpCodeType == "N")
		{}
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return tmpEndLocInRead;
}	

int getStartLocInReadOfSpecificJumpCode(
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int lastJumpCodeEndLocInRead;
	if(jumpCodeIndex == 0)
		lastJumpCodeEndLocInRead = 0;
	else
		lastJumpCodeEndLocInRead = getEndLocInReadOfSpecificJumpCode(
			cigarStringJumpCodeVec, jumpCodeIndex - 1);
	int tmpStartLocInRead = lastJumpCodeEndLocInRead + 1;
	return tmpStartLocInRead;
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

int getStartPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int lastJumpCodeEndPosInChr;
	if(jumpCodeIndex == 0)
		return startPos;
	else
	{
		lastJumpCodeEndPosInChr = getEndPosOfSpecificJumpCode(startPos,
			cigarStringJumpCodeVec, jumpCodeIndex - 1);
	}
	int tmpStartPosInChr = lastJumpCodeEndPosInChr + 1;
	return tmpStartPosInChr;
}



void generateNonBackSpliceSegmentVec(
	int mapPosInChr, string& cigarStringStr,
	vector< pair<int,int> >& nonBackSpliceSegmentLocInReadVec, 
	vector<int>& nonBackSpliceSegmentMapPosVec,
	vector< vector<Jump_Code> >& nonBackSpliceSegmentJumpCodeVecVec)
{
	vector<Jump_Code> totalJumpCodeVec;
	cigarString2jumpCodeVec(cigarStringStr, totalJumpCodeVec);
	int totalJumpCodeVecSize = totalJumpCodeVec.size();

	vector<int> backSpliceJumpCodeIndexVec;
	for(int tmp = 0; tmp < totalJumpCodeVecSize; tmp++)
	{
		int tmpJumpCodeLength = totalJumpCodeVec[tmp].len;
		string tmpJumpCodeType = totalJumpCodeVec[tmp].type;
		if((tmpJumpCodeType == "N")&&(tmpJumpCodeLength < 0))
			backSpliceJumpCodeIndexVec.push_back(tmp);
	}

	//cout << "backSpliceJumpCodeIndexVecSize: " << backSpliceJumpCodeIndexVec.size() << endl;

	vector< pair<int,int> > segmentJumpCodeIndexPairVec;
	int tmpSegmentJumpCodeIndex_start = 0;
	int tmpSegmentJumpCodeIndex_end;// = totalJumpCodeVecSize - 1;
	for(int tmpIndex_inBackSpliceJumpCodeIndexVec = 0;
		tmpIndex_inBackSpliceJumpCodeIndexVec < backSpliceJumpCodeIndexVec.size();
		tmpIndex_inBackSpliceJumpCodeIndexVec ++)
	{
		int	tmpBackSpliceJumpCodeIndex = backSpliceJumpCodeIndexVec[tmpIndex_inBackSpliceJumpCodeIndexVec];
		tmpSegmentJumpCodeIndex_end = tmpBackSpliceJumpCodeIndex - 1;
		segmentJumpCodeIndexPairVec.push_back(pair<int,int>(
			tmpSegmentJumpCodeIndex_start, tmpSegmentJumpCodeIndex_end));
		tmpSegmentJumpCodeIndex_start = tmpBackSpliceJumpCodeIndex + 1;
	}
	segmentJumpCodeIndexPairVec.push_back(pair<int,int>(tmpSegmentJumpCodeIndex_start, 
		totalJumpCodeVecSize-1));
	//cout << "segmentJumpCodeIndexPairVec.size(): " << segmentJumpCodeIndexPairVec.size() << endl;
	for(int tmpSegmentIndex = 0; tmpSegmentIndex < segmentJumpCodeIndexPairVec.size();
		tmpSegmentIndex ++)
	{
		int tmpSegmentJumpCodeIndex_start = segmentJumpCodeIndexPairVec[tmpSegmentIndex].first;
		int tmpSegmentJumpCodeIndex_end = segmentJumpCodeIndexPairVec[tmpSegmentIndex].second;
		int tmpSegmentStartLocInRead = getStartLocInReadOfSpecificJumpCode(
			totalJumpCodeVec, tmpSegmentJumpCodeIndex_start);
		int tmpSegmentEndLocInRead = getEndLocInReadOfSpecificJumpCode(
			totalJumpCodeVec, tmpSegmentJumpCodeIndex_end);
		int tmpSegmentStartPosInChr = getStartPosOfSpecificJumpCode(mapPosInChr,
			totalJumpCodeVec, tmpSegmentJumpCodeIndex_start);
		nonBackSpliceSegmentLocInReadVec.push_back(pair<int,int>(tmpSegmentStartLocInRead, tmpSegmentEndLocInRead));
		nonBackSpliceSegmentMapPosVec.push_back(tmpSegmentStartPosInChr);
		vector<Jump_Code> tmpSegmentJumpCodeVec;
		for(int tmpJumpCodeIndex = tmpSegmentJumpCodeIndex_start;
			tmpJumpCodeIndex <= tmpSegmentJumpCodeIndex_end; tmpJumpCodeIndex ++)
			tmpSegmentJumpCodeVec.push_back(totalJumpCodeVec[tmpJumpCodeIndex]);
		nonBackSpliceSegmentJumpCodeVecVec.push_back(tmpSegmentJumpCodeVec);
	}
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputOriginalSAMfile outputFolder" << endl;
		exit(1);
	}

	string inputOriginalSAMfile = argv[1];
	ifstream oriSAM_ifs(inputOriginalSAMfile.c_str());

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	

	string finalSAMfile = outputFolderStr + "final.sam";
	string reformedAllSAMfile = outputFolderStr + "reformed_all.sam";
	string reformedBackSpliceSAMfile = outputFolderStr + "reformed_backSplice.sam";
	string keptNonBackSpliceSAMfile = outputFolderStr + "kept_nonBackSplice.sam";
	string headerFile = outputFolderStr + "headerSection.txt";

	ofstream header_ofs(headerFile.c_str());
	ofstream reformedAllSAM_ofs(reformedAllSAMfile.c_str());
	ofstream reformedBackSpliceSAM_ofs(reformedBackSpliceSAMfile.c_str());
	ofstream keptNonBackSpliceSAM_ofs(keptNonBackSpliceSAMfile.c_str());

	while(!(oriSAM_ifs.eof()))
	{
		string samStr;
		getline(oriSAM_ifs, samStr);
		 if(oriSAM_ifs.eof()||(samStr == ""))
		 	break;
		if(samStr.at(0) == '@')		
		{
			header_ofs << samStr << endl;
			continue;
		}
		//cout << "samStr: " << endl << samStr << endl;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec.push_back(samStr.substr(startLoc));

		string cigarStringStr = samFieldVec[5];
		int minusSignFoundLoc = cigarStringStr.find("-");
		//cout << "minusSignFoundLoc: " << minusSignFoundLoc << endl;
		if(minusSignFoundLoc == string::npos) // nonBackSplice SAM
		{
			reformedAllSAM_ofs << samStr << endl;
			keptNonBackSpliceSAM_ofs << samStr << endl;
		}
		else // backSplice SAM
		{
			string readNameStr = samFieldVec[0];
			string flagStr = samFieldVec[1];
			string chrNameStr = samFieldVec[2];
			string chrPosStr = samFieldVec[3];
			string mapQualStr = samFieldVec[4];
			string rnextStr = samFieldVec[6];
			string pnextStr = samFieldVec[7];
			string tlenStr = samFieldVec[8];
			string readSeqStr = samFieldVec[9];
			string qualSeqStr = samFieldVec[10];
			string otherStr = samFieldVec[11];

			int chrPosInt = atoi(chrPosStr.c_str());
			vector< pair<int,int> > nonBackSpliceSegmentLocInReadVec;
			vector< int > nonBackSpliceSegmentMapPosVec;
			vector< vector<Jump_Code> > nonBackSpliceSegmentJumpCodeVecVec; 
			generateNonBackSpliceSegmentVec(chrPosInt, cigarStringStr,
				nonBackSpliceSegmentLocInReadVec, nonBackSpliceSegmentMapPosVec,
				nonBackSpliceSegmentJumpCodeVecVec);
			//cout << "segmentNum: " << nonBackSpliceSegmentLocInReadVec.size() << endl;
			for(int tmpSegment = 0; tmpSegment < nonBackSpliceSegmentLocInReadVec.size();
				tmpSegment ++)
			{
				string tmpSegmentReadName = readNameStr + "_" + int_to_str(tmpSegment+1);
				int tmpSegmentStartLocInRead = nonBackSpliceSegmentLocInReadVec[tmpSegment].first;
				int tmpSegmentEndLocInRead = nonBackSpliceSegmentLocInReadVec[tmpSegment].second;
				int tmpSegmentLengthInRead = tmpSegmentEndLocInRead - tmpSegmentStartLocInRead + 1;
				int tmpSegmentChrMapPos = nonBackSpliceSegmentMapPosVec[tmpSegment];
				string tmpSegmentCigarString = jumpCodeVec2Str(nonBackSpliceSegmentJumpCodeVecVec[tmpSegment]);
				string tmpSegmentReadSeq = readSeqStr.substr(tmpSegmentStartLocInRead - 1,
					tmpSegmentLengthInRead);
				string tmpSegmentQualSeq = "*";
				if(qualSeqStr.length() > 1)
					tmpSegmentQualSeq = qualSeqStr.substr(tmpSegmentStartLocInRead - 1,
						tmpSegmentLengthInRead);

				reformedAllSAM_ofs 
					<< tmpSegmentReadName << "\t" << flagStr << "\t"
					<< chrNameStr << "\t" << tmpSegmentChrMapPos << "\t"
					<< mapQualStr << "\t" << tmpSegmentCigarString << "\t"
					<< rnextStr << "\t" << pnextStr << "\t"
					<< tlenStr << "\t" << tmpSegmentReadSeq << "\t"
					<< tmpSegmentQualSeq << "\t" << otherStr << endl;
				reformedBackSpliceSAM_ofs
					<< tmpSegmentReadName << "\t" << flagStr << "\t"
					<< chrNameStr << "\t" << tmpSegmentChrMapPos << "\t"
					<< mapQualStr << "\t" << tmpSegmentCigarString << "\t"
					<< rnextStr << "\t" << pnextStr << "\t"
					<< tlenStr << "\t" << tmpSegmentReadSeq << "\t"
					<< tmpSegmentQualSeq << "\t" << otherStr << endl;
			}

		}
	}

	header_ofs.close();
	reformedAllSAM_ofs.close();
	reformedBackSpliceSAM_ofs.close();
	keptNonBackSpliceSAM_ofs.close();
	

	string cat2finalSAM_cmd = "cat " + headerFile + " " 
		+ reformedAllSAMfile + " > " + finalSAMfile;
	system(cat2finalSAM_cmd.c_str());

	return 0;
}