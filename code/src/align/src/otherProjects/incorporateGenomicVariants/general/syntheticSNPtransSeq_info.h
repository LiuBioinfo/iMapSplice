// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SYNTHETICSNPTRANSSEQ_INFO_H
#define SYNTHETICSNPTRANSSEQ_INFO_H

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

class SyntheticSNPtransSeq_Info
{
private:
	string chrName;
	int chrNameInt;
	int startPosInChr;
	vector<Jump_Code> cigarStringJumpCodeVec;
	int syntheticReadLength;
	int SNPlocInSyntheticSeq;
	int SNPposInChr;
	string transcriptName;

	vector<int> exon_start_pos;
	vector<int> exon_end_pos;
	int exon_num;
	vector< pair<int,int> > spliceSiteInSNPseq; //50M1000N50M : splice site is 51	
public:
	SyntheticSNPtransSeq_Info()
	{
		exon_num = 0;
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
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "D")
			{
			}
			else if(tmpJumpCodeType == "N")
			{
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return tmpEndLocInRead;
	}	

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
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

	void generate_exonInfo_spliceSiteInfo()
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
			if(tmpJumpCodeType == "M")
			{
				exon_num ++;
				int tmpExonEndPos = this->getEndPosOfSpecificJumpCode(startPosInChr, cigarStringJumpCodeVec, tmp);
				int tmpExonStartPos = tmpExonEndPos - tmpJumpCodeLength + 1;
				exon_start_pos.push_back(tmpExonStartPos);
				exon_end_pos.push_back(tmpExonEndPos);
			}
			else if(tmpJumpCodeType == "N")
			{
				int tmpSpliceSiteLoc = this->getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmp - 1) + 1;
				spliceSiteInSNPseq.push_back(pair<int,int>(tmpSpliceSiteLoc, tmpJumpCodeLength));
			}
			else
			{}
		}		
	}

	string return_chrName()
	{
		return chrName;
	}

	void initiateWithSyntheticSNPtransSeqName(string& tmpSyntheticSNPtransSeqName, Index_Info* genomeIndexInfo)
	{
		vector<string> tmpFieldVec;
		int startLoc = 0;
		for(int  tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpSyntheticSNPtransSeqName.find(":", startLoc);
			string tmpFieldStr = tmpSyntheticSNPtransSeqName.substr(startLoc, tabLoc - startLoc);
			tmpFieldVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}

		chrName = tmpFieldVec[0];
		
		chrNameInt = genomeIndexInfo->convertStringToInt(chrName);
		
		string startPosInChrStr = tmpFieldVec[1];
		startPosInChr = atoi(startPosInChrStr.c_str());
		
		string tmpCigarString = tmpFieldVec[2];
		this->cigarString2jumpCodeVec(tmpCigarString, cigarStringJumpCodeVec);
		
		string syntheticReadLengthStr = tmpFieldVec[4];
		syntheticReadLength = atoi(syntheticReadLengthStr.c_str());

		string tmpSNPlocInSyntheticSNPtransSeqStr = tmpFieldVec[3];
		SNPlocInSyntheticSeq = atoi(tmpSNPlocInSyntheticSNPtransSeqStr.c_str());

		string tmpSNPposInChrStr = tmpFieldVec[5];
		SNPposInChr = atoi(tmpSNPposInChrStr.c_str());

		transcriptName = tmpFieldVec[6];

		this->generate_exonInfo_spliceSiteInfo();
	}

	int returnMapPosInChrWithMapPosInSNPseq(int mapPosInSNPseq) 
	{
		// mapPosInTranscript starts from 0, exon_start_pos and exon_end_pos start from 1
		if(exon_num == 1)
		{
			return mapPosInSNPseq + exon_start_pos[0] - 1;
		}
		for(int tmp = 0; tmp < spliceSiteInSNPseq.size(); tmp++)
		{
			int tmp_splice_site = spliceSiteInSNPseq[tmp].first;
			if(mapPosInSNPseq < tmp_splice_site)
			{
				int tmp_dist2exonEnd = tmp_splice_site - 1 - mapPosInSNPseq;
				int tmp_exon_end_pos = exon_end_pos[tmp];
				return tmp_exon_end_pos - tmp_dist2exonEnd;
			}
		}
		int last_splice_site = spliceSiteInSNPseq[exon_num - 2].first;
		int last_exon_start_pos_chrom = exon_start_pos[exon_num - 1];
		return (last_exon_start_pos_chrom + mapPosInSNPseq - last_splice_site);
	}

	void getSJsiteAndSize(int tmpStartPos, int tmpEndPos, vector< pair<int, int> >& SJsiteAndSize)
	{
		for(int tmp = 0; tmp < spliceSiteInSNPseq.size(); tmp++)
		{
			int tmp_splice_site = spliceSiteInSNPseq[tmp].first;
			int tmp_splice_size = spliceSiteInSNPseq[tmp].second;
			if((tmp_splice_site > tmpStartPos)&&(tmp_splice_site <= tmpEndPos))
			{    
				SJsiteAndSize.push_back(pair<int,int>(tmp_splice_site, tmp_splice_size));
			}
		}
	}
};
#endif