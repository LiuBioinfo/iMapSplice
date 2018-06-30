// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef DECOYALIGN_INFO_H
#define DECOYALIGN_INFO_H

#include "stdio.h"
#include "stdlib.h"

#include "../../../general/splice_info.h"

using namespace std;

class DecoyAlignElement_Info
{
private:
	int chrNameInt_gt;
	int chrPos_gt;
	string cigarString_gt;

	int SNPlocInRead;
	int simulatedReadLength;
	int SNPposInChr;
	string transcriptID;
	int SNPposInTranscript;
	
	int simulatedReadAlignment_flag;
	int simulatedReadAlignment_chrNameInt;
	int simulatedReadAlignment_chrPos;
	string simulatedReadAlignment_cigarString;

	string oriSamStr;
public:
	DecoyAlignElement_Info()
	{}

	int return_chrNameInt_gt()
	{
		return chrNameInt_gt;
	}

	int return_SNPposInChr()
	{
		return SNPposInChr;
	}

	string return_decoyAlignElementStr()
	{
		return oriSamStr;
	}

	void initiate(string& tmpSamStr,
		int tmpChrNameInt_gt,
		int tmpChrPos_gt,
		string& tmpCigarString_gt,
		int tmpSNPlocInRead,
		int tmpSimulatedReadLength,
		int tmpSNPposInChr,
		string& tmpTranscriptID,
		int tmpSNPposInTranscript,
		int tmpSimulatedReadAlignment_flag,
		int tmpSimulatedReadAlignment_chrNameInt,
		int tmpSimulatedReadAlignment_chrPos,
		string& tmpSimulatedReadAlignment_cigarString)
	{
		chrNameInt_gt = tmpChrNameInt_gt;
		chrPos_gt = tmpChrPos_gt;
		cigarString_gt = tmpCigarString_gt;
		SNPlocInRead = tmpSNPlocInRead;
		simulatedReadLength = tmpSimulatedReadLength;
		SNPposInChr = tmpSNPposInChr;
		transcriptID = tmpTranscriptID;
		SNPposInTranscript = tmpSNPposInTranscript;
		simulatedReadAlignment_flag = tmpSimulatedReadAlignment_flag;
		simulatedReadAlignment_chrNameInt = tmpSimulatedReadAlignment_chrNameInt;
		simulatedReadAlignment_chrPos = tmpSimulatedReadAlignment_chrPos;
		simulatedReadAlignment_cigarString = tmpSimulatedReadAlignment_cigarString;

		oriSamStr = tmpSamStr;
	}
};

class DecoyAlign_Info
{
private:
	int decoyAlign_chrNameInt;
	int decoyAlign_chrPos;
	vector<DecoyAlignElement_Info> decoyAlignElementVec;
public:
	DecoyAlign_Info()
	{}

	string returnSNPinfoVecStr(Index_Info* indexInfo)
	{
		//cout << "start to do returnSNPinfoVecStr ..." << endl;
		vector<int> snpChrNameIntVec;
		vector<int> snpChrPosVec;
		int elementSize = decoyAlignElementVec.size();
		//cout << "elementSize: " << elementSize << endl;
		for(int tmp = 0; tmp < elementSize; tmp++)
		{
			int tmpSNPchrNameInt = decoyAlignElementVec[tmp].return_chrNameInt_gt();
			int tmpSNPposInChr = decoyAlignElementVec[tmp].return_SNPposInChr();
			//cout << "tmpSNPchrNameInt: " << tmpSNPchrNameInt << endl;
			//cout << "tmpSNPposInChr: " << tmpSNPposInChr << endl;
			int tmpCurrentAddedSnpInfoVecSize = snpChrNameIntVec.size();
			bool tmpSNPinfoExistsBool = false;
			for(int tmp2 = 0; tmp2 < tmpCurrentAddedSnpInfoVecSize; tmp2++)
			{
				int tmpExistingSnpChrNameInt = snpChrNameIntVec[tmp2];
				int tmpExistingSnpChrPos = snpChrPosVec[tmp2];
				if((tmpSNPchrNameInt == tmpExistingSnpChrNameInt)
					&&(tmpSNPposInChr == tmpExistingSnpChrPos))
				{
					tmpSNPinfoExistsBool = true;
					break;
				}
			}
			//cout << "tmpSNPinfoExistsBool: " << tmpSNPinfoExistsBool << endl;
			if(tmpSNPinfoExistsBool)
			{}
			else
			{
				snpChrNameIntVec.push_back(tmpSNPchrNameInt);
				snpChrPosVec.push_back(tmpSNPposInChr);
			}
		}
		string tmpStr;
		int snpNum = snpChrNameIntVec.size();
		//cout << "snpNum: " << snpNum << endl;
		string snpNumStr = int_to_str(snpNum);
		tmpStr += snpNumStr;
		//tmpStr += "\t";
		for(int tmp = 0; tmp < snpNum; tmp++)
		{
			int tmpSNPchrNameInt = snpChrNameIntVec[tmp];
			int tmpSNPchrPos = snpChrPosVec[tmp];
			string tmpSNPchrName = indexInfo->returnChrNameStr(tmpSNPchrNameInt);
			string tmpSNPchrPosStr = int_to_str(tmpSNPchrPos);
			//cout << "tmpSNPchrName: " << tmpSNPchrName << endl;
			//cout << "tmpSNPchrPosStr: " << tmpSNPchrPosStr << endl;
			tmpStr += "\t";
			tmpStr += tmpSNPchrName;
			tmpStr += ":";
			tmpStr += tmpSNPchrPosStr;
		}
		return tmpStr;
	}

	int return_decoyAlign_chrNameInt()
	{
		return decoyAlign_chrNameInt;
	}

	int return_decoyAlign_chrPos()
	{
		return decoyAlign_chrPos;
	}

	string return_decoyAlignmentStr(Index_Info* indexInfo)
	{
		string decoyAlignmentStr;
		string decoyAlign_chrName = indexInfo->returnChrNameStr(decoyAlign_chrNameInt);
		string decoyAlign_chrPosStr = int_to_str(decoyAlign_chrPos); 
		decoyAlignmentStr = decoyAlign_chrName + "\t" + decoyAlign_chrPosStr 
				+ "\t" + int_to_str(decoyAlignElementVec.size());
		int decoyAlignElementVecSize = decoyAlignElementVec.size();
		for(int tmp = 0; tmp < decoyAlignElementVecSize; tmp++)
		{
			string tmpDecoyAlignElementStr 
				= decoyAlignElementVec[tmp].return_decoyAlignElementStr();
			decoyAlignmentStr += "\n";
			decoyAlignmentStr += tmpDecoyAlignElementStr;
		}
		return decoyAlignmentStr;
	}

	void initiateAndAddNewDecoyAlignElement(string& tmpSamStr,
		int tmpSNPcorrespondingChrNameIntInSimulatedReadAlignment, 
		int tmpSNPcorrespondingPosInSimulatedReadAlignment,
		int tmpChrNameInt_gt, int tmpChrPos_gt, 
		string& tmpCigarString_gt, int tmpSNPlocInRead, 
		int tmpSimulatedReadLength, 
		int tmpSNPposInChr, string& tmpTranscriptID, 
		int tmpSNPlocInTranscript, 
		int tmpSimulatedReadAlignment_flag, 
		int tmpSimulatedReadAlignment_chrNameInt, 
		int tmpSimulatedReadAlignment_chrPos, 
		string& tmpSimulatedReadAlignment_cigarString)
	{
		decoyAlign_chrNameInt 
			= tmpSNPcorrespondingChrNameIntInSimulatedReadAlignment;
		decoyAlign_chrPos
			= tmpSNPcorrespondingPosInSimulatedReadAlignment;
		this->addNewDecoyAlignElement(tmpSamStr,
			tmpChrNameInt_gt, tmpChrPos_gt, 
			tmpCigarString_gt, tmpSNPlocInRead, 
			tmpSimulatedReadLength, tmpSNPposInChr, 
			tmpTranscriptID, tmpSNPlocInTranscript, 
			tmpSimulatedReadAlignment_flag, 
			tmpSimulatedReadAlignment_chrNameInt, 
			tmpSimulatedReadAlignment_chrPos, 
			tmpSimulatedReadAlignment_cigarString);
	}

	void addNewDecoyAlignElement(string& tmpSamStr,
		int tmpChrNameInt_gt, int tmpChrPos_gt, 
		string& tmpCigarString_gt, int tmpSNPlocInRead, 
		int tmpSimulatedReadLength, 
		int tmpSNPposInChr, string& tmpTranscriptID, 
		int tmpSNPlocInTranscript, 
		int tmpSimulatedReadAlignment_flag, 
		int tmpSimulatedReadAlignment_chrNameInt, 
		int tmpSimulatedReadAlignment_chrPos, 
		string& tmpSimulatedReadAlignment_cigarString)
	{
		DecoyAlignElement_Info tmpDecoyAlignElementInfo;
		tmpDecoyAlignElementInfo.initiate(tmpSamStr,
			tmpChrNameInt_gt, tmpChrPos_gt, 
			tmpCigarString_gt, tmpSNPlocInRead, 
			tmpSimulatedReadLength, tmpSNPposInChr, 
			tmpTranscriptID, tmpSNPlocInTranscript, 
			tmpSimulatedReadAlignment_flag, 
			tmpSimulatedReadAlignment_chrNameInt, 
			tmpSimulatedReadAlignment_chrPos, 
			tmpSimulatedReadAlignment_cigarString);
		decoyAlignElementVec.push_back(tmpDecoyAlignElementInfo);
	}
};
#endif