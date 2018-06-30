// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SNP_INFO_H
#define SNP_INFO_H

#include "stdio.h"
#include "stdlib.h"

using namespace std;

class SNP_Info
{
private:
	int chrNameInt;
	int pos;
	string referBase;
	string alterBase;
public:
	SNP_Info()
	{}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnPos()
	{
		return pos;
	}

	string returnReferBase()
	{
		return referBase;
	}

	void initiate(int tmpChrNameInt, int tmpPos,
		string& tmpAlterBase, Index_Info* indexInfo)
	{
		chrNameInt = tmpChrNameInt;//indexInfo->convertStringToInt(tmpChrNameStr);
		pos = tmpPos;		
		referBase = indexInfo->returnOneBaseStrInGenome(chrNameInt, tmpPos);
		alterBase = tmpAlterBase;
	}

	string returnAlterBase()
	{
		return alterBase;
	}

	string returnSNPinfoStr(Index_Info* indexInfo)
	{
		string tmpSNPstr;
		string tmpChrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string tmpPosStr = int_to_str(pos);
		tmpSNPstr = tmpChrNameStr + "\t" + tmpPosStr + "\t" + referBase + "\t" + alterBase;
		return tmpSNPstr;
	}

	string simulatedRead(Index_Info* indexInfo)
	{
		int simulatedReadSeqLength = 51;
		int halfLength = (simulatedReadSeqLength - 1) / 2;
		int SNPlocInReadSeq = halfLength + 1;
		string chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		int startPos = pos - halfLength;
		string startPosStr = int_to_str(startPos);
		string cigarString = int_to_str(simulatedReadSeqLength) + "M";
		string SNPlocInReadSeqStr = int_to_str(SNPlocInReadSeq);
		string simulatedRead_name = ">" + chrNameStr + "_" + startPosStr + "_" + cigarString + "_" + SNPlocInReadSeqStr;
		string simulatedRead_seq_1stHalf = indexInfo->returnChromStrSubstr(chrNameInt, startPos, halfLength);
		string simulatedRead_seq_2ndHalf = indexInfo->returnChromStrSubstr(chrNameInt, pos + 1, halfLength);
		string simulatedRead_seq = simulatedRead_seq_1stHalf + alterBase + simulatedRead_seq_2ndHalf;
		string simulatedRead = simulatedRead_name;
		simulatedRead += "\n";
		simulatedRead += simulatedRead_seq;
		return simulatedRead;
	}
};


#endif