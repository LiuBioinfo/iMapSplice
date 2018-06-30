// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <string>
#include <string.h>
//#include "splice_info.h"

using namespace std;

class JumpCode_Info
{
public:
	vector< pair< int, int > > SJforCompareVec;
	vector<Jump_Code> cigarStringJumpCode;

	JumpCode_Info()
	{}

	void jumpCodeStr2jumpCodeVec(const string& jumpCodeStr)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNID";
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
				cigarStringJumpCode.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	void jumpCodeVec2SJforCompareVec(int alignChromPos) 
	{
		int tmpPosInRead = 0; 
		int tmpDonerPosInChr = alignChromPos - 1;
		int tmpAcceptorPosInChr = 0;
		//cout << "start to get sjVec ..." << endl;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			//cout << "JumpCode[" << tmp << "]: " << cigarStringJumpCode[tmp].len << " " << cigarStringJumpCode[tmp].type << endl;
			if(cigarStringJumpCode[tmp].type == "S")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "M")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;
			}
			else if (cigarStringJumpCode[tmp].type == "N")
			{
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;

				tmpAcceptorPosInChr = tmpDonerPosInChr + cigarStringJumpCode[tmp].len + 1;
				//cout << "tmpAcceptorPosInChr: " << tmpAcceptorPosInChr << endl;
				//SpliceJunction_Alignment* tmpSJ = 
				//	new SpliceJunction_Alignment(alignChromName, tmpDonerPosInChr, tmpAcceptorPosInChr, indexInfo);
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				//tmpSJ
				SJforCompareVec.push_back(pair <int, int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "I")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "D")
			{
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else
			{
				cout << "other jumpCodeType" << endl;
			}
		}
	}

};