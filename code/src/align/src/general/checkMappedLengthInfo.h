// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef CHECKMAPPEDLENGTHINFO_H
#define CHECKMAPPEDLENGTHINFO_H

#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <map>

using namespace std;

class CheckMappedLengthInfo
{
private:
	int readLength;
	int mappedLength;
public:	
	int returnMappedLength()
	{
		return mappedLength;
	}

	CheckMappedLengthInfo()
	{}

	// return 0 -- OK: mapped, primary; 
	// return 1 -- unmap (primary); 
	// return 2 -- not primary;
	int checkMappedLength(const string& samStr)
	{
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

		string mapChrNameStr = samFieldVec[2];
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];

		if((mapChrNameStr == "*")||(cigarString == "*")||(mapChrPos == 0))
		{
			return 1;
		}

		string tagStr = samFieldVec[11];
		//cout << "tagStr: " << tagStr << endl;
		int HIpos = tagStr.find("HI");
		if(HIpos == tagStr.npos)
		{
			cout << "no HI found in tag sequence" << endl;
			exit(1);
		}
		int HIintNextTabPos = tagStr.find("\t", HIpos);
		string HIintStr = tagStr.substr(HIpos+5, HIintNextTabPos-HIpos-5);
		int HIint = atoi(HIintStr.c_str());
		if(HIint != 1)
		{
			return 2;
		}

		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);	

		string readSeqStr = samFieldVec[9];
		int readSeqLength = readSeqStr.length();

		int tmpMappedLength = 0;
		int tmpReadLength = 0;;
		int tmpSoftClippingLength = 0;
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			if(cigarStringJumpCodeVec[tmp].type == "S")
			{
				tmpReadLength += cigarStringJumpCodeVec[tmp].len;
				tmpSoftClippingLength += cigarStringJumpCodeVec[tmp].len;
			}
			else if(cigarStringJumpCodeVec[tmp].type == "M")
			{
				tmpMappedLength += cigarStringJumpCodeVec[tmp].len;
				tmpReadLength += cigarStringJumpCodeVec[tmp].len;
			}
			else if(cigarStringJumpCodeVec[tmp].type == "N")
			{

			}
			else if(cigarStringJumpCodeVec[tmp].type == "I")
			{
				tmpMappedLength += cigarStringJumpCodeVec[tmp].len;
				tmpReadLength += cigarStringJumpCodeVec[tmp].len;
			}
			else if(cigarStringJumpCodeVec[tmp].type == "D")
			{

			}
			else
			{
				cout << "other invalid cigarString type" << endl;
				exit(1);
			}									
		}

		if((tmpSoftClippingLength + tmpMappedLength == tmpReadLength)
			&&
			(tmpReadLength == readSeqLength))
		{
			mappedLength = tmpMappedLength;
			readLength = tmpReadLength;
			return 0;
		}
		else
		{
			cout << "tmpSoftClippingLength + tmpMappedLength != tmpReadLength OR" << endl;
			cout << "tmpReadLength == readSeqLength" << endl;
			exit(1);
		}
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
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
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}	
};
#endif
