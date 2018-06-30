// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef HEADERSECTION_INFO_H
#define HEADERSECTION_INFO_H

#include <stdlib.h>
#include <string>
#include <string.h>

using namespace std;

class HeaderSection_Info
{
private:
	string SQ_str;
	string headerSectionInfoStr;
public:
	string returnHeaderSectionInfoStr()
	{
		headerSectionInfoStr = SQ_str;
		return headerSectionInfoStr;
	}

	HeaderSection_Info()
	{
		SQ_str = "@SQ\tSN:Not_avaliable";
		//headerSectionInfoStr = 
	}

	HeaderSection_Info(Index_Info* indexInfo//, Option_Info* optionInfo
		)
	{
		int chromNum = indexInfo->returnChrNameStrSize();
		SQ_str = "@SQ\tSN:" + indexInfo->returnChrNameStr(0)
			+ "\tLN:" + int_to_str(indexInfo->returnChromLength(0)) + "\n";
		for(int tmp = 1; tmp < (
			chromNum - 1
			); tmp++)
		{
			SQ_str = SQ_str
				+ "@SQ\tSN:" + indexInfo->returnChrNameStr(tmp)
				+ "\tLN:" + int_to_str(indexInfo->returnChromLength(tmp)) + "\n"; 
		}
		SQ_str = SQ_str
			+ "@SQ\tSN:" + indexInfo->returnChrNameStr(chromNum - 1)
			+ "\tLN:" + int_to_str(indexInfo->returnChromLength(chromNum - 1)-1); 
	}

};

#endif