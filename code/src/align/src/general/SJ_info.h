// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SJ_INFO_H
#define SJ_INFO_H

class Splice_Junction
{
//public:
private:
	string SJchrName;
	int SJdonerEnd;
	int SJacceptorStart;
public:
	string returnSJchrName()
	{
		return SJchrName;
	}
	int returnSJdonerEnd()
	{
		return SJdonerEnd;
	}
	int returnSJacceptorStart()
	{
		return SJacceptorStart;
	}
	Splice_Junction(string chrName, int donerPos, int acceptorPos)
	{
		SJchrName = chrName;
		SJdonerEnd = donerPos;
		SJacceptorStart = acceptorPos;		
	}

	void getSpliceJunctionFromRecord(string entryString)
	{
		int tabLocation1, tabLocation2, tabLocation3;
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);		
		
		SJchrName = entryString.substr(0, tabLocation1);	

		string spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		string spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);

		SJdonerEnd = atoi(spliceStartPosString.c_str());
		SJacceptorStart = atoi(spliceEndPosString.c_str());
	}
};

class SpliceJunction_Alignment
{
private:
	//int SJposInRead;  // value = 11, if cigarString == 10M100N90M 
	string SJchrName;
	int SJdonerEnd;
	int SJacceptorStart;
	string flankString;
public:
	string returnSJchrName()
	{
		return SJchrName;
	}
	int returnSJdonerEnd()
	{
		return SJdonerEnd;
	}
	int returnSJacceptorStart()
	{
		return SJacceptorStart;
	}
	string returnFlankString()
	{
		return flankString;
	}
	SpliceJunction_Alignment(Splice_Junction& newSJ, Index_Info* indexInfo)
	{
		SJchrName = newSJ.returnSJchrName();
		SJdonerEnd = newSJ.returnSJdonerEnd();
		SJacceptorStart = newSJ.returnSJacceptorStart();
		//cout << SJdonerEnd << " " << SJacceptorStart << endl;
		this->getFlankString(indexInfo);		
	}

	SpliceJunction_Alignment(string chrName, int donerPos, int acceptorPos, Index_Info* indexInfo)
	{
		//cout << "chrName: " << chrName.length() << " " << chrName << endl;
		SJchrName = chrName;
		SJdonerEnd = donerPos;
		SJacceptorStart = acceptorPos;
		//cout << SJdonerEnd << " " << SJacceptorStart << endl;
		this->getFlankString(indexInfo);
	}

	void getFlankString(Index_Info* indexInfo)
	{
		int chromNO = indexInfo->convertStringToInt(SJchrName);
		//cout << "chromNO: " << chromNO << endl;
		//flankString = (indexInfo->chromStr)[chromNO].substr(SJdonerEnd, 2) + (indexInfo->chromStr)[chromNO].substr(SJacceptorStart-3, 2);
		flankString = indexInfo->returnChromStrSubstr(chromNO, SJdonerEnd + 1, 2) + indexInfo->returnChromStrSubstr(chromNO, SJacceptorStart-2, 2);
		//cout << "flankString: " << flankString << endl;
	}
};

#endif