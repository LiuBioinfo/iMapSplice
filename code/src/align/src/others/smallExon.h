// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <string>
#include <string.h>
#include "splice_info.h"

using namespace std;

class Small_Exon
{
public:
	//string readSeqWithDirection;
	//string mapChrName;
	int firstSegLength;
	int secondSegLength;
	int firstSegEndPosInRead;
	int secondSegStartPosInRead;
	unsigned int firstSegEndMapPos;
	unsigned int secondSegStartMapPos;

	vector< pair<int, int> > possiSJposInReadGTAG; // < small Exon Start Pos in Read, small Exon End Pos in Read >
	vector< pair<int, int> > possiSJposInReadCTAC;

	vector< pair< int, pair<int, int> > > finalSJVec; // < smallExonPosInChr, <smallExonStartPosInRead, smallExonEndPosInRead> >
	//vector< pair<int, int> > SJposInReadGTAG;
	//vector< pair<int, int> > SJposInReadCTAC;
	Small_Exon()
	{}

	Small_Exon(//const string& currentChrName, 
		int currentFirstSegLength, int currentSecondSegLength,
		int currentFirstSegEndPosInRead, int currentSecondSegStartPosInRead,
		unsigned int currentFirstSegEndMapPos, unsigned int currentSecondSegStartMapPos,
		const string& genomeStr, const string& readSeqWithDirection)
	{
		//mapChrName = currentChrName;
		firstSegLength = currentFirstSegLength;
		secondSegLength = currentSecondSegLength;
		firstSegEndPosInRead = currentFirstSegEndPosInRead; 
		secondSegStartPosInRead = currentSecondSegStartPosInRead;
		firstSegEndMapPos = currentFirstSegEndMapPos;
		secondSegStartMapPos = currentSecondSegStartMapPos;
	
		this->getPossibleSJpos(genomeStr, readSeqWithDirection);
	}

	void getPossibleSJpos(const string& chrStr, const string& readSeqWithDirection)
	{
		int leftBuffer = 4;
		int rightBuffer = 4;
		if(firstSegLength <= leftBuffer)
		{
			leftBuffer = firstSegLength - 1;
		}
		if (secondSegLength <= rightBuffer)
		{
			rightBuffer = secondSegLength - 1;
		}

		const string SJstartGTAG = "GT";
		const string SJendGTAG = "AG";

		const string SJstartCTAC = "CT";
		const string SJendCTAC = "AC";

		int pendingSeqLength = 
			leftBuffer + secondSegStartPosInRead - firstSegEndPosInRead - 1 + rightBuffer;

		//cout << "pendingReadSeq from: " << firstSegEndMapPos - leftBuffer << " Len: " << pendingSeqLength << endl;

		string pendingReadSeq = readSeqWithDirection.substr(
			firstSegEndPosInRead - leftBuffer, pendingSeqLength); 

		//cout << "pendingReadSeq: " << pendingReadSeq << endl;

		//cout << "pendingChroSeq_1 from: " << firstSegEndMapPos - leftBuffer << " Len: " << pendingSeqLength << endl;		
		string pendingChroSeq_1 = 
			chrStr.substr(firstSegEndMapPos - leftBuffer, pendingSeqLength);

		//cout << "pendingChroSeq_1: " << pendingChroSeq_1 << endl;

		//cout << "pendingChroSeq_1 from: " << secondSegStartMapPos + rightBuffer - 1 - pendingSeqLength << " Len: " << pendingSeqLength << endl;		
		string pendingChroSeq_2 = 
			chrStr.substr(secondSegStartMapPos + rightBuffer - 1 - pendingSeqLength, pendingSeqLength);

		//cout << "pendingChroSeq_2: " << pendingChroSeq_2 << endl;

		vector<int> tmpGTpos;
		vector<int> tmpAGpos;
		vector<int> tmpCTpos;
		vector<int> tmpACpos;

		int smallExonStartPosInRead;
		int smallExonEndPosInRead;

		// search for possible GTAG SJ pos
		int startSearchPos = 0;
		int foundPos = 0;
		while(1) // searching for "GT"
		{
			foundPos = pendingChroSeq_1.find(SJstartGTAG, startSearchPos);
			if(foundPos == pendingChroSeq_1.npos)
				break;
			tmpGTpos.push_back(firstSegEndPosInRead-leftBuffer+1+foundPos);
			startSearchPos = foundPos + 1;
		}

		startSearchPos = 0;
		foundPos = 0;
		while(1) // searching for "AG"
		{
			foundPos = pendingChroSeq_2.find(SJendGTAG, startSearchPos);
			if(foundPos == pendingChroSeq_2.npos)
				break;
			tmpAGpos.push_back(firstSegEndPosInRead-leftBuffer+1+foundPos);
			startSearchPos = foundPos + 1;
		}

		for(int tmp = 0; tmp < tmpGTpos.size(); tmp++)
		{
			for(int tmp2 = 0; tmp2 < tmpAGpos.size(); tmp2++)
			{
				smallExonStartPosInRead = tmpGTpos[tmp];
				smallExonEndPosInRead = tmpAGpos[tmp2] + 1;
				if(smallExonStartPosInRead < smallExonEndPosInRead)
				{
					if(this->checkMatchOrNotForTwoEndAnchor(smallExonStartPosInRead, smallExonEndPosInRead, 
							chrStr, readSeqWithDirection))
					{
						possiSJposInReadGTAG.push_back(pair<int, int> (smallExonStartPosInRead, smallExonEndPosInRead));
					}
				}
			}
		}

		// search for possible CTAC SJ pos
		startSearchPos = 0;
		foundPos = 0;
		while(1) // searching for "CT"
		{
			foundPos = pendingChroSeq_1.find(SJstartCTAC, startSearchPos);
			if(foundPos == pendingChroSeq_1.npos)
				break;
			tmpCTpos.push_back(firstSegEndPosInRead-leftBuffer+1+foundPos);
			startSearchPos = foundPos + 1;
		}

		startSearchPos = 0;
		foundPos = 0;
		while(1) // searching for "AC"
		{
			foundPos = pendingChroSeq_2.find(SJendCTAC, startSearchPos);
			if(foundPos == pendingChroSeq_2.npos)
				break;
			tmpACpos.push_back(firstSegEndPosInRead-leftBuffer+1+foundPos);
			startSearchPos = foundPos + 1;
		}

		for(int tmp = 0; tmp < tmpCTpos.size(); tmp++)
		{
			for(int tmp2 = 0; tmp2 < tmpACpos.size(); tmp2++)
			{
				smallExonStartPosInRead = tmpCTpos[tmp];
				smallExonEndPosInRead = tmpACpos[tmp2] + 1;
				if(smallExonStartPosInRead < smallExonEndPosInRead)
				{
					if(this->checkMatchOrNotForTwoEndAnchor(smallExonStartPosInRead, smallExonEndPosInRead, 
							chrStr, readSeqWithDirection))
					{
						possiSJposInReadCTAC.push_back(pair<int, int> (smallExonStartPosInRead, smallExonEndPosInRead));
					}
				}
			}
		}
	}

	bool checkMatchOrNotForTwoEndAnchor(int tmpSmallExonStartPosInRead, int tmpSmallExonEndPosInRead, 
		const string& genomeStr, const string& readSeqWithDirection)
	{
		bool firstAnchorMatch = false;
		bool secondAnchorMatch = false;

		if(tmpSmallExonStartPosInRead <= firstSegEndPosInRead + 1)
		{
			firstAnchorMatch = true;
		}
		else
		{
			size_t max_append_mismatch = (tmpSmallExonStartPosInRead - firstSegEndPosInRead - 1)/10 + 1;
			size_t mismatch_bits = 0;
			string pendingReadSeq = readSeqWithDirection.substr(firstSegEndPosInRead, tmpSmallExonStartPosInRead - firstSegEndPosInRead - 1);
			string pendingChroSeq = genomeStr.substr(firstSegEndMapPos, tmpSmallExonStartPosInRead - firstSegEndPosInRead - 1);
			firstAnchorMatch = score_string(pendingReadSeq, pendingChroSeq, max_append_mismatch, mismatch_bits);
		}

		if(tmpSmallExonEndPosInRead >= secondSegStartPosInRead - 1)
		{
			secondAnchorMatch = true;
		}
		else
		{
			size_t max_append_mismatch = (secondSegStartPosInRead - 1 - tmpSmallExonEndPosInRead)/10 + 1;
			size_t mismatch_bits = 0;
			string pendingReadSeq = 
				readSeqWithDirection.substr(tmpSmallExonEndPosInRead, secondSegStartPosInRead - 1 - tmpSmallExonEndPosInRead);
			string pendingChroSeq = 
				genomeStr.substr(secondSegStartMapPos + tmpSmallExonEndPosInRead -  secondSegStartPosInRead, secondSegStartPosInRead - 1 - tmpSmallExonEndPosInRead);
			secondAnchorMatch = score_string(pendingReadSeq, pendingChroSeq, max_append_mismatch, mismatch_bits);
		}

		return (firstAnchorMatch && secondAnchorMatch);
	}



};