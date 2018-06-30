// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef UNFIXEDHEAD_H
#define UNFIXEDHEAD_H

#include <string>
#include <string.h>
//#include "splice_info.h"

using namespace std;

class Unfixed_Head
{
//public:
private:
	string readName;
	string alignDirection; 
	unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedHeadLength;
	int midPartLength;
	string readSeqOriginal;
	//int readLength;

	int midPartMapPosInChr; // map pos of 1st base in midPart
	string midPartMapChrName;
	int midPartMapChrInt;
public:
	vector<int> possiSJposInRead; // start from 1, end at unfixedHeadLength + bufferLength

	vector<int> possiGTAGpos; // only AG end checked
	vector<int> possiCTACpos; // only AC end checked

	vector<int> possiGTAGpos_mismatch;
	vector<int> possiCTACpos_mismatch;
	vector< vector<int> > possiGTAGpos_mismatchPos;
	vector< vector<int> > possiCTACpos_mismatchPos;
	vector< vector<char> > possiGTAGpos_mismatchChar;
	vector< vector<char> > possiCTACpos_mismatchChar;

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< int > GTAGsjPos_mismatch; // sequence around the SJ has been checked including short anchor and midpart
	vector< int > CTACsjPos_mismatch; // sequence around the SJ has been checked 
	vector< vector<int> > GTAGsjPos_mismatchPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< vector<int> > CTACsjPos_mismatchPos; // sequence around the SJ has been checked 
	vector< vector<char> > GTAGsjPos_mismatchChar; // sequence around the SJ has been checked including short anchor and midpart
	vector< vector<char> > CTACsjPos_mismatchChar; // sequence around the SJ has been checked 

	vector< pair<int, int> > SJposFromRemappingVec;
	vector< int > SJposFromRemappingVec_mismatch;
	vector< vector<int> > SJposFromRemappingVec_mismatchPosVec;
	vector< vector<char> > SJposFromRemappingVec_mismatchCharVec;

	vector< pair<int, int> > SJposFromRemappingVec_candi;
	vector< int > SJposFromRemappingVec_candi_mismatch;	
	vector< vector<int> > SJposFromRemappingVec_candi_mismatchPosVec;
	vector< vector<char> > SJposFromRemappingVec_candi_mismatchCharVec;


	// vector< int > extendedHeadAcceptorSideMatchLenVec;
	// vector< int > extendedHeadSpliceJunctionSizeVec; 
	// vector< vector<Jump_Code> > extendedHeadBackwardJumpCodeVec_withAlignInfer;
	// vector< vector<int> > extendedHeadMismatchPosVec_withAlignInfer;
	// vector< vector<char> > extendedHeadMismatchCharVec_withAlignInfer;


	unsigned int returnMidPartMapPosInWholeGenome()
	{
		return midPartMapPosInWholeGenome;
	}
	int returnUnfixedHeadLength()
	{
		return unfixedHeadLength;
	}
	int returnMidPartMapPosInChr()
	{
		return midPartMapPosInChr;
	}
	int returnMidPartMapChrInt()
	{
		return midPartMapChrInt;
	}
	Unfixed_Head()
	{//readLength = 100;
	}

	Unfixed_Head(int head_length, int mapPos, string mapChrName)
	{
		unfixedHeadLength = head_length;
		midPartMapPosInChr = midPartMapPosInChr;
		midPartMapChrName = mapChrName;
	}

	~Unfixed_Head()
	{}

	int countMismatchNumInBufferSeq(int buffer, const string& readSeqWithDirection, Index_Info* indexInfo)
	{
		int countNum = 0;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		for(int tmp = 0; tmp < buffer + 1; tmp ++)
		{
			if(readSeqWithDirection.at(unfixedHeadLength + 1 + tmp - 1) != 
				//indexInfo->chromStr[chromNameInt].at(midPartMapPosInChr + tmp - 1) 
				indexInfo->returnOneBaseCharInGenome(chromNameInt, midPartMapPosInChr + tmp)
				)
			{
				countNum ++;
			}
			else
			{}
		}
		return countNum;
	}

	void getUnfixedHeadInfoFromRecord(PE_Read_Info& readInfo, bool end1, Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		if(end1)
		{
			readName = readInfo.returnReadName_1();
			readSeqOriginal = readInfo.returnReadSeq_1();
		}
		else
		{
			readName = readInfo.returnReadName_2();
			readSeqOriginal = readInfo.returnReadSeq_2();
		}
		alignDirection = alignInfo->returnAlignDirection();

		midPartMapChrName = alignInfo->returnAlignChromName();
		midPartMapPosInChr = alignInfo->returnAlignChromPos();

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			(unsigned int)midPartMapChrInt, (unsigned int)midPartMapPosInChr);
	}

	void getUnfixedHeadInfoFromRecordWithAlignInfoType(
		PE_Read_Info& readInfo, int alignInfoType, 
		Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		if((alignInfoType == 1) || (alignInfoType == 2))
		{
			readName = readInfo.returnReadName_1();
			readSeqOriginal = readInfo.returnReadSeq_1();
		}
		else
		{
			readName = readInfo.returnReadName_2();
			readSeqOriginal = readInfo.returnReadSeq_2();
		}
		alignDirection = alignInfo->returnAlignDirection();

		midPartMapChrName = alignInfo->returnAlignChromName();
		midPartMapPosInChr = alignInfo->returnAlignChromPos();

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			(unsigned int)midPartMapChrInt, (unsigned int)midPartMapPosInChr);
	}

	void getUnfixedHeadInfoFromRecordWithAlignInfoType_new(
		PE_Read_Info& readInfo, int alignInfoType, 
		//int index_peAlignInfo,
		Alignment_Info* alignInfo, 
		Index_Info* indexInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();

		midPartMapChrName = alignInfo->returnAlignChromName();
		midPartMapPosInChr = alignInfo->returnAlignChromPos();

		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);
		otherCigarString = alignInfo->otherJumpCodeVec2Str();
		unfixedHeadLength = (alignInfo->cigarStringJumpCode)[0].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[1].len;

		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			(unsigned int)midPartMapChrInt, (unsigned int)midPartMapPosInChr);
	}

	/*
	void getUnfixedHeadInfoFromRecordWithAlignInfoType_new_new(
		PE_Read_Info& readInfo, int alignInfoType, 
		int index_peAlignInfo,
		//Alignment_Info* alignInfo, 
		Index_Info* indexInfo)
	{
		if(alignInfoType == 1)
		{
			this->getUnfixedHeadInfoFromRecordWithAlignInfoType_new(
				readInfo, 1, norAlignmentInfo_PE_1[index_peAlignInfo], indexInfo);
		}
		else if(alignInfoType == 2)
		{
			this->getUnfixedHeadInfoFromRecordWithAlignInfoType_new(
				readInfo, 2, rcmAlignmentInfo_PE_1[index_peAlignInfo], indexInfo);
		}
		else if(alignInfoType == 3)
		{
			this->getUnfixedHeadInfoFromRecordWithAlignInfoType_new(
				readInfo, 3, norAlignmentInfo_PE_2[index_peAlignInfo], indexInfo);
		}
		else if(alignInfoType == 4)
		{
			this->getUnfixedHeadInfoFromRecordWithAlignInfoType_new(
				readInfo, 4, rcmAlignmentInfo_PE_2[index_peAlignInfo], indexInfo);
		}
		else
		{
			cout << "alignInfoType invalid, in unfixedHead ..." << endl;
			exit(1);
		}
	}*/

	void fromRecordStringToClass(char* unfixedHeadRecordString, int readLength)
	{
		char readNameCharForUnfixedHead[100]; 
		char alignDirectionCharForUnfixedHead;
		unsigned int midPartMapPosForUnfixedHead;
		int firstJumpCodeLengthForUnfixedHead; // type = S
		int secondJumpCodeLengthForUnfixedHead;
		char otherCigarStringCharForUnfixedHead[50]; 
		char readSeqCharForUnfixedHead[600];
		
		sscanf(unfixedHeadRecordString, "%s\t%s\t%d\t%s\t%d\t%d\t%s",
			readNameCharForUnfixedHead, &alignDirectionCharForUnfixedHead,
			&midPartMapPosForUnfixedHead, otherCigarStringCharForUnfixedHead,
			&firstJumpCodeLengthForUnfixedHead, &secondJumpCodeLengthForUnfixedHead,
			readSeqCharForUnfixedHead);

		readName = readNameCharForUnfixedHead;
		alignDirection = alignDirectionCharForUnfixedHead;
		midPartMapPosInWholeGenome = midPartMapPosForUnfixedHead;// + firstJumpCodeLengthForUnfixedHead;
		otherCigarString = otherCigarStringCharForUnfixedHead;
		unfixedHeadLength = firstJumpCodeLengthForUnfixedHead;
		midPartLength = secondJumpCodeLengthForUnfixedHead;
		readSeqOriginal = readSeqCharForUnfixedHead;
		//int readLength = 100;
		readSeqOriginal = readSeqOriginal.substr(0, readLength);
	}

	void getPossibleSJpos(const string& readSeqWithDirection,/* const string& chromSeq,*/ Index_Info* indexInfo)
	{

		int bufferLength = SINGLE_ANCHOR_TARGETMAPPING_BUFFER;
		if(bufferLength > midPartLength)
		{
			bufferLength = midPartLength;
		}

		string pendingReadSeq = readSeqWithDirection.substr(0, unfixedHeadLength+bufferLength);
		//string pendingChroSeq = chromSeq.substr(midPartMapPosInWholeGenome-unfixedHeadLength-1, unfixedHeadLength+bufferLength);		
		int pendingChroSeqStartPos = midPartMapPosInChr - unfixedHeadLength-1;
		if((pendingChroSeqStartPos < 0)||
			(pendingChroSeqStartPos + unfixedHeadLength+bufferLength > 
				//(indexInfo->chromLength[midPartMapChrInt])
				indexInfo->returnChromLength(midPartMapChrInt)
				)
			)
		{
			return;
		}	
		string pendingChroSeq = indexInfo->returnChromStrSubstr(indexInfo->convertStringToInt(midPartMapChrName),
			midPartMapPosInChr - unfixedHeadLength, unfixedHeadLength+bufferLength);

		string SJendStrAG = "AG";
		string SJendStrAC = "AC";

		//search for "AG" 
		int startSearchPos = 0;
		int foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAG, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(unfixedHeadLength-foundPos-2 > 0)
			{			
				int max_mismatch = (unfixedHeadLength-foundPos-2)/10;				
				// cout << endl << "foundPos: " << foundPos << endl;
				// cout << "max_mismatch: " << max_mismatch << endl;
				// cout << "pendingReadSeq: " << endl << pendingReadSeq << endl;
				// cout << "pendingChroSeq: " << endl << pendingChroSeq << endl;
				// cout << "1st base of MidPart: " << foundPos + 3 << endl;
				// cout << "toFixMatch Seq: " << unfixedHeadLength - foundPos - 2 << endl;
				FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
				bool matchBool = fixMatchInfo->fixMatch(
					pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
					pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
					max_mismatch, foundPos+3);
				if(matchBool)
				{
					possiGTAGpos.push_back(foundPos+3);
					possiGTAGpos_mismatch.push_back(fixMatchInfo->returnMismatchNum());

					vector<int> tmpMismatchPosVec;
					this->getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
					possiGTAGpos_mismatchPos.push_back(tmpMismatchPosVec);

					vector<char> tmpMismatchCharVec;
					this->getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
					possiGTAGpos_mismatchChar.push_back(tmpMismatchCharVec);	
				}
				delete fixMatchInfo;
			}
			else
			{
				possiGTAGpos.push_back(foundPos+3);	
				possiGTAGpos_mismatch.push_back(0);
				vector<int> tmpMismatchPosVec;
				possiGTAGpos_mismatchPos.push_back(tmpMismatchPosVec);
				vector<char> tmpMismatchCharVec;
				possiGTAGpos_mismatchChar.push_back(tmpMismatchCharVec);
			}
			startSearchPos = foundPos+1;
		}

		//cout << "searching for ACs" << endl;
		//search for "AC" 
		startSearchPos = 0;
		foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAC, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;			
			if(unfixedHeadLength-foundPos-2 > 0)
			{
				int max_mismatch = (unfixedHeadLength-foundPos-2)/10;				
				FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
				bool matchBool = fixMatchInfo->fixMatch(
					pendingReadSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
					pendingChroSeq.substr(foundPos+2, unfixedHeadLength-foundPos-2), 
					max_mismatch, foundPos+3);
				if(matchBool)
				{
					possiCTACpos.push_back(foundPos+3);
					possiCTACpos_mismatch.push_back(fixMatchInfo->returnMismatchNum());
					//if(STORE_MISMATCH_POS)
					//{
						vector<int> tmpMismatchPosVec;
						this-> getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
						possiCTACpos_mismatchPos.push_back(tmpMismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							vector<char> tmpMismatchCharVec;
							this->getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
							possiCTACpos_mismatchChar.push_back(tmpMismatchCharVec);
						//}
					//}		
				}
				delete fixMatchInfo;
			}
			else
			{
				possiCTACpos.push_back(foundPos+3);	
				possiCTACpos_mismatch.push_back(0);
				//if(STORE_MISMATCH_POS)
				//{
					vector<int> tmpMismatchPosVec;
					//this-> getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
					possiCTACpos_mismatchPos.push_back(tmpMismatchPosVec);
					//if(STORE_MISMATCH_CHA)
					//{
						vector<char> tmpMismatchCharVec;
						//this-> getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
						possiCTACpos_mismatchChar.push_back(tmpMismatchCharVec);
					//}
				//}	
			}
			startSearchPos = foundPos+1;
			//cout << "found at foundPos " << endl;
		}
		// cout << endl << "***********************************" << endl 
		// 	<< "end of getPossibleSJpos, fix unfixed heads ...." << endl
		// 	<< "************************************" << endl;

		return;
	}



	void printSJ()
	{
		cout << "...... GTAG sj ......" << endl;
		for(int tmp = 0; tmp < GTAGsjPos.size(); tmp++)
		{
			cout << "SJ:" << tmp << " headLength: " << GTAGsjPos[tmp].first - 1 
			<< " SJdistance: " <<  GTAGsjPos[tmp].second << endl;
		}
		cout << "...... CTAC sj ......" << endl;
		for(int tmp = 0; tmp < CTACsjPos.size(); tmp++)
		{
			cout << "SJ:" << tmp << " headLength: " << CTACsjPos[tmp].first - 1 
			<< " SJdistance: " <<  CTACsjPos[tmp].second << endl;
		}
	}

	void printUnfixedHeadInfo()
	{
		cout << readName << endl << alignDirection << endl << midPartMapPosInWholeGenome << endl << unfixedHeadLength << readSeqOriginal << endl; 
	}

	bool SJsearchInSJhash_areaStringHash(SJhash_Info* SJinfo, 
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize)
	{
		bool SJfoundInSJhash = false;
		int buffer = 4;
		if(buffer > midPartLength - 1)
			buffer = midPartLength - 1;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		int areaNOmin = (int)((midPartMapPosInChr - unfixedHeadLength + 1)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr + buffer)/areaSize);
		vector<int> SJacceptorSiteVec;
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			//cout << "tmpArea: " << tmpArea << endl;
			SJareaHashIter tmpSJareaHashIter 
				= ((SJinfo->SJendPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJendPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJacceptorPos = (*intSetIter);
					//cout << "tmpSJacceptorPos: " << tmpSJacceptorPos << endl;
					if( (tmpSJacceptorPos >= midPartMapPosInChr - unfixedHeadLength + 1) 
						&& (tmpSJacceptorPos <= midPartMapPosInChr + buffer) )
						SJacceptorSiteVec.push_back(tmpSJacceptorPos);
				}
			}
			else
			{}
		}
		//cout << "SJacceptorSiteVec.size(): " << SJacceptorSiteVec.size() << endl;
		string readPendingStr = readSeqWithDirection.substr(
			0, unfixedHeadLength + buffer + 1);
		if(SJacceptorSiteVec.size() > 0)
		{
			for(int tmp = 0; tmp < SJacceptorSiteVec.size(); tmp++)
			{
				vector<int> tmpDonerSiteVec;
				int tmpSJacceptorSite = SJacceptorSiteVec[tmp];
				int tmpHeadLength = unfixedHeadLength + tmpSJacceptorSite - midPartMapPosInChr;
				//cout << "tmpSJacceptorSite: " << tmpSJacceptorSite << endl;
				//cout << "tmpHeadLength: " << tmpHeadLength << endl; 
				SplicePosHashIter tmpPosHashIter 
					= (SJinfo->spliceJunctionReverse)[chromNameInt].find(tmpSJacceptorSite);
				if(tmpPosHashIter != ((SJinfo->spliceJunctionReverse)[chromNameInt]).end())
				{
					if( tmpHeadLength >= (SJinfo->anchorStringLength) )
					{
						string tmpAnchorString = readSeqWithDirection.substr(
							tmpHeadLength - (SJinfo->anchorStringLength),  (SJinfo->anchorStringLength));
						//cout << "tmpAnchorString: " << tmpAnchorString << endl;
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							//cout << "found tmpAnchorString " << endl;
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
						}
						else
						{
							//cout << "not found tmpAnchorString" << endl;
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
									tmpDonerSiteVec.push_back(*tmpIntSetIter);
							}
							else
								cout << "error! * should be found in endStrHash" << endl; 						
						}
					}
					else
					{
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
						}
						else
							cout << "error! * should be found in endStrHash" << endl; 
					}
				}
				else
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJacceptorSite should be found in hash ! " << endl;

				//cout << "tmpDonerSiteVec.size(): " << tmpDonerSiteVec.size() << endl;
				for(int tmpVecNO = 0; tmpVecNO < tmpDonerSiteVec.size(); tmpVecNO ++)
				{
					int tmpSJdonerEndPosInRead = tmpHeadLength;
					int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;

					int tmpSJdonerEndPosInChr = tmpDonerSiteVec[tmpVecNO];
					int tmpSJacceptorStartPosInChr = tmpSJacceptorSite;

					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					bool matchBool;
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();

					if(tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead >= 0)
					{
						chromDonerEndStr = indexInfo->returnChromStrSubstr(chromNameInt, 
							tmpSJdonerEndPosInChr - tmpSJdonerEndPosInRead + 1, tmpSJdonerEndPosInRead);
						chromAcceptorStartStr = indexInfo->returnChromStrSubstr(chromNameInt,
							tmpSJacceptorStartPosInChr, unfixedHeadLength + 1 - tmpSJacceptorStartPosInRead + buffer + 1);

						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;
						int max_mismatch = (unfixedHeadLength)/MATCH_BASE_PER_MISMATCH_BASE;
						matchBool = fixMatchInfo->fixMatch(readPendingStr, chromPendingStr, max_mismatch, 1);
					}
					else
						matchBool = false;

					if(matchBool)
					{
						SJfoundInSJhash = true;
						//cout << "SJ found: at " << tmpSJdonerEndPosInRead << " SJsize: " << tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1 << endl;
						SJposFromRemappingVec_candi.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
						SJposFromRemappingVec_candi_mismatch.push_back(fixMatchInfo->returnMismatchNum());

						vector<int> tmpMismatchPosVec;
						this->getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
						SJposFromRemappingVec_candi_mismatchPosVec.push_back(tmpMismatchPosVec);

						vector<char> tmpMismatchCharVec;
						this->getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
						SJposFromRemappingVec_candi_mismatchCharVec.push_back(tmpMismatchCharVec);
					}
					delete fixMatchInfo;
				}
			}
		}
		else
		{}	

		if(SJfoundInSJhash)
		{
			//int mismatchNumInBufferSeq = this->countMismatchNumInBufferSeq(buffer, readSeqWithDirection, indexInfo);
			int mismatchNumInBufferSeq = 0;
			this->filterSJposFromRemapping_candi(mismatchNumInBufferSeq);
		}
		
		return SJfoundInSJhash;
	}
	
	void generateCandiSJspliceSiteVec(SJhash_Info* SJinfo,
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize, vector< pair<int, vector<int> > >& candiSJspliceSiteVec)
	{
		int buffer = 4;
		if(buffer > midPartLength - 1)
			buffer = midPartLength - 1;

		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);

		int areaNOmin = (int)((midPartMapPosInChr - unfixedHeadLength + 1)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr + buffer)/areaSize);

		vector<int> SJacceptorSiteVec;

		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SJareaHashIter tmpSJareaHashIter 
				= ((SJinfo->SJendPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJendPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJacceptorPos = (*intSetIter);
					if( (tmpSJacceptorPos >= midPartMapPosInChr - unfixedHeadLength + 1) 
						&& (tmpSJacceptorPos <= midPartMapPosInChr + buffer) )
					{
						SJacceptorSiteVec.push_back(tmpSJacceptorPos);
					}
				}
			}
			else
			{}
		}

		if(SJacceptorSiteVec.size() > 0)
		{
			//string readPendingStr = readSeqWithDirection.substr(
			//	0, unfixedHeadLength + buffer + 1);
			for(int tmp = 0; tmp < SJacceptorSiteVec.size(); tmp++)
			{
				vector<int> tmpDonerSiteVec;

				int tmpSJacceptorSite = SJacceptorSiteVec[tmp];

				//cout << "tmpSJacceptorSite: " << tmpSJacceptorSite << endl;

				int tmpHeadLength = unfixedHeadLength + tmpSJacceptorSite - midPartMapPosInChr;
				
				//cout << "tmpHeadLength: " << tmpHeadLength << endl;

				SplicePosHashIter tmpPosHashIter 
					= (SJinfo->spliceJunctionReverse)[chromNameInt].find(tmpSJacceptorSite);
				
				if(tmpPosHashIter != ((SJinfo->spliceJunctionReverse)[chromNameInt]).end())
				{
					if( tmpHeadLength >= (SJinfo->anchorStringLength) )
					{
						string tmpAnchorString = readSeqWithDirection.substr(tmpHeadLength - (SJinfo->anchorStringLength) , (SJinfo->anchorStringLength));

						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
							{
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
								//cout << "exact anchorStr found donerSite: " << (*tmpIntSetIter) << endl;
							}
						}
						else
						{
							//cout << "error! * should be found in endStrHash" << endl;
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
								{
									tmpDonerSiteVec.push_back(*tmpIntSetIter);
									//cout << "exact anchorStr no found donerSite: " << (*tmpIntSetIter) << endl;
								}
							}
							else
							{
								cout << "error! * should be found in endStrHash" << endl; 
							}							
						}
					}
					else
					{
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin(); 
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++)
							{
								tmpDonerSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							cout << "error! * should be found in endStrHash" << endl; 
						}
					}
				}
				else
				{
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJacceptorSite should be found in hash ! " << endl;
				}

				if(tmpDonerSiteVec.size() > 0)
				{
					candiSJspliceSiteVec.push_back(pair<int, vector<int> >(tmpSJacceptorSite, tmpDonerSiteVec));
				}
			}
		}
		else
		{}			
	}
	
	bool generateInferedPathVec(vector< pair<int,int> >& spliceSitePairVec, 
		vector< vector<Jump_Code> >& inferedPathJumpCodeVecVec,
		vector< vector<int> >& inferedPathMismatchPosVecVec,
		vector< vector<char> >& inferedPathMismatchCharVecVec,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo, 
		Index_Info* indexInfo, const string& readSeqWithDirection)
	{
		//cout << "start tp do generateInferedPathVec .... in  unfixedHead.h" << endl;

		bool inferedPathGeneratedBool = false;
		int tmpChromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "spliceSitePairVecSize: " << spliceSitePairVec.size() << endl;
		for(int tmp = 0; tmp < spliceSitePairVec.size(); tmp++)
		{
			vector<Jump_Code> tmpInferedPathJumpCodeVec;
			vector<int> tmpInferedPathMismatchPosVec;
			vector<char> tmpInferedPathMismatchCharVec;
			int tmpSpliceDonerEndPosInChr = spliceSitePairVec[tmp].first;
			int tmpSpliceAcceptorStartPosInChr = spliceSitePairVec[tmp].second;
			int tmpUnfixedHeadLength 
				= unfixedHeadLength + tmpSpliceAcceptorStartPosInChr - midPartMapPosInChr;
			int tmpSpliceDonerEndPosInRead = tmpUnfixedHeadLength;
			int tmpSpliceAcceptorStartPosInRead = tmpUnfixedHeadLength + 1;

			// cout << "tmpSpliceDonerEndPosInChr: " << tmpSpliceDonerEndPosInChr << endl;
			// cout << "tmpSpliceAcceptorStartPosInChr: " << tmpSpliceAcceptorStartPosInChr << endl;
			// cout << "tmpSpliceDonerEndPosInRead: " << tmpSpliceDonerEndPosInRead << endl;
			// cout << "tmpSpliceAcceptorStartPosInRead: " << tmpSpliceAcceptorStartPosInRead << endl;
			// cout << "*********************************" << endl
			// 	<< "start to do alignInferJunctionHashInfo->generateInferedUnfixedHeadPath_alignInferHash ..." << endl;

			bool tmpInferedPathFoundBool 
				= alignInferJunctionHashInfo->generateInferedUnfixedHeadPath_alignInferHash(
					tmpInferedPathJumpCodeVec,
					tmpInferedPathMismatchPosVec,
					tmpInferedPathMismatchCharVec,
					indexInfo, tmpChromNameInt,
					unfixedHeadLength, readSeqWithDirection,
					tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
					tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);

			//cout << "********************\ntmpInferedPathFoundBool: " << tmpInferedPathFoundBool << endl;
			if(tmpInferedPathFoundBool)
			{
				inferedPathGeneratedBool = true;
				
				//cout << "tmpInferedPathJumpCodeVec: " << tmpInferedPathJumpCodeVec.size() << endl;
				// for(int tmp = 0; tmp < tmpInferedPathJumpCodeVec.size(); tmp++)
				// {
				// 	cout << "tmpJumpCodeType: " << (tmpInferedPathJumpCodeVec)[tmp].type;
				// 	cout << "  tmpJumpCodeLength: " << (tmpInferedPathJumpCodeVec)[tmp].len << endl;
				// }

				inferedPathJumpCodeVecVec.push_back(tmpInferedPathJumpCodeVec);
				inferedPathMismatchPosVecVec.push_back(tmpInferedPathMismatchPosVec);
				inferedPathMismatchCharVecVec.push_back(tmpInferedPathMismatchCharVec);
			}
		}

		//cout << "end of generateInferedPathVec .... in  unfixedHead.h" << endl;
		return inferedPathGeneratedBool;
	}

	
	bool SJsearchInSJhash_areaStringHash_withAlignInferJuncHash(SJhash_Info* SJinfo,
	 	AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize,
		vector< vector<Jump_Code> >& inferedPathJumpCodeVecVec,
		vector< vector<int> >& inferedPathMismatchPosVecVec,
		vector< vector<char> >& inferedPathMismatchCharVecVec)
	{
		//cout << "SJsearchInSJhash_areaStringHash_withAlignInferJuncHash starts for unfixed head " << endl;
		vector< pair<int, vector<int> > > candiSJspliceSiteVec_reverse; // <acceptorPos, vector<donerPos> >

		this->generateCandiSJspliceSiteVec(SJinfo, readSeqWithDirection, 
			indexInfo, areaSize, candiSJspliceSiteVec_reverse);

		bool SJfoundInSJhash = (candiSJspliceSiteVec_reverse.size() > 0);

		if(!SJfoundInSJhash)
		{
			return false;
		}
		
		vector< pair<int,int> > normalSJspliceSitePairVec; // <donerEndPos, acceptorStartPos>
		for(int tmp = 0; tmp < candiSJspliceSiteVec_reverse.size(); tmp++)
		{
			int tmpAcceptorStartPos = candiSJspliceSiteVec_reverse[tmp].first;
			for(int tmp2 = 0; tmp2 < (candiSJspliceSiteVec_reverse[tmp].second).size(); tmp2++)
			{
				int tmpDonerEndPos = (candiSJspliceSiteVec_reverse[tmp].second)[tmp2];
				normalSJspliceSitePairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				//cout << "tmpDonerEndPos: " << tmpDonerEndPos << "  tmpAcceptorStartPos: " << tmpAcceptorStartPos << endl;
			}
		}

		//vector< vector<Jump_Code> > inferedPathJumpCodeVecVec;
		//vector< vector<int> > inferedPathMismatchPosVec;
		//vector< vector<char> > inferedPathMismatchCharVec;
		//cout << "start to do generateInferedPathVec ...." << endl;
		bool inferedPathGeneratedBool = this->generateInferedPathVec(
			normalSJspliceSitePairVec, 
			inferedPathJumpCodeVecVec, inferedPathMismatchPosVecVec,
			inferedPathMismatchCharVecVec, alignInferJunctionHashInfo, 
			indexInfo, readSeqWithDirection);
		//cout << "inferedPathGeneratedBool: " << inferedPathGeneratedBool << endl;
		return inferedPathGeneratedBool;// SJfoundInSJhash;
	}


	void getTmpMismatchPosVec(FixDoubleAnchor_Match_Info* fixMatchInfo, vector<int>& targetMismatchPosVec)
	{
		for(int tmp = 0; tmp < fixMatchInfo->returnMismatchPosVecSize(); tmp++)
		{
			targetMismatchPosVec.push_back(fixMatchInfo->returnMismatchPosInRead(tmp));
		}
	}
	void getTmpMismatchCharVec(FixDoubleAnchor_Match_Info* fixMatchInfo, vector<char>& targetMismatchCharVec)
	{
		for(int tmp = 0; tmp < fixMatchInfo->returnMismatchCharVecSize(); tmp++)
		{
			targetMismatchCharVec.push_back(fixMatchInfo->returnMismatchChar(tmp));
		}
	}

	void filterSJposFromRemapping_candi(int mismatchNumInBufferSeq)
	{
		int currentBestSJposNO = 0;
		int currentBestSJposNO_SJdistance = (SJposFromRemappingVec_candi[0]).second;
		int currentBestSJposNO_mismatch = SJposFromRemappingVec_candi_mismatch[0];
		for(int tmpSJposNO = 1; tmpSJposNO < SJposFromRemappingVec_candi.size(); tmpSJposNO ++)
		{
			int tmpSJpos_mismatch = SJposFromRemappingVec_candi_mismatch[tmpSJposNO]; 
			int tmpSJpos_SJdistance = (SJposFromRemappingVec_candi[tmpSJposNO]).second;

			if(tmpSJpos_mismatch < currentBestSJposNO_mismatch)
			{
				currentBestSJposNO = tmpSJposNO;
				currentBestSJposNO_mismatch = tmpSJpos_mismatch;
				currentBestSJposNO_SJdistance = tmpSJpos_SJdistance;
			}
			else if(tmpSJpos_mismatch == currentBestSJposNO_mismatch)
			{
				if(tmpSJpos_SJdistance < currentBestSJposNO_SJdistance)
				{
					currentBestSJposNO = tmpSJposNO;
					currentBestSJposNO_mismatch = tmpSJpos_mismatch;
					currentBestSJposNO_SJdistance = tmpSJpos_SJdistance;
				}
				else
				{}
			}
			else
			{}
		}

		int SJposFromRemapping_first = (SJposFromRemappingVec_candi[currentBestSJposNO]).first;
		int SJposFromRemapping_second = (SJposFromRemappingVec_candi[currentBestSJposNO]).second;
		int SJposFromRemapping_mismatch = SJposFromRemappingVec_candi_mismatch[currentBestSJposNO];

		SJposFromRemappingVec.push_back(pair<int, int> (SJposFromRemapping_first, SJposFromRemapping_second) );
		SJposFromRemappingVec_mismatch.push_back(SJposFromRemapping_mismatch - mismatchNumInBufferSeq);
		//if(STORE_MISMATCH_POS)
		//{
			SJposFromRemappingVec_mismatchPosVec.push_back(SJposFromRemappingVec_candi_mismatchPosVec[currentBestSJposNO]);
			//if(STORE_MISMATCH_CHA)
			//{
				SJposFromRemappingVec_mismatchCharVec.push_back(SJposFromRemappingVec_candi_mismatchCharVec[currentBestSJposNO]);
			//}
		//}
	}

};

#endif