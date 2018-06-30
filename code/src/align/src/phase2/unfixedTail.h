// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef UNFIXEDTAIL_H
#define UNFIXEDTAIL_H

#include <string>
#include <string.h>
#include "../general/fixDoubleAnchorMatch_info.h"
//#include "splice_info.h"

using namespace std;

class Unfixed_Tail
{
//public:
private:
	string readName;
	string alignDirection; 
	unsigned int midPartMapPosInWholeGenome;
	string otherCigarString;
	int unfixedTailLength;
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

	vector<int> possiGTAGpos_mismatch; // only AG end checked
	vector<int> possiCTACpos_mismatch; // only AC end checked
	vector< vector<int> > possiGTAGpos_mismatchPos; // only AG end checked
	vector< vector<int> > possiCTACpos_mismatchPos; // only AC end checked
	vector< vector<char> > possiGTAGpos_mismatchChar; // only AG end checked
	vector< vector<char> > possiCTACpos_mismatchChar; // only AC end checked

	vector< pair<int,int> > GTAGsjPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< pair<int,int> > CTACsjPos; // sequence around the SJ has been checked 

	vector< int > GTAGsjPos_mismatch; // sequence around the SJ has been checked including short anchor and midpart
	vector< int > CTACsjPos_mismatch; // sequence around the SJ has been checked 
	vector< vector<int> > GTAGsjPos_mismatchPos; // sequence around the SJ has been checked including short anchor and midpart
	vector< vector<int> > CTACsjPos_mismatchPos; // sequence around the SJ has been checked 
	vector< vector<char> > GTAGsjPos_mismatchChar; // sequence around the SJ has been checked including short anchor and midpart
	vector< vector<char> > CTACsjPos_mismatchChar; // sequence around the SJ has been checked 

	vector< pair<int,int> > SJposFromRemappingVec; // // <posInRead, SJsize>
	vector< int > SJposFromRemappingVec_mismatch;
	vector< vector<int> > SJposFromRemappingVec_mismatchPosVec;
	vector< vector<char> > SJposFromRemappingVec_mismatchCharVec;

	vector< pair<int,int> > SJposFromRemappingVec_candi;
	vector< int > SJposFromRemappingVec_candi_mismatch;   
	vector< vector<int> > SJposFromRemappingVec_candi_mismatchPosVec;
	vector< vector<char> > SJposFromRemappingVec_candi_mismatchCharVec;

	int returnUnfixedTailLength()
	{
		return unfixedTailLength;
	}
	unsigned int returnMidPartMapPosInWholeGenome()
	{
		return midPartMapPosInWholeGenome;
	}
	int returnMidPartMapChrInt()
	{
		return midPartMapChrInt;
	}
	int returnMidPartMapPosInChr()
	{
		return midPartMapPosInChr;
	}

	Unfixed_Tail()
	{//readLength = 100;
	}

	Unfixed_Tail(int tail_length, int mapPos, string mapChrName)
	{
		unfixedTailLength = tail_length;
		midPartMapPosInChr = midPartMapPosInChr;
		midPartMapChrName = mapChrName;
	}

	~Unfixed_Tail()
	{}

	int countMismatchNumInBufferSeq(int buffer, const string& readSeqWithDirection, Index_Info* indexInfo)
	{
		int countNum = 0;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		int readLength = readSeqWithDirection.length();
		for(int tmp = 0; tmp < buffer + 1; tmp ++)
		{
			if(readSeqWithDirection.at(readLength - unfixedTailLength - buffer - 1 + tmp ) 
				//!= (indexInfo->chromStr[chromNameInt]).at(midPartMapPosInChr - buffer - 1 + tmp) 
				!= indexInfo->returnOneBaseCharInGenome(chromNameInt, (midPartMapPosInChr - buffer + tmp))
				)
			{
				countNum ++;
			}
			else
			{}
		}
		return countNum;
	}

	void getUnfixedTailInfoFromRecord(PE_Read_Info& readInfo, bool end1, Alignment_Info* alignInfo, Index_Info* indexInfo)
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
		//cout << "midpartMapChrPos: " << 
		midPartMapPosInChr = alignInfo->getEndMatchedPosInChr();
		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();

		int jumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-1].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-2].len;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			indexInfo->convertStringToInt(midPartMapChrName), midPartMapPosInChr);
	}

	void getUnfixedTailInfoFromRecordWithAlignInfoType(PE_Read_Info& readInfo, int alignInfoType, Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		//cout << "getUnfixedTailInfoFromRecordWithAlignInfoType  starts ...." << endl;
		//cout << "alignInfoType: " << alignInfoType << endl;
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
		//cout << "midpartMapChrPos: " << 
		midPartMapPosInChr = alignInfo->getEndMatchedPosInChr();
		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		//cout << "midPartMapChrInt: " << midPartMapChrInt << endl;
		//cout << "midPartMapPosInChr: " << midPartMapPosInChr << endl;

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();
		//cout << "otherCigarString: " << otherCigarString << endl;
		int jumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		//cout << "jumpCodeVecSize: " << jumpCodeVecSize << endl;
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-1].len;
		//cout << "unfixedTailLength: " << unfixedTailLength << endl;
		midPartLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-2].len;
		//cout << "midPartLength: " << midPartLength << endl;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			indexInfo->convertStringToInt(midPartMapChrName), midPartMapPosInChr);
	}

	void getUnfixedTailInfoFromRecordWithAlignInfoType_new(PE_Read_Info& readInfo, int alignInfoType, Alignment_Info* alignInfo, Index_Info* indexInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();

		midPartMapChrName = alignInfo->returnAlignChromName();
		//cout << "midpartMapChrPos: " << 
		midPartMapPosInChr = alignInfo->getEndMatchedPosInChr();
		midPartMapChrInt = indexInfo->convertStringToInt(midPartMapChrName);

		otherCigarString = alignInfo->otherJumpCodeVec2StrForTail();

		int jumpCodeVecSize = (alignInfo->cigarStringJumpCode).size();
		unfixedTailLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-1].len;
		midPartLength = (alignInfo->cigarStringJumpCode)[jumpCodeVecSize-2].len;
		midPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			indexInfo->convertStringToInt(midPartMapChrName), midPartMapPosInChr);
	}

	void getPossibleSJpos(const string& readSeqWithDirection, const string& chromSeq, Index_Info* indexInfo)
	{
		// cout << endl << "***********************************" << endl 
		// 	<< "start to getPossibleSJpos ... fix unfixed tails ...." << endl
		// 	<< "************************************" << endl;

		int bufferLength = SINGLE_ANCHOR_TARGETMAPPING_BUFFER;
		if(bufferLength > midPartLength)
		{
			bufferLength = midPartLength;
		}
		//cout << "readLength: " << readLength << endl;
		//cout << "bufferlength: " << bufferLength << endl;
		int readLength = readSeqWithDirection.length();
		string pendingReadSeq = readSeqWithDirection.substr(
			readLength - unfixedTailLength - bufferLength, unfixedTailLength + bufferLength);
		//string pendingChroSeq = chromSeq.substr(
		//	midPartMapPosInWholeGenome - bufferLength, unfixedTailLength + bufferLength);

		//cout << "midPartChrMapPos: " << midPartMapPosInChr << endl;
		int midPartMapChrNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "midPartMapChrNameInt: " << midPartMapChrNameInt << endl;

		//cout << "midPartMapPosInChr-bufferLength" << midPartMapPosInChr-bufferLength << endl;

		//cout << "midPartMapPosInChr + unfixedTailLength " <<midPartMapPosInChr + unfixedTailLength << endl;
		//cout << "(indexInfo->chromLength)[midPartMapChrNameInt]: " << (indexInfo->chromLength)[midPartMapChrNameInt] << endl; 
 		if((midPartMapPosInChr-bufferLength < 0)||
			(midPartMapPosInChr-bufferLength + unfixedTailLength + bufferLength > 
				//(indexInfo->chromLength)[midPartMapChrNameInt] 
				indexInfo->returnChromLength(midPartMapChrNameInt)
				))
		{
			return;
		}	

		//string pendingChroSeq = indexInfo->chromStr[midPartMapChrNameInt].substr(
		//							midPartMapPosInChr-bufferLength, unfixedTailLength + bufferLength);
		string pendingChroSeq= indexInfo->returnChromStrSubstr(midPartMapChrNameInt, midPartMapPosInChr-bufferLength+1, unfixedTailLength + bufferLength);
		//cout << "pendingReadSeq: " << pendingReadSeq << endl;
		//cout << "pendingChroSeq: " << pendingChroSeq << endl;
		const string SJendStrAG = "GT";
		const string SJendStrAC = "CT";	
		//search for "AG"
		int startSearchPos = 0;
		int foundPos = 0;

		int otherPartInReadLength = readLength - unfixedTailLength;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAG, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(foundPos > bufferLength)
			{
				/*
				size_t max_append_mismatch = (foundPos - bufferLength)/10 + 1;
				size_t mismatch_bits = 0;
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				*/
				int max_mismatch = (foundPos - bufferLength)/10;
				FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
				bool matchBool = fixMatchInfo->fixMatch(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
					pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
					max_mismatch, readLength - unfixedTailLength + 1);
				if(matchBool)
				{
					possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
					//possiGTAGpos_mismatch.push_back(mismatch_bits);
					possiGTAGpos_mismatch.push_back(fixMatchInfo->returnMismatchNum());
					//if(STORE_MISMATCH_POS)
					//{
						vector<int> tmpMismatchPosVec;
						this->getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
						possiGTAGpos_mismatchPos.push_back(tmpMismatchPosVec);
						//if(STORE_MISMATCH_CHA)
						//{
							vector<char> tmpMismatchCharVec;
							this-> getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
							possiGTAGpos_mismatchChar.push_back(tmpMismatchCharVec);
						//}
					//}
				}
				delete fixMatchInfo;
			}
			else
			{
				possiGTAGpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
				possiGTAGpos_mismatch.push_back(0);
				//if(STORE_MISMATCH_POS)
				//{
					vector<int> tmpMismatchPosVec;
					//this-> getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
					possiGTAGpos_mismatchPos.push_back(tmpMismatchPosVec);
					//if(STORE_MISMATCH_CHA)
					//{
						vector<char> tmpMismatchCharVec;
						//this-> getTmpMismatchCharVec(fixMatchInfo, tmpMismatchCharVec);
						possiGTAGpos_mismatchChar.push_back(tmpMismatchCharVec);
					//}
				//}
			}
			startSearchPos = foundPos + 1;
		}

		//search for "AC"
		startSearchPos = 0;
		foundPos = 0;
		for(int tmp = 0; ; tmp++)
		{
			foundPos = pendingChroSeq.find(SJendStrAC, startSearchPos);
			if(foundPos == pendingChroSeq.npos)
				break;
			if(foundPos > bufferLength)
			{
				/*
				size_t max_append_mismatch = (foundPos - bufferLength)/10 + 1;
				size_t mismatch_bits = 0;
				size_t comb_bits = 0;
				bool matchBool = score_string(pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
												pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
												max_append_mismatch, mismatch_bits, comb_bits);//append first
				*/
				int max_mismatch = (foundPos - bufferLength)/10;
				FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
				bool matchBool = fixMatchInfo->fixMatch(
					pendingReadSeq.substr(bufferLength, foundPos - bufferLength),
					pendingChroSeq.substr(bufferLength, foundPos - bufferLength),
					max_mismatch, readLength - unfixedTailLength + 1);
				if(matchBool)
				{
					possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
					//possiCTACpos_mismatch.push_back(mismatch_bits);
					possiCTACpos_mismatch.push_back(fixMatchInfo->returnMismatchNum());
					//if(STORE_MISMATCH_POS)
					//{
						vector<int> tmpMismatchPosVec;
						this->getTmpMismatchPosVec(fixMatchInfo, tmpMismatchPosVec);
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
				possiCTACpos.push_back(foundPos + otherPartInReadLength - bufferLength + 1);
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
			startSearchPos = foundPos + 1;
		}

		// cout << "output GTAG positions: " << endl;
		// for(int tmp = 0; tmp < possiGTAGpos.size(); tmp++)
		// {
		// 	cout << possiGTAGpos[tmp] << endl;
		// }
		// cout << "output CTAC positions: " << endl;
		// for(int tmp = 0; tmp < possiGTAGpos.size(); tmp++)
		// {
		// 	cout << possiCTACpos[tmp] << endl;
		// }		
		// cout << endl << "***********************************" << endl 
		// 	<< "end of getPossibleSJpos, fix unfixed tails ...." << endl
		// 	<< "************************************" << endl;		
	}

	bool SJsearchInSJhash_areaStringHash(SJhash_Info* SJinfo, 
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize)
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();

		int buffer = 4;

		if(buffer > midPartLength - 1)
			buffer = midPartLength - 1;

		/*if(readLength - unfixedTailLength - buffer - 1 < 0)
		{
			buffer = readLength - unfixedTailLength - 1;
		}*/
		//cout << "unfixedTailLength: " << unfixedTailLength << endl;
		//cout << "startLocInRead for readPendingStr: " << readLength - unfixedTailLength - buffer << endl;
		//cout << "seqLength for readPendingStr: " << unfixedTailLength + buffer + 1 << endl;
		string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
		//cout << "readPendingStr: " << endl << readPendingStr << endl;
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		//cout << "chromNameInt: " << chromNameInt << endl;

		int areaNOmin = (int)((midPartMapPosInChr-buffer)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr+unfixedTailLength-1)/areaSize);

		//cout << "areaNOmin: " << areaNOmin << endl;
		//cout << "areaNOmax: " << areaNOmax << endl;

		vector<int> SJdonerSiteVec;
		// search SJdonerSite in Hash
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			//cout << "tmpArea: " << tmpArea << endl;
			SJareaHashIter tmpSJareaHashIter
				= ((SJinfo->SJstartPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJstartPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJdonerPos = (*intSetIter);
					//cout << "tmpSJdonerPos: " << tmpSJdonerPos << endl;
					if( (tmpSJdonerPos >= (midPartMapPosInChr-buffer)) && (tmpSJdonerPos <= (midPartMapPosInChr+unfixedTailLength-1)) )
					{
						//cout << "push_back tmpSJdonerPos" << endl;
						SJdonerSiteVec.push_back(tmpSJdonerPos);
					}
				}
			}
			else
			{
			}
		}

		//cout << "SJdonerSiteVec.size(): " << SJdonerSiteVec.size() << endl;

		if(SJdonerSiteVec.size() > 0)
		{
			//string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
			//cout << "readPendingStr: " << endl << readPendingStr << endl;
			//search anchorString in Hash
			for(int tmp = 0; tmp < SJdonerSiteVec.size(); tmp ++)
			{
				//cout << "tmpIndex in SJdonerSiteVec: " << tmp << endl;
				vector<int> tmpAcceptorSiteVec;
				int tmpSJdonerSite = SJdonerSiteVec[tmp];
				int tmpTailLength = unfixedTailLength + midPartMapPosInChr - tmpSJdonerSite;
				//cout << "tmpSJdonerSite: " << tmpSJdonerSite << endl;
				//cout << "tmpTailLength: " << tmpTailLength << endl;
				SplicePosHashIter tmpPosHashIter
					= (SJinfo->spliceJunctionNormal)[chromNameInt].find(tmpSJdonerSite);
				if(tmpPosHashIter != (SJinfo->spliceJunctionNormal)[chromNameInt].end())
				{
					if(tmpTailLength >= (SJinfo->anchorStringLength))
					{
						//cout << "startLocInRead for tmpAnchorString: " << readLength - tmpTailLength << endl;
						string tmpAnchorString = readSeqWithDirection.substr(readLength - tmpTailLength, (SJinfo->anchorStringLength));
						//cout << "tmpAnchorString: " << tmpAnchorString << endl; 
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
							{
								//cout << "tmpAnchorString found" << endl;
								//cout << "push_back tmpAcceptorSite: " << (*tmpIntSetIter) << endl;
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							//cout << "tmpAnchorString not found" << endl;
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
								{
									//cout << "tmpAnchor == * found" << endl;
									//cout << "push_back tmpAcceptorSite: " << (*tmpIntSetIter) << endl;
									tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
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
						//cout << "tmpTailLength too short" << endl;
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++ )
							{
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
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
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJdonerSite should be found in hash ! " << endl;

				}

				// check all possible short anchor SJ candidates 
				for(int tmpVecNO = 0; tmpVecNO < tmpAcceptorSiteVec.size(); tmpVecNO ++)
				{
					int tmpSJdonerEndPosInRead = readLength - tmpTailLength;
					int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;

					int tmpSJdonerEndPosInChr = tmpSJdonerSite;
					int tmpSJacceptorStartPosInChr = tmpAcceptorSiteVec[tmpVecNO];

					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					//size_t max_append_mismatch;
					//size_t mismatch_bits;
					//size_t comb_bits;
					bool matchBool;
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					if(
						(midPartMapPosInChr + tmpSJdonerEndPosInRead - readLength + unfixedTailLength 
							<= (indexInfo->returnChromLength(midPartMapChrInt) - 1))
						&&
						(tmpSJacceptorStartPosInChr  + readLength - tmpSJacceptorStartPosInRead 
							<= (indexInfo->returnChromLength(midPartMapChrInt) - 1))
						)
					{
						//chromDonerEndStr = (indexInfo->chromStr[chromNameInt]).substr(midPartMapPosInChr - buffer - 1,
						//	tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);
						chromDonerEndStr = indexInfo->returnChromStrSubstr(chromNameInt, midPartMapPosInChr - buffer,
							tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);

						//chromAcceptorStartStr = (indexInfo->chromStr[chromNameInt]).substr(
						//	tmpSJacceptorStartPosInChr - 1, readLength - tmpSJacceptorStartPosInRead + 1);
						chromAcceptorStartStr = indexInfo->returnChromStrSubstr(chromNameInt,
							tmpSJacceptorStartPosInChr, readLength - tmpSJacceptorStartPosInRead + 1);
						//cout << "tmpChromAcceptorStartStr: " << chromAcceptorStartStr << endl;
						//cout << "chromNameInt: " << chromNameInt << endl;
						//cout << "tmpSJacceptorStartPosInChr: " << tmpSJacceptorStartPosInChr << endl;
						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;
						int max_mismatch = (unfixedTailLength)/MATCH_BASE_PER_MISMATCH_BASE;
						matchBool = fixMatchInfo->fixMatch(readPendingStr, chromPendingStr, max_mismatch, readLength - unfixedTailLength - buffer);
					}
					else
					{
						matchBool = false;
					}

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

		//cout << "SJfoundInSJhash: " << SJfoundInSJhash << endl;

		if(SJfoundInSJhash)
		{
			//int mismatchNumInBufferSeq = this->countMismatchNumInBufferSeq(buffer, readSeqWithDirection, indexInfo);
			int mismatchNumInBufferSeq = 0;
			this->filterSJposFromRemapping_candi(mismatchNumInBufferSeq);
		}
		
		return SJfoundInSJhash;
	}


	void generateCandiSJspliceSiteVec(
		SJhash_Info* SJinfo,
		const string& readSeqWithDirection, 
		Index_Info* indexInfo, int areaSize, 
		vector< pair<int, vector<int> > >& candiSJspliceSiteVec)
	{
		int readLength = readSeqWithDirection.length();

		int buffer = 4;

		if(buffer > midPartLength - 1)
			buffer = midPartLength - 1;

		/*if(readLength - unfixedTailLength - buffer - 1 < 0)
		{
			buffer = readLength - unfixedTailLength - 1;
		}*/

		string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
		int chromNameInt = indexInfo->convertStringToInt(midPartMapChrName);

		int areaNOmin = (int)((midPartMapPosInChr-buffer)/areaSize);
		int areaNOmax = (int)((midPartMapPosInChr+unfixedTailLength-1)/areaSize);

		vector<int> SJdonerSiteVec;

		
		// search SJdonerSite in Hash
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SJareaHashIter tmpSJareaHashIter
				= ((SJinfo->SJstartPosAreaHash)[chromNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJinfo->SJstartPosAreaHash)[chromNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJdonerPos = (*intSetIter);
					if( (tmpSJdonerPos >= (midPartMapPosInChr-buffer)) && (tmpSJdonerPos <= (midPartMapPosInChr+unfixedTailLength-1)) )
					{
						SJdonerSiteVec.push_back(tmpSJdonerPos);
					}
				}
			}
			else
			{}
		}

		if(SJdonerSiteVec.size() > 0)
		{
			string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
			//search anchorString in Hash
			for(int tmp = 0; tmp < SJdonerSiteVec.size(); tmp ++)
			{
				vector<int> tmpAcceptorSiteVec;

				int tmpSJdonerSite = SJdonerSiteVec[tmp];
				int tmpTailLength = unfixedTailLength + midPartMapPosInChr - tmpSJdonerSite;

				SplicePosHashIter tmpPosHashIter
					= (SJinfo->spliceJunctionNormal)[chromNameInt].find(tmpSJdonerSite);

				if(tmpPosHashIter != (SJinfo->spliceJunctionNormal)[chromNameInt].end())
				{
					if(tmpTailLength >= (SJinfo->anchorStringLength))
					{
						string tmpAnchorString = readSeqWithDirection.substr(readLength - tmpTailLength, (SJinfo->anchorStringLength));
					
						SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find(tmpAnchorString);
						if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
						{
							for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
							{
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
							}
						}
						else
						{
							tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
							if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
							{
								for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
									tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
								{
									tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
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
								tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter ++ )
							{
								tmpAcceptorSiteVec.push_back(*tmpIntSetIter);
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
					cout << "error in SJsearchInAreaAndStringHash ! tmpSJdonerSite should be found in hash ! " << endl;

				}

				if(tmpAcceptorSiteVec.size() > 0)
				{
					candiSJspliceSiteVec.push_back(pair<int, vector<int> >(
						tmpSJdonerSite, tmpAcceptorSiteVec));
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
		int tmpReadLengthWithDirection = readSeqWithDirection.length();

		bool inferedPathGeneratedBool = false;
		int tmpChromNameInt = indexInfo->convertStringToInt(midPartMapChrName);
		for(int tmp = 0; tmp < spliceSitePairVec.size(); tmp++)
		{
			vector<Jump_Code> tmpInferedPathJumpCodeVec;
			vector<int> tmpInferedPathMismatchPosVec;
			vector<char> tmpInferedPathMismatchCharVec;
			int tmpSpliceDonerEndPosInChr = spliceSitePairVec[tmp].first;
			int tmpSpliceAcceptorStartPosInChr = spliceSitePairVec[tmp].second;
			int tmpUnfixedTailLength 
				= unfixedTailLength - (tmpSpliceDonerEndPosInChr - midPartMapPosInChr);
			int tmpSpliceDonerEndPosInRead = tmpReadLengthWithDirection - tmpUnfixedTailLength;
			int tmpSpliceAcceptorStartPosInRead = tmpSpliceDonerEndPosInRead + 1;
			// cout << "tmpSpliceDonerEndPosInChr: " << tmpSpliceDonerEndPosInChr << endl;
			// cout << "tmpSpliceAcceptorStartPosInChr: " << tmpSpliceAcceptorStartPosInChr << endl;
			// cout << "tmpUnfixedTailLength: " << tmpUnfixedTailLength << endl;
			// cout << "tmpSpliceDonerEndPosInRead: " << tmpSpliceDonerEndPosInRead << endl;
			// cout << "tmpSpliceAcceptorStartPosInRead: " << tmpSpliceAcceptorStartPosInRead << endl;

			bool tmpInferedPathFoundBool 
				= alignInferJunctionHashInfo->generateInferedUnfixedTailPath_alignInferHash(
					tmpInferedPathJumpCodeVec,
					tmpInferedPathMismatchPosVec,
					tmpInferedPathMismatchCharVec,
					indexInfo, tmpChromNameInt,
					unfixedTailLength, readSeqWithDirection,
					tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
					tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
			
			//cout << "tmpInferedPathFoundBool: " << tmpInferedPathFoundBool << endl;
			if(tmpInferedPathFoundBool)
			{
				inferedPathGeneratedBool = true;
				inferedPathJumpCodeVecVec.push_back(tmpInferedPathJumpCodeVec);
				inferedPathMismatchPosVecVec.push_back(tmpInferedPathMismatchPosVec);
				inferedPathMismatchCharVecVec.push_back(tmpInferedPathMismatchCharVec);
			}
		}
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
		vector< pair< int, vector<int> > > candiSJspliceSiteVec; // <donerPos, vector<acceptorPos> >

		this->generateCandiSJspliceSiteVec(SJinfo, readSeqWithDirection, 
			indexInfo, areaSize, candiSJspliceSiteVec);
		//cout << "candiSJspliceSiteVec: " << candiSJspliceSiteVec.size() << endl;
		bool SJfoundInSJhash = (candiSJspliceSiteVec.size() > 0);

		if(!SJfoundInSJhash)
		{
			return false;
		}
		
		vector< pair<int,int> > normalSJspliceSitePairVec; // <donerEndPos, acceptorStartPos>
		for(int tmp = 0; tmp < candiSJspliceSiteVec.size(); tmp++)
		{
			int tmpDonerEndPos = candiSJspliceSiteVec[tmp].first;
			for(int tmp2 = 0; tmp2 < (candiSJspliceSiteVec[tmp].second).size(); tmp2++)
			{
				int tmpAcceptorStartPos = (candiSJspliceSiteVec[tmp].second)[tmp2];
				//cout << "tmpDonerEndPos: " << tmpDonerEndPos << "  tmpAcceptorStartPos: " << tmpAcceptorStartPos << endl; 
				normalSJspliceSitePairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
			}
		}

		//vector< vector<Jump_Code> > inferedPathJumpCodeVecVec;
		//vector< vector<int> > inferedPathMismatchPosVec;
		//vector< vector<char> > inferedPathMismatchCharVec;

		bool inferedPathGeneratedBool = this->generateInferedPathVec(
			normalSJspliceSitePairVec, 
			inferedPathJumpCodeVecVec, inferedPathMismatchPosVecVec,
			inferedPathMismatchCharVecVec, alignInferJunctionHashInfo, 
			indexInfo, readSeqWithDirection);

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
	void filterSJposFromRemapping_candi(int mismatchNumInBufferSeq)//, Index_Info* indexInfo)
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
	/*
	bool SJsearchInSJhash(SJhash_Info* SJinfo, const string& readSeqWithDirection, Index_Info* indexInfo)
	//without use of areaHash
	{
		bool SJfoundInSJhash = false;

		int readLength = readSeqWithDirection.length();

		int buffer = 4;

		if(readLength - unfixedTailLength - buffer - 1 < 0)
		{
			buffer = readLength - unfixedTailLength - 1;
		}

		string readPendingStr = readSeqWithDirection.substr(readLength - unfixedTailLength - buffer - 1, unfixedTailLength + buffer + 1);
		
		for(int tmp = -buffer; tmp < unfixedTailLength; tmp++)
		{
			int tmpSJdonerEndPosInRead = readLength - unfixedTailLength + tmp;
			int tmpSJacceptorStartPosInRead = tmpSJdonerEndPosInRead + 1;
			int tmpSJdonerEndPosInChr = midPartMapPosInChr + tmp;

			if(tmpSJdonerEndPosInChr > (indexInfo->returnChromLength(midPartMapChrInt) - 2))
			{
				continue;
			}

			if( ((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr)
				== ((SJinfo->SJintHashNormal)[midPartMapChrInt]).end() )
			{}
			else
			{
				for(set<int>::iterator tmp2 = ((((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr))->second).begin(); 
					tmp2 != ((((SJinfo->SJintHashNormal)[midPartMapChrInt]).find(tmpSJdonerEndPosInChr))->second).end(); 
					tmp2++)
				{
					int tmpSJacceptorStartPosInChr 
						= *tmp2;
					//string readPendingStr = readSeqWithDirection.substr(tmpSJdonerEndPosInRead + buffer - 1, readLength - (tmpSJdonerEndPosInRead + buffer) + 1)
					string chromDonerEndStr;
					string chromAcceptorStartStr;
					string chromPendingStr;
					size_t max_append_mismatch;
					size_t mismatch_bits;
					bool matchBool;

					if((midPartMapPosInChr + tmpSJdonerEndPosInRead - readLength + unfixedTailLength > 
						//((indexInfo->chromLength)[midPartMapChrInt] - 1)
						(indexInfo->returnChromLength(midPartMapChrInt) - 1)
						)
						||
						(tmpSJacceptorStartPosInChr  + readLength - tmpSJacceptorStartPosInRead > 
							(
								//(indexInfo->chromLength)[midPartMapChrInt] - 1
								(indexInfo->returnChromLength(midPartMapChrInt))- 1
								)))
					{
						matchBool = false;
					}	
					else
					{
						//chromDonerEndStr = indexInfo->chromStr[midPartMapChrInt].substr(
						//				midPartMapPosInChr - buffer - 1, tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);
						chromDonerEndStr = indexInfo->returnChromStrSubstr(midPartMapChrInt, 
							midPartMapPosInChr - buffer, tmpSJdonerEndPosInRead - readLength + unfixedTailLength + buffer + 1);
						
						//chromAcceptorStartStr = indexInfo->chromStr[midPartMapChrInt].substr(
						//				tmpSJacceptorStartPosInChr - 1, readLength - tmpSJacceptorStartPosInRead + 1);
						chromAcceptorStartStr = indexInfo->returnChromStrSubstr(midPartMapChrInt,
							tmpSJacceptorStartPosInChr, readLength - tmpSJacceptorStartPosInRead + 1);
						chromPendingStr = chromDonerEndStr + chromAcceptorStartStr;

						
						max_append_mismatch = (unfixedTailLength + buffer + 1)/10 + 1;
						
						mismatch_bits = 0;
						
						matchBool = score_string(readPendingStr, chromPendingStr,
													max_append_mismatch, mismatch_bits);//append first
						//cout << "readPendingStr: "  << endl << readPendingStr << endl; 
						//cout << "chroPendingStr: "  << endl << chromPendingStr << endl;
					}
					
					if(matchBool)
					{
						SJfoundInSJhash = true;
						//cout << "SJ found: at " << tmpSJdonerEndPosInRead << " SJsize: " << tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1 << endl;
						SJposFromRemappingVec.push_back(pair <int, int > (tmpSJdonerEndPosInRead, tmpSJacceptorStartPosInChr - tmpSJdonerEndPosInChr - 1));
					}
				} 
			}

		}


		return SJfoundInSJhash;
	}*/
};

#endif