// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNEVENT_INFO_H
#define ALIGNEVENT_INFO_H

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

#include "splice_info.h"

using namespace std;

class AlignEvent_Info
{
private:
	string mapChrStr;
	int mapChrPos;

	vector< Jump_Code > cigarStringJumpCode;  // <jumpCode, segment_startPosInChr> 
	// for insertion, segmentStartInChr is the next seg's startMapPos-1
	vector< int > jumpCodeIndexInSegVec; 
	// if jumpCode.type == M or I or S, seg_index starts from 0, else seg_index = -1;

	vector< pair<int, int> > segVec; 
	// < <segmentStartInRead, segmentStartInChr>, segment_length >, 
	// for insertion, segmentStartInChr is the next seg's startMapPos-1
	vector< int > segIndexInCigarStringJumpCode;

	vector<int> mismatchPosInRead;
	vector<char> mismatchPosCharInRead;
	vector<char> mismatchPosQualityInRead;

public:
	AlignEvent_Info()
	{}

	void generateSJposVec(vector< pair<int,int> >& SJposVec, 
		vector< vector<Jump_Code> >& SJfollowingJumpCodeVec,
		vector< vector< pair<int, pair<char,char> > > >& SJfollowingJumpCodeMismatchPosVec,
		bool recordFollowingJumpCode_bool, bool recordFollowingSegsMismatch_bool, bool fasta_fastq_bool)
	{
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCode[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				if((cigarStringJumpCode[tmp-1].type == "M") && (cigarStringJumpCode[tmp+1].type == "M"))
				{
					int lastJumpCodeIndexInSegVec = jumpCodeIndexInSegVec[tmp-1];
					int nextJumpCodeIndexInSegVec = jumpCodeIndexInSegVec[tmp+1];
					int donerEndPos = segVec[lastJumpCodeIndexInSegVec].second + cigarStringJumpCode[tmp-1].len - 1;
					int acceptorStartPos = segVec[nextJumpCodeIndexInSegVec].second;
					SJposVec.push_back(pair<int,int> (donerEndPos, acceptorStartPos));
					if(recordFollowingJumpCode_bool)
					{	
						vector<Jump_Code> newTmpJumpCodeVec;
						for(int tmp2 = tmp+1; tmp2 < cigarStringJumpCode.size(); tmp2++)
						{
							newTmpJumpCodeVec.push_back(cigarStringJumpCode[tmp2]);
						}
						SJfollowingJumpCodeVec.push_back(newTmpJumpCodeVec);
					}	
					if(recordFollowingSegsMismatch_bool)
					{
						int nextSegStartLocInRead = segVec[nextJumpCodeIndexInSegVec].first;
						vector< pair<int, pair<char,char> > > newTmpMismatchPosCharQualityVec;
						for(int tmpMismatchPosIndex = 0; tmpMismatchPosIndex < mismatchPosInRead.size(); tmpMismatchPosIndex++)
						{
							if(mismatchPosInRead[tmpMismatchPosIndex] >= nextSegStartLocInRead)
							{
								int tmpMismatchLocInRead = mismatchPosInRead[tmpMismatchPosIndex];
								char tmpMismatchCharInRead = mismatchPosCharInRead[tmpMismatchPosIndex];
								char tmpMismatchQualityInRead = mismatchPosQualityInRead[tmpMismatchPosIndex];
								newTmpMismatchPosCharQualityVec.push_back(pair<int, pair<char,char> >(tmpMismatchLocInRead,
									pair<char,char>(tmpMismatchCharInRead, tmpMismatchQualityInRead)));
							}
						}
						SJfollowingJumpCodeMismatchPosVec.push_back(newTmpMismatchPosCharQualityVec);
					}
				}
				else
				{
					cout << "Error: jumpCode next to N is not M !" << endl;
					exit(1);
				}
			}
		}
	}

	string returnMapChrStr()
	{
		return mapChrStr;
	}

	string returnCigarStringJumpCodeStr()
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			//cout << "jumpCodeIndexInSegVec.size(): " << jumpCodeIndexInSegVec.size() << endl;
			//cout << "jumpCodeIndexInSegVec[tmp]: " << jumpCodeIndexInSegVec[tmp] << endl;
			str = str + int_to_str(cigarStringJumpCode[tmp].len) 
				+ cigarStringJumpCode[tmp].type + ",";
		}
		return str;
	}

	string returnSegVecStr()
	{
		string str;
		for(int tmp = 0; tmp < segVec.size(); tmp++)
		{
			str = str + int_to_str(segVec[tmp].first) + ":" + int_to_str(segVec[tmp].second)
				+ "-" + int_to_str(segIndexInCigarStringJumpCode[tmp]) + ",";
		}
		return str;
	}

	void readAlignInfoFromSAM(const string& samStr, bool recordFollowingSegsMismatch_bool,
		bool fasta_fastq_bool)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 12; tmp++)
		{	
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}	
		//cout << "readName: " << samFieldVec[0] << endl;
		// mapChrStr
		mapChrStr = samFieldVec[2];
		//cout << "mapChrStr: " << mapChrStr << endl;
		// mapChrPos
		string mapChrPosStr = samFieldVec[3];
		mapChrPos = atoi(mapChrPosStr.c_str());
		//cout << "mapChrPos: " << mapChrPos << endl;
		// cigarStringJumpCodeVec
		string cigarString = samFieldVec[5];
		this->cigarString2jumpCodeVec(cigarString);
		//cout << "cigarStr: " << this->returnCigarStringJumpCodeStr() << endl;
		// exonVec
		this->jumpCode2exonVec();
		//cout << "segVecStr: " << this->returnSegVecStr() << endl << endl;
		// mismatchPosInRead		 
		
		if(recordFollowingSegsMismatch_bool)
		{	
			string readSeq = samFieldVec[9];
			string qualitySeq = samFieldVec[10]; 
			if(samFieldVec.size() >= 12)
			{
				string flagStr = samFieldVec[11];
				this->getMismatchPosInRead(flagStr, fasta_fastq_bool, readSeq, qualitySeq);
			}
			else
			{}
		}
	}

	void cigarString2jumpCodeVec(const string& jumpCodeStr)
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

	void jumpCode2exonVec()
	{
		int tmpSegIndex = 0;
		int tmpPosInChr = mapChrPos;
		int tmpPosInRead = 1;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			int tmpJumpCodeLength = cigarStringJumpCode[tmp].len;
			string tmpJumpCodeType = cigarStringJumpCode[tmp].type;
			if(tmpJumpCodeType == "S")
			{
				// update jumpCodeIndexInSegVec
				jumpCodeIndexInSegVec.push_back(tmpSegIndex);
				tmpSegIndex ++;
				// update segVec;
				int segmentStartInRead = tmpPosInRead;
				int segmentStartInChr = -1;
				segVec.push_back(pair<int,int>(segmentStartInRead, segmentStartInChr));
				// update segIndexInCigarStringJumpCode
				segIndexInCigarStringJumpCode.push_back(tmp);


				// update tmpPosInChr, tmpPosInRead
				tmpPosInRead += tmpJumpCodeLength;
			}
			else if (tmpJumpCodeType == "M")
			{
				// update jumpCodeIndexInSegVec
				jumpCodeIndexInSegVec.push_back(tmpSegIndex);
				tmpSegIndex ++;
				// update segVec;
				int segmentStartInRead = tmpPosInRead;
				int segmentStartInChr = tmpPosInChr;
				segVec.push_back(pair<int,int>(segmentStartInRead, segmentStartInChr));
				// update segIndexInCigarStringJumpCode
				segIndexInCigarStringJumpCode.push_back(tmp);		


				// update tmpPosInChr, tmpPosInRead
				tmpPosInChr += tmpJumpCodeLength;	
				tmpPosInRead += tmpJumpCodeLength;
			}
			else if (tmpJumpCodeType == "I")
			{
				// update jumpCodeIndexInSegVec
				jumpCodeIndexInSegVec.push_back(tmpSegIndex);
				tmpSegIndex ++;
				// update segVec;
				int segmentStartInRead = tmpPosInRead;
				int segmentStartInChr = -1;
				segVec.push_back(pair<int,int>(segmentStartInRead, segmentStartInChr));
				// update segIndexInCigarStringJumpCode
				segIndexInCigarStringJumpCode.push_back(tmp);	


				// update tmpPosInChr, tmpPosInRead		
				tmpPosInRead += tmpJumpCodeLength;	
			}
			else if (tmpJumpCodeType == "D")
			{
				// update jumpCodeIndexInSegVec
				jumpCodeIndexInSegVec.push_back(-1);
				// update segVec;
				// update segIndexInCigarStringJumpCode	


				// update tmpPosInChr, tmpPosInRead
				tmpPosInChr += tmpJumpCodeLength;
			}
			else if (tmpJumpCodeType == "N")
			{
				// update jumpCodeIndexInSegVec
				jumpCodeIndexInSegVec.push_back(-1);
				// update segVec;
				// update segIndexInCigarStringJumpCode	


				// update tmpPosInChr, tmpPosInRead
				tmpPosInChr += tmpJumpCodeLength;
			}						
			else
			{
				cout << "error! in jumpCode2exonVec (AlignEvent_Info)" << endl;
				exit(1);
			}
		}
	}

	void getMismatchPosInRead(const string& flagStr, bool fasta_fastq_bool,
		const string& readSeq, const string& qualitySeq)
	{
		int MDloc = flagStr.find("MD",0);
		int notIntLoc;
		if(MDloc == flagStr.npos)
		{}
		else
		{
			notIntLoc = flagStr.find_first_not_of("0123456789", MDloc+5);
		}
		if(notIntLoc != string::npos)
		{	
			string MDstr = flagStr.substr(MDloc+5, notIntLoc-MDloc-5);
			int startSearchPos = 0;
			while(1)
			{	
				int commaLoc = MDstr.find(",", startSearchPos);
				if(commaLoc == MDstr.npos)
				{
					return;
				}
				else
				{
					string tmpMismatchPosStr = MDstr.substr(startSearchPos, commaLoc - startSearchPos);
					int tmpMismatchPos = atoi(tmpMismatchPosStr.c_str());
					mismatchPosInRead.push_back(tmpMismatchPos);
					mismatchPosCharInRead.push_back(readSeq.at(tmpMismatchPos-1));
					if(!fasta_fastq_bool)
					{
						mismatchPosQualityInRead.push_back(qualitySeq.at(tmpMismatchPos-1));
					}
					else
					{
						mismatchPosQualityInRead.push_back('*');
					}
				}
				startSearchPos = commaLoc+1;
			}
		}
	}
};

#endif