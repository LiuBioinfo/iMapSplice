// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//#include "splice_info.h"
//#define SPLICE_MIN_LENGTH 50
//#define FIX_SPLICE_BUFFER 5

/*    clock_t read_file_begin, read_file_end, read_file_end2, read_file_end3, align_begin, align_end, align_cost = 0, first_align_begin, first_align_end, first_align_cost = 0,
    	overall_begin, overall_end, overall_cost = 0, other_align_begin, other_align_end, other_align_cost = 0,
    	copy_readString_begin, copy_readString_end, copy_readString_cost = 0,
    	copy_readStringRC_begin, copy_readStringRC_end, copy_readStringRC_cost = 0,
    	map_main_begin, map_main_end, map_main_cost = 0, det_begin, det_end, det_cost = 0,
    	getRcRead_begin, getRcRead_end, getRcRead_cost = 0, read_process_start, read_process_end, read_process_cost = 0, 
    	load_junction_begin, load_junction_end, remapping_begin, remapping_end,
    	searchSplice_begin, searchSplice_end, searchSplice_cost = 0,
    	checkSplice_begin, checkSplice_end, checkSplice_cost = 0,
    	getHead_begin, getHead_end, getHead_cost = 0,
    	searchSpliceOther_begin, searchSpliceOther_end, searchSpliceOther_cost = 0,
    	search_pos_begin, search_pos_end, search_pos_cost = 0,
    	search_char_begin, search_char_end, search_char_cost = 0;
*/
//#include "constantDefinitions.h"

string chromStr[30];
//string readString;
string chromString;

bool shortHeadFixedBool = false;
bool shortTailFixedBool = false;


typedef map<char, vector<int> > SpliceEndStringHash;
typedef map<int, SpliceEndStringHash> SpliceStartPosHash;

typedef SpliceStartPosHash::iterator SpliceStartPosHashIter;
typedef SpliceEndStringHash::iterator SpliceEndStringHashIter;

SpliceStartPosHash spliceJunctionNormal[22];
SpliceStartPosHash spliceJunctionReverse[22];

SpliceStartPosHashIter startPosIter;
SpliceEndStringHashIter endStringIter;

int deletionNum = 0;
int junctionNum = 0;

bool checkNoSoftClipping(const string& cigarString)
{
	return (cigarString.find("S") == cigarString.npos);
}

int getMin(int a, int b)
{
	if(a < b)
		return a;
	else
		return b;
}

int extendBackInGenome(unsigned int chrMapPos, int locInRead, const string& readString)// number of bases can be extended back
{
	int tmp;
	for(tmp = locInRead-1; tmp > 0; tmp--)
	{
		if(chromString[chrMapPos - locInRead + tmp -1] != readString[tmp-1])
			break;
	}
	return locInRead-tmp-1;
}

bool checkShortHeadSoftClipping(const string& cigarString)
{
	return (cigarString.substr(1,1) == "S");
}

bool checkShortTailSoftClipping(const string& cigarString)
{
	int cigarStringSize = cigarString.size();

	return (cigarString.substr(cigarStringSize-1, 1) == "S");
}


/*
int covertStringToInt(const string& chrName)
{
	if(chrName == "chr1")
		return 0;
	else if(chrName == "chr2")
		return 1;
	else if(chrName == "chr3")
		return 2;
	else if(chrName == "chr4")
		return 3;
	else if(chrName == "chr5")
		return 4;
	else if(chrName == "chr6")
		return 5;
	else if(chrName == "chr7")
		return 6;
	else if(chrName == "chr8")
		return 7;
	else if(chrName == "chr9")
		return 8;
	else if(chrName == "chr10")
		return 9;
	else if(chrName == "chr11")
		return 10;
	else if(chrName == "chr12")
		return 11;
	else if(chrName == "chr13")
		return 12;
	else if(chrName == "chr14")
		return 13;
	else if(chrName == "chr15")
		return 14;
	else if(chrName == "chr16")
		return 15;
	else if(chrName == "chr17")
		return 16;
	else if(chrName == "chr18")
		return 17;
	else if(chrName == "chr19")
		return 18;
	else if(chrName == "chrX")
		return 19;
	else if(chrName == "chrY")
		return 20;
	else if(chrName == "chrM")
		return 21;
}
*/


bool insertSpliceJunction2ReverseHashForHead(SpliceStartPosHash* spliceJunction,
	int chrInt, int spliceStartPos, int spliceEndPos)
{
	//////////////////////////////////insert to spliceJunction Hash///////////////////////////////////////////
	if(((spliceEndPos - spliceStartPos) < SPLICE_MIN_LENGTH)) //&& ((spliceStartPos - spliceEndPos) < SPLICE_MIN_LENGTH))
	{
		deletionNum++;
		return true;
	}
	//cout << "chrInt = " << chrInt << endl;
	//cout << "spliceStartPos = " << spliceStartPos << endl;
	//cout << "spliceEndPos = " << spliceEndPos << endl;
	char spliceEndChar;
	//string spliceEndString[6];
	SpliceStartPosHashIter iter = spliceJunction[chrInt].find(spliceEndPos); // to find spliceStartPos in Hash
	if(iter == spliceJunction[chrInt].end()) 
	{
		SpliceEndStringHash* newSpliceEndStringHash = new SpliceEndStringHash;
		vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
		spliceEndPosVec.push_back(1);
		spliceEndPosVec.push_back(spliceStartPos);
		//cout << "stop2" << endl;
		//string spliceEndString[5];
		//cout << "stop2.5" << endl;
		spliceEndChar = chromStr[chrInt].at(spliceStartPos-1);
		//spliceEndString[0] = "**" + chromStr[chrInt].substr(spliceStartPos-1, 1);
		//spliceEndString[1] = "*" + chromStr[chrInt].substr(spliceStartPos-2, 2);
		//spliceEndString[2] = chromStr[chrInt].substr(spliceStartPos-3, 3);
		//spliceEndString[3] = chromStr[chrInt].substr(spliceStartPos-6, 3) + "$";
		//spliceEndString[4] = chromStr[chrInt].substr(spliceStartPos-9, 3) + "$$";
		//spliceEndString[5] = "!";
		//cout << "stop3" << endl;
		(*newSpliceEndStringHash).insert(pair<char, vector<int> > (spliceEndChar, spliceEndPosVec)); 
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[0], spliceEndPosVec));  
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[1], spliceEndPosVec));  
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[2], spliceEndPosVec));  
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[3], spliceEndPosVec));  
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[4], spliceEndPosVec)); 
		//(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[5], spliceEndPosVec)); 		 
		//cout << "stop4" << endl;		
		spliceJunction[chrInt].insert(pair<int, SpliceEndStringHash> (spliceEndPos, (*newSpliceEndStringHash)));
	}
	else
	{
		//SpliceEndStringHash (foundSpliceEndStringHash) = iter->second;
		
		//string spliceEndString[5];
		/*spliceEndString[0] = "**" + chromStr[chrInt].substr(spliceStartPos-1, 1);
		spliceEndString[1] = "*" + chromStr[chrInt].substr(spliceStartPos-2, 2);
		spliceEndString[2] = chromStr[chrInt].substr(spliceStartPos-3, 3);
		spliceEndString[3] = chromStr[chrInt].substr(spliceStartPos-6, 3) + "$";
		spliceEndString[4] = chromStr[chrInt].substr(spliceStartPos-9, 3) + "$$";
		spliceEndString[5] = "!";*/
		spliceEndChar = chromStr[chrInt].at(spliceStartPos-1);
		/*for(int tmp = 0; tmp <= 5; tmp++)
		{
			endStringIter = (iter->second).find(spliceEndString[tmp]);
			if(endStringIter == (iter->second).end())
			{
				vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
				spliceEndPosVec.push_back(1);
				spliceEndPosVec.push_back(spliceStartPos);	
				(iter->second).insert(pair<string, vector<int> > (spliceEndString[tmp], spliceEndPosVec)); 			
			}
			else
			{
				(endStringIter->second)[0]++;
				(endStringIter->second).push_back(spliceStartPos);
			}
		}*/
		endStringIter = (iter->second).find(spliceEndChar);
		if(endStringIter == (iter->second).end())
		{
			vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
			spliceEndPosVec.push_back(1);
			spliceEndPosVec.push_back(spliceStartPos);	
			(iter->second).insert(pair<char, vector<int> > (spliceEndChar, spliceEndPosVec));			
		}
		else
		{
			(endStringIter->second)[0]++;
			(endStringIter->second).push_back(spliceStartPos);
		}
	
	}
	/*
	#ifdef DEBUG
	if((spliceEndPos == 114677905)&&(chrInt==9))
	{
		for(int tmpInt = 0; tmpInt <= 4; tmpInt++)
		{
			cout << "string in hash is: " << spliceEndString[tmpInt] << endl;
		}
	}
	#endif
	*/
	return true;
}

bool insertSpliceJunction2NormalHashForTail(SpliceStartPosHash* spliceJunction,
	int chrInt, int spliceStartPos, int spliceEndPos)
{
	//////////////////////////////////insert to spliceJunction Hash///////////////////////////////////////////
	if(((spliceEndPos - spliceStartPos) < SPLICE_MIN_LENGTH))// && ((spliceStartPos - spliceEndPos) < SPLICE_MIN_LENGTH))
	{
		deletionNum++;
		return true;
	}
	SpliceStartPosHashIter iter = spliceJunction[chrInt].find(spliceStartPos); // to find spliceStartPos in Hash
	if(iter == spliceJunction[chrInt].end()) 
	{
		SpliceEndStringHash* newSpliceEndStringHash = new SpliceEndStringHash;
		vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
		spliceEndPosVec.push_back(1);
		spliceEndPosVec.push_back(spliceEndPos);

		/*string spliceEndString[6];
		spliceEndString[0] = chromStr[chrInt].substr(spliceEndPos-1, 1) + "**";
		spliceEndString[1] = chromStr[chrInt].substr(spliceEndPos-1, 2) + "*";
		spliceEndString[2] = chromStr[chrInt].substr(spliceEndPos-1, 3);
		spliceEndString[3] = "$"+chromStr[chrInt].substr(spliceEndPos+2, 3);
		spliceEndString[4] = "$$"+chromStr[chrInt].substr(spliceEndPos+5, 3);;
		spliceEndString[5] = "!";*/
		char spliceEndChar = chromStr[chrInt].at(spliceEndPos-1);
		(*newSpliceEndStringHash).insert(pair<char, vector<int> > (spliceEndChar, spliceEndPosVec)); 
		/*(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[0], spliceEndPosVec));  
		(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[1], spliceEndPosVec));  
		(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[2], spliceEndPosVec));  
		(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[3], spliceEndPosVec));  
		(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[4], spliceEndPosVec));  
		(*newSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[5], spliceEndPosVec));*/  

		spliceJunction[chrInt].insert(pair<int, SpliceEndStringHash> (spliceStartPos, (*newSpliceEndStringHash)));
	}
	else
	{
		SpliceEndStringHash foundSpliceEndStringHash = iter->second;
		
		/*string spliceEndString[6];
		spliceEndString[0] = chromStr[chrInt].substr(spliceEndPos-1, 1) + "**";
		spliceEndString[1] = chromStr[chrInt].substr(spliceEndPos-1, 2) + "*";
		spliceEndString[2] = chromStr[chrInt].substr(spliceEndPos-1, 3);
		spliceEndString[3] = "$"+chromStr[chrInt].substr(spliceEndPos+2, 3);
		spliceEndString[4] = "$$"+chromStr[chrInt].substr(spliceEndPos+5, 3);
		spliceEndString[5] = "!";*/
		char spliceEndChar = chromStr[chrInt].at(spliceEndPos-1);

		/*for(int tmp = 0; tmp <= 5; tmp++)
		{
			endStringIter = (foundSpliceEndStringHash).find(spliceEndString[tmp]);
			if(endStringIter == (foundSpliceEndStringHash).end())
			{
				vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
				spliceEndPosVec.push_back(1);
				spliceEndPosVec.push_back(spliceEndPos);	
				(foundSpliceEndStringHash).insert(pair<string, vector<int> > (spliceEndString[tmp], spliceEndPosVec)); 			
			}
			else
			{
				(endStringIter->second)[0]++;
				(endStringIter->second).push_back(spliceEndPos);
			}
		}*/
		endStringIter = (foundSpliceEndStringHash).find(spliceEndChar);
		if(endStringIter == (foundSpliceEndStringHash).end())
		{
			vector<int> spliceEndPosVec ;//= new spliceEndPosVec;
			spliceEndPosVec.push_back(1);
			spliceEndPosVec.push_back(spliceEndPos);	
			(foundSpliceEndStringHash).insert(pair<char, vector<int> > (spliceEndChar, spliceEndPosVec)); 			
		}
		else
		{
			(endStringIter->second)[0]++;
			(endStringIter->second).push_back(spliceEndPos);
		}		
	}
	return true;
}


/*
bool searchForShortHeadSplice(
	unsigned int mapPos, unsigned int* chrNameInt, 
	unsigned int* chrMapPosInt,
	int leftBound, int rightBound, 
	vector<int> &spliceBufferDistance, 
	vector<int> &spliceDistance, int headLength, const string& readSeqStr)
{
	#ifdef DEBUG
	cout << "start to search for short head splice " << endl;
	cout << "headLength: " << headLength << endl;
	cout << "mapPos: " << mapPos << endl;
	#endif
	bool searchResult = false;
	getChrLocation(mapPos, chrNameInt, chrMapPosInt);
	SpliceStartPosHash spliceJunctionHash = spliceJunctionReverse[(*chrNameInt)];
	SpliceEndStringHash spliceJunctionStringHash;
	SpliceStartPosHashIter iterSpliceHashFound;
	int searchSpliceLeftBound = *chrMapPosInt - leftBound;
	int searchSpliceRightBound = *chrMapPosInt + rightBound;	
	//string tmpSearchString;
	char tmpSearchChar;
	for(int tmp = searchSpliceLeftBound; tmp <= searchSpliceRightBound; tmp++)
	{	
		int tmpSpliceBufferDistance = tmp - *chrMapPosInt;
		int tmpHeadLength = headLength + tmpSpliceBufferDistance;
		#ifdef DEBUG
		cout << endl << "tmpHeadLength: " << tmpHeadLength << endl;
		cout << "tmpSearchSpliceStartPos: " << tmp << endl;
		#endif
		iterSpliceHashFound = spliceJunctionHash.find(tmp);
		if(iterSpliceHashFound == spliceJunctionHash.end())
		{
			#ifdef DEBUG
			cout << "splice start pos not found in Hash" << endl;
			#endif
			continue;
		}
		else
		{
			#ifdef DEBUG
			cout << "splice start pos found in Hash" << endl;
			//cout << "tmpHeadLength = " << tmpHeadLength << endl;
			#endif
			if(tmpHeadLength <= 0)//= tailLength - tmpSpliceBufferDistance
			{
				continue;
			}
			else
			{
				tmpSearchChar = readSeqStr.at(tmpHeadLength-1);
				SpliceEndStringHash tmpSpliceJunctionStringHash = iterSpliceHashFound->second;
				SpliceEndStringHashIter iterSpliceStringHashFound = tmpSpliceJunctionStringHash.find(tmpSearchChar);
				if(iterSpliceStringHashFound != tmpSpliceJunctionStringHash.end())
				{
					searchResult = true;
					int tmpSpliceDistance;
					for(int tmpVec = 1; tmpVec <= (iterSpliceStringHashFound->second)[0]; tmpVec++)
					{
						tmpSpliceDistance = (iterSpliceStringHashFound->second)[tmpVec] - tmp;
						spliceBufferDistance.push_back(tmpSpliceBufferDistance);
						spliceDistance.push_back(tmpSpliceDistance);
					}
				}
				else
				{
					continue;
				}				
			}
		}
	}
	return searchResult;
}*/
/*
bool searchForShortHeadSplice(
	unsigned int mapPos, unsigned int* chrNameInt, 
	unsigned int* chrMapPosInt,
	int leftBound, int rightBound, 
	vector<int> &spliceBufferDistance, 
	vector<int> &spliceDistance, int headLength, const string& readSeqStr)
{
	//searchSpliceOther_begin = clock();
	#ifdef DEBUG
	cout << "start to search for short head splice " << endl;
	cout << "headLength: " << headLength << endl;
	cout << "mapPos: " << mapPos << endl;
	#endif
	bool searchResult = false;
	getChrLocation(mapPos, chrNameInt, chrMapPosInt);

	//SpliceStartPosHash *spliceJunctionHash) = spliceJunctionReverse[(*chrNameInt)];

	//SpliceEndStringHash spliceJunctionStringHash;
	SpliceStartPosHashIter iterSpliceHashFound;
	int searchSpliceLeftBound = *chrMapPosInt - leftBound;
	int searchSpliceRightBound = *chrMapPosInt + rightBound;	
	//string tmpSearchString;
	char tmpSearchChar;
	//searchSpliceOther_end = clock();	
	//searchSpliceOther_cost = searchSpliceOther_cost + searchSpliceOther_end - searchSpliceOther_begin;

	for(int tmp = searchSpliceLeftBound; tmp <= searchSpliceRightBound; tmp++)
	{	
		//search_pos_begin = clock();
		int tmpSpliceBufferDistance = tmp - *chrMapPosInt;
		int tmpHeadLength = headLength + tmpSpliceBufferDistance;
		#ifdef DEBUG
		cout << endl << "tmpHeadLength: " << tmpHeadLength << endl;
		cout << "tmpSearchSpliceStartPos: " << tmp << endl;
		#endif

		
		iterSpliceHashFound = (spliceJunctionReverse[(*chrNameInt)]).find(tmp);
		//search_pos_end = clock();
		//search_pos_cost = search_pos_cost + search_pos_end - search_pos_begin;

		//search_char_begin = clock();
		if(iterSpliceHashFound == (spliceJunctionReverse[(*chrNameInt)]).end())
		{
			#ifdef DEBUG
			cout << "splice start pos not found in Hash" << endl;
			#endif
			continue;
		}
		else
		{
			#ifdef DEBUG
			cout << "splice start pos found in Hash" << endl;
			//cout << "tmpHeadLength = " << tmpHeadLength << endl;
			#endif
			if(tmpHeadLength <= 0)//= tailLength - tmpSpliceBufferDistance
			{
				continue;
			}
			else
			{
				SpliceEndStringHash tmpSpliceJunctionStringHash = iterSpliceHashFound->second;
				for(SpliceEndStringHashIter iterSpliceStringHashFound = tmpSpliceJunctionStringHash.begin(); 
					iterSpliceStringHashFound != tmpSpliceJunctionStringHash.end(); iterSpliceStringHashFound++)
				{
					searchResult = true;
					int tmpSpliceDistance;
					for(int tmpVec = 1; tmpVec <= (iterSpliceStringHashFound->second)[0]; tmpVec++)
					{
						tmpSpliceDistance = (iterSpliceStringHashFound->second)[tmpVec] - tmp;
						spliceBufferDistance.push_back(tmpSpliceBufferDistance);
						spliceDistance.push_back(tmpSpliceDistance);
					}					
				}

			}
		}
		//search_char_end = clock();
		//search_char_cost = search_char_cost + search_char_end - search_char_begin;
	}
	return searchResult;
}*/

/*string covertCharToReverseComplement(const string& Ori_Char)
{
	if(Ori_Char == "A")
	{
		return "T";
	}
	else if(Ori_Char == "T")
	{
		return "A";
	}
	else if(Ori_Char == "G")
	{
		return "C";
	}
	else if(Ori_Char == "C")
	{
		return "G";
	}
	else if(Ori_Char == "N")
	{
		return "N";
	}	
}*/

/*string covertStringToReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}*/

