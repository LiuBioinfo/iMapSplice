// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INDEX_INFO_H
#define INDEX_INFO_H

#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

using namespace std;

int INDEX_KMER_LENGTH = 14;

int baseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 100,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};

int baseCharCount2intArray[14][26] = {0};


class Index_Info
{
private:
	string chromString;
	unsigned int genomeLength;

	int chromNum;

	vector<string> chrNameStr; // size = chromNum
	vector<int> chromLength; // size = chromNum
	vector<string> chromStr;
	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	//map<string, int>::iterator chrNameMapIter;

	int secondLevelIndexNormalSize;// = 3000000;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	set<int> invalidSecondLevelIndexNOset;
	//vector<int> secondLevelIndexLengthVec;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 

	vector<int> chrNameIndexArray;
	int chrNameIndexIntervalSize;
	map<int, set<int> > chrNameIndexArrayMap; // index_chrNameIndexArray, set<chrNameInt>
public:
	
	void initiateChrNameIndexArray(int intervalSize)
	{
		chrNameIndexIntervalSize = intervalSize;
		int chrNameIndexArraySize = (chrEndPosInGenome[chromNum-1])/chrNameIndexIntervalSize + 1;
		for(int tmp = 0; tmp < chrNameIndexArraySize; tmp++)
		{
			chrNameIndexArray.push_back(-2);
		}

		int indexArray_start_1 = 1/chrNameIndexIntervalSize;
		int indexArray_end_1 = chrEndPosInGenome[0]/chrNameIndexIntervalSize;
		for(int tmp = indexArray_start_1; tmp <= indexArray_end_1; tmp++)
		{
			chrNameIndexArray[tmp] = 0;
		}

		for(int tmp = 1; tmp < chromNum; tmp++)
		{
			int indexArray_start_tmp = (chrEndPosInGenome[tmp-1]+1)/chrNameIndexIntervalSize;
			int indexArray_end_tmp = chrEndPosInGenome[tmp]/chrNameIndexIntervalSize;
			for(int tmp2 = indexArray_start_tmp; tmp2 <= indexArray_end_tmp; tmp2++)
			{
				if(chrNameIndexArray[tmp2] == -2)
				{
					chrNameIndexArray[tmp2] = tmp;
				}
				else if(chrNameIndexArray[tmp2] == -1)
				{
					map<int, set<int> >::iterator mapIter_found = chrNameIndexArrayMap.find(tmp2);
					(mapIter_found->second).insert(tmp);
				}
				else
				{
					set<int> tmpIntSet;
					tmpIntSet.insert(tmp-1);
					tmpIntSet.insert(tmp);
					chrNameIndexArrayMap.insert(pair<int, set<int> >(tmp2, tmpIntSet));
					chrNameIndexArray[tmp2] = -1;
				}
			}
		}
	}

	/*
	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		#ifdef CAL_TIME
		getChrLocation_begin = clock();
		#endif
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			(*chr_name_int) = 0;
			*chr_local_location = locationInWholeGenome;
		}
		else
		{
			for(int tmp = 1; tmp < chrEndPosInGenome.size(); tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					*chr_name_int = tmp;
					*chr_local_location = locationInWholeGenome - chrEndPosInGenome[tmp-1] - 2;
				}
				else
				{
					continue;
				}
			}
		}
		#ifdef CAL_TIME
		getChrLocation_end = clock();
		getChrLocation_cost = getChrLocation_cost + getChrLocation_end - getChrLocation_begin;
		#endif		
	}*/

	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		int chrNameInt = this->getChr(locationInWholeGenome);
		if(chrNameInt == 0)
		{
			(*chr_name_int) = 0;
			(*chr_local_location) = locationInWholeGenome;
		}
		else
		{
			(*chr_name_int) = chrNameInt;
			(*chr_local_location) = locationInWholeGenome - chrEndPosInGenome[chrNameInt - 1] -2;
		}
	}

	/*int getChr(unsigned int locationInWholeGenome)
	{
		#ifdef CAL_TIME
		getChrLocation_begin = clock();
		#endif		
		int chrInt;
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			chrInt = 0;
		}
		else
		{
			for(int tmp = 1; tmp < chromNum; tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					chrInt = tmp;
					break;
				}
				else
				{
					continue;
				}				
			}
		}
		#ifdef CAL_TIME
		getChrLocation_end = clock();
		getChrLocation_cost = getChrLocation_cost + getChrLocation_end - getChrLocation_begin;
		#endif				
		return chrInt;
	}*/

	int getChr(unsigned int locationInWholeGenome)
	{
		int chrNameIndexArray_index = locationInWholeGenome/chrNameIndexIntervalSize;
		int index_chrNameVec = chrNameIndexArray[chrNameIndexArray_index];
		
		int getChrInt = -3;
		if(index_chrNameVec == -1)
		{
			map<int, set<int> >::iterator mapIter_found = 
				chrNameIndexArrayMap.find(chrNameIndexArray_index);

			if(mapIter_found == chrNameIndexArrayMap.end())
			{
				cout << "error in get chr -- search in chrNameIndexArrayMap" << endl;
				exit(1);
			}
			else
			{
				for(set<int>::iterator setIter = (mapIter_found->second).begin(); 
					setIter != (mapIter_found->second).end(); setIter ++) 
				{
					int tmpIndex_chrNameVec = (*setIter);
					int tmpChr_startPos, tmpChr_endPos;
					if(tmpIndex_chrNameVec == 0)
					{
						tmpChr_startPos = 1;
						tmpChr_endPos = chrEndPosInGenome[0];
					}
					else
					{
						tmpChr_startPos = chrEndPosInGenome[tmpIndex_chrNameVec-1]+2;
						tmpChr_endPos = chrEndPosInGenome[tmpIndex_chrNameVec];
					}

					if( (locationInWholeGenome >= tmpChr_startPos) && (locationInWholeGenome <= tmpChr_endPos) )
					{
						getChrInt = tmpIndex_chrNameVec;
						break;
					}
				}
				//cout << "some error in getChr -- 1" << endl;
			}
		}
		else if(index_chrNameVec == -2)
		{
			cout << "error in getChr ..." << endl;
			exit(1);
		}
		else
		{
			getChrInt = index_chrNameVec;
		}

		if(getChrInt == -3)
		{
			cout << "getChrInt == -3" << endl;
			cout << "locationInWholeGenome: " << locationInWholeGenome << endl;
			cout << "chrNameIndexArray_index: " << chrNameIndexArray_index << endl;
			cout << "index_chrNameVec: " << index_chrNameVec << endl;
			/*if(index_chrNameVec == -1)
			{
				map<int, set<int> >::iterator mapIter_found = 
					chrNameIndexArrayMap.find(chrNameIndexArray_index);			
				for(set<int>::iterator setIter = (mapIter_found->second).begin(); 
					setIter != (mapIter_found->second).begin(); setIter ++) 
				{

				}				
			}*/
			exit(1);
		}
		return getChrInt;
	}

	bool returnInvalidSecondLevelIndexNOset_find_bool(int value_to_find)
	{
		return ((invalidSecondLevelIndexNOset).find(value_to_find)
							!= (invalidSecondLevelIndexNOset).end()
							);
	}
	void insert2invalidSecondLevelIndexNOset(int value_to_insert)
	{
		invalidSecondLevelIndexNOset.insert(value_to_insert);
	} 

	int returnSecondLevelIndexNormalSize()
	{
		return secondLevelIndexNormalSize;
	}
	int returnSecondLevelIndexPartsNum(int secondLevelIndexPartsNum_index)
	{
		return secondLevelIndexPartsNum[secondLevelIndexPartsNum_index];
	}
	int returnSecondLevelIndexPartsNumSum()
	{
		return secondLevelIndexPartsNumSum;
	}
	unsigned int returnNull_num()
	{
		return null_num;
	}
	unsigned int returnIndexSize()
	{
		return indexSize;
	}
	int returnChrNameStrSize()
	{
		return chrNameStr.size();
	}
	const string& returnChrNameStr(int chrName_index)
	{
		return chrNameStr[chrName_index];
	}
	int returnChromLength(int chromLength_index)
	{
		return chromLength[chromLength_index];
	}
	const string& returnChromStr(int chromStr_index)
	{		
		return chromStr[chromStr_index];
	}
	//const string& 
	string returnChromStrSubstr(int chromStr_index, int start_pos, int seq_length)
	{
		return chromStr[chromStr_index].substr(start_pos - 1, seq_length);
	}
	unsigned int returnChrEndPosInGenome(int chrEndPosInGenome_index)
	{
		return chrEndPosInGenome[chrEndPosInGenome_index];
	}
	void initiate()
	{
		//(indexInfo->chromStr).push_back((indexInfo->returnChromString()).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
		(chromStr).push_back(this->returnChromStringSubstr(1, (chrEndPosInGenome)[0]+1));
		(chromLength).push_back(((chrEndPosInGenome)[0]+1));
		for(int tmp = 1; tmp < this->returnChromNum(); tmp++)
		{
			//chromStr[tmp] = 
			//(indexInfo->chromStr).push_back((indexInfo->returnChromString()).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			//	(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
			(chromStr).push_back(
				this->returnChromStringSubstr((chrEndPosInGenome)[tmp-1]+3,
				(chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1));// (indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1);

			(chromLength).push_back(((chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1));
		}		
	}

	unsigned int returnChromStringLength()
	{
		return chromString.length();
	}
	void readGenome(char* chrom)
	{
		chromString = chrom;
	}
	const string& returnChromString()
	{
		return chromString;
	}

	string returnChromStringSubstr(unsigned int start_pos, unsigned int string_length)
	{
		return chromString.substr(start_pos - 1, string_length);
	}

	unsigned int returnGenomeLength()
	{
		return genomeLength;
	}
	int returnChromNum()
	{
		return chromNum;
	}

	string returnFlankString(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		return (chromStr[chrNameInt]).substr(chromPos_doner, 2) 
			+ (chromStr[chrNameInt]).substr(chromPos_acceptor-3, 2);
	}

	char returnOneBaseCharInGenome(int chrNameInt, int chromPos)
	{
		return (chromStr[chrNameInt]).at(chromPos-1);
	}

	string returnGenomeSubstr(int chrNameInt, int startPosInChrom, int substrLength)
	{
		return (chromStr[chrNameInt]).substr(startPosInChrom - 1, substrLength);
	}

	string getReferenceGenomeSubstr(const string& chrName, int startPos, int endPos) // can be used in debugging on some special case
	{
		int chrNameInt = this->convertStringToInt(chrName);
		return (this->chromStr)[chrNameInt].substr(startPos-1, endPos - startPos + 1);
	}

	string getInvalidSecondLevelIndexNOstr()
	{
		string tmpStr = "invalidSecondLevelIndexNO: \n";
		for(set<int>::iterator setIter = invalidSecondLevelIndexNOset.begin(); setIter != invalidSecondLevelIndexNOset.end(); setIter ++)
		{
			tmpStr += int_to_str(*setIter);
			tmpStr += ",";
		} 
		tmpStr += "\n";
		return tmpStr;
	}

	Index_Info()
	{}

	Index_Info(ifstream& inputIndexInfoFile, ofstream& outputIndexInfoFile)
	{

		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	Index_Info(ifstream& inputIndexInfoFile)
	{

		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		//outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		//outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		//outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		//outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		//outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		//outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		//outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
		}

	}

	int getSecondLevelIndexFromChrAndPos(int chrNameInt, int chrMapPos)
	{
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}
		return	(partsTimeBase + tmpTimes + 1); 
	}

	int getSecondLevelIndexFromChrStrAndPos(string chrNameStr, int chrMapPos)
	{
		int chrNameInt = this->convertStringToInt(chrNameStr);
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}
		return	(partsTimeBase + tmpTimes + 1); 
	} 

	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}

		int tmpSecondLevelIndexNO = secondLevelIndexNum - partsTimeBase;
		return ( (tmpSecondLevelIndexNO-1) * secondLevelIndexNormalSize + secondLevelIndexPos);	
	}

	int getChrNameIntFromSecondLevelIndexNO(int secondLevelIndexNum)
	{
		int tmpIndexNOsum = 0;
		int tmp;
		for(tmp = 0; tmp < chrNameStr.size(); tmp++)
		{
			tmpIndexNOsum += secondLevelIndexPartsNum[tmp];
			if(tmpIndexNOsum >= secondLevelIndexNum)
				return tmp;
		}
	}


	unsigned int getWholeGenomeLocation(unsigned int chr_name_int, unsigned int locationInWholeGenome)
	{
		//cout << "in function chr_name_int: " << chr_name_int << endl;
		unsigned int chr_local_location;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			//cout << "< chromNum chr_name_int: " << chr_name_int << endl;
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2; 
			//cout << "chr_local_location: " << chr_local_location << endl;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	unsigned int getWholeGenomeLocation(const string& chromNameStr, unsigned int locationInWholeGenome)
	{
		//cout << "in function chr_name_int: " << chr_name_int << endl;
		
		int chr_name_int = this->convertStringToInt(chromNameStr);
		
		unsigned int chr_local_location;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			//cout << "< chromNum chr_name_int: " << chr_name_int << endl;
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2; 
			//cout << "chr_local_location: " << chr_local_location << endl;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	int convertStringToInt(const string& chrName)
	{
		map<string, int>::iterator chrNameMapIter;
		int chrNameInt = 1000;
		chrNameMapIter = chrNameMap.find(chrName);
		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			chrNameInt = -1;
			//cout << "...... chrom name error! ...... " << endl;
		}
		return chrNameInt;
	}

};

#endif