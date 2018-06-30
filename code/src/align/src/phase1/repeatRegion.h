// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REPEATREGION_H
#define REPEATREGION_H

#include <stdlib.h>
#include <stdio.h>
//#include "enhanced_suffix_array_info.h"

using namespace std;

class RepeatRegion_Info
{
public:
	int index_repeatRegionVec_currentSize;

	map < pair<unsigned int, unsigned int>, int> repeatRegionMap; // <SA_interval_begin, SA_interval_end>, index_repeatRegionVec

	vector< pair<int, int> > repeatRegionVec;


	RepeatRegion_Info()
	{
		index_repeatRegionVec_currentSize = 0;
	}

	int push2RepeatRegionInfo(unsigned int interval_begin, unsigned int interval_end, unsigned int* sa, Index_Info* indexInfo)
	{
		//cout << "interval_begin: " << interval_begin << " interval_end: " << interval_end << endl;
		map<pair<unsigned int, unsigned int>, int>::iterator repeatRegionMapIter;
		//cout << "start to search " << endl;
		//cout << "map size: " << repeatRegionMap.size() << endl;
		repeatRegionMapIter = repeatRegionMap.find(pair<unsigned int, unsigned int> (interval_begin, interval_end) );
		//cout << "finish searching " << endl;
		if(repeatRegionMapIter != repeatRegionMap.end())
		{
			//cout << "found in repeatRegionMap !" << endl;
			return repeatRegionMapIter->second;
		}
		else
		{
			//cout << "not found in repeatRegionMap !" << endl; 
			repeatRegionMap.insert(pair< pair<unsigned int, unsigned int>, int> 
				(pair<unsigned int, unsigned int>(interval_begin, interval_end), index_repeatRegionVec_currentSize) );
			//cout << "insert finished!" << endl;


			unsigned int posInWholeGenome = sa[interval_begin] + 1;
			//cout << "posInWholeGenome " << posInWholeGenome << endl;

			unsigned int tmpChrNameInt;
			unsigned int tmpChrPosInt;
			indexInfo->getChrLocation(posInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);


			int chrNameInt = (tmpChrNameInt);
			int chrPosInt = (tmpChrPosInt);  

			//cout << "chrNameInt " << chrNameInt << " chrPosInt" << chrPosInt << endl;

			repeatRegionVec.push_back(pair<int,int> (chrNameInt, chrPosInt));



			index_repeatRegionVec_currentSize ++;
			return (index_repeatRegionVec_currentSize - 1); 
		}
	}

	int push2RepeatRegionInfo(Enhanced_Suffix_Array_Info* esaInfo,
		unsigned int interval_begin, unsigned int interval_end, Index_Info* indexInfo)
	{
		//cout << "interval_begin: " << interval_begin << " interval_end: " << interval_end << endl;
		map<pair<unsigned int, unsigned int>, int>::iterator repeatRegionMapIter;
		//cout << "start to search " << endl;
		//cout << "map size: " << repeatRegionMap.size() << endl;
		repeatRegionMapIter = repeatRegionMap.find(pair<unsigned int, unsigned int> (interval_begin, interval_end) );
		//cout << "finish searching " << endl;
		if(repeatRegionMapIter != repeatRegionMap.end())
		{
			//cout << "found in repeatRegionMap !" << endl;
			return repeatRegionMapIter->second;
		}
		else
		{
			//cout << "not found in repeatRegionMap !" << endl; 
			repeatRegionMap.insert(pair< pair<unsigned int, unsigned int>, int> 
				(pair<unsigned int, unsigned int>(interval_begin, interval_end), index_repeatRegionVec_currentSize) );
			//cout << "insert finished!" << endl;


			unsigned int posInWholeGenome = esaInfo->returnSA(interval_begin) + 1;
			//cout << "posInWholeGenome " << posInWholeGenome << endl;

			unsigned int tmpChrNameInt;
			unsigned int tmpChrPosInt;
			indexInfo->getChrLocation(posInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);


			int chrNameInt = (tmpChrNameInt);
			int chrPosInt = (tmpChrPosInt);  

			//cout << "chrNameInt " << chrNameInt << " chrPosInt" << chrPosInt << endl;

			repeatRegionVec.push_back(pair<int,int> (chrNameInt, chrPosInt));



			index_repeatRegionVec_currentSize ++;
			return (index_repeatRegionVec_currentSize - 1); 
		}
	}

	string getRepeatRegionStr(int repeatRegionInfo_tmpNum, Index_Info* indexInfo, 
		unsigned int* sa, int max_outputIntervalSize)
	{
		int repeatRegionMapSize = repeatRegionMap.size();
		vector< pair<unsigned int, unsigned int> > repeatRegion_SAintervalVec(repeatRegionMapSize);

		for(map < pair<unsigned int, unsigned int>, int>::iterator tmpIter = repeatRegionMap.begin();
			tmpIter != repeatRegionMap.end(); tmpIter ++)
		{
			int tmpIndex = tmpIter->second;
			unsigned int tmpSAinterval_begin = (tmpIter->first).first;
			unsigned int tmpSAinterval_end = (tmpIter->first).second;
			repeatRegion_SAintervalVec[tmpIndex].first = tmpSAinterval_begin;
			repeatRegion_SAintervalVec[tmpIndex].second = tmpSAinterval_end;
		}
	
		string repeatRegionStr;// = "Repeat Region Info " + int_to_str(repeatRegionInfo_tmpNum + 1) + " :\n";

		for(int tmp = 0; tmp < repeatRegionMapSize; tmp++)
		{
			repeatRegionStr = repeatRegionStr + int_to_str(repeatRegionInfo_tmpNum + 1) + " -- RP:" + int_to_str(tmp);
			for(int tmp_SA_index = repeatRegion_SAintervalVec[tmp].first; tmp_SA_index <= repeatRegion_SAintervalVec[tmp].second;
				tmp_SA_index ++)
			{
				int tmpNO = tmp_SA_index - repeatRegion_SAintervalVec[tmp].first + 1;

				if(tmpNO > max_outputIntervalSize)
				{
					break;
				}

				unsigned int tmpPosInWholeGenome = sa[tmp_SA_index] + 1;

				unsigned int tmpChrNameInt;
				unsigned int tmpChrPosInt;
				indexInfo->getChrLocation(tmpPosInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);

				int chrNameInt = (tmpChrNameInt);
				int chrPosInt = (tmpChrPosInt);  

				repeatRegionStr = repeatRegionStr + " " + int_to_str(tmpNO) + "-" + indexInfo->returnChrNameStr(chrNameInt) + "," + int_to_str(chrPosInt);
			}  
		}
	}

	string getRepeatRegionStr(Enhanced_Suffix_Array_Info* esaInfo,
		int repeatRegionInfo_tmpNum, Index_Info* indexInfo, 
		//unsigned int* sa, 
		int max_outputIntervalSize)
	{
		int repeatRegionMapSize = repeatRegionMap.size();
		vector< pair<unsigned int, unsigned int> > repeatRegion_SAintervalVec(repeatRegionMapSize);

		for(map < pair<unsigned int, unsigned int>, int>::iterator tmpIter = repeatRegionMap.begin();
			tmpIter != repeatRegionMap.end(); tmpIter ++)
		{
			int tmpIndex = tmpIter->second;
			unsigned int tmpSAinterval_begin = (tmpIter->first).first;
			unsigned int tmpSAinterval_end = (tmpIter->first).second;
			repeatRegion_SAintervalVec[tmpIndex].first = tmpSAinterval_begin;
			repeatRegion_SAintervalVec[tmpIndex].second = tmpSAinterval_end;
		}
	
		string repeatRegionStr;// = "Repeat Region Info " + int_to_str(repeatRegionInfo_tmpNum + 1) + " :\n";

		for(int tmp = 0; tmp < repeatRegionMapSize; tmp++)
		{
			repeatRegionStr = repeatRegionStr + int_to_str(repeatRegionInfo_tmpNum + 1) + " -- RP:" + int_to_str(tmp);
			for(int tmp_SA_index = repeatRegion_SAintervalVec[tmp].first; tmp_SA_index <= repeatRegion_SAintervalVec[tmp].second;
				tmp_SA_index ++)
			{
				int tmpNO = tmp_SA_index - repeatRegion_SAintervalVec[tmp].first + 1;

				if(tmpNO > max_outputIntervalSize)
				{
					break;
				}

				unsigned int tmpPosInWholeGenome = esaInfo->returnSA(tmp_SA_index) + 1;

				unsigned int tmpChrNameInt;
				unsigned int tmpChrPosInt;
				indexInfo->getChrLocation(tmpPosInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);

				int chrNameInt = (tmpChrNameInt);
				int chrPosInt = (tmpChrPosInt);  

				repeatRegionStr = repeatRegionStr + " " + int_to_str(tmpNO) + "-" + indexInfo->returnChrNameStr(chrNameInt) + "," + int_to_str(chrPosInt);
			}  
		}
	}

	void outputRepeatRegion(int repeatRegionInfo_tmpNum, Index_Info* indexInfo, 
		unsigned int* sa, int max_outputIntervalSize, ofstream& repeatRegion_ofs)
	{
		int repeatRegionMapSize = repeatRegionMap.size();
		vector< pair<unsigned int, unsigned int> > repeatRegion_SAintervalVec(repeatRegionMapSize);

		for(map < pair<unsigned int, unsigned int>, int>::iterator tmpIter = repeatRegionMap.begin();
			tmpIter != repeatRegionMap.end(); tmpIter ++)
		{
			int tmpIndex = tmpIter->second;
			unsigned int tmpSAinterval_begin = (tmpIter->first).first;
			unsigned int tmpSAinterval_end = (tmpIter->first).second;
			repeatRegion_SAintervalVec[tmpIndex].first = tmpSAinterval_begin;
			repeatRegion_SAintervalVec[tmpIndex].second = tmpSAinterval_end;
		}
	
		//string repeatRegionStr;// = "Repeat Region Info " + int_to_str(repeatRegionInfo_tmpNum + 1) + " :\n";

		for(int tmp = 0; tmp < repeatRegionMapSize; tmp++)
		{
			repeatRegion_ofs << repeatRegionInfo_tmpNum << " -- RP:" << tmp;// << " size: " << repeatRegion_SAintervalVec[tmp].second - repeatRegion_SAintervalVec[tmp].first + 1;
			for(int tmp_SA_index = repeatRegion_SAintervalVec[tmp].first; tmp_SA_index <= repeatRegion_SAintervalVec[tmp].second;
				tmp_SA_index ++)
			{
				int tmpNO = tmp_SA_index - repeatRegion_SAintervalVec[tmp].first + 1;

				if(tmpNO > max_outputIntervalSize)
				{
					break;
				}

				unsigned int tmpPosInWholeGenome = sa[tmp_SA_index] + 1;

				unsigned int tmpChrNameInt;
				unsigned int tmpChrPosInt;
				indexInfo->getChrLocation(tmpPosInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);

				int chrNameInt = (tmpChrNameInt);
				int chrPosInt = (tmpChrPosInt);  

				repeatRegion_ofs << "    " << tmpNO << ".." << indexInfo->returnChrNameStr(chrNameInt) << "," << int_to_str(chrPosInt);
			}  
			repeatRegion_ofs << endl;
		}
	}

	void outputRepeatRegion(Enhanced_Suffix_Array_Info* esaInfo,
		int repeatRegionInfo_tmpNum, Index_Info* indexInfo, 
		//unsigned int* sa, 
		int max_outputIntervalSize, ofstream& repeatRegion_ofs)
	{
		int repeatRegionMapSize = repeatRegionMap.size();
		vector< pair<unsigned int, unsigned int> > repeatRegion_SAintervalVec(repeatRegionMapSize);

		for(map < pair<unsigned int, unsigned int>, int>::iterator tmpIter = repeatRegionMap.begin();
			tmpIter != repeatRegionMap.end(); tmpIter ++)
		{
			int tmpIndex = tmpIter->second;
			unsigned int tmpSAinterval_begin = (tmpIter->first).first;
			unsigned int tmpSAinterval_end = (tmpIter->first).second;
			repeatRegion_SAintervalVec[tmpIndex].first = tmpSAinterval_begin;
			repeatRegion_SAintervalVec[tmpIndex].second = tmpSAinterval_end;
		}
	
		//string repeatRegionStr;// = "Repeat Region Info " + int_to_str(repeatRegionInfo_tmpNum + 1) + " :\n";

		for(int tmp = 0; tmp < repeatRegionMapSize; tmp++)
		{
			repeatRegion_ofs << repeatRegionInfo_tmpNum << " -- RP:" << tmp;// << " size: " << repeatRegion_SAintervalVec[tmp].second - repeatRegion_SAintervalVec[tmp].first + 1;
			for(int tmp_SA_index = repeatRegion_SAintervalVec[tmp].first; tmp_SA_index <= repeatRegion_SAintervalVec[tmp].second;
				tmp_SA_index ++)
			{
				int tmpNO = tmp_SA_index - repeatRegion_SAintervalVec[tmp].first + 1;

				if(tmpNO > max_outputIntervalSize)
				{
					break;
				}

				unsigned int tmpPosInWholeGenome = esaInfo->returnSA(tmp_SA_index) + 1;

				unsigned int tmpChrNameInt;
				unsigned int tmpChrPosInt;
				indexInfo->getChrLocation(tmpPosInWholeGenome, &tmpChrNameInt, &tmpChrPosInt);

				int chrNameInt = (tmpChrNameInt);
				int chrPosInt = (tmpChrPosInt);  

				repeatRegion_ofs << "    " << tmpNO << ".." << indexInfo->returnChrNameStr(chrNameInt) << "," << int_to_str(chrPosInt);
			}  
			repeatRegion_ofs << endl;
		}
	}

};

#endif