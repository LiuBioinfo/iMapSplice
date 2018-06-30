// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdlib.h>
#include <stdio.h>

using namespace std;


class Frag_Info
{
public:
	int mapLabelNum;
	unsigned int segMapRangeStart[SEGMENTNUM * CANDALILOC];
	unsigned int segMapRangeEnd[SEGMENTNUM * CANDALILOC];
	unsigned int segMapLoc[SEGMENTNUM * CANDALILOC];

	Frag_Info()
	{
		mapLabelNum = 0;
	}

	int detAliLocShortSegIncluded_FragInfo( Seg_Info* segInfo, unsigned int oldMapLabel, 
		//unsigned int* oldSegMapRangeStart, unsigned int* oldSegMapRangeEnd, 
		unsigned int* oldSegMapPos, 
		Index_Info* indexInfo
		)
	{
		bool detAliLocShortSegIncluded = false;
		unsigned int segmentNum_max = SEGMENTNUM;
		unsigned int segmentAliNum_max = CANDALILOC;
		map <unsigned int, unsigned int> candAliLocMap;
		map <unsigned int, unsigned int> ::iterator candAliLoc_it;
		unsigned int mapLabel = 0; // normal
		unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal


		#ifdef DEBUG
		cout << "debug -- start second level detAliLoc for alignment..." << endl;
		#endif
		for (unsigned int tmpSegNum = 1; tmpSegNum <= segInfo->segmentNum/*segmentNum*/; tmpSegNum++)
		{
			//if (//(*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
			//	(*(segmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
			if( segInfo->norSegmentAlignNum[tmpSegNum-1] > POSSIBLE_MAP_CASES_MAX)
			{
				continue;
			}		

			for(unsigned int tmpSegAliLoc = 1; tmpSegAliLoc <= segInfo->norSegmentAlignNum[tmpSegNum-1]; tmpSegAliLoc++) 
			{				
				//unsigned int tmpAlignLoc = *(segmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
				//- *(segmentLocInRead + tmpSegNum - 1) + 1;
				
				unsigned int tmpAlignLoc = segInfo->norSegmentAlignLoc[(tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1]
					- segInfo->norSegmentLocInRead[tmpSegNum-1] + 1;

				if(tmpAlignLoc > indexInfo->indexSize)
				{
					continue;
				}
				candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

				//if ((*(segmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				if ( (segInfo->norSegmentLength[tmpSegNum-1] >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()) ) //cannot find, insert
				{				
					candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
					//candAliLocWeight[mapLabel] = *(segmentLength + tmpSegNum - 1);
					candAliLocWeight[mapLabel] = segInfo->norSegmentLength[tmpSegNum-1];
					segMapRangeStart[mapLabel] = tmpSegNum;
					segMapRangeEnd[mapLabel] = tmpSegNum;
					segMapLoc[mapLabel] = tmpAlignLoc;
					mapLabel ++;
				}
				else if ((/**(segmentLength + tmpSegNum - 1)*/segInfo->norSegmentLength[tmpSegNum-1] < minValSegLength) 
					&& (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				{	
					if(checkShortSegWithCandLoc(tmpAlignLoc, oldMapLabel, oldSegMapPos, indexInfo))
					{			
						candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
						candAliLocWeight[mapLabel] = segInfo->norSegmentLength[tmpSegNum-1]; //*(segmentLength + tmpSegNum - 1);

						segMapRangeStart[mapLabel] = tmpSegNum;
						segMapRangeEnd[mapLabel] = tmpSegNum;
						segMapLoc[mapLabel] = tmpAlignLoc;
						mapLabel ++;
					}
				}
				else if ((//*(segmentLength + tmpSegNum - 1) 
					segInfo->norSegmentLength[tmpSegNum-1] >= minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
				{

					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] 
						+ segInfo->norSegmentLength[tmpSegNum-1]; //*(segmentLength + tmpSegNum - 1);
					segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
				else //((*(segmentLength + tmpSegNum - 1) < minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
				{
					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] 
						+ segInfo->norSegmentLength[tmpSegNum-1];//*(segmentLength + tmpSegNum - 1);
					segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
			}
		}
		candAliLocMap.clear();	

		return mapLabel;
	}

	bool detAliLoc_FragInfo( Seg_Info* segInfo,
		char* read, char* chrom, 
		Index_Info* indexInfo)
	{
		#ifdef DEBUG
		cout << "debug -- start detAliLoc function..." << endl;
		#endif

		//int* finalMapLabel;
		//(*finalMapLabel) = 0;

		int finalMapLabel;
		bool detAliLoc = false;
		unsigned int segmentNum_max = SEGMENTNUM;
		unsigned int segmentAliNum_max = CANDALILOC;
		
		unsigned int norSegMapRangeStart[SEGMENTNUM * CANDALILOC] = {0};
		unsigned int norSegMapRangeEnd[SEGMENTNUM * CANDALILOC] = {0};
		
		unsigned int norSegMapLoc[SEGMENTNUM * CANDALILOC] = {0};

		unsigned int norValSegRange[2] = {50, 0};

		map <unsigned int, unsigned int> candAliLocMap;
		map <unsigned int, unsigned int> ::iterator candAliLoc_it;
		unsigned int mapLabel = 0; // normal
		unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal

		#ifdef DEBUG
		cout << "debug -- start detAliLoc for alignment..." << endl;
		#endif
		for (unsigned int tmpSegNum = 1; tmpSegNum <= segInfo->segmentNum; tmpSegNum++)
		{
			//if ((*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
			//	(*(norSegmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
			if( ( segInfo->norSegmentLength[tmpSegNum-1] < minValSegLength ) || 
				( segInfo->norSegmentAlignNum[tmpSegNum - 1] > POSSIBLE_MAP_CASES_MAX ) )
			{
				continue;
			}			

			for(unsigned int tmpSegAliLoc = 1; tmpSegAliLoc <= segInfo->norSegmentAlignNum[tmpSegNum-1]; tmpSegAliLoc++) 
			{				
				//unsigned int tmpAlignLoc = *(norSegmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
				//- *(norSegmentLocInRead + tmpSegNum - 1) + 1;
				
				unsigned int tmpAlignLoc = segInfo->norSegmentAlignLoc[(tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1]
					- segInfo->norSegmentLocInRead[tmpSegNum - 1] + 1;

				if(tmpAlignLoc > (indexInfo->indexSize))
				{
					continue;
				}
				candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

				//if ((*(norSegmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				if( ( segInfo->norSegmentLength[tmpSegNum-1] >= minValSegLength ) 
					&& ( candAliLoc_it == candAliLocMap.end() ) )
				{				
					candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
					//candAliLocWeight[mapLabel] = *(norSegmentLength + tmpSegNum - 1);
					candAliLocWeight[mapLabel] = segInfo->norSegmentLength[tmpSegNum-1];

					norSegMapRangeStart[mapLabel] = tmpSegNum;
					norSegMapRangeEnd[mapLabel] = tmpSegNum;
					norSegMapLoc[mapLabel] = tmpAlignLoc;
					mapLabel ++;
				}
				//else if (*(norSegmentLength + tmpSegNum - 1) >= minValSegLength)// find, add corresponding weight 
				else if( segInfo->norSegmentLength[tmpSegNum-1] >= minValSegLength ) // find, add corresponding weight
				{
					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + segInfo->norSegmentLength[tmpSegNum-1];//*(norSegmentLength + tmpSegNum - 1);
					norSegMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
				else
				{
				}
			}
		}
		candAliLocMap.clear();

		if((mapLabel <= 0) || (mapLabel > (POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/)))
		{
			#ifdef DEBUG
			cout << "debug -- mapLabel too small or large : " << mapLabel << endl;
			#endif
			(finalMapLabel) = 0;
			return false;
		}

		#ifdef DEBUG
		cout << "# of Label is : " << mapLabel << endl;
		for (unsigned int tmp_label = 0; tmp_label < mapLabel; tmp_label++)
		{
			unsigned int tmpSegMapChr, tmpSegMapChrPos;
			getChrLocation(norSegMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

			cout << "label = " << tmp_label << ", loc = " << //norSegMapLoc[tmp_label] 
			chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
				<< norSegMapRangeStart[tmp_label] << "~" << norSegMapRangeEnd[tmp_label] << endl;
		}
		#endif

		//////////////try second Level detAlignLoc

		(finalMapLabel) = this->detAliLocShortSegIncluded_FragInfo(segInfo,
			mapLabel, norSegMapLoc, //segMapRangeStart, segMapRangeEnd, segMapLoc, 
			indexInfo);	

		if(((finalMapLabel) <= 0) || ((finalMapLabel) > POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/))
		{
			#ifdef DEBUG
			cout << "debug -- mapLabel too small or large : " << (finalMapLabel) << endl;
			#endif
			(finalMapLabel) = 0;
			return false;
		}	
		
		#ifdef DEBUG
		cout << "# of Label is : " << (finalMapLabel) << endl;
		for (unsigned int tmp_label = 0; tmp_label < (finalMapLabel); tmp_label++)
		{
			unsigned int tmpSegMapChr, tmpSegMapChrPos;
			getChrLocation(segMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

			cout << "label = " << tmp_label << ", loc = " << //segMapLoc[tmp_label] 
			chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
				<< segMapRangeStart[tmp_label] << "~" << segMapRangeEnd[tmp_label] << endl;
		}
		#endif

		mapLabelNum = (finalMapLabel);
		detAliLoc = true;
		return detAliLoc;
	}

	/*void convertPathInfo2FragInfo(Path_Info* pathInfo, Seg_Info* segInfo)
	{
		map < unsigned int, set<int> >  mapLabelMap;
		typedef map < unsigned int, set<int> >::iterator mapLabelMapIter; 
		for(int tmpPath = 0; tmpPath < (pathInfo->PathVec_seg).size(); tmpPath++)
		{
			for(int tmpPathElement = 0; tmpPathElement < (pathInfo->PathVec_seg)[tmpPath].size(); tmpPathElement ++)
			{
				int tmpSegGroupNO = ((pathInfo->PathVec_seg[tmpPath])[tmpPathElement]).first;
				int tmpSegCandiNO = ((pathInfo->PathVec_seg[tmpPath])[tmpPathElement]).second;
				unsigned int tmp1stBaseAlignLoc = *(segInfo->norSegmentAlignLoc + tmpSegGroupNO * SEGMENTNUM + tmpSegCandiNO) 
					- segInfo->norSegmentLocInRead[tmpSegGroupNO] + 1;
				mapLabelMapIter tmpIter = mapLabelMap.find(tmp1stBaseAlignLoc);
				if(tmpIter == mapLabelMap.end())
				{
					set<int> newMapLabelSet;
					newMapLabelSet.insert(tmpSegGroupNO);
					mapLabelMap.insert(pair< unsigned int, set<int> > 
						(tmp1stBaseAlignLoc, newMapLabelSet));
				} 
				else
				{
					(tmpIter->second).insert(tmpSegGroupNO);
				}
			}
		}

		for(mapLabelMapIter tmpMapIter = mapLabelMap.first(); 
			tmpMapIter < mapLabelMap.end(); tmpMapIter ++)
		{
			segMapRangeStart[mapLabelNum] = (*((tmpMapIter->second).first()));
			segMapRangeEnd[mapLabelNum] = (*((tmpMapIter->second).end()));
			segMapLoc[mapLabelNum] = tmpMapIter->first;
			mapLabelNum++;
		}
	}*/

	string fragInfoStr()
	{
		string fragInfoStr;

		return fragInfoStr;
	}

};


class FinalMapLabel_Info
{
public:
	int finalNorMapLabel, finalRcmMapLabel;
	
	unsigned int finalNorSegMapRangeStart[SEGMENTNUM * CANDALILOC];
	unsigned int finalNorSegMapRangeEnd[SEGMENTNUM * CANDALILOC];
	unsigned int finalNorSegMapLoc[SEGMENTNUM * CANDALILOC];

	unsigned int finalRcmSegMapRangeStart[SEGMENTNUM * CANDALILOC];
	unsigned int finalRcmSegMapRangeEnd[SEGMENTNUM * CANDALILOC];
	unsigned int finalRcmSegMapLoc[SEGMENTNUM * CANDALILOC];

	FinalMapLabel_Info()
	{}

	bool filterMapLabelCasesWithPairEndInformation_FinalMapLabelInfo(
		int mapLabel_1, unsigned int* segMapRangeStart_1, unsigned int* segMapRangeEnd_1, unsigned int* segMapLoc_1,
		int mapLabel_2, unsigned int* segMapRangeStart_2, unsigned int* segMapRangeEnd_2, unsigned int* segMapLoc_2,
		int* finalMapLabel_1, unsigned int* finalSegMapRangeStart_1, unsigned int* finalSegMapRangeEnd_1, unsigned int* finalSegMapLoc_1,
		int* finalMapLabel_2, unsigned int* finalSegMapRangeStart_2, unsigned int* finalSegMapRangeEnd_2, unsigned int* finalSegMapLoc_2,
		Index_Info* indexInfo/*,
		int* filterOutMapLabel_1, unsigned int* filterOutSegMapRangeStart_1, unsigned int* filterOutMapRangeEnd_1, unsigned int* filterOutSegMapLoc_1,
		int* filterOutMapLabel_2, unsigned int* filterOutSegMapRangeStart_2, unsigned int* filterOutMapRangeEnd_2, unsigned int* filterOutSegMapLoc_2*/)
	{
		typedef set<int> MapLabelSet; // <firstMapLabel, Vec<secondMapLabel> >
		MapLabelSet mapLabelSet_1, mapLabelSet_2;
		unsigned int tmpChrInt_1, tmpChrPosInt_1;
		unsigned int tmpChrInt_2, tmpChrPosInt_2;

		for(int tmpLabel_1 = 0; tmpLabel_1 < mapLabel_1; tmpLabel_1++)
		{
			for(int tmpLabel_2 = 0; tmpLabel_2 < mapLabel_2; tmpLabel_2++)
			{
				indexInfo->getChrLocation(segMapLoc_1[tmpLabel_1], &tmpChrInt_1, &tmpChrPosInt_1);
				indexInfo->getChrLocation(segMapLoc_2[tmpLabel_2], &tmpChrInt_2, &tmpChrPosInt_2);
				if(tmpChrInt_1 == tmpChrInt_2)
				{
					if(checkTwoIntInRange((int)tmpChrPosInt_1, (int)tmpChrPosInt_2, PE_MAP_DISTANCE))
					{
						mapLabelSet_1.insert(tmpLabel_1);
						mapLabelSet_2.insert(tmpLabel_2);
					}
					else
					{
						continue;
					}
				}
				else
				{
					continue;
				}	
			}
		}

		int tmpFinalMapLabel = 0;
		for(MapLabelSet::iterator tmpIter = mapLabelSet_1.begin(); tmpIter != mapLabelSet_1.end(); tmpIter++)
		{
			int tmpMapLabel = (*tmpIter);
			finalSegMapRangeStart_1[tmpFinalMapLabel] = segMapRangeStart_1[tmpMapLabel];
			finalSegMapRangeEnd_1[tmpFinalMapLabel] = segMapRangeEnd_1[tmpMapLabel];		
			finalSegMapLoc_1[tmpFinalMapLabel] = segMapLoc_1[tmpMapLabel];
			tmpFinalMapLabel ++;
		}

		tmpFinalMapLabel = 0;
		for(MapLabelSet::iterator tmpIter = mapLabelSet_2.begin(); tmpIter != mapLabelSet_2.end(); tmpIter++)
		{
			int tmpMapLabel = (*tmpIter);
			finalSegMapRangeStart_2[tmpFinalMapLabel] = segMapRangeStart_2[tmpMapLabel];
			finalSegMapRangeEnd_2[tmpFinalMapLabel] = segMapRangeEnd_2[tmpMapLabel];		
			finalSegMapLoc_2[tmpFinalMapLabel] = segMapLoc_2[tmpMapLabel];
			tmpFinalMapLabel ++;
		}

		(*finalMapLabel_1) = mapLabelSet_1.size();

		(*finalMapLabel_2) = mapLabelSet_2.size();

		return true;
	}
};