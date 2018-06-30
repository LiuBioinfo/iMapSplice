// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef MISSINGLONGSEG_INFO_H
#define MISSINGLONGSEG_INFO_H

#include <vector>

class MissingLongSeg_Info
{
private:
	int missingSegGroupNum;
	vector<int> missingSegGroupIndexInOriSegInfo;// index in Seg_Info 
	vector<unsigned int> missingSegLength;
	vector<unsigned int> missingSegAlignNum;
	vector< vector<unsigned int> > missingSegAlignLoc;

public:
	MissingLongSeg_Info()
	{
		missingSegGroupNum = 0;
	}

	string segInfoStr(Index_Info* indexInfo)
	{
		string segInfoStr;
		segInfoStr += "\nsegment Info: \n";
		for(int tmpSegGroupNO = 0; tmpSegGroupNO < missingSegGroupNum; tmpSegGroupNO++)
		{
			segInfoStr = segInfoStr + "... segment " + int_to_str(tmpSegGroupNO+1) + " len: " 
				+ int_to_str(missingSegLength[tmpSegGroupNO]) + " index_oriSegInfo: " 
				+ int_to_str(missingSegGroupIndexInOriSegInfo[tmpSegGroupNO]) + " Num: "
				+ int_to_str(missingSegAlignNum[tmpSegGroupNO]) + "\n";
			for(int tmpAlignNO = 0; tmpAlignNO < missingSegAlignNum[tmpSegGroupNO]; tmpAlignNO++)
			{
				unsigned int align_chr, align_chr_location;
				indexInfo->getChrLocation((((missingSegAlignLoc[tmpSegGroupNO])[tmpAlignNO])+1), 
					&align_chr, &align_chr_location);
				segInfoStr = segInfoStr + "\t" + int_to_str(tmpAlignNO+1)
					+ ". " + (indexInfo->returnChrNameStr(align_chr))
					+ " " + int_to_str(align_chr_location) + "\n";
			}
		}
		return segInfoStr;
	}

	int returnMissingSegGroupNum()
	{
		return missingSegGroupNum;
	}

	int returnMissingSegGroupIndexInOriSegInfo(int tmp)
	{
		return missingSegGroupIndexInOriSegInfo[tmp];
	}

	unsigned int returnMissingSegLength(int tmp)
	{
		return missingSegLength[tmp];
	}

	unsigned int returnMissingSegAlignNum(int tmp)
	{
		return missingSegAlignNum[tmp];
	}

	unsigned int returnMissingSegAlignLoc(int tmpGroupNO, int tmpIndexInGroup)
	{
		return (missingSegAlignLoc[tmpGroupNO])[tmpIndexInGroup];
	}
	void pushBackNewMissingLongSegGroup(
		unsigned int interval_begin_ori, 
		unsigned int interval_end_ori, 
		unsigned int currentSegLength,
		int index_oriSegInfo, 
		unsigned int* sa, Index_Info* indexInfo)
	{
		missingSegGroupNum ++;
		missingSegGroupIndexInOriSegInfo.push_back(index_oriSegInfo);
		missingSegLength.push_back(currentSegLength);
		missingSegAlignNum.push_back(interval_end_ori - interval_begin_ori + 1);
		vector<unsigned int> tmpAlignLocVec;
		for(unsigned int tmp = interval_begin_ori; tmp <= interval_end_ori; tmp++)
		{
			unsigned int tmpPosInWholeGenome = sa[tmp] + 1;
			//unsigned int tmpChr, tmpChrPos;
			//indexInfo->getChrLocation(tmpPosInWholeGenome, &tmpChr, &tmpChrPos);	
			//string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
			//cout << "Pos." << tmp - interval_begin_ori + 1 << ": " << tmpChrName << " " << tmpChrPos << endl;
			tmpAlignLocVec.push_back(tmpPosInWholeGenome);
		}		
		missingSegAlignLoc.push_back(tmpAlignLocVec);
	}

	bool checkSeg_missingLong_orNot(
		unsigned int interval_begin_ori, 
		unsigned int interval_end_ori, 
		unsigned int interval_begin, 
		unsigned int interval_end, 
		int currentSegmentLength)
	{	
		if(interval_end < interval_begin)
		{
			return false;
		}

		if(
			((interval_end - interval_begin) < (interval_end_ori - interval_begin_ori))
			&& (currentSegmentLength > 30)
			&& ((interval_end_ori - interval_begin_ori + 1) <= CANDALILOC)
			)
		{
			return true;
		}
		else
			return false;
	}
					
};

#endif