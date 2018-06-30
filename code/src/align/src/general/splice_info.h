// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SPLICE_INFO_H
#define SPLICE_INFO_H

#include <string>
#include <string.h>
#include "mismatch.h"
using namespace std;

class Bwtmap_Info;


class Jump_Code
{
public:
	int len; // length fo current feature
	string type; //feature type, M/N/I/D/S

	Jump_Code(int _len, string _type)
	{
		len=_len;
		type=_type;
	}
	
	~Jump_Code()
	{}
	
	string toString()
	{
		return int_to_str(len) + type;
	}
	
	string toString(int add_value)
	{
		return int_to_str(len + add_value) + type;
	}
	
	string toStringConvertDeletion(int max_del_len)
	{
		if(type == "N" && len <= max_del_len)
		{
			return int_to_str(len) + "D";
		}
		return int_to_str(len) + type;
	}
};

class Splice_Info
{
//private:
public:
	vector<Jump_Code> jump_code; 
	vector<Jump_Code> final_jump_code;
	vector<string> junc_flank_seq;
	string chrom;
	string strand;
	int start_pos;
	int end_pos;
	int buffer_len;
	int mapped_len;
	int start_seg_no;
	int end_seg_no;
	size_t start_contig;
	size_t end_contig;

public:
	Splice_Info()
	{ 
		mapped_len = -1;	
	}

	Splice_Info(const Splice_Info& copy_info)
	{
		start_pos = copy_info.start_pos;
		end_pos = copy_info.end_pos;
		start_contig = copy_info.start_contig;
		end_contig = copy_info.end_contig;
		chrom = copy_info.chrom;
		strand = copy_info.strand;
		buffer_len = copy_info.buffer_len;
		mapped_len = copy_info.mapped_len;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		for(size_t i = 0; i< copy_info.jump_code.size(); i++)
		{
			jump_code.push_back(copy_info.jump_code[i]);
		}
		for(size_t i = 0; i< copy_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(copy_info.junc_flank_seq[i]);
		}
	}

	void addOriginalJumpCode_incompleteHead(vector<Jump_Code>& oriAlignInfoJumpCodeVec)
	{
		for(int tmp = 2; tmp < oriAlignInfoJumpCodeVec.size(); tmp++) // 1st Jump_code: Soft clipping; 2nd jump_code: midPartSeg Match
		{
			final_jump_code.push_back(oriAlignInfoJumpCodeVec[tmp]);
		}
	}

	void addOriginalJumpCode_incompleteTail(vector<Jump_Code>& oriAlignInfoJumpCodeVec)
	{
		vector<Jump_Code> tmpFinalJumpCodeVec;
		tmpFinalJumpCodeVec.clear();
		for(int tmp = 0; tmp < final_jump_code.size(); tmp++) // 1st Jump_code: Soft clipping; 2nd jump_code: midPartSeg Match
		{
			tmpFinalJumpCodeVec.push_back(final_jump_code[tmp]);
		}
		final_jump_code.clear();
		for(int tmp = 0; tmp < oriAlignInfoJumpCodeVec.size()-2; tmp++)
		{
			final_jump_code.push_back(oriAlignInfoJumpCodeVec[tmp]);
		}
		for(int tmp = 0; tmp < tmpFinalJumpCodeVec.size(); tmp++)
		{
			final_jump_code.push_back(tmpFinalJumpCodeVec[tmp]);
		}
	}

	void copy(const Splice_Info& copy_info)
	{
		start_pos = copy_info.start_pos;
		end_pos = copy_info.end_pos;
		start_contig = copy_info.start_contig;
		end_contig = copy_info.end_contig;
		chrom = copy_info.chrom;
		strand = copy_info.strand;
		buffer_len = copy_info.buffer_len;
		mapped_len = copy_info.mapped_len;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		for(size_t i = 0; i< copy_info.jump_code.size(); i++)
		{
			jump_code.push_back(copy_info.jump_code[i]);
		}
		for(size_t i = 0; i< copy_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(copy_info.junc_flank_seq[i]);
		}
	}

	void cleanAndCopy(const Splice_Info& copy_info)
	{
		start_pos = copy_info.start_pos;
		end_pos = copy_info.end_pos;
		start_contig = copy_info.start_contig;
		end_contig = copy_info.end_contig;
		chrom = copy_info.chrom;
		strand = copy_info.strand;
		buffer_len = copy_info.buffer_len;
		mapped_len = copy_info.mapped_len;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		jump_code.clear();
		junc_flank_seq.clear();
		for(size_t i = 0; i< copy_info.jump_code.size(); i++)
		{
			jump_code.push_back(copy_info.jump_code[i]);
		}
		for(size_t i = 0; i< copy_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(copy_info.junc_flank_seq[i]);
		}
	}
	
	void Get_Extra_Info(string& chrom_seq)
	{
		mapped_len = 0;
		int	end_tmp = start_pos - 1;
		for(size_t i = 0; i < jump_code.size(); i++)
		{
			if(jump_code[i].type == "M")
			{
				end_tmp += jump_code[i].len;
				mapped_len += jump_code[i].len;
			}
			else if(jump_code[i].type == "N")
			{
				int junc_start = end_tmp + 1;
				end_tmp += jump_code[i].len;
				int junc_end =  end_tmp;
				string flank_seq = chrom_seq.substr(junc_start, 2) + chrom_seq.substr(junc_end - 1, 2);
				junc_flank_seq.push_back(flank_seq);
			}
			else   // "I"
			{
				mapped_len += jump_code[i].len;
			}
		}	
	}
	
	void set_value(string line, bool pairend)
	{
	}

	void set_value_single(string line)
	{
	}

	void set_value_pairend(string line)
	{
	}

	~Splice_Info()
	{
	}

	int GetMapLen() // get total mapped(M) length
	{
	  if(mapped_len == -1)
	  {
	  	mapped_len = 0;
	  	for(size_t i = 0; i < jump_code.size(); i++)
	  	{
	  		if(jump_code[i].type == "M" || jump_code[i].type == "I")
	  			mapped_len += jump_code[i].len;
	  	}
		}
	  return mapped_len;
	}

	int countFinalMapLength()
	{
		int mapped_len = 0;
	  	for(size_t i = 0; i < final_jump_code.size(); i++)
	  	{
	  		if(final_jump_code[i].type == "M" || final_jump_code[i].type == "I" || final_jump_code[i].type == "S" ) 
	  			mapped_len += final_jump_code[i].len;
	  	}
	  	return mapped_len;
	}

	int Count_Gap() //get total number of gaps
	{
		int num_gap = 0;
		for(size_t i = 0; i < jump_code.size(); i++)
	  {
	  	if(jump_code[i].type != "M")
	  		num_gap ++;
	  }
	  return num_gap;
	}

	string printJumpCode()
	{
		//cout << "start printJumpCode" << endl;
		string jumpCodeStr;
		string tmpStr;
		//cout << jump_code.size() << endl;
		for(size_t i=0; i < jump_code.size(); i++)
		{
			tmpStr = jump_code[i].toString();
			jumpCodeStr = jumpCodeStr + tmpStr;
		}
		//cout << "finish printJumpCode" << endl;
		return jumpCodeStr;
	}

	string printFinalJumpCode()
	{
		string jumpCodeStr;
		string tmpStr;
		for(size_t i=0; i < final_jump_code.size(); i++)
		{
			tmpStr = final_jump_code[i].toString();
			jumpCodeStr = jumpCodeStr + tmpStr;
		}
		return jumpCodeStr;
	}

	string printOtherCigarStringForShortHead() // except the first and second jumpCode
	{
		string jumpCodeStr;
		string tmpStr;
		if(final_jump_code.size() == 2) 
			return "*";
		for(size_t i=2; i < final_jump_code.size(); i++)
		{
			tmpStr = final_jump_code[i].toString();
			jumpCodeStr = jumpCodeStr + tmpStr;
		}
		return jumpCodeStr;		
	}

	string printOtherCigarStringForShortTail()
	{
		string jumpCodeStr;
		string tmpStr;
		if(final_jump_code.size() == 2)
			return "*";
		for(size_t i = 0; i < final_jump_code.size()-2; i++)
		{
			tmpStr = final_jump_code[i].toString();
			jumpCodeStr = jumpCodeStr + tmpStr;
		}
		return jumpCodeStr;
	}

	void appendJumpCode(Splice_Info* anotherSpliceInfo)
	{
		//cout << "here 1" << endl;
		//cout << jump_code.size() << endl;		
		for(int tmp = 0; tmp < anotherSpliceInfo->jump_code.size(); tmp++)
		{
			//cout << "tmp = " << tmp << endl;
			jump_code.push_back(anotherSpliceInfo->jump_code[tmp]);
		}
		//cout << "here 2" << endl;
		//cout << jump_code.size() << endl;
	}

	void appendJumpCode(Splice_Info* anotherSpliceInfo1, Splice_Info* anotherSpliceInfo2)
	{
		//cout << "anotherSpliceInfo1->jump_code.size() = " 
		//<< anotherSpliceInfo1->jump_code.size() << endl;
		for(int tmp = 0; tmp < anotherSpliceInfo1->jump_code.size(); tmp++)
		{
			//cout << "insert!" << endl;
			jump_code.push_back(anotherSpliceInfo1->jump_code[tmp]);
		}
		//cout << "to here"<< endl;
		//cout << "anotherSpliceInfo2->jump_code.size() = "
		//<< anotherSpliceInfo2->jump_code.size();		
		for(int tmp = 0; tmp < anotherSpliceInfo2->jump_code.size(); tmp++)
		{
			jump_code.push_back(anotherSpliceInfo2->jump_code[tmp]);
		}
	}

	void appendJumpCode(Splice_Info* anotherSpliceInfo1, Splice_Info* anotherSpliceInfo2, Splice_Info* anotherSpliceInfo3)
	{
		for(int tmp = 0; tmp < anotherSpliceInfo1->jump_code.size(); tmp++)
		{
			jump_code.push_back(anotherSpliceInfo1->jump_code[tmp]);
		}
		for(int tmp = 0; tmp < anotherSpliceInfo2->jump_code.size(); tmp++)
		{
			jump_code.push_back(anotherSpliceInfo2->jump_code[tmp]);
		}
		for(int tmp = 0; tmp < anotherSpliceInfo3->jump_code.size(); tmp++)
		{
			jump_code.push_back(anotherSpliceInfo3->jump_code[tmp]);
		}
	}
	void getFinalJumpCode()
	{
		int jumpCodeSize = jump_code.size();
		if(jumpCodeSize == 0)
			return;

		final_jump_code.push_back(jump_code[0]);
		string currentJumpCodeType = jump_code[0].type;
		//cout << "tmpJumpCode: " << jump_code[0].len << "," << jump_code[0].type << endl;
		for(int tmp = 1; tmp < jump_code.size(); tmp++)
		{
			//cout << "tmpJumpCode: " << jump_code[tmp].len << "," << jump_code[tmp].type << endl;			
			string tmpJumpCodeType = jump_code[tmp].type;
			int final_jump_code_size = final_jump_code.size();
			currentJumpCodeType = final_jump_code[final_jump_code_size - 1].type;
			if(tmpJumpCodeType == currentJumpCodeType)
			{	
				int currentJumpCodeLength = final_jump_code[final_jump_code_size - 1].len;
				int newJumpCodeLength = currentJumpCodeLength + jump_code[tmp].len;
				final_jump_code[final_jump_code_size - 1].len = newJumpCodeLength;
			}
			else
			{
				final_jump_code.push_back(jump_code[tmp]);
			}
		}
		return;
	}

	void getFinalJumpCode_old()
	{
		//Jump_Code currentJumpCode;
		int jumpCodeSize = jump_code.size();
		//cout << "jumpCodeSize" << jumpCodeSize << endl;
		if (jumpCodeSize == 0)
		{
			return;
		}
		else if(jumpCodeSize == 1)
		{
			if((jump_code[0].len <= 0)||(jump_code[0].type != "M"))
			{
				cout << "invalid jump code " << endl;
				return;
			}
			else
			{	
				final_jump_code.push_back(jump_code[0]);
				return;
			}
		}
		else
		{	
			Jump_Code currentJumpCode = jump_code[0];
			Jump_Code nextJumpCode = jump_code[1];
			for(int tmp = 0; tmp < jump_code.size() - 1; tmp++)
			{
				int tmp2 = tmp+1;
				nextJumpCode = jump_code[tmp2];
				if(currentJumpCode.type != nextJumpCode.type)
				{	
					final_jump_code.push_back(currentJumpCode);
					currentJumpCode = nextJumpCode;
				}
				else
				{
					currentJumpCode.len = currentJumpCode.len + nextJumpCode.len;
				}
			}
			final_jump_code.push_back(currentJumpCode);
		}
	}

	bool shortHeadSoftClippingOrNot()
	{
		//cout << "shortHeadSoftClippingOrNot" << endl;
		if((final_jump_code[0].len < 10) && 
			(final_jump_code[0].type == "S"))
		{
			//cout << "true" << endl;
			return true;
		}
		else 
		{	
			//cout << "false" << endl;
			return false;
		}
	}

	bool shortTailSoftClippingOrNot()
	{
		//cout << "shortTailSoftClippingOrNot" << endl;
		int jumpCodeSize = final_jump_code.size();
		if((final_jump_code[jumpCodeSize-1].len < 10) 
			&& (final_jump_code[jumpCodeSize-1].type == "S"))
		{
			//cout << "true" << endl;
			return true;
		}
		else 
		{
			//cout << "false" << endl;
			return false;
		}
	}

	bool longHeadSoftClippingOrNot()
	{
		if((final_jump_code[0].len >= 10) && 
			(final_jump_code[0].type == "S"))
		{
			//cout << "true" << endl;
			return true;
		}
		else 
		{	
			//cout << "false" << endl;
			return false;
		}		
	}

	bool longTailSoftClippingOrNot()
	{
		int jumpCodeSize = final_jump_code.size();
		if((final_jump_code[jumpCodeSize-1].len >= 10) 
			&& (final_jump_code[jumpCodeSize-1].type == "S"))
		{
			//cout << "true" << endl;
			return true;
		}
		else 
		{
			//cout << "false" << endl;
			return false;
		}
	}

	string getFirstJumpCodeLengthStr()
	{
		return int_to_str(final_jump_code[0].len);
	}

	string getSecondJumpCodeLengthStr()
	{
		return int_to_str(final_jump_code[1].len);
	}

	string getPenultimateJumpCodeLengthStr()
	{
		int jumpCodeSize = final_jump_code.size();
		return int_to_str(final_jump_code[jumpCodeSize-2].len);
	}

	string getEndJumpCodeLengthStr()
	{
		int jumpCodeSize = final_jump_code.size();
		return int_to_str(final_jump_code[jumpCodeSize-1].len);
	}

	unsigned int getLastMatchSegCorreAlignPosForShortTail(unsigned int mapPos, int readLength) // assume the first matched base of read's map pos is 1
	{
		unsigned int pos = 0;
		int finalJumpCodeSize = final_jump_code.size();
		for(int tmp = 0; tmp < finalJumpCodeSize-1; tmp++)
		{
			if ((final_jump_code[tmp].type == "S")||(final_jump_code[tmp].type == "I"))
			{
				continue;
			}
			else if ((final_jump_code[tmp].type == "M")||(final_jump_code[tmp].type == "N")||(final_jump_code[tmp].type == "D"))
			{
				pos += final_jump_code[tmp].len;
			}
			else //if (final_jump_code[tmp].type == "N")
			{
				cout << "error in getLastMatchSegCorreAlignPosForShortTail function of splice_info.h: unknown jumpcode type!!!" << endl;
			}
		}
		pos = mapPos + pos - 1;
		int readLengthExceptTail = readLength - final_jump_code[finalJumpCodeSize-1].len;
		pos = pos - readLengthExceptTail + 1;

		return pos;
	}

	int getEndBaseMapPos_jump_code(int startMapPos)
	{
		int tmpMapPos = 0;
		int jumpCodeSize = jump_code.size();
		for(int tmp = 0; tmp < jumpCodeSize; tmp++)
		{
			if((jump_code[tmp].type == "S")||(jump_code[tmp].type == "I"))
			{
				continue;
			}
			else if((jump_code[tmp].type == "M")||(jump_code[tmp].type == "N")||(jump_code[tmp].type == "D"))
			{
				tmpMapPos += jump_code[tmp].len;
			}
			else
			{
				cout << "invalid jumpCode type, error in splice_info.h: getEndBaseMapPos_jump_code !" << endl;
			}
		}
		int endBaseMapPos = startMapPos + tmpMapPos - 1;
		return endBaseMapPos;
	}




	string printSecondLevelFinalCigarStringForLongHead(int secondJumpCodeLength, const string& otherCigarString)
	{
		int jumpCodeNum = final_jump_code.size();
		final_jump_code[jumpCodeNum-1].len += secondJumpCodeLength;

		string jumpCodeStr;
		string tmpStr;
		for(size_t i=0; i < final_jump_code.size(); i++)
		{
			tmpStr = final_jump_code[i].toString();
			jumpCodeStr = jumpCodeStr + tmpStr;
		}
		string finalCigarString;
		if(otherCigarString == "*")
		{	
			finalCigarString = jumpCodeStr;
		}
		else
		{
			finalCigarString = jumpCodeStr + otherCigarString;
		}
		return finalCigarString; 
	}

	bool noNullNegativeMatchJumpCode(int secondJumpCodeLength)
	{
		int jumpCodeNum = final_jump_code.size();
		for(int tmp=0; tmp<jumpCodeNum-1; tmp++)
		{
			if(final_jump_code[tmp].len <= 0)
				return false;
		}

		int lastMatchJumpCodeLengthInHead = final_jump_code[jumpCodeNum-1].len + secondJumpCodeLength;
		if(lastMatchJumpCodeLengthInHead > 0)
			return true;
		else 
			return false;
	}

	bool allFinalJumpCodeValid()
	{
		#ifdef DETECT_CIRCULAR_RNA
		bool allJumpCodeValidBool = true;
		for(int tmp = 0; tmp < final_jump_code.size(); tmp ++)
		{
			if(final_jump_code[tmp].len <= 0)
			{
				if(final_jump_code[tmp].type != "N")
				{
					allJumpCodeValidBool = false;
				}
			}
		}
		return allJumpCodeValidBool;
		#else
		bool allJumpCodeValidBool = true;
		for(int tmp = 0; tmp < final_jump_code.size(); tmp ++)
		{
			if(final_jump_code[tmp].len <= 0)
			{
				allJumpCodeValidBool = false;
			}
		}
		return allJumpCodeValidBool;
		#endif
	}

	bool allFinalJumpCodeValid_cirRNA()
	{
		bool allJumpCodeValidBool = true;
		for(int tmp = 0; tmp < final_jump_code.size(); tmp ++)
		{
			if(final_jump_code[tmp].len <= 0)
			{
				if(final_jump_code[tmp].type != "N")
				{
					allJumpCodeValidBool = false;
				}
			}
		}

		return allJumpCodeValidBool;
	}

};


class Fusion_Splice
{
public:
	Splice_Info first_splice;
	Splice_Info second_splice;
	string flank_seq;

	Fusion_Splice()
	{}
	
	Fusion_Splice (Splice_Info& my_first_splice, Splice_Info& my_second_splice, string& my_flank_seq)
	{
	 	first_splice = my_first_splice;
		second_splice = my_second_splice;
		flank_seq = my_flank_seq; 
	}
	
	Fusion_Splice (const Fusion_Splice& fusion_splice_info)
	{
		first_splice = fusion_splice_info.first_splice;
		second_splice = fusion_splice_info.second_splice;
		flank_seq = fusion_splice_info.flank_seq;
	}
	
	~Fusion_Splice ()
	{}
};


#endif

