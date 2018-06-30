// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef BWTMAP_INFO_H
#define BWTMAP_INFO_H

#include <string>
#include "mismatch.h"
#include "splice_info.h"

using namespace std;

class Bwtmap_Info
{
public:
	Bwtmap_Info()
	{ 
			head = -1;
			tail = -1;
			incomplete_head = false;
			incomplete_tail = false;
			num_mis = 0;
	}

	~Bwtmap_Info()
	{
	}

	bool check_head(int check_seg_no) // check if this segment has been fixed with a previous segment
	{
		if(head==-1)
			return true;
		if(strand=="+" && check_seg_no>=head) // 
			return true;
		else if(strand=="-" && check_seg_no<=head)
			return true;
		return false;
	}

	bool check_tail(int check_seg_no) // check if this segment has been fixed with a later segment
	{
		if(tail==-1)
			return true;
		if(strand=="+" && check_seg_no<=tail)
			return true;
		else if(strand=="-" && check_seg_no>=tail)
			return true;
		return false;
	}

	void copy(const Bwtmap_Info& copy_info)
	{
		read_id = copy_info.read_id;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		strand = copy_info.strand;
		chrom = copy_info.chrom;
		start = copy_info.start;
		end = copy_info.end;
		other_map = copy_info.other_map;
		head = copy_info.head;
		tail = copy_info.tail;
		incomplete_head = copy_info.incomplete_head;
		incomplete_tail = copy_info.incomplete_tail;
		pair_no = copy_info.pair_no;
		splice_head.clear();
		for(size_t i = 0; i < copy_info.splice_head.size(); i++)
			splice_head.push_back(copy_info.splice_head[i]);
		splice_internal.clear();
		for(size_t i = 0; i < copy_info.splice_internal.size(); i++)
			splice_internal.push_back(copy_info.splice_internal[i]);
		splice_tail.clear();
		for(size_t i = 0; i < copy_info.splice_tail.size(); i++)
			splice_tail.push_back(copy_info.splice_tail[i]);
		num_mis = copy_info.num_mis;
	}

	
	void to_splice_info(Splice_Info& my_splice, int segment_length) //convert to Splice_Info
	{
		my_splice.start_pos = start;
		my_splice.end_pos = end;
		my_splice.chrom = chrom;
		my_splice.strand = strand;
		my_splice.start_seg_no = start_seg_no;
		my_splice.end_seg_no = end_seg_no;
		Jump_Code new_junmp_code(segment_length, "M");
		my_splice.jump_code.push_back(new_junmp_code);
	}

	void simple_assign(string readId, int startSegNo, int endSegNo, string strandId, 
		string chromId, int startPos, int endPos) //X: write for assign the class for TEST
	{
		read_id = readId;
		start_seg_no = startSegNo;
		end_seg_no = endSegNo;
		strand = strandId;
		chrom = chromId;
		start = startPos;
		end = endPos;
	}

	string read_id;
	int start_seg_no; //x:?
	int end_seg_no; //x: ?
	string strand; //xinan "+" or "-"
	string chrom;
	int start;  //Xinan: map position - start
	int end;    //Xinan: map position - end
	int other_map;
	vector<Splice_Info> splice_internal;
	vector<Splice_Info> splice_head;
	vector<Splice_Info> splice_tail;
	int num_mis;
	
	int head;
	int tail;
	int pair_no;
	bool incomplete_head;
	bool incomplete_tail;
};

#endif


