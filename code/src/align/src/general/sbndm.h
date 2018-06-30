// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include <math.h>
#include <stdint.h>
#include "bwtmap_info.h"

//typedef unit16_t unsigned int; 

using namespace std;

const int MASK_BIT = 0x0707;  
//for big endian
/*const int AA = ((65 << 8) | 65) & MASK_BIT;
const int AG = ((65 << 8) | 71) & MASK_BIT;
const int AT = ((65 << 8) | 84) & MASK_BIT;
const int AC = ((65 << 8) | 67) & MASK_BIT;
const int AN = ((65 << 8) | 78) & MASK_BIT;
const int GA = ((71 << 8) | 65) & MASK_BIT;
const int GG = ((71 << 8) | 71) & MASK_BIT;
const int GT = ((71 << 8) | 84) & MASK_BIT;
const int GC = ((71 << 8) | 67) & MASK_BIT;
const int GN = ((71 << 8) | 78) & MASK_BIT;
const int TA = ((84 << 8) | 65) & MASK_BIT;
const int TG = ((84 << 8) | 71) & MASK_BIT;
const int TT = ((84 << 8) | 84) & MASK_BIT;
const int TC = ((84 << 8) | 67) & MASK_BIT;
const int TN = ((84 << 8) | 78) & MASK_BIT;
const int CA = ((67 << 8) | 65) & MASK_BIT;
const int CG = ((67 << 8) | 71) & MASK_BIT;
const int CT = ((67 << 8) | 84) & MASK_BIT;
const int CC = ((67 << 8) | 67) & MASK_BIT;
const int CN = ((67 << 8) | 78) & MASK_BIT;
const int NA = ((78 << 8) | 65) & MASK_BIT;
const int NG = ((78 << 8) | 71) & MASK_BIT;
const int NT = ((78 << 8) | 84) & MASK_BIT;
const int NC = ((78 << 8) | 67) & MASK_BIT;
const int NN = ((78 << 8) | 78) & MASK_BIT;
*/


//for little endian
const int AA = ((65 << 8) | 65);
const int AG = ((71 << 8) | 65);
const int AT = ((84 << 8) | 65);
const int AC = ((67 << 8) | 65);
const int AN = ((78 << 8) | 65);
const int GA = ((65 << 8) | 71);
const int GG = ((71 << 8) | 71);
const int GT = ((84 << 8) | 71);
const int GC = ((67 << 8) | 71);
const int GN = ((78 << 8) | 71);
const int TA = ((65 << 8) | 84);
const int TG = ((71 << 8) | 84);
const int TT = ((84 << 8) | 84);
const int TC = ((67 << 8) | 84);
const int TN = ((78 << 8) | 84);
const int CA = ((65 << 8) | 67);
const int CG = ((71 << 8) | 67);
const int CT = ((84 << 8) | 67);
const int CC = ((67 << 8) | 67);
const int CN = ((78 << 8) | 67);
const int NA = ((65 << 8) | 78);
const int NG = ((71 << 8) | 78);
const int NT = ((84 << 8) | 78);
const int NC = ((67 << 8) | 78);
const int NN = ((78 << 8) | 78);


inline bool check_pattern(string& pattern) // check if pattern(anchor) contains N
{
	if(pattern.find("N") != string::npos)
		return false;
	else
		return true;
}


class SBNDM_FORWARD // class of pattern(anchor) searching in forward direction
{
public:
	SBNDM_FORWARD()
	{
		alphabet_size = 5;
		q_size = 4;
	}

	~SBNDM_FORWARD()
	{
	}

	inline void compute_bit_vector(string& pattern, unsigned int* B) 
	{
		B['A'] = 0;		
		B['C'] = 0;	
		B['G'] = 0;	
		B['T'] = 0;	
		B['N'] = 0;	
		unsigned int s = 1;
		for (int i = pattern_size - 1; i >= 0; i--)
		{
			B[ pattern[i] ] |= s;
			s <<= 1;
		}
	}

	inline void compute_shift(string& pattern, unsigned int* B, int& s0)
	{
		unsigned int S = B[ pattern[pattern_size - 1] ];
		s0 = pattern_size;
		for(int i = pattern_size - 2; i >= 0; i--)
		{
			if((S & check_highest_bit) != 0)
				s0 = i + 1;
			S = (S << 1) & B[ pattern[i] ];
		}
	}

	inline void compute_2gram_table(unsigned int* two_gram, unsigned int *B)
	{
		two_gram[AA] = B['A'] & (B['A'] << 1);
		two_gram[AG] = B['A'] & (B['G'] << 1);
		two_gram[AT] = B['A'] & (B['T'] << 1);
		two_gram[AC] = B['A'] & (B['C'] << 1);
		two_gram[AN] = B['A'] & (B['N'] << 1);
		two_gram[GA] = B['G'] & (B['A'] << 1);
		two_gram[GG] = B['G'] & (B['G'] << 1);
		two_gram[GT] = B['G'] & (B['T'] << 1);
		two_gram[GC] = B['G'] & (B['C'] << 1);
		two_gram[GN] = B['G'] & (B['N'] << 1);
		two_gram[TA] = B['T'] & (B['A'] << 1);
		two_gram[TG] = B['T'] & (B['G'] << 1);
		two_gram[TT] = B['T'] & (B['T'] << 1);
		two_gram[TC] = B['T'] & (B['C'] << 1);
		two_gram[TN] = B['T'] & (B['N'] << 1);
		two_gram[CA] = B['C'] & (B['A'] << 1);
		two_gram[CG] = B['C'] & (B['G'] << 1);
		two_gram[CT] = B['C'] & (B['T'] << 1);
		two_gram[CC] = B['C'] & (B['C'] << 1);
		two_gram[CN] = B['C'] & (B['N'] << 1);
		two_gram[NA] = B['N'] & (B['A'] << 1);
		two_gram[NG] = B['N'] & (B['G'] << 1);
		two_gram[NT] = B['N'] & (B['T'] << 1);
		two_gram[NC] = B['N'] & (B['C'] << 1);
		two_gram[NN] = B['N'] & (B['N'] << 1);
	}
	
	inline unsigned int q_gram(char* text_ptr, unsigned int* two_gram, int index)
	{
		uint16_t* text_ptr16 = (uint16_t*)(text_ptr + index); // X: unit16_t was not declared 08/06
		//unsigned int * text_ptr16 = (unsigned int*)(text_ptr + index);
		return two_gram[ text_ptr16[0] ] & ((two_gram[ text_ptr16[1] ] << 2));
	}

	// search the next pattern(anchor) match in text(genome), worker function called by "search_next_simple()"
	bool search_next_SBNDM(string& text, int& start_pos, int& end_pos, unsigned int* B, unsigned int* two_gram, int& pos) 
	{
		/* Searching phase */
		int i = start_pos + q_point, j = 0;
		unsigned int d = 0;
		while (i <= end_pos - q_size + 1)
		{
			d = q_gram(text_ptr, two_gram, i);
			if(d != 0)
			{
				j = i - q_point - 1;
				do
				{
					i--;
					d = (d << 1) & B[ text[i] ];
				}while(d != 0);
				if(j == i)
				{
					pos = j + 1;
					return true;
				}
			}
			i += q_point + 1;
		}
		return false;
	}

	/*bool search_next(string& text, int& start_pos_forward, int& start_pos_reverse, int& end_pos, int &pos, bool& map_fw)
	{
		if(search_fw && found_fw && update_fw)
		{
			found_fw = search_SBNDM(text, start_pos_forward, end_pos, B_forward, two_gram_forward, current_forward_anchor);
			start_pos_forward = current_forward_anchor + s0_forward;
			update_fw = false;
		}
		if(search_rc && found_rc && update_rc)
		{
			found_rc = search_SBNDM(text, start_pos_reverse, end_pos, B_reverse, two_gram_reverse, current_reverse_anchor);
			start_pos_reverse = current_reverse_anchor + s0_reverse;
			update_rc = false;
		}
		int report = 0;
		if(found_fw && found_rc)
		{
			if(current_forward_anchor <= current_reverse_anchor)	
				report = 1;
			else
				report = -1;
		}
		else if(found_fw)
			report = 1;
		else if(found_rc)
			report = -1;
		if(report == 1)
		{
			pos = current_forward_anchor;
			map_fw = true;
			update_fw = true;
			return true;
		}
		else if(report == -1)	
		{
			pos = current_reverse_anchor;
			map_fw = false;
			update_rc = true;
			return true;
		}	
		else
			return false;
	}*/
	
	// actual anchor searching function called by MapSplice, search both forward and reverse strand, find the next match
	bool search_next_simple(string& text, int& start_pos_forward, int& start_pos_reverse, int& end_pos, int &pos, bool& map_fw) 
	{  
		if(search_fw)
		{
			found_fw = search_next_SBNDM(text, start_pos_forward, end_pos, B_forward, two_gram_forward, current_forward_anchor);
			start_pos_forward = current_forward_anchor + s0_forward;
			pos = current_forward_anchor;
			cur_hits ++;
			return (found_fw && cur_hits <= max_hits);	
		}
		else if(search_rc)
		{
			found_rc = search_next_SBNDM(text, start_pos_reverse, end_pos, B_reverse, two_gram_reverse, current_reverse_anchor);
			start_pos_reverse = current_reverse_anchor + s0_reverse;
			pos = current_reverse_anchor;
			cur_hits ++;
			return (found_rc && cur_hits <= max_hits);
		}
		else
			return false;
	}
	
	// find all the pattern matches in a given range of text, worker function called by "find_all_match()"
	bool search_all_SBNDM(string& text, int& start_pos, int& end_pos, unsigned int* B, unsigned int* two_gram, int& s0, vector<Bwtmap_Info>& bwt_vector, string read_id, string chrom, string strand, int pair_no, int seg_no)
	{
		/* Searching phase */
		int i = start_pos + q_point, j = 0;
		unsigned int d = 0;
		while (i <= end_pos - q_size + 1)
		{
			d = q_gram(text_ptr, two_gram, i);
			if(d != 0)
			{
				j = i - q_point - 1;
				do
				{
					i--;
					d = (d << 1) & B[ text[i] ];
				}while(d != 0);
				if(j == i)
				{
					Bwtmap_Info new_info;
					new_info.start = j + 1;
					new_info.end = new_info.start + pattern_size - 1;
					new_info.read_id = read_id;
					new_info.chrom = chrom;
					new_info.strand = strand;
					new_info.start_seg_no = seg_no;
					new_info.end_seg_no = seg_no;
					new_info.pair_no = pair_no;
					bwt_vector.push_back(new_info);
					i += s0;
				}
			}
			i += q_point + 1;
		}
		return false;
	}
	
	// actual function called by MapSplice to find all the pattern matches in a given range of text
	void find_all_match(string& text, int& start_pos, int& end_pos, vector<Bwtmap_Info>& bwt_vector, string read_id, string chrom, int pair_no, int seg_no)
	{
		if(search_fw)
			search_all_SBNDM(text, start_pos, end_pos, B_forward, two_gram_forward, s0_forward, bwt_vector, read_id, chrom, "+", pair_no, seg_no);
		else if(search_rc)
			search_all_SBNDM(text, start_pos, end_pos, B_reverse, two_gram_reverse, s0_reverse, bwt_vector, read_id, chrom, "-", pair_no, seg_no);
	}
	
	bool set(string& pattern, string& text, bool _search_fw, bool _search_rc, int _max_hits) //preprocessing of pattern
	{
		/* Pre processing */
		if(!check_pattern(pattern))
			return false;
		pattern_size = (int)(pattern.length());
		check_highest_bit = 1 << (pattern_size - 1);  
		q_point = pattern_size - q_size;
		text_ptr = const_cast<char*> (text.data());
		search_fw = _search_fw;
		search_rc = _search_rc;
		found_fw = search_fw;
		found_rc = search_rc;
		update_fw = true;
		update_rc = true;
		current_forward_anchor = -1;
		current_reverse_anchor = -1;
		max_hits = _max_hits;
		cur_hits = 0;
		if(search_fw)
		{
			anchor_forward = pattern;
			compute_bit_vector(anchor_forward, B_forward);
			compute_shift(anchor_forward, B_forward, s0_forward);
			compute_2gram_table(two_gram_forward, B_forward);
		}
 		if(search_rc)
 		{
 			anchor_reverse = revcomp(pattern);
			compute_bit_vector(anchor_reverse, B_reverse);
			compute_shift(anchor_reverse, B_reverse, s0_reverse);
			compute_2gram_table(two_gram_reverse, B_reverse);
 		}
		return true;
	}

private:
	char* text_ptr;
	int pattern_size;
	int text_size;
	int alphabet_size;
	int q_size;
	int q_point;
	unsigned int check_highest_bit;
	string anchor_forward;
	string anchor_reverse;
	unsigned int B_forward[256];
	unsigned int B_reverse[256];
	int s0_forward;
	int s0_reverse;
	unsigned int two_gram_forward[21589];
	unsigned int two_gram_reverse[21589];
	int current_forward_anchor;
	int current_reverse_anchor;
	bool search_fw;
	bool search_rc;
	bool found_fw;
	bool found_rc;
	bool update_fw;
	bool update_rc;
	int max_hits;
	int cur_hits;
};





class SBNDM_BACKWARD  // class of pattern(anchor) searching in backward direction
{
public:
	SBNDM_BACKWARD()
	{
		alphabet_size = 5;
		q_size = 4;
	}

	~SBNDM_BACKWARD()
	{
	}

	inline void compute_bit_vector(string& pattern, unsigned int* B)
	{
		B['A'] = 0;		
		B['C'] = 0;	
		B['G'] = 0;	
		B['T'] = 0;	
		B['N'] = 0;	
		unsigned int s = 1;
		for (int i = pattern_size - 1; i >= 0; i--)
		{
			B[ pattern[i] ] |= s;
			s <<= 1;
		}
	}

	inline void compute_shift(string& pattern, unsigned int* B, int& s0)
	{
		unsigned int S = B[ pattern[0] ];
		s0 = pattern_size;
		for(int i = 1; i < pattern_size; i++)
		{
			if((S & check_lowest_bit) != 0)
				s0 = pattern_size - i;
			S = (S >> 1) & B[ pattern[i] ];
		}
	}

	inline void compute_2gram_table(unsigned int* two_gram, unsigned int *B)
	{
		two_gram[AA] = (B['A'] >> 1) & B['A'];
		two_gram[AG] = (B['A'] >> 1) & B['G'];
		two_gram[AT] = (B['A'] >> 1) & B['T'];
		two_gram[AC] = (B['A'] >> 1) & B['C'];
		two_gram[AN] = (B['A'] >> 1) & B['N'];
		two_gram[GA] = (B['G'] >> 1) & B['A'];
		two_gram[GG] = (B['G'] >> 1) & B['G'];
		two_gram[GT] = (B['G'] >> 1) & B['T'];
		two_gram[GC] = (B['G'] >> 1) & B['C'];
		two_gram[GN] = (B['G'] >> 1) & B['N'];
		two_gram[TA] = (B['T'] >> 1) & B['A'];
		two_gram[TG] = (B['T'] >> 1) & B['G'];
		two_gram[TT] = (B['T'] >> 1) & B['T'];
		two_gram[TC] = (B['T'] >> 1) & B['C'];
		two_gram[TN] = (B['T'] >> 1) & B['N'];
		two_gram[CA] = (B['C'] >> 1) & B['A'];
		two_gram[CG] = (B['C'] >> 1) & B['G'];
		two_gram[CT] = (B['C'] >> 1) & B['T'];
		two_gram[CC] = (B['C'] >> 1) & B['C'];
		two_gram[CN] = (B['C'] >> 1) & B['N'];
		two_gram[NA] = (B['N'] >> 1) & B['A'];
		two_gram[NG] = (B['N'] >> 1) & B['G'];
		two_gram[NT] = (B['N'] >> 1) & B['T'];
		two_gram[NC] = (B['N'] >> 1) & B['C'];
		two_gram[NN] = (B['N'] >> 1) & B['N'];
	}
	
	inline unsigned int q_gram(char* text_ptr, unsigned int* two_gram, int index)
	{
		uint16_t* text_ptr16 = (uint16_t*)(text_ptr + index - q_size + 1); 
		return two_gram[ text_ptr16[1] ] & ((two_gram[ text_ptr16[0] ] >> 2));
	}

// search the next pattern(anchor) match in text(genome), worker function called by "search_next_simple()"
	bool search_next_SBNDM(string& text, int& start_pos, int& end_pos, unsigned int* B, unsigned int* two_gram, int& pos)
	{
		/* Searching phase */

		int i = end_pos - q_point, j = 0;
		unsigned int d = 0;
		while (i >= start_pos + q_size - 1)
		{   
			d = q_gram(text_ptr, two_gram, i);
			if(d != 0)
			{
				j = i + q_point + 1;
				do
				{
					i++;
					d = (d >> 1) & B[ text[i] ];
				}while(d != 0);
				if(j == i)
				{
					pos = j - 1;
					return true;
				}
			}
			i -= q_point + 1;
		}
		return false;
	}

	/*bool search_next(string& text, int& start_pos_forward, int& start_pos_reverse, int& end_pos, int &pos, bool& map_fw)
	{
		if(search_fw && found_fw && update_fw)
		{
			found_fw = search_SBNDM(text, start_pos_forward, end_pos, B_forward, two_gram_forward, current_forward_anchor);
			start_pos_forward = current_forward_anchor + s0_forward;
			update_fw = false;
		}
		if(search_rc && found_rc && update_rc)
		{
			found_rc = search_SBNDM(text, start_pos_reverse, end_pos, B_reverse, two_gram_reverse, current_reverse_anchor);
			start_pos_reverse = current_reverse_anchor + s0_reverse;
			update_rc = false;
		}
		int report = 0;
		if(found_fw && found_rc)
		{
			if(current_forward_anchor <= current_reverse_anchor)	
				report = 1;
			else
				report = -1;
		}
		else if(found_fw)
			report = 1;
		else if(found_rc)
			report = -1;
		if(report == 1)
		{
			pos = current_forward_anchor;
			map_fw = true;
			update_fw = true;
			return true;
		}
		else if(report == -1)	
		{
			pos = current_reverse_anchor;
			map_fw = false;
			update_rc = true;
			return true;
		}	
		else
			return false;
	}*/
	
	// actual anchor searching function called by MapSplice, search both forward and reverse strand, find the next match
	bool search_next_simple(string& text, int& start_pos, int& end_pos_forward, int& end_pos_reverse, int &pos, bool& map_fw)
	{

		if(search_fw)
		{
			found_fw = search_next_SBNDM(text, start_pos, end_pos_forward, B_forward, two_gram_forward, current_forward_anchor);

			end_pos_forward = current_forward_anchor - s0_forward;
			pos = current_forward_anchor - pattern_size + 1;
			cur_hits ++;
			return (found_fw && cur_hits <= max_hits);
		}
		else if(search_rc)
		{
			found_rc = search_next_SBNDM(text, start_pos, end_pos_reverse, B_reverse, two_gram_reverse, current_reverse_anchor);
			end_pos_reverse = current_reverse_anchor - s0_reverse;
			pos = current_reverse_anchor - pattern_size + 1;
			cur_hits ++;
			return (found_rc && cur_hits <= max_hits);
		}
		else
			return false;
	}
	
	// find all the pattern matches in a given range of text, worker function called by "find_all_match()"
	bool search_all_SBNDM(string& text, int& start_pos, int& end_pos, unsigned int* B, unsigned int* two_gram, int& s0, vector<Bwtmap_Info>& bwt_vector, string read_id, string chrom, string strand, int pair_no, int seg_no)
	{
		/* Searching phase */
		int i = end_pos - q_point, j = 0;
		unsigned int d = 0;
		while (i >= start_pos + q_size - 1)
		{
			d = q_gram(text_ptr, two_gram, i);
			if(d != 0)
			{
				j = i + q_point + 1;
				do
				{
					i++;
					d = (d >> 1) & B[ text[i] ];
				}while(d != 0);
				if(j == i)
				{
					Bwtmap_Info new_info;
					new_info.start = j - pattern_size;
					new_info.end = j - 1;
					new_info.read_id = read_id;
					new_info.chrom = chrom;
					new_info.strand = strand;
					new_info.start_seg_no = seg_no;
					new_info.end_seg_no = seg_no;
					new_info.head = -1;
					new_info.tail = -1;
					new_info.pair_no = pair_no;
					bwt_vector.push_back(new_info);
					i -= s0;
				}
			}
			i -= q_point + 1;
		}
		return false;
	}
	
	// actual function called by MapSplice to find all the pattern matches in a given range of text	
	void find_all_match(string& text, int& start_pos, int& end_pos, vector<Bwtmap_Info>& bwt_vector, string read_id, string chrom, int pair_no, int seg_no)
	{
		if(search_fw)
			search_all_SBNDM(text, start_pos, end_pos, B_forward, two_gram_forward, s0_forward, bwt_vector, read_id, chrom, "+", pair_no, seg_no);
		else if(search_rc)
			search_all_SBNDM(text, start_pos, end_pos, B_reverse, two_gram_reverse, s0_reverse, bwt_vector, read_id, chrom, "-", pair_no, seg_no);
	}
	
	bool set(string& pattern, string& text, bool _search_fw, bool _search_rc, int _max_hits)
	{
		/* Pre processing */
		if(!check_pattern(pattern))
			return false;
		pattern_size = (int)(pattern.length());
		check_lowest_bit = 1;  
		q_point = pattern_size - q_size;
		text_ptr = const_cast<char*> (text.data());
		search_fw = _search_fw;
		search_rc = _search_rc;
		found_fw = search_fw;
		found_rc = search_rc;
		update_fw = true;
		update_rc = true;
		current_forward_anchor = -1;
		current_reverse_anchor = -1;
		max_hits = _max_hits;
		cur_hits = 0;
		if(search_fw)
		{
			anchor_forward = pattern;
			compute_bit_vector(anchor_forward, B_forward);
			compute_shift(anchor_forward, B_forward, s0_forward);
			compute_2gram_table(two_gram_forward, B_forward);
		}
 		if(search_rc)
 		{
 			anchor_reverse = revcomp(pattern);
			compute_bit_vector(anchor_reverse, B_reverse);
			compute_shift(anchor_reverse, B_reverse, s0_reverse);
			compute_2gram_table(two_gram_reverse, B_reverse);
 		}
		return true;
	}
	
private:
	char* text_ptr;
	int pattern_size;
	int text_size;
	int alphabet_size;
	int q_size;
	int q_point;
	unsigned int check_lowest_bit;
	string anchor_forward;
	string anchor_reverse;
	unsigned int B_forward[256];
	unsigned int B_reverse[256];
	int s0_forward;
	int s0_reverse;
	unsigned int two_gram_forward[21589];
	unsigned int two_gram_reverse[21589];
	int current_forward_anchor;
	int current_reverse_anchor;
	bool search_fw;
	bool search_rc;
	bool found_fw;
	bool found_rc;
	bool update_fw;
	bool update_rc;
	int max_hits;
	int cur_hits;
};


