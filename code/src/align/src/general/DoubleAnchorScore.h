// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef DOUBLE_ANCHOR
#define DOUBLE_ANCHOR

#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <ext/hash_map> //g++ only
#include <ext/hash_set>
#define _CRT_SECURE_NO_WARNINGS

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;

namespace __gnu_cxx
{
	template<typename Traits, typename Allocator>
	struct hash<std::basic_string<char, Traits, Allocator> >
	{
		size_t operator()(const std::basic_string<char, Traits, Allocator>& __s) const
		{
			return __stl_hash_string(__s.c_str());
		}
	};
}

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
//#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>
#include "splice_info.h"


using namespace std;

const size_t THIRTY_TWO = 32;
const size_t ALL_BITS_ON = static_cast<size_t>(-1);
const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

// This can be changed, if you have the memory and want to use it://XINAN: changed from 31 to 50
static const size_t MAX_SEED_WIDTH = 52;

static const size_t LEAST_SIG_BIT = static_cast<size_t>(1);
static const size_t SECOND_LSB = (static_cast<size_t>(1) << 1);
//static const size_t ALL_BITS_ON = static_cast<size_t>(-1);

// Set up for 64-bit ONLY!!
static const size_t SIXTY_FOUR = 64;
static const size_t MOST_SIG_BIT = static_cast<size_t>(0x8000000000000000);


static const size_t bit_GT_upper = 3;
static const size_t bit_GT_lower = 1;

static const size_t bit_TG_upper = 3;
static const size_t bit_TG_lower = 2;

static const size_t bit_AG_upper = 1;
static const size_t bit_AG_lower = 0;

static const size_t bit_GA_upper = 2;
static const size_t bit_GA_lower = 0;

static const size_t bit_GC_upper = 2;
static const size_t bit_GC_lower = 1;

static const size_t bit_CG_upper = 1;
static const size_t bit_CG_lower = 2;

static const size_t bit_AT_upper = 1;
static const size_t bit_AT_lower = 1;

static const size_t bit_TA_upper = 2;
static const size_t bit_TA_lower = 2;

static const size_t bit_AC_upper = 0;
static const size_t bit_AC_lower = 1;

static const size_t bit_CA_upper = 0;
static const size_t bit_CA_lower = 2;

static const size_t bit_CT_upper = 1;
static const size_t bit_CT_lower = 3;

static const size_t bit_TC_upper = 2;
static const size_t bit_TC_lower = 3;

//CTAC

//CTGC

//GTAT

static const size_t bit_GTAG_upper = bit_GT_upper << 2 | bit_AG_upper;
static const size_t bit_GTAG_lower = bit_GT_lower << 2 | bit_AG_lower;
static const size_t bit_GTAG = bit_GTAG_upper << 4 | bit_GTAG_lower;

static const size_t bit_GATG_upper = bit_GA_upper << 2 | bit_TG_upper;
static const size_t bit_GATG_lower = bit_GA_lower << 2 | bit_TG_lower;
static const size_t bit_GATG = bit_GATG_upper << 4 | bit_GATG_lower;

static const size_t bit_GCAG_upper = bit_GC_upper << 2 | bit_AG_upper;
static const size_t bit_GCAG_lower = bit_GC_lower << 2 | bit_AG_lower;
static const size_t bit_GCAG = bit_GCAG_upper << 4 | bit_GCAG_lower;

static const size_t bit_GACG_upper = bit_GA_upper << 2 | bit_CG_upper;
static const size_t bit_GACG_lower = bit_GA_lower << 2 | bit_CG_lower;
static const size_t bit_GACG = bit_GACG_upper << 4 | bit_GACG_lower;

static const size_t bit_ATAC_upper = bit_AT_upper << 2 | bit_AC_upper;
static const size_t bit_ATAC_lower = bit_AT_lower << 2 | bit_AC_lower;
static const size_t bit_ATAC = bit_ATAC_upper << 4 | bit_ATAC_lower;

static const size_t bit_CATA_upper = bit_CA_upper << 2 | bit_TA_upper;
static const size_t bit_CATA_lower = bit_CA_lower << 2 | bit_TA_lower;
static const size_t bit_CATA = bit_CATA_upper << 4 | bit_CATA_lower;

static const size_t bit_CTAC_upper = bit_CT_upper << 2 | bit_AC_upper;
static const size_t bit_CTAC_lower = bit_CT_lower << 2 | bit_AC_lower;
static const size_t bit_CTAC = bit_CTAC_upper << 4 | bit_CTAC_lower;

static const size_t bit_CATC_upper = bit_CA_upper << 2 | bit_TC_upper;
static const size_t bit_CATC_lower = bit_CA_lower << 2 | bit_TC_lower;
static const size_t bit_CATC = bit_CATC_upper << 4 | bit_CATC_lower;

static const size_t bit_CTGC_upper = bit_CT_upper << 2 | bit_GC_upper;
static const size_t bit_CTGC_lower = bit_CT_lower << 2 | bit_GC_lower;
static const size_t bit_CTGC = bit_CTGC_upper << 4 | bit_CTGC_lower;

static const size_t bit_CGTC_upper = bit_CG_upper << 2 | bit_TC_upper;
static const size_t bit_CGTC_lower = bit_CG_lower << 2 | bit_TC_lower;
static const size_t bit_CGTC = bit_CGTC_upper << 4 | bit_CGTC_lower;

static const size_t bit_GTAT_upper = bit_GT_upper << 2 | bit_AT_upper;
static const size_t bit_GTAT_lower = bit_GT_lower << 2 | bit_AT_lower;
static const size_t bit_GTAT = bit_GTAT_upper << 4 | bit_GTAT_lower;

static const size_t bit_TATG_upper = bit_TA_upper << 2 | bit_TG_upper;
static const size_t bit_TATG_lower = bit_TA_lower << 2 | bit_TG_lower;
static const size_t bit_TATG = bit_TATG_upper << 4 | bit_TATG_lower;

static const size_t LAST_FOUR_BIT = 0xf;
static const size_t LAST_TWO_BIT = 3;
static const size_t LAST_THIRD_FOUTH = LAST_TWO_BIT << 2;

// Assumes 4 nucleotide DNA alphabet
static const size_t alphabet_size = 4;

struct Kmer{
	bool bad;
	unsigned kmer;
	Kmer(bool bd, unsigned km) : bad(bd), kmer(km) {}
};

inline size_t
base2int(char c) {
	switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 2;
  case 't' : return 3; 
	}
	return 4;
}

inline bool
isvalid(char c) {
	return (base2int(c) != 4);
}

void
bit2misinfo(string read, string chrom_seq, size_t mis_bit, size_t pre_map_len, vector<Mismatch>& mis_pos)
{
	size_t merged_len = read.length();
	size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);
	const char* read_cstr = read.c_str();
	const char* chrom_cstr = chrom_seq.c_str();
	for (size_t i = 0; i < merged_len; ++i)
	{
		if ((selector_bit >> i) & mis_bit)
		{
			Mismatch new_mismatch(0, (int)(pre_map_len + i), chrom_cstr[i], read_cstr[i]);
			mis_pos.push_back(new_mismatch);
		}
	}
}

void
bit2misinfo(string read, string chrom_seq, size_t mis_bit, size_t pre_map_len, size_t prefix_len, vector<Mismatch>& mis_pos1, vector<Mismatch>& mis_pos2)
{
	size_t merged_len = read.length();
	size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);
	const char* read_cstr = read.c_str();
	const char* chrom_cstr = chrom_seq.c_str();
	for (size_t i = 0; i < merged_len; ++i)
	{
		if ((selector_bit >> i) & mis_bit)
		{
			Mismatch new_mismatch(0, (int)(pre_map_len + i), chrom_cstr[i], read_cstr[i]);
			if(i < prefix_len)
				mis_pos1.push_back(new_mismatch);
			else
				mis_pos2.push_back(new_mismatch);
		}
	}
}

/*void
bit2misinfo(string read, string chrom_seq, size_t mis_bit, size_t pre_map_len, vector<Mismatch>& mis_pos)
{
	size_t merged_len = read.length();
	size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);
	const char* read_cstr = read.c_str();
	const char* chrom_cstr = chrom_seq.c_str();
	for (size_t i = merged_len; i > 0; --i)
	{
		if ((selector_bit >> (i - 1)) & mis_bit)
		{
			Mismatch new_mismatch(1, (int)(pre_map_len + (merged_len - i)), chrom_cstr[merged_len - i], read_cstr[merged_len - i]);
			mis_pos.push_back(new_mismatch);
		}
	}
}*/

struct WordPair {
	WordPair(const string &s);
	WordPair() : upper(0), lower(0), bads(0) {}

	string tostring(size_t mask) const;
	string tostring2(size_t mask) const;
	string tostring3(size_t mask, size_t bitsnum) const;
	char get_char(size_t mask, size_t pos) const;

	size_t score(const WordPair &other, size_t mask) const;

	size_t score(const WordPair &other, size_t mask, size_t& rbit) const;

	size_t score_hmer(const WordPair &other, size_t mask) const;
	size_t ps_score(const WordPair &other,  const WordPair & wp, const size_t mid_width,  const size_t seed_length, const size_t mask, size_t & loc) ;
	void update_key(const size_t update_bit, const size_t mask, size_t &key) const;

	void shift(const size_t i);
	void shift_reserve(const WordPair &other);

	void shift_reserve(const WordPair &other, const size_t reserve_bit);
	void combine(const WordPair &other, size_t shift, WordPair &wp) const;
	void right_shift (const size_t i);
	void left_shift(const size_t i);
	void clear()
	{
		upper = 0;
		lower = 0;
		bads = 0; 
	}
	void ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp);

	void shift_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const size_t leftshift, const size_t rightshift, WordPair suffix_wp, WordPair &wp);
	//the following function is used to Duplicate middle part in a word
	void duplicate(const size_t prefix_mask, const size_t suffix_mask, const size_t mid_buff_mask, const size_t leftshift, const size_t rightshift, WordPair &wp);

	void duplicate_self(const size_t leftshift, WordPair &wp)
	{
		wp.upper = ((upper << leftshift)  + upper) ;
		wp.lower = ((lower << leftshift)  + lower); 

		//why?
		wp.bads =  ((bads << leftshift) + bads);  
	}

	Kmer get_kmer(size_t st, size_t end)
	{
		unsigned kmer = 0;

		bool good_hit = true;

		int st_1 = (int)st - 1;

		for (int i = (int)end - 1; i >= st_1; --i)
		{
			size_t cur_bit = LEAST_SIG_BIT << i;
			if ((cur_bit) && bads)
			{
				good_hit = false;
				break;
			}
			else
			{
				kmer = (kmer << 2) + (((cur_bit & upper) != 0) << 1) + ((cur_bit & lower) != 0);
			}				
		}

		return Kmer(good_hit, kmer);
	}

	void get_prefix(const size_t shift, WordPair &prefix_wp);

	void get_suffix(const size_t mask, WordPair &suffix_wp);

	size_t upper;
	size_t lower;
	size_t bads;

	static string bits2string(size_t mask, size_t bits);
	static inline size_t get_upper(const size_t i) {return i > 1;}
	static inline size_t get_lower(const size_t i) {return (i % 2);}
	static inline size_t get_bads(char c) {return !isvalid(c);}
};

inline string 
bits2string2(size_t mask, size_t bits) {
	string s;
	size_t selector = MOST_SIG_BIT;
	for (size_t i = 0; i < SIXTY_FOUR; ++i) {
		s += (selector & bits & mask) ? '1' : '0';
		selector >>= 1;
	}
	return s;
}

inline void
WordPair::combine(const WordPair &other, 
				  const size_t shift, WordPair &wp) const {
					  wp.upper = upper >> shift;
					  wp.lower = lower >> shift;
					  wp.bads  = bads >> shift;
					  if (shift != 0) {
						  const size_t other_shift = (SIXTY_FOUR - shift);
						  wp.upper |= (other.upper << other_shift);
						  wp.lower |= (other.lower << other_shift);
						  wp.bads  |= (other.bads  << other_shift);
					  }  
}

inline void
WordPair::shift(const size_t i) {
	//remove one base on the left and add a new base on the right. 
	upper = ((upper << 1) + (i > 1));
	lower = ((lower << 1) + (i % 2));
}

inline void
WordPair::shift_reserve(const WordPair &other) {
	//the word left shift one bit
	//add another bit which is the most_sig_bit of the other
	//itself and other creates and 64 + 64 buffer on the genome, which is used to align reads within it 

	upper = (upper << 1) + ((other.upper & MOST_SIG_BIT) != 0);
	lower = (lower << 1) + ((other.lower & MOST_SIG_BIT) != 0);
}

inline void
WordPair::shift_reserve(const WordPair &other, const size_t reserve_bit) {
	//the word left shift one bit
	//add another bit which is the most_sig_bit of the other
	//itself and other creates and 64 + 64 buffer on the genome, which is used to align reads within it 

	upper = (upper << 1) + ((other.upper & reserve_bit) != 0);
	lower = (lower << 1) + ((other.lower & reserve_bit) != 0);
}

//return mismatches between two word pairs. 
inline size_t
WordPair::score(const WordPair &other, size_t mask) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
WordPair::score(const WordPair &other, size_t mask, size_t& rbit) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	rbit = bits;

	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}


bool
score_string_old(const string& s1, const string& s2, size_t max_mismatch, size_t& num_mismatch, size_t& comb_bits)
{
	if (s1.length() != s2.length())
	{
		//cout << "different length of two strings in score_string from Double_anchored_score.h"<<endl;// <<s1 <<endl << s2 <<endl ;
		return false;
	}

	size_t mask = ALL_BITS_ON >> (SIXTY_FOUR - s2.length());

	WordPair w1(s1), w2(s2);

	register size_t bits = ((w1.upper ^ w2.upper) | 
		(w1.lower ^ w2.lower) | w1.bads | w2.bads) & mask;

	comb_bits = bits;
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word

	num_mismatch = ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF); 
	//cout << "num_mismatch = " << num_mismatch << endl;
	//cout << "max_mismatch = " << max_mismatch << endl;
	if(num_mismatch <= max_mismatch)
		return true;
	else
		return false;
}
bool
score_string(const string& s1, const string& s2, size_t max_mismatch, size_t& num_mismatch, size_t& comb_bits)
{
	/*#ifdef CAL_TIME
	score_string_begin = clock();
	#endif*/
	bool match_bool = true;
	if(s1.length() != s2.length())
	{
		cout << "s1.length != s2.length in score_string DoubleAnchorScore.h";
		match_bool = false;
	}
	int mismatchNum = 0;
	for(int tmp = 0; tmp < s1.length(); tmp++)
	{
		if(s1[tmp] != s2[tmp])
		{	
			mismatchNum++;
			if(mismatchNum > max_mismatch)
			{
				match_bool = false; // not match
				break;
			}
		}
	}
	num_mismatch = mismatchNum;
	/*#ifdef CAL_TIME
	score_string_end = clock();
	score_string_cost = score_string_cost + score_string_end - score_string_begin;
	#endif*/	
	return match_bool;
}

bool
score_string(const string& s1, const string& s2, size_t max_mismatch, size_t& comb_bits)
{
	size_t num_mismatch = 0;
	return score_string(s1, s2, max_mismatch, num_mismatch, comb_bits);
}

bool
score_string_max_extension_forward(int read_length, int total_mismatch, size_t mis_bit,  int max_mismatch, double max_mismatch_percent, int& max_extension)
{
	size_t selector_bit = LEAST_SIG_BIT;
	int current_extend_len = 0;
	for (int i = 0; i < read_length; ++i)
	{
		if ((selector_bit << i) & mis_bit)
		{
			if(total_mismatch <= max_mismatch && (1 / (double)(current_extend_len + 1)) <= max_mismatch_percent)
			{
				max_extension = read_length - i + current_extend_len;
				return max_extension > 0;
			}
			else
			{
				current_extend_len = 0;
			  total_mismatch --;	
			}
		}
		else
		{
			current_extend_len ++;
		}
	}
	max_extension = current_extend_len;
	return max_extension > 0;
}

bool
score_string_max_extension_backward(int read_length, int total_mismatch, size_t mis_bit, int max_mismatch, double max_mismatch_percent, int& max_extension)
{
	size_t selector_bit = LEAST_SIG_BIT << (read_length - 1);
	int current_extend_len = 0;
	for (int i = 0; i < read_length; ++i)
	{
		if ((selector_bit >> i) & mis_bit)
		{
			if(total_mismatch <= max_mismatch && (1 / (double)(current_extend_len + 1)) <= max_mismatch_percent)
			{
				max_extension = read_length - i + current_extend_len;
				return max_extension > 0;
			}
			else
			{
				current_extend_len = 0;
				total_mismatch --;
			}
		}
		else
		{
			current_extend_len ++;
		}
	}
	max_extension = current_extend_len;
	return max_extension > 0;
}


/*bool
score_string_max_extension_forward(int read_length, size_t mis_bit, int max_mismatch, double max_mismatch_percent, int& max_extension)
{
	size_t selector_bit = LEAST_SIG_BIT << (read_length - 1);
	int current_extend_len = 0;
	int current_mismatch = 0;
	int total_mismatch = 0;
	bool last_mismatch = true;
	for (int i = 0; i < read_length; i++)
	{
		if ((selector_bit >> i) & mis_bit)
		{
			if(!last_mismatch)
			{
				if(current_mismatch == 0)
				{
					current_extend_len = i;
				}
				else
				{
					double mismatch_match_ratio = (double)(current_mismatch) / (i - current_extend_len - current_mismatch);
					current_mismatch = 0;
					if(mismatch_match_ratio <= max_mismatch_percent)
						current_extend_len = i;
					else
						break;
				}
			}
			current_mismatch ++;
			total_mismatch ++;
			last_mismatch = true;
			if(total_mismatch > max_mismatch)
				break;
		}
		else
		{
			last_mismatch = false;
			if(i == read_length - 1)
			{
				double mismatch_match_ratio = (double)(current_mismatch) / (i + 1 - current_extend_len - current_mismatch);
				if(mismatch_match_ratio <= max_mismatch_percent)
					current_extend_len = i + 1;
			}
		}
	}
	max_extension = current_extend_len;
	return (max_extension > 0);
}



bool
score_string_max_extension_backward(int read_length, size_t mis_bit, int max_mismatch, double max_mismatch_percent, int& max_extension)
{
	size_t selector_bit = LEAST_SIG_BIT;
	int current_extend_len = 0;
	int current_mismatch = 0;
	int total_mismatch = 0;
	bool last_mismatch = true;
	for (int i = 0; i < read_length; i++)
	{
		if ((selector_bit << i) & mis_bit)
		{
			if(!last_mismatch)
			{
				if(current_mismatch == 0)
				{
					current_extend_len = i;
				}
				else
				{
					double mismatch_match_ratio = (double)(current_mismatch) / (i - current_extend_len - current_mismatch);
					current_mismatch = 0;
					if(mismatch_match_ratio <= max_mismatch_percent)
						current_extend_len = i;
					else
						break;
				}
			}
			current_mismatch ++;
			total_mismatch ++;
			last_mismatch = true;
			if(total_mismatch > max_mismatch)
				break;
		}
		else
		{
			last_mismatch = false;
			if(i == read_length - 1)
			{
				double mismatch_match_ratio = (double)(current_mismatch) / (i + 1 - current_extend_len - current_mismatch);
				if(mismatch_match_ratio <= max_mismatch_percent)
					current_extend_len = i + 1;
			}
		}
	}
	max_extension = current_extend_len;
	return (max_extension > 0);
}*/

inline size_t
WordPair::score_hmer(const WordPair &other, size_t mask) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	//assume anchor width is less than 8
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	return bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	//// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	//return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
socreBits(register size_t bits, size_t anchor_width) 
{
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
WordPair::ps_score(const WordPair &other, const WordPair & wp,  const size_t mid_width, const size_t seed_width, const size_t mask_all, size_t& loc) {

	size_t mask_midright = (ALL_BITS_ON >> (SIXTY_FOUR - mid_width))  << seed_width; 
	size_t mask_ps = mask_all - mask_midright; 
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = wp.score (other, mask_ps);

	bool debug = true;
	if (debug) {
		cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, mask_midright) << endl;
		cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, mask_ps) << endl;
		cerr << "suffix_wp" << endl << other.tostring(ALL_BITS_ON) << endl;
		cerr << "prefix_wp" << endl << wp.tostring(ALL_BITS_ON) << endl;
	}
	//   	cerr >> "mask_midright" >> endl >> bits2string2(ALL_BITS_ON, mask_midright);

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((other.upper ^ wp.upper) | 
		(other.lower ^ wp.lower) | other.bads | wp.bads) & mask_all;
	if (debug) {
		cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	}
	//the following are two pointers each point to the bit to be turned on and turned off
	size_t selector1 = LEAST_SIG_BIT << seed_width ; //to be turned on 
	size_t selector2 = selector1 << mid_width;//to be turned off

	loc = mid_width; // loc is the index to the left
	size_t mins = s;

	if (debug) cerr << "score :" << mins << endl;
	//check out the accumulation score of mismatches when the selectors are moving...
	for (size_t i = 1; i <= mid_width; ++i) {

		if (debug )
		{
			cerr << "loop " << i << "score " << s << endl;
			cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
			cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
		}



		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;     	
		if (mins > s){
			mins = s;
			loc = mid_width - i;
		}
		selector1 <<= 1;
		selector2 <<= 1;

	}
	if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;

	return mins;
}

// this basically merges upper and lower together into one value. 
// the original Key was shifted left 2 bits. 

inline void
WordPair::update_key(const size_t bit, const size_t mask, size_t &key) const {
	key = (((((key << 1) + ((upper & bit) != 0)) << 1) + 
		((lower & bit) != 0)) & mask);
}

ostream& 
operator<<(ostream& s, const WordPair& wp) {
	return s << wp.tostring2(static_cast<size_t>(-1));
}

char
WordPair::get_char(size_t mask, size_t pos) const {
	// 00 -> A, 01 -> C, 10 -> G, 11 -> T
	const size_t selector = (LEAST_SIG_BIT << (pos - 1));

	if ((mask & bads) & selector)
		return 'N';

	const bool upper_bit = ((mask & upper) & selector);
	const bool lower_bit = ((mask & lower) & selector);
	if (upper_bit) return (lower_bit) ? 'T' : 'G';
	else return (lower_bit) ? 'C' : 'A';
}

WordPair::WordPair(const string &s) : upper(0), lower(0), bads(0) {
	string::const_iterator i = s.begin();
	const string::const_iterator limit = s.end();
	while (i != limit) {
		const char c = base2int(*i) & static_cast<size_t>(3);
		upper = ((upper << 1) + get_upper(c));
		lower = ((lower << 1) + get_lower(c));
		bads  = ((bads << 1) + get_bads(*i));
		++i;
	}
}

void
WordPair::get_prefix(const size_t shift, WordPair &prefix_wp)
{
	prefix_wp.upper = upper>>shift;
	prefix_wp.lower = lower>>shift;
	prefix_wp.bads = bads>>shift;
}

void
WordPair::get_suffix(const size_t mask, WordPair &suffix_wp)
{
	suffix_wp.upper = upper & mask;
	suffix_wp.lower = lower & mask;
	suffix_wp.bads = bads & mask;
}

//right most bits corresponds to the combined word

void WordPair::left_shift(const size_t i)
{
	upper = upper << i; //only the left seed_width + mid_buff_size is useful though
	lower = lower << i;
	bads =  bads << i; 
}

void WordPair::right_shift (const size_t i){
	upper = upper >> i; //only the left seed_width + mid_buff_size is useful though
	lower = lower >> i;
	bads =  bads >> i; 
}

inline
void WordPair::ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp)
{
	wp.upper = ((upper &  prefix_mask) + (suffix_wp.upper & suffix_mask) ) & big_buff_mask;
	wp.lower = ((lower &  prefix_mask) + (suffix_wp.lower & suffix_mask) ) & big_buff_mask;
	wp.bads = ((bads &  prefix_mask) + (suffix_wp.bads & suffix_mask) ) & big_buff_mask;
}

inline
void WordPair::duplicate(const size_t prefix_mask, const size_t suffix_mask, const size_t mid_buff_mask, const size_t leftshift, const size_t rightshift, WordPair &wp)
{
	wp.upper = (((upper << leftshift) & prefix_mask)  + ((upper >> rightshift) & suffix_mask)) & mid_buff_mask;
	wp.lower = (((lower << leftshift) & prefix_mask)  + ((lower >> rightshift) & suffix_mask)) & mid_buff_mask; 

	//why?
	wp.bads = (((bads << leftshift) & prefix_mask)  + ((bads >> rightshift) & suffix_mask)) & mid_buff_mask;  
}

string
WordPair::bits2string(size_t mask, size_t bits) {
	string s;
	size_t selector = MOST_SIG_BIT;
	for (size_t i = 0; i < SIXTY_FOUR; ++i) {
		s += (selector & bits & mask) ? '1' : '0';
		selector >>= 1;
	}
	return s;
}

string
WordPair::tostring(size_t mask) const {
	const string s(bits2string(mask, upper) + "\n" +
		bits2string(mask, lower) + "\n" + 
		bits2string(mask, bads) + "\n");
	string seq;
	for (size_t i = SIXTY_FOUR; i > 0; --i)
		seq += get_char(mask, i);
	return s + seq;
}

string
WordPair::tostring2(size_t mask) const {
	string seq;
	for (size_t i = SIXTY_FOUR; i > 0; --i)
		seq += get_char(mask, i);
	return seq;
}

string
WordPair::tostring3(size_t mask, size_t bitsnum) const {
	string seq;
	for (size_t i = bitsnum; i > 0; --i)
		seq += get_char(mask, i);
	return seq;
}

struct Masks{
	Masks(/*const size_t read_width, */const size_t seed_width, const size_t max_mismatch, const size_t anchor_width, const size_t num_anchor, /*const size_t num_seg, */const size_t seg_width, const size_t extend_bits);

	void Set(const size_t read_width, const size_t seed_width);

	Masks (size_t seg_len)
	{
		//duplicate reads
		comp_buff_width = 2 * seg_len;

		comb_seg_bits_on = (ALL_BITS_ON >> (SIXTY_FOUR - comp_buff_width)); 

		suffix_seg_bits_on = (ALL_BITS_ON >> (SIXTY_FOUR - seg_len)); ;

		prefix_seg_bits_on = suffix_seg_bits_on << seg_len;

		comb_seg_first_selector_rt = LEAST_SIG_BIT;

		comb_seg_first_selector_lt = LEAST_SIG_BIT << seg_len;

		score_seg_buf_width = seg_len;
	}

	size_t mid_buff_width;

	// marks bits corresponding to the seed_key
	size_t small_mask;

	size_t kmer_mask;

	// marks  rightmost bits numbering the width of a read
	size_t big_mask;

	size_t suffix_mask;

	size_t prefix_mask;

	size_t big_buff_mask;

	size_t mask_ps;

	size_t mask_midright;

	// bit-vector indicating the location of the hit key in the frame
	size_t bad_base_mask;

	size_t bad_base_maskII;

	size_t bad_kmer_mask;

	vector<size_t> bad_kmer_masks;

	vector<size_t> bad_kmer_masks_upper;

	vector<size_t> bad_kmer_masks_lower;

	// marks bit where the key starts in the chromosome frame, the seed_width from the left
	size_t key_update_bit;

	size_t key_update_bitII;

	vector<size_t> kmer_update_bits;

	vector<size_t> kmer_update_bits_upper;

	vector<size_t> kmer_update_bits_lower;

	//vector<size_t> imer_update_bits;

	size_t first_selector_rt;

	size_t first_selector_lt;

	size_t mid_buff_width_ext;

	size_t first_half_mask;

	size_t second_half_mask;

	size_t first_half;

	size_t second_half;

	//////

	size_t comp_buff_width;

	size_t comp_left_shift_width;

	size_t comp_right_shift_width;

	size_t comp_first_half_mask;

	size_t comp_second_half_mask;

	size_t comp_big_buff_mask;

	size_t score_buff_width;

	size_t score_first_half_mask;

	size_t score_second_half_mask;

	size_t score_big_buff_mask;

	size_t score_first_selector_rt;

	size_t score_first_selector_lt;

	size_t mis_first_selector_lt;

	size_t comp_first_selector_rt;

	size_t comp_first_selector_lt;

	size_t comp_flankstr_left_shift_width;

	size_t comp_flankstr_right_shift_width;

	///fix hole

	//duplicate reads
	size_t comb_seg_bits_on;

	size_t suffix_seg_bits_on;

	size_t prefix_seg_bits_on;

	size_t right_shift_seg_width;

	size_t left_shift_seg_width;

	size_t comb_seg_first_selector_rt;

	size_t comb_seg_first_selector_lt;

	size_t score_seg_buf_width;

	//fix tail
	size_t ft_comb_prefix_half_mask;

	size_t ft_comb_suffix_half_mask;

	size_t ft_combined_mask;

	//duplicated reads extend
	size_t comb_seg_bits_on_ext;

	size_t suffix_seg_bits_on_ext;

	size_t prefix_seg_bits_on_ext;

	size_t right_shift_seg_width_ext;

	size_t left_shift_seg_width_ext;

	size_t comb_seg_first_selector_rt_ext;

	size_t comb_seg_first_selector_lt_ext;

	size_t score_seg_buf_width_ext;

	size_t mis_selector_lt;

	size_t mis_selector_rt;

	size_t prefix_ext_mask;

	size_t suffix_ext_mask;

	size_t prefix_append_mask;

	size_t prefix_append_shifted_mask;

	size_t prefix_suffix_append_mask;

	size_t reserve_bit;

};

Masks::Masks(const size_t seed_width, const size_t max_mismatch, const size_t anchor_width, const size_t num_anchor, 
	const size_t seg_width, const size_t extend_bits) : 
	bad_kmer_masks(num_anchor), bad_kmer_masks_upper(num_anchor), bad_kmer_masks_lower(num_anchor), kmer_update_bits(num_anchor), kmer_update_bits_upper(num_anchor), kmer_update_bits_lower(num_anchor)/*, imer_update_bits(num_anchor)*/
{
	comp_left_shift_width = comp_buff_width;// - anchor_width;

	comp_flankstr_left_shift_width = comp_buff_width - (anchor_width * 2);

	comp_flankstr_right_shift_width = anchor_width;

	comp_right_shift_width = 0;//anchor_width;

	comp_second_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - comp_buff_width)); 

	comp_first_half_mask = comp_second_half_mask << comp_buff_width;

	comp_big_buff_mask = comp_second_half_mask | comp_first_half_mask;

	score_buff_width = seg_width - anchor_width;// read_width - anchor_width;//read_width / 2 - anchor_width;

	score_second_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - score_buff_width)) << anchor_width;

	score_first_half_mask = score_second_half_mask << score_buff_width << anchor_width;

	score_big_buff_mask = score_second_half_mask | score_first_half_mask;

	comp_first_selector_rt = LEAST_SIG_BIT;

	comp_first_selector_lt = comp_first_selector_rt << comp_buff_width;

	score_first_selector_rt = LEAST_SIG_BIT << anchor_width;

	score_first_selector_lt = score_first_selector_rt << seg_width;//score_buff_width << anchor_width;

	mis_first_selector_lt = LEAST_SIG_BIT << (2 * seg_width - 1);

	small_mask = (LEAST_SIG_BIT << 2*seed_width) - 1;

	kmer_mask = (LEAST_SIG_BIT << 2*anchor_width) - 1;

	bad_kmer_mask = ((LEAST_SIG_BIT << anchor_width) - 1);

	for (size_t i = 0; i < num_anchor; ++i)
	{
		bad_kmer_masks[i] = bad_kmer_mask << (i * anchor_width);

		bad_kmer_masks_upper[i] = bad_kmer_mask << (i * anchor_width) << 2;

		bad_kmer_masks_lower[i] = bad_kmer_mask << (seg_width - ((i + 1) * anchor_width)) << 2;
	}

	for (size_t i = 0; i < num_anchor; ++i)
	{
		kmer_update_bits[i] = LEAST_SIG_BIT << (i * anchor_width);

		kmer_update_bits_upper[i] = LEAST_SIG_BIT << (i * anchor_width) << 2;

		kmer_update_bits_lower[i] = LEAST_SIG_BIT << (seg_width - ((i + 1) * anchor_width)) << 2;

	}

	mask_midright = (ALL_BITS_ON >> (SIXTY_FOUR - mid_buff_width - max_mismatch))  << (seed_width - max_mismatch);//- max_mismatch

	mask_ps = big_buff_mask - mask_midright;

	mid_buff_width_ext = mid_buff_width + max_mismatch;

	first_selector_rt = LEAST_SIG_BIT << (seed_width - max_mismatch);

	first_selector_lt = first_selector_rt << mid_buff_width_ext;//to be turned off

	comb_seg_bits_on = ALL_BITS_ON >> (SIXTY_FOUR - (2 * seg_width));

	suffix_seg_bits_on = ALL_BITS_ON >> (SIXTY_FOUR - seg_width);

	prefix_seg_bits_on = suffix_seg_bits_on << seg_width;

	right_shift_seg_width = 0;

	left_shift_seg_width = seg_width;

	comb_seg_first_selector_rt = LEAST_SIG_BIT;

	comb_seg_first_selector_lt = LEAST_SIG_BIT << seg_width;

	score_seg_buf_width = seg_width;

	//fix hole extend
	comb_seg_bits_on_ext = ALL_BITS_ON >> (SIXTY_FOUR - (2 * (seg_width + (2 * extend_bits))));

	suffix_seg_bits_on_ext = ALL_BITS_ON >> (SIXTY_FOUR - seg_width - (2 * extend_bits));

	prefix_seg_bits_on_ext = suffix_seg_bits_on_ext << (seg_width + (2 * extend_bits));

	right_shift_seg_width_ext = 0;

	left_shift_seg_width_ext = seg_width + (2 * extend_bits);

	comb_seg_first_selector_rt_ext = LEAST_SIG_BIT;

	comb_seg_first_selector_lt_ext = LEAST_SIG_BIT << (seg_width + (2 * extend_bits));

	score_seg_buf_width_ext = seg_width + (2 * extend_bits);

    mis_selector_lt = LEAST_SIG_BIT << ((score_seg_buf_width_ext * 2) - 1);

	mis_selector_rt = LEAST_SIG_BIT;

	suffix_ext_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits);

	prefix_ext_mask = suffix_ext_mask << seg_width;

	prefix_append_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits - seg_width);

	prefix_append_shifted_mask = prefix_append_mask << extend_bits;

	prefix_suffix_append_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits - extend_bits - seg_width);

	//fix tail
	ft_comb_suffix_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - seg_width));

	ft_comb_prefix_half_mask = ft_comb_suffix_half_mask << seg_width;

	ft_combined_mask = ft_comb_suffix_half_mask | ft_comb_prefix_half_mask;

	//fix hmer
	
	reserve_bit = LEAST_SIG_BIT << (seg_width - 1);
}


class GenomeScan{
public:
	bool Double_anchored_score(string read, string chrom_prefix, string chrom_suffix, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical);

	bool Double_anchored_score(string read, string chrom_prefix, string chrom_suffix, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq);
	
	bool Double_anchored_score(string read, string chrom_prefix, string chrom_suffix, size_t& prefix_length, size_t max_mismatch);
	
	bool Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq, size_t& msimatch_num);	

	bool Double_anchored_score_least_mis(string read, string chrom_prefix, string chrom_suffix, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical);

	bool Double_anchored_score_least_mis(string read, string chrom_prefix, string chrom_suffix, size_t& prefix_length, size_t max_mismatch);

	bool Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, size_t& mismatch_num);
	
	bool Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits);

	bool Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits, size_t& num_mismatch);

private:
	bool FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	size_t Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	string FlankString(size_t matched_flank, bool bad);

	WordPair m_five_prim_suffix, m_three_prim_prefix;

	hash_map<size_t, Masks> m_hash_map_masks;

	size_t m_matched_flank, m_matched_bads;
};

inline bool
GenomeScan::FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = score_buf_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);
	
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		//combine_upper = ((five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.upper & three_prim_mask)) >> i;
		//combine_lower = ((five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask)) >> i;
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i; 
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			if (s < mins)
			{
				mins = s;
				loc = score_buf_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = score_buf_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_ATAC:
			{
				if ( s < mins)
				{
					prim = 1;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				//if ( s < mins)
				//{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//}
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				if ( s < mins)
				{
					prim = 3;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				if ( s < mins)
				{
					prim = 4;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				if ( s < mins)
				{
					prim = 2;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;

		switch (combine_words)
		{
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		}
	}
	return false;
}

inline size_t 
GenomeScan::Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (prim < 5)
			{
				if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}


inline size_t 
GenomeScan::Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	//size_t five_prim_mask = LAST_THIRD_FOUTH;
	//size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;

	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch)// && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
	{
		mins = s;

		loc = mask_ptr->score_seg_buf_width;

		if (!mins) 
			return mins;
	}

	//five_prim_mask <<= 1;
	//three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off

	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on

	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;

	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (mins > s)
		{
			mins = s;

			loc = mask_ptr->score_seg_buf_width - i;
			if (!mins)
				return mins;
		}

		
		//if (s <= left_mismatch)
		//{
		//	if (!prim)
		//	{
		//		if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}
		//	else if (prim < 5)
		//	{
		//		if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}
		//	else if (s < mins)
		//	{
		//		if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}

		//}
		selector1 <<= 1;
		selector2 <<= 1;
		//five_prim_mask <<= 1;
		//three_prim_mask <<= 1;
	}

	return mins;
}


inline string
GenomeScan::FlankString(size_t matched_flank, bool bad)
{
	string flankstr = "";
	if (bad)
		flankstr = "BADS";
	else
	{
		for (size_t i = 7; i >= 4; --i)
		{
			if (matched_flank & (LEAST_SIG_BIT << i))
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "T";
				else
					flankstr += "G";
			}
			else
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "C";
				else
					flankstr += "A";
			}
		}
	}

	return flankstr;
}

bool 
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch)
{
	size_t comb_bits = 0;
	return Double_anchored_score(tobe_fixed_str, doner_str, acceptor_str, prefix_length, max_mismatch, comb_bits, false);	
}

bool      // select longer prefix 
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical)
{
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	//if(do_noncanonical)
		
	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	//else
		
	score = Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
			
	if(!do_noncanonical)
	{
		string flank_string = FlankString(m_matched_flank, m_matched_bads > 0);
		if(flank_string != "GTAG" && flank_string != "CTAC" && flank_string != "ATAC" && flank_string != "GTAT" && flank_string != "CTGC" && flank_string != "GCAG")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}

bool      // select longer prefix 
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq)
{
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	//if(do_noncanonical)
		
	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	//else
		
	score = Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	flank_seq = FlankString(m_matched_flank, m_matched_bads > 0);
					
	if(!do_noncanonical)
	{
		if(flank_seq != "GTAG" && flank_seq != "CTAC" && flank_seq != "ATAC" && flank_seq != "GTAT" && flank_seq != "CTGC" && flank_seq != "GCAG")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}

bool
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq, size_t& msimatch_num)
{
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	//if(do_noncanonical)
		
	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	//else
		
	score = Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	flank_seq = FlankString(m_matched_flank, m_matched_bads > 0);
					
	if(!do_noncanonical)
	{
		if(flank_seq != "GTAG" && flank_seq != "CTAC" && flank_seq != "ATAC" && flank_seq != "GTAT" && flank_seq != "CTGC" && flank_seq != "GCAG")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		msimatch_num = score;
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}

bool 
GenomeScan::Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch)
{
	size_t comb_bits = 0;
	return Double_anchored_score_least_mis(tobe_fixed_str, doner_str, acceptor_str, prefix_length, max_mismatch, comb_bits, false);	
}

bool      // select longer prefix 
GenomeScan::Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical)
{
	
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	//if(do_noncanonical)
		
	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	//else
		
	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	if(!do_noncanonical)
	{
		string flank_string = FlankString(m_matched_flank, m_matched_bads > 0);
		if(flank_string != "GTAG" && flank_string != "CTAC" && flank_string != "ATAC" && flank_string != "GTAT" && flank_string != "CTGC" && flank_string != "GCAG")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}

bool      // select longer prefix 
GenomeScan::Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, size_t& mismatch_num)
{
	
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	//if(do_noncanonical)
		
	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	//else
		
	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	if(!do_noncanonical)
	{
		string flank_string = FlankString(m_matched_flank, m_matched_bads > 0);
		if(flank_string != "GTAG" && flank_string != "CTAC" && flank_string != "ATAC" && flank_string != "GTAT" && flank_string != "CTGC" && flank_string != "GCAG")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			mismatch_num = score;
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}

/*bool    //select shorter prefix
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical)
{
	tobe_fixed_str = revcomp(tobe_fixed_str);

	string temp_doner = doner_str, temp_acceptor = acceptor_str;

	doner_str = revcomp(temp_acceptor);

	acceptor_str = revcomp(temp_doner);

	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);
	
	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;
	
	if(do_noncanonical)
		
		score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
	
	else
		
		score = Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);
			
	if(!do_noncanonical)
	{
		string flank_string = FlankString(m_matched_flank, m_matched_bads > 0);
		if(flank_string != "GTAG" && flank_string != "CTAC")
			return false;
	}
	
	if(score <= max_mismatch)
	{
		prefix_length = tobe_fixed_len - max_loc;;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}		
		return true;
	}
	else
		return false;
}*/



bool 
GenomeScan::Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits)
{
	//cout << "stop3" << endl;
	size_t tobe_fixed_len = chrom_seq.length();

	Masks mask(tobe_fixed_len);
	//cout << "stop4" << endl;
	string prefix_read = tobe_fixed_str.substr(0, tobe_fixed_len);

	string suffix_read = tobe_fixed_str.substr(tobe_fixed_str.length() - tobe_fixed_len, tobe_fixed_len);

	string comb_read = prefix_read + suffix_read;
	//cout << "stop5" << endl;
	WordPair comb_read_wp(comb_read);

	WordPair dup_chrom_wp, chrom_wp(chrom_seq);

	chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);
	//cout << "stop6" << endl;
	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = 2;

	size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatches, rbits, &mask);

	vector<pair<size_t, pair<char, char> > > mis_pos;

	if(/*max_loc != 0 && max_loc != chrom_seq.length() && */score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_read_str = comb_read.substr(0, max_loc) + comb_read.substr(comb_read.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}

		return true;
	}
	else
		return false;
}

bool 
GenomeScan::Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits, size_t& num_mismatch)
{
	size_t tobe_fixed_len = chrom_seq.length();

	Masks mask(tobe_fixed_len);

	string prefix_read = tobe_fixed_str.substr(0, tobe_fixed_len);

	string suffix_read = tobe_fixed_str.substr(tobe_fixed_str.length() - tobe_fixed_len, tobe_fixed_len);

	string comb_read = prefix_read + suffix_read;

	WordPair comb_read_wp(comb_read);

	WordPair dup_chrom_wp, chrom_wp(chrom_seq);

	chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = 2;

	size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatches, rbits, &mask);

	vector<pair<size_t, pair<char, char> > > mis_pos;

	if(/*max_loc != 0 && max_loc != chrom_seq.length() && */score <= max_mismatch)
	{
		//cout << "score: " << score << endl;
		(num_mismatch) = score;
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_read_str = comb_read.substr(0, max_loc) + comb_read.substr(comb_read.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}

		return true;
	}
	else
		return false;
}


#endif

