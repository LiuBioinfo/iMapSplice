// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PARAMETER_INFO_H
#define PARAMETER_INFO_H

#include <string>
#include <string.h>

using namespace std;

class Parameter_Info
{
private:
	// int SEGMENTNUM;
	// int POSSIBLE_MAP_CASES_MAX;
	// int MAX_SPLICE_LENGTH;

	// int START_READ_NUM;
	// int SPLICEDISTANCE;

	// unsigned int NULL_NUM;
	// int MAX_INSERTION_LENGTH;
	// int MAX_DELETION_LENGTH;
	// int MAPCASES_MAX;

	// int ANCHOR_LENGTH;
	// int MIN_LENGTH_TO_STITCH;
	// int FIX_SPLICE_BUFFER;

	// int SHORT_EXON_BUFFER;
	// int FIX_MULTI_SPLICE_BUFFER;
	// int SHORT_EXON_MIN;

	// int SPLICE_MIN_LENGTH; 

	// int CANDALILOC;

	// int mindValSegLength;

	// int SPLICE_MIN_LENGTH;

	// int PE_MAP_DISTANCE;

	int pair_output_multi_alignment_num_max;
	int pair_mapped_length_min;
	int read_mismatch_num_max;
public:
	Parameter_Info()
	{

	}

	void set_pair_output_multi_alignment_num_max(int m)
	{
		pair_output_multi_alignment_num_max = m;
	}

	void set_pair_mapped_pair_length_min(int l)
	{
		pair_mapped_pair_length_min = l;
	}

	void set_read_mismatch_num_max(int m)
	{
		read_mismatch_num_max = m;
	}

	int return_pair_output_multi_alignment_num_max()
	{
		return pair_output_multi_alignment_num_max;
	}

	int return_pair_mapped_length_min()
	{
		return pair_mapped_length_min;
	}

	int return_mismatch_num_max()
	{
		
	}
};


#endif