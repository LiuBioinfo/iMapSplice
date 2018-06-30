// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef MISMATCH_H
#define MISMATCH_H

class Mismatch
{
public:
	int seg_no; //segment number
	int pos; // mismatch position in the segment
	char src_char; // source base sequence
	char dst_char;

	Mismatch()
	{
	}

	Mismatch(int s, int p, char src, char dst)
	{
		seg_no=s;
		pos=p;
		src_char=src;
		dst_char=dst;
	}

	~Mismatch()
	{
	}
};

#endif


