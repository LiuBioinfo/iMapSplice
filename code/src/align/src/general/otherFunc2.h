// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef OTHERFUNC_2_H
#define OTHERFUNC_2_H

//#include "index_info.h"
#include <string>
#include <string.h>

using namespace std;

bool tooManyValSegs_bool_phase1(Seg_Info* segInfo, int readLength, bool oriMap_bool)
{
	if(!oriMap_bool)
		return false;

	bool tooManySegs_or_not = segInfo->returnTooManyValSegsOrNot(
		readLength, LONG_SEG_LENGTH_THRESHOLD_PHASE1, 5);
	
	return tooManySegs_or_not;
}


#endif