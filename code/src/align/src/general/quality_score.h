// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
/*    
 *    quality_score.h		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng
 *
 *    Authors: Zeng Zheng
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef QUALITY_SCORE_H
#define QUALITY_SCORE_H

#include <string>
#include <math.h>
using namespace std;
//#pragma warning(disable:4996)


double CalcPSanger(char inchar)
{  
       int intchar =  inchar;
       if (intchar >= 64 && intchar <= 126)
               intchar = intchar - 64;
       return pow(double(10), double(-intchar)/double(10));
}

double CalcPSanger_v1_8(char inchar)
{
       int intchar =  inchar;
       if (intchar >= 33 && intchar <= 126)
               intchar = intchar - 33;
       return pow(double(10), double(-intchar)/double(10));
}

double CalcPSolexa(char inchar)
{
       int intchar =  inchar;
       if (intchar >= 64 && intchar <= 126)
               intchar = intchar - 64;
       double ppow = pow(double(10), double(-intchar)/double(10));
       return ppow / (double(1) + ppow);
}

int GetQualityScore(const vector<Mismatch>& mismatch_vec, const string quality_str, string quality_scale)
{
	double qual_score = 0;
	vector<bool> mismatch_score(quality_str.length(), false);
	vector<Mismatch>::const_iterator mis_iter;
	const   char   *quality_str_c=quality_str.c_str(); 
	for (mis_iter = mismatch_vec.begin(); mis_iter != mismatch_vec.end(); ++mis_iter)
	{
		mismatch_score[mis_iter->pos] = true;
		if(quality_scale == "phred64")
			qual_score += log(CalcPSanger(quality_str_c[mis_iter->pos]) /(double(1) - double(0.25)));
		else if(quality_scale == "phred33")
			qual_score += log(CalcPSanger_v1_8(quality_str_c[mis_iter->pos]) /(double(1) - double(0.25)));
		else if(quality_scale == "solexa64")
			qual_score += log(CalcPSolexa(quality_str_c[mis_iter->pos]) /(double(1) - double(0.25)));
	}
	for (size_t i = 0; i < mismatch_score.size(); ++i)
	{
		if (!mismatch_score[i])
		{
			if(quality_scale == "phred64")
				qual_score += log((double(1) - CalcPSanger(quality_str_c[i])) / double(0.25));
			else if(quality_scale == "phred33")
				qual_score += log((double(1) - CalcPSanger_v1_8(quality_str_c[i])) / double(0.25));
			else if(quality_scale == "solexa64")
				qual_score += log((double(1) - CalcPSolexa(quality_str_c[i])) / double(0.25));
		}
	}
	if(qual_score<0)
		qual_score=0;
	if(qual_score>255)
		qual_score=255;
	return (int)qual_score;
}

int GetQualityScore_new(const vector<int>& mismatchPosVec, const string quality_str, string quality_scale)
{
	double qual_score = 0;
	vector<bool> mismatch_score(quality_str.length(), false);
	vector<int>::const_iterator mis_iter;
	const char *quality_str_c=quality_str.c_str(); 
	//for (mis_iter = mismatchPosVec.begin(); mis_iter != mismatchPosVec.end(); ++mis_iter)
	for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
	{
		int tmpMismatchPos = mismatchPosVec[tmp]-1; // correct? mismatchPos here is 0-coordinate ?
		mismatch_score[tmpMismatchPos] = true;
		if(quality_scale == "phred64")
			qual_score += log(CalcPSanger(quality_str_c[tmpMismatchPos]) /(double(1) - double(0.25)));
		else if(quality_scale == "phred33")
			qual_score += log(CalcPSanger_v1_8(quality_str_c[tmpMismatchPos]) /(double(1) - double(0.25)));
		else if(quality_scale == "solexa64")
			qual_score += log(CalcPSolexa(quality_str_c[tmpMismatchPos]) /(double(1) - double(0.25)));
	}
	for (size_t i = 0; i < mismatch_score.size(); ++i)
	{
		if (!mismatch_score[i])
		{
			if(quality_scale == "phred64")
				qual_score += log((double(1) - CalcPSanger(quality_str_c[i])) / double(0.25));
			else if(quality_scale == "phred33")
				qual_score += log((double(1) - CalcPSanger_v1_8(quality_str_c[i])) / double(0.25));
			else if(quality_scale == "solexa64")
				qual_score += log((double(1) - CalcPSolexa(quality_str_c[i])) / double(0.25));
		}
	}
	if(qual_score<0)
		qual_score=0;
	if(qual_score>255)
		qual_score=255;
	return (int)qual_score;
}

#endif

