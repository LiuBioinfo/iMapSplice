// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef CONSTANTDEFINITIONS_H
#define CONSTANTDEFINITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

using namespace std;

#define Debug_SNPmap_Bool false

#define MIN_SEGMENTLENGTH_TOEXTEND_DURINGSEGMENTMAPPING 30

#define INDEL_BESIDE_SJ_LEN_MAX 3

//parameters for parallel processes
#define ReadNumInReadArray_Phase1 1000000
#define InputTimeWeight_Phase1 1000000
#define OutputTimeWeight_Phase1 1

#define ReadNumInReadArray_FixOneEndUnmapped 1000000		
#define InputTimeWeight_FixOneEndUnmapped 1000000
#define OutputTimeWeight_FixOneEndUnmapped 1

#define ReadNumInReadArray_FixHeadTail 1000000	
#define InputTimeWeight_FixHeadTail 1000000
#define OutputTimeWeight_FixHeadTail 1



//some stricts in alignment pairs # to avoid too much time cost in some specific repeat-region alignments
#define MAX_INTER_ALIGN_PAIR_NUM_BEFORE_FILTER_OVERLAP_ONES 100

#define ToProcessInterAlignPairMax 20
#define ToOutputFinalAlignPairMax 5

//penalyies in choosing Better alignment
#define DISTANT_INSERT_SIZE_MIN 100000
#define DISTANT_INSERT_SIZE_PENALTY 8

#define MISMATCH_PENALTY 2
#define DELETION_PENALTY 3
#define INSERTION_PENALTY 4
#define SEMICANONICALSJ_PENALTY 8
#define NONCANONICALSJ_PENALTY 10

//#define ALIGNMENT_SCORE_MIN_OUTPUT 160
#define ALIGNMENT_SCORE_MIN_OUTPUT_COMPLEMENT_PerHundredBases 25
#define Min_Map_Score_keptAlignmentPair 50
#define Min_Map_Score_keptAlignmentPair_PerHundredBases 50
//penalties in choosing best splice site when do fixDoubleAnchor splice
#define CANONICAL_SJ_PENALTY 0
#define SEMICANONICAL_SJ_PENALTY 2
#define NONCANONICAL_SJ_PENALTY 3

#define STORE_MISMATCH_POS true//required to be true; Do not change it !
#define STORE_MISMATCH_CHA true//false
#define Do_detectNovelSJ_bool true
#define FIX_INDEL_AROUND_SJ_BOOL true

#define MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH 5
#define MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH 8

#define MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH_WITH_ANNOTATED_SJ 2
#define MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH_WITH_ANNOTATED_SJ 5

#define MATCH_BASE_PER_MISMATCH_BASE 5

#define LONG_SEG_LENGTH_THRESHOLD_PHASE1 20//default: 20
#define LONG_SEG_LENGTH_THRESHOLD_PHASE2  15// default: 15

#define BEST_ALIGNMENT_DIFFERENCE_MAX 0.0001

#define MIN_MAP_LEN_PAIR 75

#ifdef SPEED_MODE
#define SEGMENTNUM 5
#define CANDALILOC 5
#else
#define SEGMENTNUM 20
#define CANDALILOC 40
#endif
//#define POSSIBLE_MAP_CASES_MAX 40

#define READ_ALIGN_AREA_LENGTH 500000
#define START_READ_NUM 0

#define PAIR_READ_DISTANCE_MAX 500000

#define PAIR_READ_DISTANCE_SHORT_MAX 500

#define MAX_INSERTION_LENGTH 10  // required: MAX_INSERTION_LENGTH == MAX_DELETION_LENGTH
#define MAX_DELETION_LENGTH 10
#define MAPCASES_MAX 200
//#define POSSIBLE_MAP_CASES_MAX 40
#define ANCHOR_LENGTH 10
#define MIN_LENGTH_TO_STITCH 20
#define FIX_SPLICE_BUFFER 5

#define SHORT_EXON_BUFFER 2
#define FIX_MULTI_SPLICE_BUFFER 3
#define SHORT_EXON_MIN 4

#define minValSegLength 20

typedef unsigned char BYTE;

#define min_anchor_length 8

#define SPLICE_MIN_LENGTH 50

#define ToRecordIntermediateAlignInfo_MAX 40

#define PE_MAP_DISTANCE 300000

#define READ_LENGTH_MAX 200

//#define SHORT_LONG_END_THRESHOLD 10

#define CONFIDENT_SEG_LENGTH_FIX_LONG_END 12

#define DETECT_NONCANONICAL_SJ false

//#define MAX_SPLICE_LENGTH 300000
//#define SPLICEDISTANCE 300000
#define MAX_SPLICE_LENGTH 300000
#define MAX_SPLICE_DISTANCE_PHASE1 300000
#define MAX_SPLICE_LENGTH_PHASE1 300000  // required: MAX_SPLICE_LENGTH_PHASE1 = MAX_SPLICE_DISTANCE_PHASE1
#define SPLICE_MIN_LENGTH 50
#define MAX_SPLICE_DISTANCE_PHASE2 30000
#define MAX_SHORT_SPLICE_DISTANCE 10000

#define MAX_SPLICE_DISTANCE_TARGETMAPPING_LONGSINGLEANCHOR 50000
#define MAX_SPLICE_DISTANCE_TARGETMAPPING_SHORTSINGLEANCHOR 10000

#define NO_SPLICE_JUNCTION_CASE 0
#define SPLICE_JUNCTION_CANONICAL_ONLY 1
#define SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL 2 // (no noncanonical splice junctions exist, semi-canonical SJ exists)
#define SPLICE_JUNCTION_NONCANONICAL 3 // noncanonical splice junction exists

#define FIX_SPLICE_BUFFER_MAX 20
#define PreIndexSize 268435456

#define SINGLE_ANCHOR_TARGETMAPPING_BUFFER 6

#define ALIGNINFER_BASE_NUM_MAX 30

#define LengthOfSeqPerMismatchAllowed_REMAPPING 5
#define LengthOfSeqPerMismatchAllowed_INDEL_READ_END 8

#define FUSIONBREAKPOINTSEARCH_ENCOMPASSING_AREA_MAX 500000

#define GlobalMapForFusionDetectionUnfixedHeadTailBuffer 4
#define InsertedBaseNumMaxBetweenFusedTranscript 1
#define AllowedMismatchNumMaxPerBaseNum 5

#define IGNORED_SEG_LENGTH_WHEN_UPDATE_SEGINFO_WITH_TARGETMAP2SNPSEQINDEX_MAX 8
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

#endif