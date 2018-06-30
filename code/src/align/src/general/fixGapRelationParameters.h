// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE

bool DO_REMAPPING_DuringFirstMapping= false;
bool BUILD_SPLICE_HASH_FROM_FILE = true;
bool BUILD_SPLICE_HASH_FROM_FIRST_MAPPING = false;
bool PRINT_JUNC = false;

typedef map<unsigned int, unsigned int> SpliceWeightMapForPrint; 
typedef map <unsigned int, SpliceWeightMapForPrint > SpliceJunctionHashForPrint; 
typedef SpliceJunctionHashForPrint::iterator SpliceJunctionHashIterForPrint;
typedef SpliceWeightMapForPrint::iterator SpliceWeightMapIterForPrint;

SpliceJunctionHashForPrint spliceJunctionForPrint[22];
SpliceJunctionHashIterForPrint spliceJuncHashIterForPrint;
SpliceWeightMapIterForPrint weightMapIterForPrint;

bool DO_NONCANONICAL = false;

bool DO_NONCANONICAL_SHORT_ANCHOR = false;

const int FIX_NO_RELATIONSHIP = 0;
const int FIX_TOO_CLOSE = 1;
const int FIX_DOUBLE_ANCHOR = 2;
const int FIX_TOO_FAR = 3;
const int FIX_INSERTION_GAP = 4;
const int FIX_INSERTION_NEIGHBOUR = 5;
const int FIX_DELETION_GAP = 6;
const int FIX_DELETION_NEIGHBOUR = 7;
const int FIX_SPLICE_GAP = 8;
const int FIX_SPLICE_NEIGHBOUR = 9;
const int FIX_MATCH = 10;
const int FIX_SOFTCLIPPING = 11;
const int FIX_REMAPPING_SHORT_HEAD = 12;
const int FIX_REMAPPING_SHORT_TAIL = 13;
const int FIX_CIRCULAR_RNA = 14;
const int FIX_CIRCULAR_RNA_GAP = 15;
const int FIX_CIRCULAR_RNA_NEIGHBOUR = 16;