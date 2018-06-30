// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
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
#include <omp.h>
#include <time.h>

#include "read_block_test.h"
//#include "bwtmap_info.h"
//#include "DoubleAnchorScore.h"
//#include "sbndm.h"
//#include "otherFunc.h"
#include "index_info.h"

int main(int argc, char**argv)
{
    if(argc < 4)
	{
		//cout << "Executable <InputReads> <InputReads_PE> <OutputSAM> <threads_num> <Fasta_or_Fastq> <HumanOrMouse>" << endl;
		cout << "Executable <indexFilePrefix> <chrName> <pos_start> <pos_end>" << endl;

		exit(0); 
	}

	bool load2ndLevelIndexBool = true;
	bool detectExactFusionSiteBool = true;
	bool load2ndLevelIndexBool_compressedSize = true;

	//string inputRecordsStr = argv[1];

	string indexStr = argv[1];

	string chrNameStr = argv[2];

	string posStartStr = argv[3];

	string posEndStr = argv[4];

	int pos_start = atoi(posStartStr.c_str());
	int pos_end = atoi(posEndStr.c_str());

    //log_ofs << "indexStr: " << indexStr << endl;
   	
   	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);


	//log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	//log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	//log_ofs << "chromSize = " <<indexInfo->returnChromStringLength() << endl;
	//log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	

	cout << "end of loading chromsomes" << endl;

	string chromSubstr = indexInfo->getReferenceGenomeSubstr(chrNameStr, pos_start, pos_end);

	cout << "chromSubstr: " << chrNameStr << " " << pos_start << " ~ " << pos_end << " length: " << pos_end - pos_start + 1 << endl;
	cout << chromSubstr << endl;
	return 0;
}