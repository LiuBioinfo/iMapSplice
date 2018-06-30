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
#include "otherFunc.h"

using namespace std;  


int main(int argc, char**argv)
{
    if(argc < 2)
	{
		//cout << "Executable <InputReads> <InputReads_PE> <OutputSAM> <threads_num> <Fasta_or_Fastq> <HumanOrMouse>" << endl;
		cout << "Executable <outputReadFilePrefix> <ReadNum>" << endl;
		exit(0);
	}
	string outputReadFilePrefix = argv[1];

	string outputReadFilePrefix_1 = outputReadFilePrefix + ".mapped.fa.1";
	string outputReadFilePrefix_2 = outputReadFilePrefix + ".mapped.fa.2";

	string outputReadFilePrefix_1_unmap = outputReadFilePrefix + ".unmapped.fa.1";
	string outputReadFilePrefix_2_unmap = outputReadFilePrefix + ".unmapped.fa.2";

	ofstream outputReadFilePrefix_1_ofs(outputReadFilePrefix_1.c_str());
	ofstream outputReadFilePrefix_2_ofs(outputReadFilePrefix_2.c_str());

	ofstream outputReadFilePrefix_1_unmap_ofs(outputReadFilePrefix_1_unmap.c_str());
	ofstream outputReadFilePrefix_2_unmap_ofs(outputReadFilePrefix_2_unmap.c_str());

	//string forward_or_reverseComplement = argv[2];

	string readNumStr = argv[2];

	int readNum = atoi(readNumStr.c_str());

	string readSeq_1_nor_mapped = "AGCAGATATTAAGGCTGGGAAGTATTTGGAGCATGGGGAGTACGAAGGAAACCTGTATGGAACCAAGATCGACTCTATTCTCGAGGTTGTCCAGACTGGG";
	string readSeq_2_rcm_mapped = "CCCGCATCCACCACAGCCTTGTGCATGGCACGCAAGGTCTCTAGCTCTGGAGCAACAATAAACACCACATAAGGCATGAATTCAGAAGTCCTTAACACTT";
	string readSeq_1_rcm_unmapped = "CCCGCATCCACCACAGCCTTGTGCATGGCACGCAAGGTCTCTAGCTCTGGAGCAACAATAAACACCACATAAGGCATGAATTCAGAAGTCCTTAACACTT";
	//convertStringToReverseComplement(readSeq_2_rcm_mapped);
	string readSeq_2_nor_unmapped = "AGCAGATATTAAGGCTGGGAAGTATTTGGAGCATGGGGAGTACGAAGGAAACCTGTATGGAACCAAGATCGACTCTATTCTCGAGGTTGTCCAGACTGGG";


	for(int tmp = 0; tmp < readNum; 
		tmp++)
	{
		outputReadFilePrefix_1_ofs << ">seq." << tmp+1 << "/1\n";
		outputReadFilePrefix_1_ofs << readSeq_1_nor_mapped << endl;
		outputReadFilePrefix_2_ofs << ">seq." << tmp+1 << "/2\n";
		outputReadFilePrefix_2_ofs << readSeq_2_rcm_mapped << endl;
		outputReadFilePrefix_1_unmap_ofs << ">seq." << tmp+1 << "/1\n";
		outputReadFilePrefix_1_unmap_ofs << readSeq_1_rcm_unmapped << endl; 
		outputReadFilePrefix_2_unmap_ofs << ">seq." << tmp+1 << "/2\n";
		outputReadFilePrefix_2_unmap_ofs << readSeq_2_nor_unmapped << endl;
	}

	outputReadFilePrefix_1_ofs.close();
	outputReadFilePrefix_2_ofs.close();

	outputReadFilePrefix_1_unmap_ofs.close();
	outputReadFilePrefix_2_unmap_ofs.close();

	return 0;
}