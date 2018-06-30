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

#include "buildIndex_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		cout << "Executable <InputChromIndexPrefix> <CompressedIndexFilePrefix> <IndexSize>" << endl;
		exit(0);
	}

	string chromIndexPrefixFileStr = argv[1];
	string compressedIndexFileStr = argv[2];
	string IndexSizeStr = argv[3]; //stringLength + 1; stringLength is used in build Original Size index

	unsigned int indexSize = strtoul((IndexSizeStr.c_str()), NULL, 10);//(IndexSizeStr.c_str());

	cout << "indexSize: " << indexSize << endl;

	string lcp_file = chromIndexPrefixFileStr; lcp_file.append("_lcp"); ifstream lcp_file_ifs(lcp_file.c_str(), ios::binary);

	string lcpCompress_file = chromIndexPrefixFileStr; lcpCompress_file.append("_lcpCompress"); ofstream lcpCompress_file_ofs(lcpCompress_file.c_str(),ios::binary);


	cout << "start to load lcp" << endl;
    unsigned int *lcp;
    lcp = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
    lcp_file_ifs.read((char*)lcp, (indexSize * sizeof(unsigned int)));

	BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();

	BYTE *lcpCompress = (BYTE*)malloc(indexSize * sizeof(BYTE));
	tmpBuildIndexInfo->compressLcp2Lcpcompress(lcp, 
		lcpCompress, indexSize);

	lcpCompress_file_ofs.write((const char*) lcpCompress, indexSize * sizeof(BYTE));	
	
	free(lcpCompress);
	free(lcp);
	lcp_file_ifs.close();
	lcpCompress_file_ofs.close();

	return 0;
}