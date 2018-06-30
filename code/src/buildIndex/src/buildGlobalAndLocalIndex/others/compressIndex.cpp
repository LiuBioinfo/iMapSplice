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
		cout << "Executable <InputChromIndexPrefix> <CompressedIndexFile> <IndexSize>" << endl;
		exit(0);
	}

	string chromIndexPrefixFileStr = argv[1];
	string compressedIndexFileStr = argv[2];
	string IndexSizeStr = argv[3]; //stringLength + 1; stringLength is used in build Original Size index; == MAX

	cout << "InputChromIndexPrefix: " << chromIndexPrefixFileStr << endl;
	cout << "compressedIndexPrefix: " << compressedIndexFileStr << endl;
	//cout << "indexSize: " << IndexSizeStr << endl;

	unsigned int indexSize = strtoul((IndexSizeStr.c_str()), NULL, 10);//(IndexSizeStr.c_str());

	cout << "indexSize: " << indexSize << endl;

	string up_file = chromIndexPrefixFileStr; up_file.append("_up"); ifstream up_file_ifs(up_file.c_str(), ios::binary);
	string down_file = chromIndexPrefixFileStr; down_file.append("_down"); ifstream down_file_ifs(down_file.c_str(), ios::binary);
	string next_file = chromIndexPrefixFileStr; next_file.append("_next"); ifstream next_file_ifs(next_file.c_str(), ios::binary);
	
	string childTab_file = compressedIndexFileStr; childTab_file.append("_childTab"); ofstream childTab_file_ofs(childTab_file.c_str(),ios::binary);
	string detChild_file = compressedIndexFileStr; detChild_file.append("_detChild"); ofstream detChild_file_ofs(detChild_file.c_str(), ios::binary);

	cout << "start to load up" << endl;
    unsigned int *up;
    up = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
    up_file_ifs.read((char*)up, (indexSize * sizeof(unsigned int)));

	cout << "start to load down" << endl;
    unsigned int *down;
    down = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
    down_file_ifs.read((char*)down, (indexSize * sizeof(unsigned int)));

	cout << "start to load next" << endl;
    unsigned int *next;
    next = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
    next_file_ifs.read((char*)next, (indexSize * sizeof(unsigned int)));

    cout << "All index files loaded ..." << endl;

    unsigned int *childTab;
    childTab = (unsigned int*)malloc(indexSize * sizeof(unsigned int));

   	BYTE *verifyChild;
	verifyChild = (BYTE*)malloc(indexSize * sizeof(BYTE));


	BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();

	tmpBuildIndexInfo->compressUpDownNext2ChildtabVerifyChild(
		up, down, next, childTab, verifyChild, indexSize);

	childTab_file_ofs.write((const char*) childTab, indexSize * sizeof(unsigned int));
	detChild_file_ofs.write((const char*) verifyChild, indexSize * sizeof(BYTE));	

	childTab_file_ofs.close();
	detChild_file_ofs.close();

	cout << "finish compressing index" << endl; 
	free(up); free(down); free(next); free(childTab); free(verifyChild);
	return 0;
}