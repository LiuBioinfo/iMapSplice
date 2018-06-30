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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <InputSJ> <OutputFolderPath>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;

	string inputJuncFile = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	
	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert SJ into SJmap" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfo->insertJuncFromJuncFile_withAlterSpliceSiteAnchorSimilarity(
		inputJuncFile, indexInfo);
	
	string outputSJstr = outputFolderStr + "/";
	string outputSJpath_kept 
		= outputSJstr + "keptSJalignInferHash";
	string outputSJpath_filterOut 
		= outputSJstr + "filterOutSJalignInferHash";
	string outputSJpath_filterOut_extension 
		= outputSJstr + "filterOutSJalignInferHash_extension";
	string outputSJpath_filterOut_alterSpliceSite 
		= outputSJstr + "filterOutSJalignInferHash_alterSpliceSite";
	string outputSJpath_filterOut_lowSupportNonCanonical 
		= outputSJstr + "filterOutSJalignInferHash_lowSupportNonCanonical";
	string outputSJpath_filterOut_lowSupportNonCanonical_canBeRefinedWithKeptSJspliceSite
		= outputSJstr + "filterOutSJalignInferHash_lowSupportNonCanonical_canBeRefinedWithKeptSJspliceSite";
	string outputSJpath_filterOut_lowSupportNonCanonical_canNotBeRefined
		= outputSJstr + "filterOutSJalignInferHash_lowSupportNonCanonical_canNotBeRefined";		
	//alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, outputSJstr);
	cout << "start to output SJ" << endl;
	alignInferJunctionHashInfo->outputSJ_kept(outputSJpath_kept, indexInfo);
	alignInferJunctionHashInfo->outputSJ_filterOut(
		outputSJpath_filterOut, 
		outputSJpath_filterOut_extension,
		outputSJpath_filterOut_alterSpliceSite,
		outputSJpath_filterOut_lowSupportNonCanonical,
		indexInfo);
	cout << "end of output SJ" << endl;
	cout << "start to output SJs that can be refined " << endl;
	alignInferJunctionHashInfo->outputSJ_filterOut_canBeRefinedWithKeptSJspliceSite(
		outputSJpath_filterOut_lowSupportNonCanonical_canBeRefinedWithKeptSJspliceSite,
		outputSJpath_filterOut_lowSupportNonCanonical_canNotBeRefined, indexInfo);
	cout << "end of output SJ" << endl;
	delete alignInferJunctionHashInfo;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}