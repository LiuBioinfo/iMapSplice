// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
/*    
 *    buildAll2ndLevelIndex.cpp
 *	  MapSplice3
 *
 *    Copyright (C) 2016 University of Kentucky and
 *                       Xinan Liu
 *
 *    Authors: Xinan Liu
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
#include "build2ndLevelIndex_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		cout << "Executable <InputIndexPrefix> <Output2ndLevelIndexPrefix> " << endl;
		exit(0);
	}

	string inputIndexPrefix = argv[1];
	string Output2ndLevelIndexPrefix = argv[2];

	inputIndexPrefix += "/";
	string localIndexLogStr = inputIndexPrefix + "buildLocalIndex.log";
	ofstream log_ofs(localIndexLogStr.c_str());
	cout << "inputIndexPrefix: " << inputIndexPrefix << endl;
	log_ofs << "inputIndexPrefix: " << inputIndexPrefix << endl;
	Output2ndLevelIndexPrefix += "/";
	cout << "Output2ndLevelIndexPrefix: " << Output2ndLevelIndexPrefix << endl; 
	log_ofs << "Output2ndLevelIndexPrefix: " << Output2ndLevelIndexPrefix << endl; 
	//cout << "mkdir " << endl;
	//string mkdir_cmd_outputIndex = "mkdir -p " + Output2ndLevelIndexPrefix;
	//system(mkdir_cmd_outputIndex.c_str());

	string chrom_bit_file = inputIndexPrefix; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = inputIndexPrefix; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	SecondLevelIndex_Info* secondLevelIndexInfo = new SecondLevelIndex_Info(parameter_file_ifs, log_ofs);
	
	cout << "start to load whole genome" << endl;
	log_ofs << "start to load whole genome" << endl;
	char *chrom;

	cout << "indexSize: " << secondLevelIndexInfo->indexSize << endl;
	log_ofs << "indexSize: " << secondLevelIndexInfo->indexSize << endl;

	chrom = (char*)malloc((secondLevelIndexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (secondLevelIndexInfo->indexSize) * sizeof(char)); 

	string genomeString = chrom;

	cout << "genomeString length: " << genomeString.length() << endl;
	log_ofs << "genomeString length: " << genomeString.length() << endl;
	free(chrom);
	
	cout << endl << "start to divide Genome to chromosomes" << endl;
	log_ofs << endl << "start to divide Genome to chromosomes" << endl;	
	secondLevelIndexInfo->Genome2ChromString(genomeString);
	cout << "finish dividing Genome to chromosomes" << endl;
	log_ofs << "finish dividing Genome to chromosomes" << endl;

	cout << endl <<  "start to create subfolders" << endl;
	log_ofs << endl <<  "start to create subfolders" << endl;
	secondLevelIndexInfo->creatAllSubFolder(Output2ndLevelIndexPrefix);
	cout << "finish creating subfolders" << endl;
	log_ofs << "finish creating subfolders" << endl; 
	
	cout << endl << "start to generate 2nd level index chrom file" << endl;
	log_ofs << endl << "start to generate 2nd level index chrom file" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexChrom(Output2ndLevelIndexPrefix);
	cout << "finish generating 2nd level index chrom file" << endl;
	log_ofs << "finish generating 2nd level index chrom file" << endl;

	cout << endl <<  "start to generate 2nd level index with original size" << endl;
	log_ofs << endl <<  "start to generate 2nd level index with original size" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexOriginalSize(Output2ndLevelIndexPrefix, log_ofs);
	cout << "finish generating 2nd level index with original size" << endl;
	log_ofs << "finish generating 2nd level index with original size" << endl;

	cout << endl << "start to generate 2nd level index with Compressed size" << endl;
	log_ofs << endl << "start to generate 2nd level index with Compressed size" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexCompressedSize(Output2ndLevelIndexPrefix);
	cout << "finish generating 2nd level index with Compressed size" << endl;
	log_ofs << "finish generating 2nd level index with Compressed size" << endl;
	
	cout << endl << "finish all 2nd level index" << endl;
	log_ofs << endl << "finish all 2nd level index" << endl;
	
	chrom_bit_file_ifs.close();
	parameter_file_ifs.close();

	delete(secondLevelIndexInfo);

	//free(chrom);
	return 0;
}