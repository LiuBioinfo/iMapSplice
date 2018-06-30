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

using namespace std;

void extractChrNamePosFromRedCloverGenePosStr(string& tmpRedCoverGenePosStr, 
	string& tmpGeneChrName, int& tmpGeneStartPos, 
	int& tmpGeneEndPos, string& tmpGeneName)
{
	int tabLoc_1 = tmpRedCoverGenePosStr.find("\t");
	int tabLoc_2 = tmpRedCoverGenePosStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpRedCoverGenePosStr.find("\t", tabLoc_2 + 1);	
	int tabLoc_4 = tmpRedCoverGenePosStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpRedCoverGenePosStr.find("\t", tabLoc_4 + 1);
	int tabLoc_6 = tmpRedCoverGenePosStr.find("\t", tabLoc_5 + 1);
	int tabLoc_7 = tmpRedCoverGenePosStr.find("\t", tabLoc_6 + 1);
	int tabLoc_8 = tmpRedCoverGenePosStr.find("\t", tabLoc_7 + 1);
	tmpGeneChrName = tmpRedCoverGenePosStr.substr(0, tabLoc_1);
	string tmpGeneStartPosStr = tmpRedCoverGenePosStr.substr(
		tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	tmpGeneStartPos = atoi(tmpGeneStartPosStr.c_str());
	string tmpGeneEndPosStr = tmpRedCoverGenePosStr.substr(
		tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
	tmpGeneEndPos = atoi(tmpGeneEndPosStr.c_str());
	string tmpOtherStr = tmpRedCoverGenePosStr.substr(tabLoc_8 + 1);
	int IdLoc = tmpRedCoverGenePosStr.find("ID=", tabLoc_8 + 1);
	int tmpGeneNameStrLen = IdLoc - 1 - tabLoc_8 - 6;
	tmpGeneName = tmpRedCoverGenePosStr.substr(tabLoc_8 + 6, tmpGeneNameStrLen);
}

void extractChrNamePosFromAsmStr(
	string& tmpAsmStr, string& tmpChrName, 
	int& tmpStartPos, int& tmpEndPos, string& tmpOtherStr)
{
	int tabLoc_1 = tmpAsmStr.find("\t");
	int tabLoc_2 = tmpAsmStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpAsmStr.find("\t", tabLoc_2 + 1);
	tmpChrName = tmpAsmStr.substr(0, tabLoc_1);
	string tmpStartPosStr = tmpAsmStr.substr(
		tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEndPosStr = tmpAsmStr.substr(
		tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	tmpStartPos = atoi(tmpStartPosStr.c_str());
	tmpEndPos = atoi(tmpEndPosStr.c_str());
	tmpOtherStr = tmpAsmStr.substr(tabLoc_3 + 1);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputGenePosFile inputASMfile outputASMwithGeneNameFile" << endl;
		exit(1);
	}

	string inputGenePosFile = argv[1];
	string inputASMfile = argv[2];
	string outputASMwithGeneNameFile = argv[3];

	cout << "start to read RedCover gene file " << endl;
	vector<string> geneChrNameVec;
	vector<int> geneStartPosVec;
	vector<int> geneEndPosVec;
	vector<string> geneNameVec;
	ifstream genePos_ifs(inputGenePosFile.c_str());
	while(!genePos_ifs.eof())
	{
		string tmpGeneStr;
		getline(genePos_ifs, tmpGeneStr);
		if(tmpGeneStr == "")
			break;
		string tmpGeneChrName;
		int tmpGeneStartPos;
		int tmpGeneEndPos;
		string tmpGeneName;
		extractChrNamePosFromRedCloverGenePosStr(tmpGeneStr, tmpGeneChrName, 
			tmpGeneStartPos, tmpGeneEndPos, tmpGeneName);
		geneChrNameVec.push_back(tmpGeneChrName);
		geneStartPosVec.push_back(tmpGeneStartPos);
		geneEndPosVec.push_back(tmpGeneEndPos);
		cout << "tmpGeneName: " << tmpGeneName << endl;
		geneNameVec.push_back(tmpGeneName);
	}
	genePos_ifs.close();

	cout << "start to read ASM file" << endl;
	ifstream asm_ifs(inputASMfile.c_str());
	vector<string> asmChrNameVec;
	vector<int> asmStartPosVec;
	vector<int> asmEndPosVec;
	vector<string> asmOtherStrVec;
	string asmHeaderStr;
	getline(asm_ifs, asmHeaderStr);
	while(!asm_ifs.eof())
	{
		string tmpAsmStr;
		getline(asm_ifs, tmpAsmStr);
		if(tmpAsmStr == "")
			break;
		string tmpAsmChrName;
		int tmpAsmStartPos;
		int tmpAsmEndPos;
		string tmpAsmOtherStr;
		extractChrNamePosFromAsmStr(tmpAsmStr, tmpAsmChrName, 
			tmpAsmStartPos, tmpAsmEndPos, tmpAsmOtherStr);
		asmChrNameVec.push_back(tmpAsmChrName);
		asmStartPosVec.push_back(tmpAsmStartPos);
		asmEndPosVec.push_back(tmpAsmEndPos);
		asmOtherStrVec.push_back(tmpAsmOtherStr);
	}
	asm_ifs.close();

	ofstream ASMwithGeneName_ofs(outputASMwithGeneNameFile.c_str());
	cout << "start to assign gene names to ASMs and output " << endl;
	for(int tmpASM = 0; tmpASM < asmChrNameVec.size(); tmpASM++)
	{
		string tmpASM_chrName = asmChrNameVec[tmpASM];
		int tmpASM_startPos = asmStartPosVec[tmpASM];
		int tmpASM_endPos = asmEndPosVec[tmpASM];
		string tmpASM_otherStr = asmOtherStrVec[tmpASM];
		string tmpASM_assignedGeneName = "";
		for(int tmpGene = 0; tmpGene < geneChrNameVec.size(); tmpGene++)
		{
			string tmpGene_chrName = geneChrNameVec[tmpGene];
			int tmpGene_startPos = geneStartPosVec[tmpGene];
			int tmpGene_endPos = geneEndPosVec[tmpGene];
			if(tmpASM_chrName == tmpGene_chrName)
			{
				if(!((tmpASM_startPos > tmpGene_endPos)||(tmpASM_endPos < tmpGene_startPos)))
				{
					tmpASM_assignedGeneName += geneNameVec[tmpGene];
					tmpASM_assignedGeneName += ",";
				}
				else
				{}
			}
		}
		if(tmpASM_assignedGeneName == "")
			tmpASM_assignedGeneName = "NULL";
		ASMwithGeneName_ofs << tmpASM_chrName << "\t" << tmpASM_startPos << "\t"
			<< tmpASM_endPos << "\t" << tmpASM_otherStr << "\t" << tmpASM_assignedGeneName << endl;
	}
	ASMwithGeneName_ofs.close();
	cout << "end of running assignGeneName2ASM" << endl;
	return 0;
}