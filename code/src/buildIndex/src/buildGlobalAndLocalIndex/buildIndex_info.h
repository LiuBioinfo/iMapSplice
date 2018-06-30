// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

using namespace std;

#define LONGLCP 255

typedef unsigned char BYTE;

class BuildIndex_Info
{
public:
	
	unsigned int genomeLength;

	int chromNum;

	vector<string> chrNameStr; // size = chromNum

	vector<int> chromLength; // size = chromNum

	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	
	int secondLevelIndexNormalSize;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 
	//int NULL_NUM; // 2654911540 for mm9_noRandom genome
	//int MAX; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 

	BuildIndex_Info()
	{}

	void getIndexInfoFromParameterFile(ifstream& inputIndexInfoFile)
	{
		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNumLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromNameLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		getline(inputIndexInfoFile, s);
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);

		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 1;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
		}
	}

	void compressLcp2Lcpcompress(unsigned int *lcp, 
		BYTE *lcpCompress, unsigned int IndexLength)
	{
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if(lcp[Num] < LONGLCP+1)
				lcpCompress[Num] = lcp[Num];
			else
			{
				lcpCompress[Num] = LONGLCP;
			}
		}
		return;
	}

	void compressUpDownNext2ChildtabVerifyChild(
		unsigned int *up, unsigned int *down, unsigned int *next,
		unsigned int *child, BYTE *verifyChild, 
		unsigned int IndexLength)	
	{
		/////////////////////  build ChildTable ////////////////
		//next -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			child[Num] = next[Num];
		}
		//down -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((down[Num] > 0) && /*(child[Num] == 0)*/(down[Num] != up[next[Num]]))
			{
				if(child[Num]!=0)
				{
					cout << "error child[Num]!= 0 " << endl;
					exit;
				}
				child[Num] = down[Num];
			}
		}
		//up -> child
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((up[Num] > 0) /*&&(child[Num] > 0)*/)
			{
				if(child[Num-1] != 0)
				{
					cout << "error child[Num-1] > 0" << endl;
					exit;
				}
				child[Num-1] = up[Num]; 
			}
		}

		/////////////////////  build VerifyChild ////////////////
		for(unsigned int index = 0; index < IndexLength; index++)
		{
   	    	verifyChild[index] = 0;
    		if(up[index] > 0) //up
    		{
    			verifyChild[index-1] = 1;
    		}
    		if((down[index] > 0)&&(next[index] == 0))  //down 
    		{	
    			verifyChild[index] = 2;
    		}  	  	
    		if(next[index] > 0)                //next
    		{
    			verifyChild[index] = 3;
    		}
			if((down[index] > 0)&&(next[index] > 0)) 
    		{
    			verifyChild[index] = 4;
    		}
		}

		return;
	}

	bool compressUpDownNext2ChildtabVerifyChild_bool(
		unsigned int *up, unsigned int *down, unsigned int *next,
		unsigned int *child, BYTE *verifyChild, 
		unsigned int IndexLength)	
	{
		bool compressCorrect = true;
		/////////////////////  build ChildTable ////////////////
		//next -> child
		cout << "next->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			child[Num] = next[Num];
		}
		//down -> child
		cout << "down->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((down[Num] > 0) && ((next[Num] < 0)||((next[Num] > IndexLength-1))))
			{
				compressCorrect = false;
				continue;
			}
			if((down[Num] > 0) && /*(child[Num] == 0)*/(down[Num] != up[next[Num]]))
			{
				if(child[Num]!=0)
				{
					//cout << "error child[Num]!= 0 " << endl;
					compressCorrect = false;
					//exit;
				}
				child[Num] = down[Num];
			}
		}
		//up -> child
		cout << "up->child" << endl;
		for(unsigned int Num = 0; Num < IndexLength; Num++)
		{
			if((up[Num] > 0) /*&&(child[Num] > 0)*/)
			{
				if(Num == 0)
				{
					compressCorrect = false;
					continue;
				}
				if(child[Num-1] != 0)
				{
					compressCorrect = false;
					//cout << "error child[Num-1] > 0" << endl;
					//exit;
				}
				child[Num-1] = up[Num]; 
			}
		}
		cout << "build VerifyChild" << endl;
		/////////////////////  build VerifyChild ////////////////
		for(unsigned int index = 0; index < IndexLength; index++)
		{
   	    	verifyChild[index] = 0;
    		if(up[index] > 0) //up
    		{
				if(index == 0)
				{
					compressCorrect = false;
					//continue;
				}
				else
				{
	    			verifyChild[index-1] = 1;
    			}
    		}
    		if((down[index] > 0)&&(next[index] == 0))  //down 
    		{	
    			verifyChild[index] = 2;
    		}  	  	
    		if(next[index] > 0)                //next
    		{
    			verifyChild[index] = 3;
    		}
			if((down[index] > 0)&&(next[index] > 0)) 
    		{
    			verifyChild[index] = 4;
    		}
		}

		return compressCorrect;
	}

	int convertStringToInt(const string& chrName)
	{
		int chrNameInt;
		map<string, int>::iterator chrNameMapIter;
		chrNameMapIter = chrNameMap.find(chrName);

		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			cout << "...... chrom name error! ...... " << endl;
		}
		return chrNameInt;
	}

};