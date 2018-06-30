// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef TRANSCRIPT_SET_H
#define TRANSCRIPT_SET_H
#include "index_info.h"
#include "normal_transcript_info.h" 

class Transcript_Set
{
private:
	vector<Normal_Transcript_Info*> normal_transcript_vec;
	map <string, int> transcriptSetIndexMap;
	//map<string, int> transcript_set_map;
public:
	Transcript_Set()
	{

	}

	int returnTranscriptExonEndPosWithExonIndex(int tmpTranscriptIndex, int tmpExonIndex)
	{
		return normal_transcript_vec[tmpTranscriptIndex]->returnExonEndPos(tmpExonIndex + 1);
	}

	int returnTranscriptExonStartPosWithExonIndex(int tmpTranscriptIndex, int tmpExonIndex)
	{
		return normal_transcript_vec[tmpTranscriptIndex]->returnExonStartPos(tmpExonIndex + 1);
	}
	
	string returnTranscriptStrand(int index)
	{	
		return normal_transcript_vec[index]->returnTranscriptStrand();
	}

	// int returnTranscriptSJnum(int index)
	// {
	// 	return normal_transcript_vec[index]->returnTranscriptSJnum();
	// }

	/*void outputTranscriptInfo(ofstream& output_ofs)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			output_ofs << normal_transcript_vec[tmp]->returnNormalTranscriptInfo() << endl;
		}
	}*/
	void generateGenomicJumpCodeVecWithStartAndEndLocInTranscript(
		int index, int startLocInTranscript, 
		int endLocInTranscript, vector<Jump_Code>& generatedJumpCodeVec)
	{
		normal_transcript_vec[index]->generateGenomicJumpCodeVecWithStartAndEndLocInTranscript(
			startLocInTranscript, endLocInTranscript, generatedJumpCodeVec);
	}

	int returnGenomicPosWithLocInTranscript(int index, int locInTranscript)
	{
		return normal_transcript_vec[index]->returnGenomicPosWithLocInTranscript(
			locInTranscript);
	}

	int returnMapPosInChrWithMapPosInTranscript(int index, int mapPosInTranscript)
	{
		return normal_transcript_vec[index]->returnMapPosInChrWithMapPosInTranscript(
			mapPosInTranscript);
	}

	int returnTranscriptExonNum(int index)
	{
		return normal_transcript_vec[index]->returnExonNum();
	}

	string returnTranscriptID(int index)
	{
		return normal_transcript_vec[index]->returnTranscriptName();
	}

	string returnTranscriptChromName(int index)
	{
		return normal_transcript_vec[index]->returnChromName();
	}

	string returnChromName(int index)
	{
		return normal_transcript_vec[index]->returnChromName();
	}

	void copyExonPos2anotherVec(int index,
		vector<int>& tmpAnotherExonStartPosVec, vector<int>& tmpAnotherExonEndPosVec)
	{
		normal_transcript_vec[index]->copyExonPos2anotherVec(tmpAnotherExonStartPosVec, 
			tmpAnotherExonEndPosVec);
	}

	int returnTranscriptStartPosInChr(int index)
	{
		return normal_transcript_vec[index]->returnTranscriptStartPosInChr();
	}

	int returnTranscriptEndPosInChr(int index)
	{
		return normal_transcript_vec[index]->returnTranscriptEndPosInChr();
	}

	void getSJsiteAndSize(int index, int alignmentStartPosInTranscript, 
		int alignmentEndPosInTranscript,
		vector< pair<int,int> >& SJsiteAndSize) // vector < SJsiteInExonArea, SJsize >
	{
		//cout << "start to getSJsiteAndSize in transcript_set" << endl;
		normal_transcript_vec[index]->generateSJsiteAndSizeInAlignment(
			alignmentStartPosInTranscript,
			alignmentEndPosInTranscript, SJsiteAndSize);
		//cout << "finish generateSJsiteAndSizeInAlignment" << endl;
	}

	int returnTranscriptSetIndex(const string& transcript_name)
	{
		map<string, int>::iterator mapIter;
		mapIter = transcriptSetIndexMap.find(transcript_name);
		if(mapIter == transcriptSetIndexMap.end())
		{
			return -1;
		}
		else
		{
			return mapIter->second;
		}
	}

	void extractTranscript(ifstream& annotation_ifs, Index_Info* indexInfo, const string& formatStr)
	{
		if(formatStr == "GAF")
		{
			this->extractTranscript_GAF(annotation_ifs, indexInfo);
		}
		else if(formatStr == "BEER")
		{
			this->extractTranscript_BEER(annotation_ifs, indexInfo);
		}
		else if(formatStr == "GTF")
		{
			this->extractTranscript_GTF(annotation_ifs, indexInfo);
		}
		else
		{
			cout << "format invalid !" << endl;
			exit(1);
		}
	}	
	void extractTranscript_GTF(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		/*while(1)
		{
			while(!annotation_ifs.eof())
			{
				string tmpAnnotationEntry;
				getline(annotation_ifs, tmpAnnotationEntry);
				if(tmpAnnotationEntry == "")
					break;
			}

		}*/
	}

	void extractTranscript_GAF(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		//cout << "start to extractTranscript from GAF file" << endl;
		int tmp_index = 0;
		while(!annotation_ifs.eof())
		{
			string tmpAnnotationEntry;
			getline(annotation_ifs, tmpAnnotationEntry);
			if(tmp_index == 0)
			{	
				if(tmpAnnotationEntry.at(0) == '#')
					continue;
			}
			if(tmpAnnotationEntry == "")
				break;
			Normal_Transcript_Info* normal_transcript_info = new Normal_Transcript_Info();
			bool parse_bool = normal_transcript_info->Normal_Transcript_Info_parseGAF(tmpAnnotationEntry, indexInfo);
			if(!parse_bool)
			{
				continue;
			}
			normal_transcript_vec.push_back(normal_transcript_info);
			transcriptSetIndexMap.insert(pair<string, int> (
				normal_transcript_info->returnTranscriptName(), tmp_index));
			tmp_index ++;
		}
		return;
	}

	void extractTranscript_BEER(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		int tmp_index = 0; 
		while(!annotation_ifs.eof())
		{
			//bool newGeneBool = false;
			if(annotation_ifs.eof())
			{
				break;
			}
			string headLineEachGeneStr;
			getline(annotation_ifs, headLineEachGeneStr);
			//cout << "headLineEachGeneStr: " << headLineEachGeneStr << " " << headLineEachGeneStr.length() << endl;
			if(headLineEachGeneStr.length() <= 0)
				return;
			if(headLineEachGeneStr == "---------")
			{
				string geneNameStr, exonStartStr, exonEndStr, strandStr, chrStr;
				getline(annotation_ifs, geneNameStr);
				//cout << "geneNameStr: " << geneNameStr << endl;
				getline(annotation_ifs, exonStartStr);
				//cout << "exonStartStr: " << exonStartStr << endl;
				getline(annotation_ifs, exonEndStr);
				//cout << "exonEndStr: " << exonEndStr << endl;
				getline(annotation_ifs, strandStr);
				//cout << "strandStr: " << strandStr << endl;
				getline(annotation_ifs, chrStr);
				//cout << "chrStr: " << chrStr << endl;
				Normal_Transcript_Info* normal_transcript_info = new Normal_Transcript_Info();
				bool parse_bool = normal_transcript_info->Normal_Transcript_Info_parseBEER(
					geneNameStr, exonStartStr, exonEndStr, strandStr, chrStr, indexInfo);
				if(!parse_bool)
				{
					continue;
				}
				normal_transcript_vec.push_back(normal_transcript_info);
				transcriptSetIndexMap.insert(pair<string, int> (
					normal_transcript_info->returnTranscriptName(), tmp_index));
				//cout << "transcript_name: " << normal_transcript_info->returnTranscriptName()
				//	<< endl << "tmp_index: " << tmp_index << endl;
				tmp_index ++;
			}
			else
			{}
			//cout << "turn ..." << endl;
		}
	}

	string returnTranscriptSeq(int index, Index_Info* indexInfo)
	{
		return normal_transcript_vec[index]->getTranscriptSeq_noDir(indexInfo);
	}

	void outputRefSeq(const string& outputRefFolderStr, Index_Info* indexInfo, int baseNumInEachLine)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			normal_transcript_vec[tmp]->outputRefSeq(outputRefFolderStr, indexInfo, baseNumInEachLine);
		}		
	}

	void outputTranscriptInfo(ofstream& transcriptInfo_ofs)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			normal_transcript_vec[tmp]->outputTranscriptInfo(transcriptInfo_ofs);
		}				
	}

	int returnTranscriptNum()
	{
		return normal_transcript_vec.size();
	}

	Normal_Transcript_Info* returnTranscript(int vec_index)
	{
		return normal_transcript_vec[vec_index];
	}

	void memoryFree()
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
			delete normal_transcript_vec[tmp];
	}
};

#endif