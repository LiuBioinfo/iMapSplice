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

	void output_GTF(ofstream& GTF_ofs, const string& annotation_name)
	{
		//cout << "normal_transcript_vec.size(): " << normal_transcript_vec.size() << endl;
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl; 
			GTF_ofs << normal_transcript_vec[tmp]->returnTranscriptInfo_GTF(annotation_name) << endl;
		}		
	}

	void output_GAF(ofstream& GAF_ofs)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			GAF_ofs << normal_transcript_vec[tmp]->returnTranscriptInfo_GAF() << endl;
		}
	}

	void output_GAF_standardID(ofstream& GAF_ofs)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			GAF_ofs << normal_transcript_vec[tmp]->returnTranscriptInfo_GAF_standardID(tmp+1) << endl;
		}
	}

	void output_reissuedTranscriptID_geneName_oriTranscriptName(string& reissuedTranscriptID_geneName_oriTranscriptName_path)
	{
		ofstream reissuedTranscriptID_geneName_oriTranscriptName_ofs(reissuedTranscriptID_geneName_oriTranscriptName_path.c_str());
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			string tmpTranscriptName = normal_transcript_vec[tmp]->returnTranscriptName();
			string tmpReissuedTranscriptID = "id_" + int_to_str(tmp+1);
			string tmpGeneName = normal_transcript_vec[tmp]->returnGeneName();
			reissuedTranscriptID_geneName_oriTranscriptName_ofs << tmpReissuedTranscriptID << "\t" << tmpGeneName << "\t" << tmpTranscriptName << endl;
		}
		reissuedTranscriptID_geneName_oriTranscriptName_ofs.close();
	}

	void memoryFree()
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			delete normal_transcript_vec[tmp];
		}
	}
	/*void outputTranscriptInfo(ofstream& output_ofs)
	{
		for(int tmp = 0; tmp < normal_transcript_vec.size(); tmp++)
		{
			output_ofs << normal_transcript_vec[tmp]->returnNormalTranscriptInfo() << endl;
		}
	}*/
	int returnMapPosInChrWithMapPosInTranscript(int index, int mapPosInTranscript)
	{
		return normal_transcript_vec[index]->returnMapPosInChrWithMapPosInTranscript(
			mapPosInTranscript);
	}

	string returnStrand(int index)
	{
		return normal_transcript_vec[index]->returnStrand();
	}

	string returnChromName(int index)
	{
		return normal_transcript_vec[index]->returnChromName();
	}

	int returnTranscriptStartPosInChr(int index)
	{
		return normal_transcript_vec[index]->returnTranscriptStartPosInChr();
	}

	int returnTranscriptEndPosInChr(int index)
	{
		return normal_transcript_vec[index]->returnTranscriptEndPosInChr();
	}

	void getSJsitePairVec(int index, vector< pair<int,int> >& SJsitePairVec)
	{
		normal_transcript_vec[index]->generateSJsitePairVec(SJsitePairVec);
	}

	void getSJsiteAndSize(int index, int alignmentStartPosInTranscript, 
		int alignmentEndPosInTranscript,
		vector< pair<int,int> >& SJsiteAndSize) // vector < SJsiteInExonArea, SJsize >
	{
		normal_transcript_vec[index]->generateSJsiteAndSizeInAlignment(
			alignmentStartPosInTranscript,
			alignmentEndPosInTranscript, SJsiteAndSize);
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

	void extractTranscript_GAF(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		int tmp_index = 0;
		string headLine;
		getline(annotation_ifs, headLine);		
		//int transcript_set_map_index = 0;
		while(!annotation_ifs.eof())
		{
			//cout << "transcript_index: " << tmp_index + 1 << endl;
			string tmpAnnotationEntry;
			getline(annotation_ifs, tmpAnnotationEntry);
			if(tmpAnnotationEntry == "")
				break;
			//cout << "tmpAnnotationEntry: " << tmpAnnotationEntry << endl;
			Normal_Transcript_Info* normal_transcript_info = new Normal_Transcript_Info();
			bool parse_bool = normal_transcript_info->Normal_Transcript_Info_parseGAF(tmpAnnotationEntry, indexInfo);
			//cout << "parse_bool: " << parse_bool << endl;
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

	string getTranscriptIDfromGTFline(const string& lineStr)
	{
		int lastTabLocation = lineStr.rfind("\t");
		return lineStr.substr(lastTabLocation+1);
	}

	void extractTranscript_GTF(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		//cout << "start to extract transcript from GTF file" << endl;
		string transcript_id_old;
		vector<string> transcriptLineVec;

		string firstLineStr;
		getline(annotation_ifs, firstLineStr);
		transcript_id_old = this->getTranscriptIDfromGTFline(firstLineStr);
		transcriptLineVec.push_back(firstLineStr);
		//string transcript_id_current;
		int tmpIndex = 0;
		//cout << "start to read each line from file" << endl;
		while(1)
		{
			if(annotation_ifs.eof())
			{
				Normal_Transcript_Info* normal_transcript_info = new Normal_Transcript_Info();
				bool parse_bool = normal_transcript_info->Normal_Transcript_Info_parseGTF(
					transcriptLineVec, indexInfo);
				if(!parse_bool)
				{
					break;
					//continue;
				}
				normal_transcript_vec.push_back(normal_transcript_info);
				transcriptSetIndexMap.insert(pair<string,int> (
					normal_transcript_info->returnTranscriptName(), tmpIndex));
				tmpIndex ++;
				break;
			}
			string tmpLineStr;
			getline(annotation_ifs, tmpLineStr);
			//cout << "tmpLineStr: " << endl << tmpLineStr << endl;
			string tmpTranscriptID = this->getTranscriptIDfromGTFline(tmpLineStr);
			//cout << "tmpTranscriptID: " << endl << tmpTranscriptID << endl;
			if(tmpTranscriptID == transcript_id_old)
			{
				transcriptLineVec.push_back(tmpLineStr);
			}
			else
			{
				//cout << "old transcript id" << endl << transcript_id_old << endl;
				//cout << "new transcript id" << endl << tmpTranscriptID << endl;
				//cout << "new transcript info " << endl;
				Normal_Transcript_Info* normal_transcript_info = new Normal_Transcript_Info();
				bool parse_bool = normal_transcript_info->Normal_Transcript_Info_parseGTF(
					transcriptLineVec, indexInfo);
				//cout << "parse_bool: " << endl << parse_bool << endl;
				if(!parse_bool)
				{
					continue;
				}
				normal_transcript_vec.push_back(normal_transcript_info);
				transcriptSetIndexMap.insert(pair<string,int> (
					normal_transcript_info->returnTranscriptName(), tmpIndex));
				transcriptLineVec.clear();
				transcriptLineVec.push_back(tmpLineStr);
				transcript_id_old = tmpTranscriptID;
				tmpIndex ++;
			}
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

	void extractTranscript_BEER_outputInvalidTranscriptInfo(Index_Info* indexInfo,
		const string& inputBeersFile, const string& invalidBeersTranscriptInfoOutputFile,
		ofstream& log_ofs)
	{
		ifstream annotation_ifs(inputBeersFile.c_str());
		ofstream invalidBeersTranscriptInfo_ofs(invalidBeersTranscriptInfoOutputFile.c_str());

		bool validBeersTranscript_bool = true;
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
					validBeersTranscript_bool = false;
					invalidBeersTranscriptInfo_ofs << headLineEachGeneStr << endl;
					invalidBeersTranscriptInfo_ofs << geneNameStr << endl;
					invalidBeersTranscriptInfo_ofs << exonStartStr << endl;
					invalidBeersTranscriptInfo_ofs << exonEndStr << endl;
					invalidBeersTranscriptInfo_ofs << strandStr << endl;
					invalidBeersTranscriptInfo_ofs << chrStr << endl;
					continue;
				}
				validBeersTranscript_bool = true;
				normal_transcript_vec.push_back(normal_transcript_info);
				transcriptSetIndexMap.insert(pair<string, int> (
					normal_transcript_info->returnTranscriptName(), tmp_index));
				//cout << "transcript_name: " << normal_transcript_info->returnTranscriptName()
				//	<< endl << "tmp_index: " << tmp_index << endl;
				tmp_index ++;
			}
			else
			{
				if(!validBeersTranscript_bool)
				{
					invalidBeersTranscriptInfo_ofs << headLineEachGeneStr << endl;
				}
			}
			//cout << "turn ..." << endl;
		}

		annotation_ifs.close();
		invalidBeersTranscriptInfo_ofs.close();
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
};

#endif