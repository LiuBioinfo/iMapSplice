// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//#include "otherFunc.h"

typedef map< pair<int, int>, pair<int, vector< pair<int, int> > > > CandidateFusionMap;
// pair<int, int> (left-end readMapPosAreaNO, right-end readMapPos2ndLevelNO), pair<int, pair<int,int> > (voting score, pair<1stEndRead_endMapPos, 2ndEndRead_startMapPos>)    
typedef map< int, vector<pair<int, pair<int, int> > > > FusionEndSet; // <int, int> (mapPosAreaNO_1, (mapPosAreaNO_2, 1stEndRead_endMapPos_max) )

typedef map< pair<int, int>, int > CandidateFusionMap2ArrayIndex;

typedef map< pair<int, int>, int > CandidateFusionMap_type;

class FusionDetection_Info
{

public:
	CandidateFusionMap fusionCandidate;//fusionCandidate_ForwForw, fusionCandidate_ForwReve, fusionCandidate_ReveForw, fusionCandidate_ReveReve,
		//fusionCandidate_Forw, fusionCandidate_Reve; // forw: end1 -- end2; reve: end2 -- end1

	FusionEndSet candidateFusionEndSet_1, candidateFusionEndSet_2; 
		//candidateFusionEndSet_Forw_1, candidateFusionEndSet_Forw_2,
		//candidateFusionEndSet_Reve_1, candidateFusionEndSet_Reve_2;

	CandidateFusionMap2ArrayIndex candidateFusionArrayIndexMap;
		//candidateFusionArrayIndexMap_Forw, candidateFusionArrayIndexMap_Reve;

	vector< vector< pair<int, int> > > candidateFusionExactSiteVec;
		//candidateFusionExactSiteVec_Forw, candidateFusionExactSiteVec_Reve;
	
	int minSupportNum;
	int mapIntervalHalfSize;// = 300000;
	//vector< pair<int, int> > candidateFusionVec;// vector<pair<int,int>(2ndLevelIndexNO_1, 2ndLevelIndexNO_2)>

	FusionDetection_Info()
	{
		minSupportNum = 20;
		mapIntervalHalfSize = 300000;
	}

	string convertStringToReverseComplement(const string& originalString)
	{
		int stringLength = originalString.size();
		string resultString = this->covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
		for (int tmp = 1; tmp < stringLength; tmp++)
		{
			resultString = resultString + this->covertCharToReverseComplement(
				originalString.substr(stringLength-1-tmp, 1));
		}
		return resultString;
	}

	string covertCharToReverseComplement(const string& Ori_Char)
	{
		if(Ori_Char == "A")
		{
			return "T";
		}
		else if(Ori_Char == "T")
		{
			return "A";
		}
		else if(Ori_Char == "G")
		{
			return "C";
		}
		else if(Ori_Char == "C")
		{
			return "G";
		}
		else if(Ori_Char == "N")
		{
			return "N";
		}	
	}

	void generateCandidateFusionMap(const string& fileName, Index_Info* indexInfo)
	{
		ifstream inputRecord_ifs(fileName.c_str()); // unpaired alignments file
		string line1, line2, line3, line4, line5, line6, line7,
			line8, line9, line10, line11; 
		getline(inputRecord_ifs, line1);
		int readNO = 0;
		while(!inputRecord_ifs.eof())
		{
			readNO ++;
			getline(inputRecord_ifs, line1);
			getline(inputRecord_ifs, line2);
			getline(inputRecord_ifs, line3);
			getline(inputRecord_ifs, line4);
			getline(inputRecord_ifs, line5);
			getline(inputRecord_ifs, line6);
			getline(inputRecord_ifs, line7);
			getline(inputRecord_ifs, line8);
			getline(inputRecord_ifs, line9);
			getline(inputRecord_ifs, line10);
			getline(inputRecord_ifs, line11);

				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
				string readNameStr_1, readNameStr_2;

				int startSearchPos = 0, foundSearchPos;
				foundSearchPos = line1.find("\t", startSearchPos);
				readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1.find("\t", startSearchPos);
				Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1.find("\t", startSearchPos);		
				Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

				startSearchPos = 0; 
				foundSearchPos = line4.find("\t", startSearchPos);
				readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4.find("\t", startSearchPos);
				Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4.find("\t", startSearchPos);		
				Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

				int readLength_1 = line2.length() - 1;
				int readLength_2 = line5.length() - 1;

				PE_Read_Info* peReadInfo = new PE_Read_Info();

				peReadInfo->getFastaFormatReadInfo(readNameStr_1, readNameStr_2,
					line2.substr(0, readLength_1), line5.substr(0, readLength_2));	

				char* read = const_cast<char*>((peReadInfo->readInfo_pe1).readSeq.c_str());
				char* read_PE = const_cast<char*>((peReadInfo->readInfo_pe2).readSeq.c_str());

 		   		char read_RC_PE[(peReadInfo->readInfo_pe1).readSeqLength], read_RC[(peReadInfo->readInfo_pe2).readSeqLength];
    			for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe1).readSeqLength; read_RC_loc++) // get read_RC
    				*(read_RC + read_RC_loc) = reverseComplement(*(read + (peReadInfo->readInfo_pe1).readSeqLength - 1 - read_RC_loc)); 			
    	 		for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe2).readSeqLength; read_RC_loc++) // get read_RC_PE
    				*(read_RC_PE + read_RC_loc) = reverseComplement(*(read_PE + (peReadInfo->readInfo_pe2).readSeqLength - 1 - read_RC_loc)); 

	    		string rcmReadSeq_1 = read_RC;
    			(peReadInfo->readInfo_pe1).rcmReadSeq = rcmReadSeq_1.substr(0, (peReadInfo->readInfo_pe1).readSeqLength);

	    		string rcmReadSeq_2 = read_RC_PE;
    			(peReadInfo->readInfo_pe2).rcmReadSeq = rcmReadSeq_2.substr(0, (peReadInfo->readInfo_pe2).readSeqLength);


				PE_Read_Alignment_Info* peAlignInfo = 
					new PE_Read_Alignment_Info(line7, line8, line9, line10,
					Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);

			//peAlignInfo->pairingAlignment();
			//peAlignInfo->chooseBestAlignment();		

			peAlignInfo->chooseBestAlignment_selectRandomOneIfMulti();

			bool pairExistsBool = peAlignInfo->finalPairExistsBool();

			if(!pairExistsBool)
			{
				vector< pair< int, pair<int, int> > > alignChromPosVec_1;
				vector< pair< int, pair<int, int> > > alignChromPosVec_2;
				
				// to fix
				if(peAlignInfo->norAlignmentInfo_PE_1.size() + peAlignInfo->rcmAlignmentInfo_PE_1.size()
					+ peAlignInfo->norAlignmentInfo_PE_2.size() + peAlignInfo->rcmAlignmentInfo_PE_2.size() > 2)
				{
					//continue;
				}	

				for(int tmp = 0; tmp < peAlignInfo->norAlignmentInfo_PE_1.size(); tmp ++)
				{
					int chromNameInt = indexInfo->convertStringToInt((peAlignInfo->norAlignmentInfo_PE_1)[tmp]->alignChromName);
					int startPos = (peAlignInfo->norAlignmentInfo_PE_1)[tmp]->alignChromPos;
					int endPos = (peAlignInfo->norAlignmentInfo_PE_1)[tmp]->endMatchedPosInChr;

					int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, startPos);

					alignChromPosVec_1.push_back( pair<int, pair<int, int> > (secondLevelIndexNO, pair<int,int>(startPos, endPos) ) );
				}

				for(int tmp = 0; tmp < peAlignInfo->rcmAlignmentInfo_PE_1.size(); tmp ++)
				{
					int chromNameInt = indexInfo->convertStringToInt((peAlignInfo->rcmAlignmentInfo_PE_1)[tmp]->alignChromName);
					int startPos = (peAlignInfo->rcmAlignmentInfo_PE_1)[tmp]->alignChromPos;
					int endPos = (peAlignInfo->rcmAlignmentInfo_PE_1)[tmp]->endMatchedPosInChr;

					int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, startPos);

					alignChromPosVec_1.push_back( pair<int, pair<int, int> > (secondLevelIndexNO, pair<int,int>(startPos, endPos) ) );
				}

				for(int tmp = 0; tmp < peAlignInfo->norAlignmentInfo_PE_2.size(); tmp ++)
				{
					int chromNameInt = indexInfo->convertStringToInt((peAlignInfo->norAlignmentInfo_PE_2)[tmp]->alignChromName);
					int startPos = (peAlignInfo->norAlignmentInfo_PE_2)[tmp]->alignChromPos;
					int endPos = (peAlignInfo->norAlignmentInfo_PE_2)[tmp]->endMatchedPosInChr;

					int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, startPos);

					alignChromPosVec_2.push_back( pair<int, pair<int, int> > (secondLevelIndexNO, pair<int,int>(startPos, endPos) ) );
				}

				for(int tmp = 0; tmp < peAlignInfo->rcmAlignmentInfo_PE_2.size(); tmp ++)
				{
					int chromNameInt = indexInfo->convertStringToInt((peAlignInfo->rcmAlignmentInfo_PE_2)[tmp]->alignChromName);
					int startPos = (peAlignInfo->rcmAlignmentInfo_PE_2)[tmp]->alignChromPos;
					int endPos = (peAlignInfo->rcmAlignmentInfo_PE_2)[tmp]->endMatchedPosInChr;

					int secondLevelIndexNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, startPos);

					alignChromPosVec_2.push_back( pair<int, pair<int, int> > (secondLevelIndexNO, pair<int,int>(startPos, endPos) ) );
				}



				for(int tmp = 0; tmp < alignChromPosVec_1.size(); tmp++)
				{
					for(int tmp2 = 0; tmp2 < alignChromPosVec_2.size(); tmp2++)
					{
						int tmp2ndLevelIndexNO_1 = alignChromPosVec_1[tmp].first;
						int tmp2ndLevelIndexNO_2 = alignChromPosVec_2[tmp2].first;

						if(tmp2ndLevelIndexNO_1 != tmp2ndLevelIndexNO_2)
						{
							/*
							this->insertCandidateFusionMap_Forw(tmp2ndLevelIndexNO_1, 
								tmp2ndLevelIndexNO_2, (alignChromPosVec_1[tmp].second).second,
								 (alignChromPosVec_2[tmp].second).first);
							this->insertCandidateFusionMap_Reve(tmp2ndLevelIndexNO_2, 
								tmp2ndLevelIndexNO_1, (alignChromPosVec_2[tmp].second).second,
								 (alignChromPosVec_1[tmp].second).first);*/
							/*if((tmp2ndLevelIndexNO_1 == 0)||(tmp2ndLevelIndexNO_2 == 0)
								||((alignChromPosVec_1[tmp].second).second <= 0)
								||((alignChromPosVec_2[tmp2].second).first <= 0)
								||((alignChromPosVec_2[tmp2].second).second <= 0)
								||((alignChromPosVec_1[tmp].second).first <= 0))
							{

								cout << line1 << endl << line4 << endl << line7 << endl << line8 << endl << line9 << endl << line10 << endl;
								cout << "tmp2ndLevelIndexNO_1 "<< tmp2ndLevelIndexNO_1 << endl 
									<< "tmp2ndLevelIndexNO_2 "<< tmp2ndLevelIndexNO_2 << endl
									<< "(alignChromPosVec_1[tmp].second).second " << (alignChromPosVec_1[tmp].second).second << endl
									<< "(alignChromPosVec_2[tmp2].second).first " << (alignChromPosVec_2[tmp2].second).first << endl
									<< "(alignChromPosVec_2[tmp2].second).second " << (alignChromPosVec_2[tmp2].second).second << endl
									<< "(alignChromPosVec_1[tmp].second).first " << (alignChromPosVec_1[tmp].second).first << endl; 

								return;
							}*/
							this->insertCandidateFusionMap(tmp2ndLevelIndexNO_1, 
								tmp2ndLevelIndexNO_2, (alignChromPosVec_1[tmp].second).second,
								 (alignChromPosVec_2[tmp2].second).first);
							this->insertCandidateFusionMap(tmp2ndLevelIndexNO_2, 
								tmp2ndLevelIndexNO_1, (alignChromPosVec_2[tmp2].second).second,
								 (alignChromPosVec_1[tmp].second).first);
						}
					}
				}
			}
			else
			{
				//cout << " pairedReads: " << readNameStr_1 << endl;
			}
				
		}
		inputRecord_ifs.close();
	}

	void insertCandidateFusionMap(int secondLevelIndexNO_1, int secondLevelIndexNO_2, int tmpChrPos_1, int tmpChrPos_2)
	{
		CandidateFusionMap::iterator fusionMapIter;
		fusionMapIter = fusionCandidate.find(pair<int,int> (secondLevelIndexNO_1, secondLevelIndexNO_2));
		if(fusionMapIter 
			!= fusionCandidate.end())
		{
			(fusionMapIter->second).first ++;
			((fusionMapIter->second).second).push_back(pair<int,int>(tmpChrPos_1, tmpChrPos_2));
		}
		else
		{
			vector< pair<int,int> > tmpVec;
			tmpVec.push_back (pair<int,int>(tmpChrPos_1, tmpChrPos_2));

			fusionCandidate.insert(pair< pair<int, int>, pair<int, vector<pair<int,int> > > > 
				(pair<int,int>(secondLevelIndexNO_1, secondLevelIndexNO_2),  pair<int, vector<pair<int,int> > >(1, tmpVec) ));
			//fusionCandidate_type.insert( pair< pair<int,int>, int > (pair<int,int>(norEnd2ndLevelIndexNO, rcmEnd2ndLevelIndexNO), type ) );
		}
	}

	string getCandidateFusionMapStr()
	{
		string fusionInfoStr;
		CandidateFusionMap::iterator tmpIter;
		int elementNo = 0;

		for(tmpIter = fusionCandidate.begin(); tmpIter != fusionCandidate.end(); tmpIter ++)
		{
			int tmp2ndLevelIndexNO_1 = (tmpIter->first).first;
			int tmp2ndLevelIndexNO_2 = (tmpIter->first).second;
			int tmpScore = (tmpIter->second).first;

			if(tmpScore < minSupportNum)
				continue;

			fusionInfoStr = fusionInfoStr + "1stSiteAreaNO: " + int_to_str(tmp2ndLevelIndexNO_1) 
				+ " 2ndSiteAreaNO: " + int_to_str(tmp2ndLevelIndexNO_2) + " FORWARD" + " vote: " + int_to_str(tmpScore) + " FusionSJposPair: ";

			for(int tmp = 0; tmp < ((tmpIter->second).second).size(); tmp ++)
			{
				int tmpMapPos_1 = ((tmpIter->second).second)[tmp].first;
				int tmpMapPos_2 = ((tmpIter->second).second)[tmp].second;

				fusionInfoStr = fusionInfoStr + int_to_str(tmpMapPos_1) + "~" + int_to_str(tmpMapPos_2) + ",";
			}						

			fusionInfoStr += "\n";
		}

		return fusionInfoStr;
	}

	void generateCandidateFusionSetAndCandidateFusionArrayIndexMap()
	{
		CandidateFusionMap::iterator tmpIter;
		FusionEndSet::iterator tmpFusionEndSetIter_1, tmpFusionEndSetIter_2;
		
		int tmpArrayIndex = 0;//_Forw = 0, tmpArrayIndex_Reve = 0;//, tmpArrayIndex_ReveForw = 0, tmpArrayIndex_ReveReve = 0;

		for(tmpIter = fusionCandidate.begin(); tmpIter != fusionCandidate.end(); tmpIter ++)
		{
			int tmp2ndLevelIndexNO_1 = (tmpIter->first).first;
			int tmp2ndLevelIndexNO_2 = (tmpIter->first).second;
			int tmpScore = (tmpIter->second).first;

			if(tmpScore >= minSupportNum)
			{
			    // generate CandidateFusionArrayIndexMap
				candidateFusionArrayIndexMap.insert(pair< pair<int,int>, int > (
					pair<int,int>(tmp2ndLevelIndexNO_1, tmp2ndLevelIndexNO_2), tmpArrayIndex));
				vector< pair<int,int> > newTmpVec;
				candidateFusionExactSiteVec.push_back(newTmpVec);
				tmpArrayIndex ++;

				//generate candidateFusionSet
				int tmpEnd1Pos_Max = ((tmpIter->second).second)[0].first;
				int tmpEnd2Pos_Min = ((tmpIter->second).second)[0].second;

				for(int tmp = 1; tmp < ((tmpIter->second).second).size(); tmp++)
				{
					int tmpEnd1Pos = ((tmpIter->second).second)[tmp].first;
					int tmpEnd2Pos = ((tmpIter->second).second)[tmp].second;

					if(tmpEnd1Pos > tmpEnd1Pos_Max)
						tmpEnd1Pos_Max = tmpEnd1Pos;
					if(tmpEnd2Pos < tmpEnd2Pos_Min)
						tmpEnd2Pos_Min = tmpEnd2Pos;
				}


				//candidateFusionVec.push_back(pair<int,int>(tmp2ndLevelIndexNO_1, tmp2ndLevelIndexNO_2));
				tmpFusionEndSetIter_1 = candidateFusionEndSet_1.find(tmp2ndLevelIndexNO_1);
				tmpFusionEndSetIter_2 = candidateFusionEndSet_2.find(tmp2ndLevelIndexNO_2);

				// generate candidateFusionEndSet_1;
				if(tmpFusionEndSetIter_1 != candidateFusionEndSet_1.end())// found;
				{
					(tmpFusionEndSetIter_1->second).push_back(
						pair<int, pair<int, int> >(tmp2ndLevelIndexNO_2, pair<int,int>(tmpEnd1Pos_Max, tmpEnd2Pos_Min) ) );
				}
				else // new fusionSplicePair 
				{
					vector<pair<int, pair<int, int> > > tmpVec;
					tmpVec.push_back(pair<int, pair<int, int> >(tmp2ndLevelIndexNO_2, pair<int,int>(tmpEnd1Pos_Max, tmpEnd2Pos_Min) ) );
					candidateFusionEndSet_1.insert(pair<int, vector<pair<int, pair<int, int> > > > (tmp2ndLevelIndexNO_1, tmpVec));
				}

				// generate candidateFusionEndSet_2;
				if(tmpFusionEndSetIter_2 != candidateFusionEndSet_2.end())// found;
				{
					(tmpFusionEndSetIter_2->second).push_back(
						pair<int, pair<int, int> >(tmp2ndLevelIndexNO_1, pair<int,int>(tmpEnd1Pos_Max, tmpEnd2Pos_Min) ) );
				}
				else
				{
					vector<pair<int, pair<int, int> > > tmpVec;
					tmpVec.push_back(pair<int, pair<int, int> >(tmp2ndLevelIndexNO_1, pair<int,int>(tmpEnd1Pos_Max, tmpEnd2Pos_Min) ) );
					candidateFusionEndSet_2.insert(pair<int, vector<pair<int, pair<int, int> > > > (tmp2ndLevelIndexNO_2, tmpVec));
				}

			}
		}	

	}

	string getFusionSetStr()
	{
		string fusionInfoStr;

		FusionEndSet::iterator tmpSetIter;
		
		fusionInfoStr = fusionInfoStr + "candidateFusionEndSet_1 info: \n"; 
		for(tmpSetIter = candidateFusionEndSet_1.begin(); tmpSetIter != candidateFusionEndSet_1.end(); tmpSetIter ++)
		{
			int tmp2ndLevelIndexNO_1 = tmpSetIter->first;
			fusionInfoStr = fusionInfoStr + "1stSiteAreaNO: " + int_to_str(tmp2ndLevelIndexNO_1) + "\n";
			for(int tmp = 0; tmp < (tmpSetIter->second).size(); tmp ++)
			{
				fusionInfoStr = fusionInfoStr + "...2ndSiteAreaNO: " +  int_to_str((tmpSetIter->second)[tmp].first) + " FORWARD_1 "
					+ " tmpPos: " + int_to_str(((tmpSetIter->second)[tmp].second).first) + "~" 
					+  int_to_str(((tmpSetIter->second)[tmp].second).second) + "\n"; 
			}
		}

		fusionInfoStr = fusionInfoStr + "\n\ncandidateFusionEndSet_2 info: \n"; 
		for(tmpSetIter = candidateFusionEndSet_2.begin(); tmpSetIter != candidateFusionEndSet_2.end(); tmpSetIter ++)
		{
			int tmp2ndLevelIndexNO_1 = tmpSetIter->first;
			fusionInfoStr = fusionInfoStr + "2ndSiteAreaNO: " + int_to_str(tmp2ndLevelIndexNO_1) + "\n";
			for(int tmp = 0; tmp < (tmpSetIter->second).size(); tmp ++)
			{
				fusionInfoStr = fusionInfoStr + "...1stSiteAreaNO: " +  int_to_str((tmpSetIter->second)[tmp].first) + " FORWARD_2 "
					+ " tmpPos: " + int_to_str(((tmpSetIter->second)[tmp].second).first) + "~" 
					+  int_to_str(((tmpSetIter->second)[tmp].second).second) + "\n"; 
			}
		}

		return fusionInfoStr;
	}

	void detectExactFusionFromPairedAlignFile(
		const string& inputPairedAlignFile,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild, Index_Info* indexInfo)
	{
		ifstream inputRecord_ifs(inputPairedAlignFile.c_str()); // unpaired alignments file
		string line1, line2, line3, line4, line5, line6, line7,
			line8, line9, line10, line11; 
		getline(inputRecord_ifs, line1);
		int readNO = 0;
		while(!inputRecord_ifs.eof())
		{
			readNO ++;
			getline(inputRecord_ifs, line1);
			getline(inputRecord_ifs, line2);
			getline(inputRecord_ifs, line3);
			getline(inputRecord_ifs, line4);
			getline(inputRecord_ifs, line5);
			getline(inputRecord_ifs, line6);
			getline(inputRecord_ifs, line7);
			getline(inputRecord_ifs, line8);
			getline(inputRecord_ifs, line9);
			getline(inputRecord_ifs, line10);
			getline(inputRecord_ifs, line11);

				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int Nor1Num = 0, Rcm1Num = 0, Nor2Num = 0, Rcm2Num = 0;
				string readNameStr_1, readNameStr_2;

				int startSearchPos = 0, foundSearchPos;
				foundSearchPos = line1.find("\t", startSearchPos);
				readNameStr_1 = line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1.find("\t", startSearchPos);
				Nor1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line1.find("\t", startSearchPos);		
				Rcm1Num = atoi((line1.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());			

				startSearchPos = 0; 
				foundSearchPos = line4.find("\t", startSearchPos);
				readNameStr_2 = line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1);
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4.find("\t", startSearchPos);
				Nor2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());		
				startSearchPos = foundSearchPos + 1;
				foundSearchPos = line4.find("\t", startSearchPos);		
				Rcm2Num = atoi((line4.substr(startSearchPos, foundSearchPos-1-startSearchPos+1)).c_str());	

				int readLength_1 = line2.length() - 1;
				int readLength_2 = line5.length() - 1;

				PE_Read_Info* peReadInfo = new PE_Read_Info();

				peReadInfo->getFastaFormatReadInfo(readNameStr_1, readNameStr_2,
					line2.substr(0, readLength_1), line5.substr(0, readLength_2));	

				char* read = const_cast<char*>((peReadInfo->readInfo_pe1).readSeq.c_str());
				char* read_PE = const_cast<char*>((peReadInfo->readInfo_pe2).readSeq.c_str());

 		   		char read_RC_PE[(peReadInfo->readInfo_pe1).readSeqLength], read_RC[(peReadInfo->readInfo_pe2).readSeqLength];
    			for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe1).readSeqLength; read_RC_loc++) // get read_RC
    				*(read_RC + read_RC_loc) = reverseComplement(*(read + (peReadInfo->readInfo_pe1).readSeqLength - 1 - read_RC_loc)); 			
    	 		for(int read_RC_loc = 0; read_RC_loc < (peReadInfo->readInfo_pe2).readSeqLength; read_RC_loc++) // get read_RC_PE
    				*(read_RC_PE + read_RC_loc) = reverseComplement(*(read_PE + (peReadInfo->readInfo_pe2).readSeqLength - 1 - read_RC_loc)); 

	    		string rcmReadSeq_1 = read_RC;
    			(peReadInfo->readInfo_pe1).rcmReadSeq = rcmReadSeq_1.substr(0, (peReadInfo->readInfo_pe1).readSeqLength);

	    		string rcmReadSeq_2 = read_RC_PE;
    			(peReadInfo->readInfo_pe2).rcmReadSeq = rcmReadSeq_2.substr(0, (peReadInfo->readInfo_pe2).readSeqLength);


				PE_Read_Alignment_Info* peAlignInfo = 
					new PE_Read_Alignment_Info(line7, line8, line9, line10,
					Nor1Num, Rcm1Num, Nor2Num, Rcm2Num);

				PE_Fusion_Info* peFusionInfo = new PE_Fusion_Info();

			//peAlignInfo->pairingAlignment();
			//peAlignInfo->chooseBestAlignment();		

			peAlignInfo->chooseBestAlignment_selectRandomOneIfMulti();

			bool pairExistsBool = peAlignInfo->finalPairExistsBool();	
			if(pairExistsBool)
			{
				this->detectExactFusionSite_pairedAlignments(peFusionInfo, peAlignInfo, peReadInfo, indexInfo, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild );
			}
			else
			{
				cout << "error in detectExactFusionFromPairedAlignFile ... " << endl;
			}

			peFusionInfo->refineFusionAlignment_all(peAlignInfo, peReadInfo, indexInfo);
		}
		inputRecord_ifs.close();
	}

	void detectExactFusionSite_pairedAlignments(
		PE_Fusion_Info* peFusionInfo, 
		PE_Read_Alignment_Info* peAlignInfo, PE_Read_Info* peReadInfo, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{
		for(int tmp = 0; tmp < peAlignInfo->finalAlignPair_Nor1Rcm2.size(); tmp ++)
		{
			int index_Nor1 = ( peAlignInfo->finalAlignPair_Nor1Rcm2[tmp]).first;
			int index_Rcm2 = ( peAlignInfo->finalAlignPair_Nor1Rcm2[tmp]).second;
			
			this->detectExactFusionSite_PairedAlignments_Nor1Rcm2(peFusionInfo, 
				index_Nor1, index_Rcm2, 
				peAlignInfo, peReadInfo, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}

		for(int tmp = 0; tmp < peAlignInfo->finalAlignPair_Nor2Rcm1.size(); tmp ++)
		{
			int index_Nor2 = peAlignInfo->finalAlignPair_Nor2Rcm1[tmp].first;
			int index_Rcm1 = peAlignInfo->finalAlignPair_Nor2Rcm1[tmp].second;
			this->detectExactFusionSite_PairedAlignments_Nor2Rcm1(peFusionInfo, 
				index_Nor2, index_Rcm1, 
				peAlignInfo, peReadInfo, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}
	}

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2(
		PE_Fusion_Info* peFusionInfo,
	 	int index_Nor1, int index_Rcm2, 
		PE_Read_Alignment_Info* peAlignInfo, PE_Read_Info* peReadInfo, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{

		bool unfixedHeadBool_Nor1 = 
			(peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->unfixedHeadExistsBool();
		bool unfixedTailBool_Rcm2 = 
			(peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->unfixedTailExistsBool();

		if( (!unfixedHeadBool_Nor1) && (!unfixedTailBool_Rcm2) )
		{
			return;
		}	

		FusionEndSet::iterator fusionEndSetIter;

		string alignChromNameStr 
			= (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromName;

		int alignChromNameInt = indexInfo->convertStringToInt(alignChromNameStr);

		if(unfixedHeadBool_Nor1)	
		{
			this->detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedHead(
				peFusionInfo, peAlignInfo, peReadInfo, 
				index_Nor1, alignChromNameInt, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}	

		if(unfixedTailBool_Rcm2)
		{
			this->detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedTail(
				peFusionInfo, peAlignInfo, peReadInfo, 
				index_Rcm2, alignChromNameInt, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}
	}

	void detectExactFusionSite_PairedAlignments_Nor2Rcm1(
		PE_Fusion_Info* peFusionInfo,
	 	int index_Nor2, int index_Rcm1, 
		PE_Read_Alignment_Info* peAlignInfo, PE_Read_Info* peReadInfo, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{

		bool unfixedHeadBool_Nor2 = 
			(peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->unfixedHeadExistsBool();
		bool unfixedTailBool_Rcm1 = 
			(peAlignInfo->rcmAlignmentInfo_PE_1)[index_Rcm1]->unfixedTailExistsBool();

		if( (!unfixedHeadBool_Nor2) && (!unfixedTailBool_Rcm1) )
		{
			return;
		}	

		FusionEndSet::iterator fusionEndSetIter;

		string alignChromNameStr 
			= (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->alignChromName;

		int alignChromNameInt = indexInfo->convertStringToInt(alignChromNameStr);

		if(unfixedHeadBool_Nor2)	
		{
			this->detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedHead(
				peFusionInfo, peAlignInfo, peReadInfo, 
				index_Nor2, alignChromNameInt, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}	

		if(unfixedTailBool_Rcm1)
		{
			this->detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedTail(
				peFusionInfo, peAlignInfo, peReadInfo, 
				index_Rcm1, alignChromNameInt, indexInfo, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild );
		}
	}

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedHead(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, 
		PE_Read_Info* peReadInfo, int index_Nor1, int chromNameInt, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{
			int alignChromPos_start 
				= (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromPos;
			int mapStartPos_areaNO 
				= indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, alignChromPos_start);
			FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_2.find(mapStartPos_areaNO);			
			if(fusionEndSetIter != candidateFusionEndSet_2.end())
			{
				int unfixedHeadLength= (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->unfixedHeadLength();
				string unfixedHeadSeq_nor =  ((peReadInfo->readInfo_pe1).readSeq).substr(0, unfixedHeadLength);
				string unfixedHeadSeq_rcm = this->convertStringToReverseComplement(unfixedHeadSeq_nor);
				for(int tmp = 0; tmp < (fusionEndSetIter->second).size(); tmp ++)
				{
					int targetRegion_areaNO = (fusionEndSetIter->second)[tmp].first;
					int targetRegion_endPos = ((fusionEndSetIter->second)[tmp].second).first;
					int mappedRegion_startPos = ((fusionEndSetIter->second)[tmp].second).second;

					
					int targetRegion_2ndLevelIndexNO = targetRegion_areaNO - 1;
					int otherFusionSJendChrNameInt = indexInfo->getChrNameIntFromSecondLevelIndexNO(targetRegion_2ndLevelIndexNO+1);
					int otherEndPosInCandidateFusionSet =  targetRegion_endPos;
					int otherEndMapIntervalStart = otherEndPosInCandidateFusionSet - mapIntervalHalfSize;
					int otherEndMapIntervalEnd = otherEndPosInCandidateFusionSet + mapIntervalHalfSize;
					if(otherEndMapIntervalStart <= 0)
					{
						otherEndMapIntervalStart = 1;
					}
					if(otherEndMapIntervalEnd >= ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedHeadLength - 1))
					{
						otherEndMapIntervalEnd = ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedHeadLength - 1);
					}

					string otherFusionSJendChrNameStr = (indexInfo->chrNameStr)[otherFusionSJendChrNameInt];
					int otherEndMapPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
						otherFusionSJendChrNameInt, (targetRegion_2ndLevelIndexNO+1), 1);

					char* unfixedHeadSeq_nor_char = const_cast<char*>(unfixedHeadSeq_nor.c_str());
					char* unfixedHeadSeq_rcm_char = const_cast<char*>(unfixedHeadSeq_rcm.c_str()); 
					// try to fix nor sequence
					//cout << "try to detect Fusion exact site ..." << endl;
					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo, 
						unfixedHeadSeq_nor_char, 
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_nor, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Nor1, true, true, true);

					// try to fix rcm sequence
					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo, unfixedHeadSeq_rcm_char,
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_rcm, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Nor1, true, true, false);
				}
			}
	}

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedTail(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, 
		PE_Read_Info* peReadInfo, int index_Rcm2, int chromNameInt, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{
		int alignChromPos_end = (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->endMatchedPosInChr;
		int mapEndPos_areaNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, alignChromPos_end);

		FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_1.find(mapEndPos_areaNO);	
		if(fusionEndSetIter != candidateFusionEndSet_1.end())
		{
			int unfixedTailLength = (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->unfixedTailLength();
			int readLength_2 = (peReadInfo->readInfo_pe2).rcmReadSeq.length();
			string unfixedTailSeq_nor =  ((peReadInfo->readInfo_pe2).rcmReadSeq).substr(readLength_2 - unfixedTailLength, unfixedTailLength);
			string unfixedTailSeq_rcm = this->convertStringToReverseComplement(unfixedTailSeq_nor);			
			
				for(int tmp = 0; tmp < (fusionEndSetIter->second).size(); tmp ++)
				{
					int targetRegion_areaNO = (fusionEndSetIter->second)[tmp].first;
					int targetRegion_startPos = ((fusionEndSetIter->second)[tmp].second).second;
					int mappedRegion_endPos = ((fusionEndSetIter->second)[tmp].second).first;

					
					int targetRegion_2ndLevelIndexNO = targetRegion_areaNO - 1;
					int otherFusionSJendChrNameInt = indexInfo->getChrNameIntFromSecondLevelIndexNO(targetRegion_2ndLevelIndexNO+1);
					int otherEndPosInCandidateFusionSet =  targetRegion_startPos;
					int otherEndMapIntervalStart = otherEndPosInCandidateFusionSet - mapIntervalHalfSize;
					int otherEndMapIntervalEnd = otherEndPosInCandidateFusionSet + mapIntervalHalfSize;
					if(otherEndMapIntervalStart <= 0)
					{
						otherEndMapIntervalStart = 1;
					}
					if(otherEndMapIntervalEnd >= ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedTailLength - 1))
					{
						otherEndMapIntervalEnd = ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedTailLength - 1);
					}

					string otherFusionSJendChrNameStr = (indexInfo->chrNameStr)[otherFusionSJendChrNameInt];
					int otherEndMapPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
						otherFusionSJendChrNameInt, (targetRegion_2ndLevelIndexNO+1), 1);

					char* unfixedTailSeq_nor_char = const_cast<char*>(unfixedTailSeq_nor.c_str());
					char* unfixedTailSeq_rcm_char = const_cast<char*>(unfixedTailSeq_rcm.c_str()); 
					// try to fix nor sequence

					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo,
						unfixedTailSeq_nor_char, 
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedTailSeq_nor, unfixedTailLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Rcm2,
						true, false, true);
					// try to fix rcm sequence

					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo,
						unfixedTailSeq_rcm_char,
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedTailSeq_rcm, unfixedTailLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Rcm2,
						true, false, false);
				}			
		}
	}

	void detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedHead(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, 
		PE_Read_Info* peReadInfo, int index_Nor2, int chromNameInt, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{
			int alignChromPos_start 
				= (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->alignChromPos;
			int mapStartPos_areaNO 
				= indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, alignChromPos_start);
			FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_2.find(mapStartPos_areaNO);			
			if(fusionEndSetIter != candidateFusionEndSet_2.end())
			{
				int unfixedHeadLength= (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->unfixedHeadLength();
				string unfixedHeadSeq_nor =  ((peReadInfo->readInfo_pe2).readSeq).substr(0, unfixedHeadLength);
				string unfixedHeadSeq_rcm = this->convertStringToReverseComplement(unfixedHeadSeq_nor);
				for(int tmp = 0; tmp < (fusionEndSetIter->second).size(); tmp ++)
				{
					int targetRegion_areaNO = (fusionEndSetIter->second)[tmp].first;
					int targetRegion_endPos = ((fusionEndSetIter->second)[tmp].second).first;
					int mappedRegion_startPos = ((fusionEndSetIter->second)[tmp].second).second;

					
					int targetRegion_2ndLevelIndexNO = targetRegion_areaNO - 1;
					int otherFusionSJendChrNameInt = indexInfo->getChrNameIntFromSecondLevelIndexNO(targetRegion_2ndLevelIndexNO+1);
					int otherEndPosInCandidateFusionSet =  targetRegion_endPos;
					int otherEndMapIntervalStart = otherEndPosInCandidateFusionSet - mapIntervalHalfSize;
					int otherEndMapIntervalEnd = otherEndPosInCandidateFusionSet + mapIntervalHalfSize;
					if(otherEndMapIntervalStart <= 0)
					{
						otherEndMapIntervalStart = 1;
					}
					if(otherEndMapIntervalEnd >= ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedHeadLength - 1))
					{
						otherEndMapIntervalEnd = ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedHeadLength - 1);
					}

					string otherFusionSJendChrNameStr = (indexInfo->chrNameStr)[otherFusionSJendChrNameInt];
					int otherEndMapPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
						otherFusionSJendChrNameInt, (targetRegion_2ndLevelIndexNO+1), 1);

					char* unfixedHeadSeq_nor_char = const_cast<char*>(unfixedHeadSeq_nor.c_str());
					char* unfixedHeadSeq_rcm_char = const_cast<char*>(unfixedHeadSeq_rcm.c_str()); 
					// try to fix nor sequence
					//cout << "try to detect Fusion exact site ..." << endl;
					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo, 
						unfixedHeadSeq_nor_char, 
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_nor, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Nor2, false, true, true);

					// try to fix rcm sequence
					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo, unfixedHeadSeq_rcm_char,
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_rcm, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Nor2, false, true, false);
				}
			}
	}

	void detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedTail(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, 
		PE_Read_Info* peReadInfo, int index_Rcm1, int chromNameInt, Index_Info* indexInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild )
	{
		int alignChromPos_end = (peAlignInfo->rcmAlignmentInfo_PE_1)[index_Rcm1]->endMatchedPosInChr;
		int mapEndPos_areaNO = indexInfo->getSecondLevelIndexFromChrAndPos(chromNameInt, alignChromPos_end);

		FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_1.find(mapEndPos_areaNO);	
		if(fusionEndSetIter != candidateFusionEndSet_1.end())
		{
			int unfixedTailLength = (peAlignInfo->rcmAlignmentInfo_PE_1)[index_Rcm1]->unfixedTailLength();
			int readLength_1 = (peReadInfo->readInfo_pe1).rcmReadSeq.length();
			string unfixedTailSeq_nor =  ((peReadInfo->readInfo_pe1).rcmReadSeq).substr(readLength_1 - unfixedTailLength, unfixedTailLength);
			string unfixedTailSeq_rcm = this->convertStringToReverseComplement(unfixedTailSeq_nor);			
			
				for(int tmp = 0; tmp < (fusionEndSetIter->second).size(); tmp ++)
				{
					int targetRegion_areaNO = (fusionEndSetIter->second)[tmp].first;
					int targetRegion_startPos = ((fusionEndSetIter->second)[tmp].second).second;
					int mappedRegion_endPos = ((fusionEndSetIter->second)[tmp].second).first;

					
					int targetRegion_2ndLevelIndexNO = targetRegion_areaNO - 1;
					int otherFusionSJendChrNameInt = indexInfo->getChrNameIntFromSecondLevelIndexNO(targetRegion_2ndLevelIndexNO+1);
					int otherEndPosInCandidateFusionSet =  targetRegion_startPos;
					int otherEndMapIntervalStart = otherEndPosInCandidateFusionSet - mapIntervalHalfSize;
					int otherEndMapIntervalEnd = otherEndPosInCandidateFusionSet + mapIntervalHalfSize;
					if(otherEndMapIntervalStart <= 0)
					{
						otherEndMapIntervalStart = 1;
					}
					if(otherEndMapIntervalEnd >= ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedTailLength - 1))
					{
						otherEndMapIntervalEnd = ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - unfixedTailLength - 1);
					}

					string otherFusionSJendChrNameStr = (indexInfo->chrNameStr)[otherFusionSJendChrNameInt];
					int otherEndMapPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
						otherFusionSJendChrNameInt, (targetRegion_2ndLevelIndexNO+1), 1);

					char* unfixedTailSeq_nor_char = const_cast<char*>(unfixedTailSeq_nor.c_str());
					char* unfixedTailSeq_rcm_char = const_cast<char*>(unfixedTailSeq_rcm.c_str()); 
					// try to fix nor sequence

					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo,
						unfixedTailSeq_nor_char, 
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedTailSeq_nor, unfixedTailLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Rcm1,
						false, false, true);
					// try to fix rcm sequence

					this->detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
						peFusionInfo, peAlignInfo,
						unfixedTailSeq_rcm_char,
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedTailSeq_rcm, unfixedTailLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						index_Rcm1,
						false, false, false);
				}			
		}
	}	

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2OrNor2Rcm1_unfixedHeadOrTail_NorOrRcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		char* unfixedHeadOrTailSeq_NorOrRcm_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedHeadOrTailSeq_NorOrRcm, int unfixedHeadOrTailLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int index_Nor1OrNor2OrRcm1OrRcm2, 
		bool Nor1Rcm2OrNor2Rcm1Bool, bool unfixedHeadOrTailBool, bool targetSeq_NorOrRcm)
	{
		// try to fix nor sequence
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
		bool unfixedHeadOrTailMappedBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									unfixedHeadOrTailSeq_NorOrRcm_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedHeadOrTailLength, indexInfo);	
		if(!unfixedHeadOrTailMappedBool)
		{
			delete(seg2ndOriInfo);
			return;
			//continue;
		}

		Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg(segInfo);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum = pathInfo->pathValidNumInt();
		if(pathValidNum > 10)
		{
			pathInfo->memoryFree();
			delete(pathInfo);
			delete(segInfo);
			delete(seg2ndOriInfo);
			return;// pathInfo_nor;
		}
	
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath(pathInfo, segInfo, 
			indexInfo, unfixedHeadOrTailSeq_NorOrRcm, unfixedHeadOrTailLength);

		this->pushBackPathInfo_PairedAlignments(
			Nor1Rcm2OrNor2Rcm1Bool, unfixedHeadOrTailBool, targetSeq_NorOrRcm,
			peFusionInfo, peAlignInfo, pathInfo, indexInfo, index_Nor1OrNor2OrRcm1OrRcm2);

		delete(gapInfo);
		pathInfo->memoryFree();
		delete(pathInfo);
		delete(segInfo);
		delete(seg2ndOriInfo);
		return;// pathInfo_nor;
	}

	void pushBackPathInfo_PairedAlignments(
		bool Nor1Rcm2OrNor2Rcm1Bool,
		bool unfixedHeadOrTailBool, bool targetSeq_NorOrRcmBool,
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, Path_Info* newPathInfo,
		Index_Info* indexInfo, int index_Nor1OrNor2OrRcm1OrRcm2)
	{
		string strandStr;
		if( (unfixedHeadOrTailBool && targetSeq_NorOrRcmBool) || ((!unfixedHeadOrTailBool) && (!targetSeq_NorOrRcmBool)) )
		{
			strandStr = "+";
		}
		else
		{
			strandStr = "-";
		}
		for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info(strandStr, indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

			cout << "......tmpPath: " << tmpPath+1 << endl << "......Nor1Rcm2OrNor2Rcm1: " << Nor1Rcm2OrNor2Rcm1Bool << endl
				<< "......unfixedHeadOrTail: " << unfixedHeadOrTailBool << endl 
				<< "......targetSeq_NorOrRcm: " << targetSeq_NorOrRcmBool << endl 
				<< "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " 
						<< (peAlignInfo->returnAlignInfoInPeAlignInfo(Nor1Rcm2OrNor2Rcm1Bool, (unfixedHeadOrTailBool), 
								index_Nor1OrNor2OrRcm1OrRcm2))->alignChromName << " -- " 
						<< (peAlignInfo->returnAlignInfoInPeAlignInfo(Nor1Rcm2OrNor2Rcm1Bool, (unfixedHeadOrTailBool), 
								index_Nor1OrNor2OrRcm1OrRcm2))->alignChromPos << " -- " 
						<< (peAlignInfo->returnAlignInfoInPeAlignInfo(Nor1Rcm2OrNor2Rcm1Bool, (unfixedHeadOrTailBool), 
								index_Nor1OrNor2OrRcm1OrRcm2))->jumpCodeVec2Str() << endl;  

			peFusionInfo->pushBackFusionTargetRegionAlignInfo(Nor1Rcm2OrNor2Rcm1Bool,
				unfixedHeadOrTailBool, targetSeq_NorOrRcmBool,
				index_Nor1OrNor2OrRcm1OrRcm2, tmpAlignmentInfo);
		}
	}


	/*int chooseBestFusionSiteInMultiAlignment(int expectedFusionSiteInCandidateSet, Path_Info* newPathInfo, Index_Info* indexInfo)
	{
		int tmpBestFusionSite = 0;
		int tmpDistanceFromExpectedSite = 1000000;
		//cout << "start to choose best fusion site ... finalPathSize: " << (newPathInfo->finalPathVec).size() << endl;
		for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
		{
			int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
			int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = 0;

			Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
				mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			int tmpEndMapPos = tmpAlignmentInfo->getEndMatchedPosInChr();
			if(abs(tmpEndMapPos - expectedFusionSiteInCandidateSet) < tmpDistanceFromExpectedSite)
			{
				tmpDistanceFromExpectedSite = abs(tmpEndMapPos - expectedFusionSiteInCandidateSet);
				tmpBestFusionSite = tmpEndMapPos;
			}
		}	
		return tmpBestFusionSite;
	}*/

	void detectExactFusionSite_PairedAlignments_Nor2Rcm1(//int index_Nor2, int index_Rcm1, 
		//PE_Read_Alignment_Info* peAlignInfo, PE_Read_Info* peReadInfo, Index_Info* indexInfo
		)
	{

	}	


};