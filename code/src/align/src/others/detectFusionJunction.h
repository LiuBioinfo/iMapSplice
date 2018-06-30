// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
typedef map< pair<int, int>, pair<int, vector< pair<int, int> > > > CandidateFusionMap;
// pair<int, int> (left-end readMapPos2ndLevelNO, right-end readMapPos2ndLevelNO), pair<int, pair<int,int> > (voting score, pair<1stEndRead_endMapPos, 2ndEndRead_startMapPos>)    
typedef map< int, vector<pair<int, pair<int, int> > > > FusionEndSet; // <int, int> (2ndLevelIndexNO_1, (2ndLevelIndexNO_2, 1stEndRead_endMapPos_max) )

typedef map< pair<int, int>, int > CandidateFusionMap2ArrayIndex;

class Fusion_Forw_Info
{
public:
	CandidateFusionMap fusionCandidate;

	CandidateFusionMap2ArrayIndex candidateFusionArrayIndexMap;

	vector< vector< pair<int, int> > > candidateFusionExactSiteVec;
	int minSupportNum;
	int mapIntervalHalfSize;// = 300000;
	//vector< pair<int, int> > candidateFusionVec;// vector<pair<int,int>(2ndLevelIndexNO_1, 2ndLevelIndexNO_2)>

	FusionEndSet candidateFusionEndSet_1;
	FusionEndSet candidateFusionEndSet_2;

	//set<int> candidateFusionEndSet_1;
	//set<int> candidateFusionEndSet_2;

	Fusion_Forw_Info()
	{
		minSupportNum = 20;
		mapIntervalHalfSize = 300000;
	}

	void detectExactFusionSite(const string& fileName, vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		vector<unsigned int*>& secondLevelLcp,
		vector<unsigned int*>& secondLevelUp,
		vector<unsigned int*>& secondLevelDown,
		vector<unsigned int*>& secondLevelNext, 
		//bool load2ndLevelIndexBool_compressedSize,
		Index_Info* indexInfo
		)
	{
		//cout << "start to detect exact fusion site" << endl << endl;
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

				//cout << "readName: " << (peReadInfo->readInfo_pe1).readName << endl;
				//cout << "readSeq: " << (peReadInfo->readInfo_pe1).readSeq << endl;

			//peAlignInfo->pairingAlignment();
			//peAlignInfo->chooseBestAlignment();		

			peAlignInfo->chooseBestAlignment_selectRandomOneIfMulti();

			bool pairExistsBool = peAlignInfo->finalPairExistsBool();

			if(pairExistsBool &&
				((peAlignInfo->finalAlignPair_Nor1Rcm2).size() + (peAlignInfo->finalAlignPair_Nor2Rcm1).size() == 1) )
			{
				//cout << "readName: " << (peReadInfo->readInfo_pe1).readName << endl;
				//cout << "readSeq: " << (peReadInfo->readInfo_pe1).readSeq << endl;

				bool Nor1Rcm2orNor2Rcm1Bool;
				int norNO, rcmNO;
				string tmpChrNameStr;
				int tmpChrMapPosIntStart;
				int tmpChrMapPosIntEnd;
				bool pairHeadUnfixedBool, pairTailUnfixedBool; 
				int pairHeadLength, pairTailLength;
				string pairHeadSeq, pairTailSeq;
				if((peAlignInfo->finalAlignPair_Nor1Rcm2).size() == 1)
				{
					Nor1Rcm2orNor2Rcm1Bool = true;
					norNO = (peAlignInfo->finalAlignPair_Nor1Rcm2)[0].first;
					rcmNO = (peAlignInfo->finalAlignPair_Nor1Rcm2)[0].second;
					
					tmpChrNameStr = (peAlignInfo->norAlignmentInfo_PE_1)[norNO]->alignChromName;
					tmpChrMapPosIntStart = (peAlignInfo->norAlignmentInfo_PE_1)[norNO]->alignChromPos;
					tmpChrMapPosIntEnd = (peAlignInfo->rcmAlignmentInfo_PE_2)[rcmNO]->endMatchedPosInChr;
					pairHeadUnfixedBool = (peAlignInfo->norAlignmentInfo_PE_1)[norNO]->unfixedHeadExistsBool();
					pairTailUnfixedBool = (peAlignInfo->rcmAlignmentInfo_PE_2)[rcmNO]->unfixedTailExistsBool();
					if(pairHeadUnfixedBool)
					{
						pairHeadLength = (peAlignInfo->norAlignmentInfo_PE_1)[norNO]->unfixedHeadLength();
						pairHeadSeq = ((peReadInfo->readInfo_pe1).readSeq).substr(0, pairHeadLength);
					}
					if(pairTailUnfixedBool)
					{
						pairTailLength = (peAlignInfo->rcmAlignmentInfo_PE_2)[rcmNO]->unfixedTailLength();
						int tmpReadLength = (peReadInfo->readInfo_pe2).readSeqLength;
						pairTailSeq = ((peReadInfo->readInfo_pe2).rcmReadSeq).substr(tmpReadLength - pairTailLength, pairTailLength);
					}
				}
				else
				{
					Nor1Rcm2orNor2Rcm1Bool = false;
					norNO = (peAlignInfo->finalAlignPair_Nor2Rcm1)[0].first;
					rcmNO = (peAlignInfo->finalAlignPair_Nor2Rcm1)[0].second;
					
					tmpChrNameStr = (peAlignInfo->norAlignmentInfo_PE_2)[norNO]->alignChromName;
					tmpChrMapPosIntStart = (peAlignInfo->norAlignmentInfo_PE_2)[norNO]->alignChromPos;
					tmpChrMapPosIntEnd = (peAlignInfo->rcmAlignmentInfo_PE_1)[rcmNO]->endMatchedPosInChr;
					pairHeadUnfixedBool = (peAlignInfo->norAlignmentInfo_PE_2)[norNO]->unfixedHeadExistsBool();
					pairTailUnfixedBool = (peAlignInfo->rcmAlignmentInfo_PE_1)[rcmNO]->unfixedTailExistsBool();
					if(pairHeadUnfixedBool)
					{
						pairHeadLength = (peAlignInfo->norAlignmentInfo_PE_2)[norNO]->unfixedHeadLength();
						pairHeadSeq = ((peReadInfo->readInfo_pe2).readSeq).substr(0, pairHeadLength);
					}
					if(pairTailUnfixedBool)
					{
						pairTailLength = (peAlignInfo->rcmAlignmentInfo_PE_1)[rcmNO]->unfixedTailLength();
						int tmpReadLength = (peReadInfo->readInfo_pe1).readSeqLength;
						pairTailSeq = ((peReadInfo->readInfo_pe1).rcmReadSeq).substr(tmpReadLength - pairTailLength, pairTailLength);
					}
				}
			
				int tmp2ndLevelIndexNO_start = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr, tmpChrMapPosIntStart);
				int tmp2ndLevelIndexNO_end = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr, tmpChrMapPosIntEnd);


				if(pairHeadUnfixedBool)
				{
					FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_2.find(tmp2ndLevelIndexNO_start);
					if(fusionEndSetIter == candidateFusionEndSet_2.end())
					{
						continue;
					}
					//cout << "tmp2ndLevelIndexNO_start: " << tmp2ndLevelIndexNO_start << endl;

					for(int tmp = 0; tmp < (fusionEndSetIter->second).size(); tmp ++)
					{
						int secondLevelIndexNO = ((fusionEndSetIter->second)[tmp]).first;
						secondLevelIndexNO --; // Xinan: to debug, '--' correct?
						
						char* pairHeadSeqChar = const_cast<char*>(pairHeadSeq.c_str());
						Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();
						if(
							((indexInfo->invalidSecondLevelIndexNOset).find(secondLevelIndexNO+1))
							!= (indexInfo->invalidSecondLevelIndexNOset).end()
							)
						{
							delete(seg2ndOriInfo);
							continue;
						}			

						bool pairHeadMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
									pairHeadSeqChar,
									secondLevelSa[secondLevelIndexNO], 
									secondLevelLcpCompress[secondLevelIndexNO],
									secondLevelChildTab[secondLevelIndexNO],
									secondLevelChrom[secondLevelIndexNO], 
									secondLevelDetChild[secondLevelIndexNO],
									pairHeadLength, indexInfo);	
						if(!pairHeadMapBool)
						{
							delete(seg2ndOriInfo);
							continue;
						}

						int otherFusionSJendChrNameInt = indexInfo->getChrNameIntFromSecondLevelIndexNO(secondLevelIndexNO+1);
						int otherEndPosInCandidateFusionSet =  (((fusionEndSetIter->second)[tmp]).second).first;
						int otherEndMapIntervalStart = otherEndPosInCandidateFusionSet - mapIntervalHalfSize;
						int otherEndMapIntervalEnd = otherEndPosInCandidateFusionSet + mapIntervalHalfSize;
						if(otherEndMapIntervalStart <= 0)
						{
							otherEndMapIntervalStart = 1;
						}
						if(otherEndMapIntervalEnd >= ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - pairHeadLength - 1))
						{
							otherEndMapIntervalEnd = ((indexInfo->chromLength)[otherFusionSJendChrNameInt] - pairHeadLength - 1);
						}

						string otherFusionSJendChrNameStr = (indexInfo->chrNameStr)[otherFusionSJendChrNameInt];
						int otherEndMapPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
							otherFusionSJendChrNameInt, (secondLevelIndexNO+1), 1);

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
							//delete(unmapEndInfo);								
							continue;
						}
	
						Gap_Info* gapInfo = new Gap_Info();
						gapInfo->fixGapInPath(pathInfo, segInfo, 
							indexInfo, pairHeadSeq, pairHeadLength);
							
						//// to delete, just to output some intermedate results here for now ....
						int tmpBestFusionSite = this->chooseBestFusionSiteInMultiAlignment(otherEndPosInCandidateFusionSet, pathInfo, indexInfo);

						if((pathInfo->finalPathVec).size() > 0)
						{
							cout << "1stIndexNO: " << secondLevelIndexNO + 1 << " 2ndIndexNO: " << tmp2ndLevelIndexNO_start << endl;
							cout << "possibleFusionSite: " << otherFusionSJendChrNameStr << " " << tmpBestFusionSite 
								<< " -- " << tmpChrNameStr << " " << tmpChrMapPosIntStart << endl << endl;
						}
						delete(gapInfo);
						pathInfo->memoryFree();
						delete(pathInfo);
						delete(segInfo);
						delete(seg2ndOriInfo);
						//delete(unmapEndInfo);
					}
				}

				if(pairTailUnfixedBool)
				{
					FusionEndSet::iterator fusionEndSetIter = candidateFusionEndSet_1.find(tmp2ndLevelIndexNO_end);
					if(fusionEndSetIter == candidateFusionEndSet_1.end())
					{
						continue;
					}

					//cout << "tmp2ndLevelIndexNO_end: " << tmp2ndLevelIndexNO_end << endl;
				}

			}
			else
			{
				continue;
			}

				
		}
		inputRecord_ifs.close();
	}

	int chooseBestFusionSiteInMultiAlignment(int expectedFusionSiteInCandidateSet, Path_Info* newPathInfo, Index_Info* indexInfo)
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
				if(((peAlignInfo->norAlignmentInfo_PE_1).size() == 1) 
					&& ((peAlignInfo->rcmAlignmentInfo_PE_2).size() == 1)
					&& ((peAlignInfo->rcmAlignmentInfo_PE_1).size() + (peAlignInfo->norAlignmentInfo_PE_2).size() == 0))
		
				{
					string tmpChrNameStr_1 = (peAlignInfo->norAlignmentInfo_PE_1)[0]->alignChromName;
					string tmpChrNameStr_2 = (peAlignInfo->rcmAlignmentInfo_PE_2)[0]->alignChromName;
					int tmpChrPos_1 = (peAlignInfo->norAlignmentInfo_PE_1)[0]->alignChromPos;
					int tmpChrPos_2 = (peAlignInfo->rcmAlignmentInfo_PE_2)[0]->alignChromPos;

					int tmpChrPos_1_end = (peAlignInfo->norAlignmentInfo_PE_1)[0]->getEndMatchedPosInChr();

					int tmp2ndLevelIndexNO_1 = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr_1, tmpChrPos_1_end);
					int tmp2ndLevelIndexNO_2 = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr_2, tmpChrPos_2);

					if(tmp2ndLevelIndexNO_1 != tmp2ndLevelIndexNO_2)
						this->insertCandidateFusionMap(tmp2ndLevelIndexNO_1, tmp2ndLevelIndexNO_2, tmpChrPos_1_end, tmpChrPos_2);
					//cout << "read to be used for detecting fusion junction: " << readNameStr_1 << " readNO: " << readNO << endl;
					//cout << "tmp2ndLevelIndexNO_1: " << tmp2ndLevelIndexNO_1 << endl << "tmp2ndLevelIndexNO_2: " 
					//	<< tmp2ndLevelIndexNO_2 << endl;
				}
				else if(((peAlignInfo->norAlignmentInfo_PE_2).size() == 1) 
					&& ((peAlignInfo->rcmAlignmentInfo_PE_1).size() == 1)
					&& ((peAlignInfo->rcmAlignmentInfo_PE_2).size() + (peAlignInfo->norAlignmentInfo_PE_1).size() == 0))
				{
					string tmpChrNameStr_1 = (peAlignInfo->norAlignmentInfo_PE_2)[0]->alignChromName;
					string tmpChrNameStr_2 = (peAlignInfo->rcmAlignmentInfo_PE_1)[0]->alignChromName;
					int tmpChrPos_1 = (peAlignInfo->norAlignmentInfo_PE_2)[0]->alignChromPos;
					int tmpChrPos_2 = (peAlignInfo->rcmAlignmentInfo_PE_1)[0]->alignChromPos;

					int tmpChrPos_1_end = (peAlignInfo->norAlignmentInfo_PE_2)[0]->getEndMatchedPosInChr();

					int tmp2ndLevelIndexNO_1 = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr_1, tmpChrPos_1_end);
					int tmp2ndLevelIndexNO_2 = indexInfo->getSecondLevelIndexFromChrStrAndPos(tmpChrNameStr_2, tmpChrPos_2);

					if(tmp2ndLevelIndexNO_1 != tmp2ndLevelIndexNO_2) 
						this->insertCandidateFusionMap(tmp2ndLevelIndexNO_1, tmp2ndLevelIndexNO_2, tmpChrPos_1_end, tmpChrPos_2);

					//cout << "read to be used for detecting fusion junction: " << readNameStr_1 << " readNO: " << readNO << endl;
					//cout << "tmp2ndLevelIndexNO_1: " << tmp2ndLevelIndexNO_1 << endl << "tmp2ndLevelIndexNO_2: " 
					//	<< tmp2ndLevelIndexNO_2 << endl;					
				}
				else
				{
					//cout << "unpaired reads but can not be used to detect fusion: " << readNameStr_1 << endl;
				}

			}
			else
			{
				//cout << " pairedReads: " << readNameStr_1 << endl;
			}
				
		}
		inputRecord_ifs.close();
	}

	void insertCandidateFusionMap(int norEnd2ndLevelIndexNO, int rcmEnd2ndLevelIndexNO, int tmpChrPos_1_end, int tmpChrPos_2)
	{
		CandidateFusionMap::iterator fusionMapIter;
		fusionMapIter = fusionCandidate.find(pair<int,int> (norEnd2ndLevelIndexNO, rcmEnd2ndLevelIndexNO));
		if(fusionMapIter 
			!= fusionCandidate.end())
		{
			(fusionMapIter->second).first ++;
			//((fusionMapIter->second).second).insert(readNO);
			((fusionMapIter->second).second).push_back(pair<int,int>(tmpChrPos_1_end, tmpChrPos_2));
		}
		else
		{
			//set<int> tmpSet;
			//tmpSet.insert(readNO);
			vector< pair<int,int> > tmpVec;
			tmpVec.push_back (pair<int,int>(tmpChrPos_1_end, tmpChrPos_2));

			fusionCandidate.insert(pair< pair<int, int>, pair<int, vector<pair<int,int> > > > 
				(pair<int,int>(norEnd2ndLevelIndexNO, rcmEnd2ndLevelIndexNO),  pair<int, vector<pair<int,int> > >(1, tmpVec) ));
		}
	}

	void generateCandidateFusionArrayIndexMap()
	{
		CandidateFusionMap::iterator tmpIter;
		int tmpArrayIndex = 0;
		for(tmpIter = fusionCandidate.begin(); tmpIter != fusionCandidate.end(); tmpIter ++)
		{
			int tmp2ndLevelIndexNO_1 = (tmpIter->first).first;
			int tmp2ndLevelIndexNO_2 = (tmpIter->first).second;
			int tmpScore = (tmpIter->second).first;
			if(tmpScore >= minSupportNum)
			{
				candidateFusionArrayIndexMap.insert(pair< pair<int,int>, int > (
					pair<int,int>(tmp2ndLevelIndexNO_1, tmp2ndLevelIndexNO_2), tmpArrayIndex));
				vector< pair<int,int> > newTmpVec;
				candidateFusionExactSiteVec.push_back(newTmpVec);
			}
			tmpArrayIndex ++;
		}
	}

	string getFusionCandidateStr()
	{
		string fusionInfoStr;
		CandidateFusionMap::iterator tmpIter;
		int elementNo = 0;
		for(tmpIter = fusionCandidate.begin(); tmpIter != fusionCandidate.end(); tmpIter ++)
		{
			//elementNo ++ ;
			//cout << "elementNo: " << elementNo << endl;
			//cout << "here 1"<< endl;
			int tmp2ndLevelIndexNO_1 = (tmpIter->first).first;
			int tmp2ndLevelIndexNO_2 = (tmpIter->first).second;
			int tmpScore = (tmpIter->second).first;

			if(tmpScore < minSupportNum)
				continue;
			//cout << "here 2" << endl;
			fusionInfoStr = fusionInfoStr + "1stMapPosIndexNO: " + int_to_str(tmp2ndLevelIndexNO_1) 
				+ " 2ndMapPosIndexNO: " + int_to_str(tmp2ndLevelIndexNO_2) + " vote: " + int_to_str(tmpScore) + " FusionSJposPair: ";
			//cout << "here 3" << endl;
			/*set<int>::iterator setIter;
			for(setIter = ((tmpIter->second).second).begin(); setIter != ((tmpIter->second).second).end(); setIter ++)
			{
				fusionInfoStr = fusionInfoStr + int_to_str((*setIter)) + ",";
			}*/
			//cout << "here 4" << endl;
			for(int tmp = 0; tmp < ((tmpIter->second).second).size(); tmp ++)
			{
				int tmpMapPos_1 = ((tmpIter->second).second)[tmp].first;
				int tmpMapPos_2 = ((tmpIter->second).second)[tmp].second;

				fusionInfoStr = fusionInfoStr + int_to_str(tmpMapPos_1) + "--" + int_to_str(tmpMapPos_2) + ",";
			}						

			fusionInfoStr += "\n";
		}	
		return fusionInfoStr;
	}


	void generateCandidateFusionSet()
	{
		CandidateFusionMap::iterator tmpIter;
		FusionEndSet::iterator tmpFusionEndSetIter_1, tmpFusionEndSetIter_2;
		for(tmpIter = fusionCandidate.begin(); tmpIter != fusionCandidate.end(); tmpIter ++)
		{
			int tmp2ndLevelIndexNO_1 = (tmpIter->first).first;
			int tmp2ndLevelIndexNO_2 = (tmpIter->first).second;
			int tmpScore = (tmpIter->second).first;

			if(tmpScore >= minSupportNum)
			{
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

				//candidateFusionEndSet_1.insert(tmp2ndLevelIndexNO_1);
				//candidateFusionEndSet_2.insert(tmp2ndLevelIndexNO_2);
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
			fusionInfoStr = fusionInfoStr + "tmp2ndLevelIndexNO_1: " + int_to_str(tmp2ndLevelIndexNO_1) + "\n";
			for(int tmp = 0; tmp < (tmpSetIter->second).size(); tmp ++)
			{
				fusionInfoStr = fusionInfoStr + "...tmp2ndLevelIndexNO_2: " +  int_to_str((tmpSetIter->second)[tmp].first) 
					+ " tmpEndPos: " + int_to_str(((tmpSetIter->second)[tmp].second).first) + "--" 
					+  int_to_str(((tmpSetIter->second)[tmp].second).second) + "\n"; 
			}
		}
		fusionInfoStr = fusionInfoStr + "candidateFusionEndSet_2 info: \n"; 
		for(tmpSetIter = candidateFusionEndSet_2.begin(); tmpSetIter != candidateFusionEndSet_2.end(); tmpSetIter ++)
		{
			int tmp2ndLevelIndexNO_1 = tmpSetIter->first;
			fusionInfoStr = fusionInfoStr + "tmp2ndLevelIndexNO_2: " + int_to_str(tmp2ndLevelIndexNO_1) + "\n";
			for(int tmp = 0; tmp < (tmpSetIter->second).size(); tmp ++)
			{
				fusionInfoStr = fusionInfoStr + "...tmp2ndLevelIndexNO_1: " +  int_to_str((tmpSetIter->second)[tmp].first) 
					+ " tmpEndPos: " + int_to_str(((tmpSetIter->second)[tmp].second).first) + "--" 
					+  int_to_str(((tmpSetIter->second)[tmp].second).second) + "\n"; 
			}
		}
		return fusionInfoStr;
	}




	/*void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedHead_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		char* unfixedHeadSeq_nor_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedHeadSeq_nor, int unfixedHeadLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int otherEndPosInCandidateFusionSet,
		int alignChromPos_start,
		int mapStartPos_areaNO,
		int chromNameInt_mappedPart, int index_Nor1)
	{
		// try to fix nor sequence
		Seg2ndOri_Info* seg2ndOriInfo_nor = new Seg2ndOri_Info();
		bool unfixedHeadMappedBool_nor = seg2ndOriInfo_nor->mapMainSecondLevel_compressedIndex(
									unfixedHeadSeq_nor_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedHeadLength, indexInfo);	
		if(!unfixedHeadMappedBool_nor)
		{
			delete(seg2ndOriInfo_nor);
			return;
			//continue;
		}

		Seg_Info* segInfo_nor = new Seg_Info(seg2ndOriInfo_nor, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo_nor = new Path_Info();
		pathInfo_nor->getPossiPathFromSeg(segInfo_nor);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum_nor = pathInfo_nor->pathValidNumInt();
		if(pathValidNum_nor > 10)
		{
			pathInfo_nor->memoryFree();
			delete(pathInfo_nor);
			delete(segInfo_nor);
			delete(seg2ndOriInfo_nor);
			return;// pathInfo_nor;
		}
	
		Gap_Info* gapInfo_nor = new Gap_Info();
		gapInfo_nor->fixGapInPath(pathInfo_nor, segInfo_nor, 
			indexInfo, unfixedHeadSeq_nor, unfixedHeadLength);
							
		//// to delete, just to output some intermedate results here for now ....
		int tmpBestFusionSite_nor = this->chooseBestFusionSiteInMultiAlignment(
			otherEndPosInCandidateFusionSet, pathInfo_nor, indexInfo);

		if((pathInfo_nor->finalPathVec).size() > 0)
		{
			cout << endl << "****************" << endl;
			cout << "unfixed Head Nor:" << endl;
			cout << "1stIndexNO: " << targetRegion_2ndLevelIndexNO + 1 
				<< " 2ndIndexNO: " << mapStartPos_areaNO << endl;
			cout << "possibleFusionSite: " << otherFusionSJendChrNameStr << " " << tmpBestFusionSite_nor 
				<< " -- " << indexInfo->chrNameStr[chromNameInt_mappedPart] << " " << alignChromPos_start << endl << endl;
		}

		this->pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedHead_nor(
			peFusionInfo, peAlignInfo, pathInfo_nor,
			indexInfo, index_Nor1);


		delete(gapInfo_nor);
		pathInfo_nor->memoryFree();
		delete(pathInfo_nor);
		delete(segInfo_nor);
		delete(seg2ndOriInfo_nor);
		return;// pathInfo_nor;
	}

	void pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedHead_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Nor1)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Nor1_head_nor).push_back(
					pair<int, Alignment_Info*> (index_Nor1, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromName << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromPos << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->jumpCodeVec2Str() << endl;  
				//norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}				
	}

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedHead_rcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		char* unfixedHeadSeq_nor_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedHeadSeq_nor, int unfixedHeadLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int otherEndPosInCandidateFusionSet,
		int alignChromPos_start,
		int mapStartPos_areaNO,
		int chromNameInt_mappedPart, int index_Nor1)
	{
		// try to fix nor sequence
		Seg2ndOri_Info* seg2ndOriInfo_nor = new Seg2ndOri_Info();
		bool unfixedHeadMappedBool_nor = seg2ndOriInfo_nor->mapMainSecondLevel_compressedIndex(
									unfixedHeadSeq_nor_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedHeadLength, indexInfo);	
		if(!unfixedHeadMappedBool_nor)
		{
			delete(seg2ndOriInfo_nor);
			return;
			//continue;
		}

		Seg_Info* segInfo_nor = new Seg_Info(seg2ndOriInfo_nor, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo_nor = new Path_Info();
		pathInfo_nor->getPossiPathFromSeg(segInfo_nor);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum_nor = pathInfo_nor->pathValidNumInt();
		if(pathValidNum_nor > 10)
		{
			pathInfo_nor->memoryFree();
			delete(pathInfo_nor);
			delete(segInfo_nor);
			delete(seg2ndOriInfo_nor);
			return;
			//delete(unmapEndInfo);								
			//continue;
		}
	
		Gap_Info* gapInfo_nor = new Gap_Info();
		gapInfo_nor->fixGapInPath(pathInfo_nor, segInfo_nor, 
			indexInfo, unfixedHeadSeq_nor, unfixedHeadLength);
							
		//// to delete, just to output some intermedate results here for now ....
		int tmpBestFusionSite_nor = this->chooseBestFusionSiteInMultiAlignment(otherEndPosInCandidateFusionSet, pathInfo_nor, indexInfo);
		
		if((pathInfo_nor->finalPathVec).size() > 0)
		{
			cout << endl << "****************" << endl;
			cout << "unfixed Head Rcm:" << endl;
			cout << "1stIndexNO: " << targetRegion_2ndLevelIndexNO + 1 
				<< " 2ndIndexNO: " << mapStartPos_areaNO << endl;
			cout << "possibleFusionSite: " << otherFusionSJendChrNameStr << " " << tmpBestFusionSite_nor 
				<< " -- " << indexInfo->chrNameStr[chromNameInt_mappedPart] << " " << alignChromPos_start << endl << endl;
		}

		this->pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedHead_rcm(
			peFusionInfo, peAlignInfo, pathInfo_nor,
			indexInfo, index_Nor1);

		delete(gapInfo_nor);
		pathInfo_nor->memoryFree();
		delete(pathInfo_nor);
		delete(segInfo_nor);
		delete(seg2ndOriInfo_nor);
		//delete(unmapEndInfo);
	}

	void pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedHead_rcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Nor1)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Nor1_head_rcm).push_back(
					pair<int, Alignment_Info*> (index_Nor1, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromName << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->alignChromPos << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_1)[index_Nor1]->jumpCodeVec2Str() << endl;  
			}				
	}

	void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedTail_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		char* unfixedTailSeq_nor_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedTailSeq_nor, int unfixedTailLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int otherEndPosInCandidateFusionSet,
		int alignChromPos_end,
		int mapEndPos_areaNO,
		int chromNameInt_mappedPart, int index_Rcm2)
	{
		// try to fix nor sequence
		//int targetRegion_2ndLevelIndexNO = 
		Seg2ndOri_Info* seg2ndOriInfo_nor = new Seg2ndOri_Info();
		bool unfixedTailMappedBool_nor = seg2ndOriInfo_nor->mapMainSecondLevel_compressedIndex(
									unfixedTailSeq_nor_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedTailLength, indexInfo);	
		if(!unfixedTailMappedBool_nor)
		{
			delete(seg2ndOriInfo_nor);
			return;
			//continue;
		}

		Seg_Info* segInfo_nor = new Seg_Info(seg2ndOriInfo_nor, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo_nor = new Path_Info();
		pathInfo_nor->getPossiPathFromSeg(segInfo_nor);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum_nor = pathInfo_nor->pathValidNumInt();
		if(pathValidNum_nor > 10)
		{
			pathInfo_nor->memoryFree();
			delete(pathInfo_nor);
			delete(segInfo_nor);
			delete(seg2ndOriInfo_nor);
			return;
			//delete(unmapEndInfo);								
			//continue;
		}
	
		Gap_Info* gapInfo_nor = new Gap_Info();
		gapInfo_nor->fixGapInPath(pathInfo_nor, segInfo_nor, 
			indexInfo, unfixedTailSeq_nor, unfixedTailLength);
							
		//// to delete, just to output some intermedate results here for now ....
		int tmpBestFusionSite_nor = this->chooseBestFusionSiteInMultiAlignment(
			otherEndPosInCandidateFusionSet, pathInfo_nor, indexInfo);

		if((pathInfo_nor->finalPathVec).size() > 0)
		{
			cout << endl << "****************" << endl;
			cout << "unfixed Tail Nor:" << endl;
			cout << "1stIndexNO: " << mapEndPos_areaNO  
				<< " 2ndIndexNO: " << targetRegion_2ndLevelIndexNO + 1 << endl;
			cout << "possibleFusionSite: " << indexInfo->chrNameStr[chromNameInt_mappedPart] << " " << alignChromPos_end 
				<< " -- " <<  otherFusionSJendChrNameStr << " " << tmpBestFusionSite_nor << endl << endl;
		}

		this->pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedTail_nor(
			peFusionInfo, peAlignInfo, pathInfo_nor,
			indexInfo, index_Rcm2);

		delete(gapInfo_nor);
		pathInfo_nor->memoryFree();
		delete(pathInfo_nor);
		delete(segInfo_nor);
		delete(seg2ndOriInfo_nor);
		//delete(unmapEndInfo);
	}


	void pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedTail_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Rcm2)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("-", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Rcm2_tail_nor).push_back(
					pair<int, Alignment_Info*> (index_Rcm2, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromName << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromPos << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->jumpCodeVec2Str() << endl;  
			}				
	}


	void detectExactFusionSite_PairedAlignments_Nor1Rcm2_unfixedTail_rcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo, 
		char* unfixedTailSeq_nor_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedTailSeq_nor, int unfixedTailLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int otherEndPosInCandidateFusionSet,
		int alignChromPos_end,
		int mapEndPos_areaNO,
		int chromNameInt_mappedPart, int index_Rcm2)
	{
		// try to fix nor sequence
		Seg2ndOri_Info* seg2ndOriInfo_nor = new Seg2ndOri_Info();
		bool unfixedTailMappedBool_nor = seg2ndOriInfo_nor->mapMainSecondLevel_compressedIndex(
									unfixedTailSeq_nor_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedTailLength, indexInfo);	
		if(!unfixedTailMappedBool_nor)
		{
			delete(seg2ndOriInfo_nor);
			return;
			//continue;
		}

		Seg_Info* segInfo_nor = new Seg_Info(seg2ndOriInfo_nor, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo_nor = new Path_Info();
		pathInfo_nor->getPossiPathFromSeg(segInfo_nor);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum_nor = pathInfo_nor->pathValidNumInt();
		if(pathValidNum_nor > 10)
		{
			pathInfo_nor->memoryFree();
			delete(pathInfo_nor);
			delete(segInfo_nor);
			delete(seg2ndOriInfo_nor);
			return;
			//delete(unmapEndInfo);								
			//continue;
		}
	
		Gap_Info* gapInfo_nor = new Gap_Info();
		gapInfo_nor->fixGapInPath(pathInfo_nor, segInfo_nor, 
			indexInfo, unfixedTailSeq_nor, unfixedTailLength);
							
		//// to delete, just to output some intermedate results here for now ....
		int tmpBestFusionSite_nor = this->chooseBestFusionSiteInMultiAlignment(otherEndPosInCandidateFusionSet, pathInfo_nor, indexInfo);
		
		if((pathInfo_nor->finalPathVec).size() > 0)
		{
			cout << endl << "****************" << endl;
			cout << "unfixed Tail Rcm:" << endl;
			cout << "1stIndexNO: " << mapEndPos_areaNO  
				<< " 2ndIndexNO: " << targetRegion_2ndLevelIndexNO + 1 << endl;
			cout << "possibleFusionSite: " << indexInfo->chrNameStr[chromNameInt_mappedPart] << " " << alignChromPos_end 
				<< " -- " <<  otherFusionSJendChrNameStr << " " << tmpBestFusionSite_nor << endl << endl;
		}

		this->pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedTail_rcm(
			peFusionInfo, peAlignInfo, pathInfo_nor,
			indexInfo, index_Rcm2);

		delete(gapInfo_nor);
		pathInfo_nor->memoryFree();
		delete(pathInfo_nor);
		delete(segInfo_nor);
		delete(seg2ndOriInfo_nor);
		//delete(unmapEndInfo);
	}


	void pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedTail_rcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Rcm2)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Rcm2_tail_rcm).push_back(
					pair<int, Alignment_Info*> (index_Rcm2, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromName << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromPos << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->jumpCodeVec2Str() << endl;  
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
					this->detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedHead_nor(
						peFusionInfo, peAlignInfo, unfixedHeadSeq_nor_char, 
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_nor, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						otherEndPosInCandidateFusionSet, alignChromPos_start, mapStartPos_areaNO,
						chromNameInt, index_Nor2);
					// try to fix rcm sequence
					this->detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedHead_rcm(
						peFusionInfo, peAlignInfo, unfixedHeadSeq_rcm_char,
						otherEndMapIntervalStart, otherEndMapIntervalEnd, 
						otherEndMapPosStartIn2ndLevelIndex, indexInfo, otherFusionSJendChrNameStr,
						unfixedHeadSeq_rcm, unfixedHeadLength, targetRegion_2ndLevelIndexNO,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						otherEndPosInCandidateFusionSet, alignChromPos_start, mapStartPos_areaNO,
						chromNameInt, index_Nor2);
				}
			}
	}

	void detectExactFusionSite_PairedAlignments_Nor2Rcm1_unfixedHead_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		char* unfixedHeadSeq_nor_char, 
		int otherEndMapIntervalStart, int otherEndMapIntervalEnd, 
		int otherEndMapPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& otherFusionSJendChrNameStr, 
		const string& unfixedHeadSeq_nor, int unfixedHeadLength, int targetRegion_2ndLevelIndexNO,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		int otherEndPosInCandidateFusionSet,
		int alignChromPos_start,
		int mapStartPos_areaNO,
		int chromNameInt_mappedPart, int index_Nor2)
	{
		// try to fix nor sequence
		Seg2ndOri_Info* seg2ndOriInfo_nor = new Seg2ndOri_Info();
		bool unfixedHeadMappedBool_nor = seg2ndOriInfo_nor->mapMainSecondLevel_compressedIndex(
									unfixedHeadSeq_nor_char,
									secondLevelSa[targetRegion_2ndLevelIndexNO], 
									secondLevelLcpCompress[targetRegion_2ndLevelIndexNO],
									secondLevelChildTab[targetRegion_2ndLevelIndexNO],
									secondLevelChrom[targetRegion_2ndLevelIndexNO], 
									secondLevelDetChild[targetRegion_2ndLevelIndexNO],
									unfixedHeadLength, indexInfo);	
		if(!unfixedHeadMappedBool_nor)
		{
			delete(seg2ndOriInfo_nor);
			return;
			//continue;
		}

		Seg_Info* segInfo_nor = new Seg_Info(seg2ndOriInfo_nor, otherEndMapIntervalStart,
				otherEndMapIntervalEnd, otherEndMapPosStartIn2ndLevelIndex,
				indexInfo, otherFusionSJendChrNameStr);

		Path_Info* pathInfo_nor = new Path_Info();
		pathInfo_nor->getPossiPathFromSeg(segInfo_nor);

		//cout << pathInfo->possiPathStr() << endl;
		int pathValidNum_nor = pathInfo_nor->pathValidNumInt();
		if(pathValidNum_nor > 10)
		{
			pathInfo_nor->memoryFree();
			delete(pathInfo_nor);
			delete(segInfo_nor);
			delete(seg2ndOriInfo_nor);
			return;// pathInfo_nor;
		}
	
		Gap_Info* gapInfo_nor = new Gap_Info();
		gapInfo_nor->fixGapInPath(pathInfo_nor, segInfo_nor, 
			indexInfo, unfixedHeadSeq_nor, unfixedHeadLength);
							
		//// to delete, just to output some intermedate results here for now ....
		int tmpBestFusionSite_nor = this->chooseBestFusionSiteInMultiAlignment(
			otherEndPosInCandidateFusionSet, pathInfo_nor, indexInfo);

		if((pathInfo_nor->finalPathVec).size() > 0)
		{
			cout << endl << "****************" << endl;
			cout << "unfixed Head Nor:" << endl;
			cout << "1stIndexNO: " << targetRegion_2ndLevelIndexNO + 1 
				<< " 2ndIndexNO: " << mapStartPos_areaNO << endl;
			cout << "possibleFusionSite: " << otherFusionSJendChrNameStr << " " << tmpBestFusionSite_nor 
				<< " -- " << indexInfo->chrNameStr[chromNameInt_mappedPart] << " " << alignChromPos_start << endl << endl;
		}

		this->pushBackPathInfo_PairedAlignments_Nor2Rcm1_unfixedHead_nor(bool Nor1Rcm2OrNor2Rcm1,
			bool unfixedHeadOrTail, bool fusionNorOrRcm,
			peFusionInfo, peAlignInfo, pathInfo_nor,
			indexInfo, index_Nor2);


		delete(gapInfo_nor);
		pathInfo_nor->memoryFree();
		delete(pathInfo_nor);
		delete(segInfo_nor);
		delete(seg2ndOriInfo_nor);
		return;// pathInfo_nor;
	}



	void pushBackPathInfo_PairedAlignments_Nor2Rcm1_unfixedHead_nor(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Nor2)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Nor2_head_nor).push_back(
					pair<int, Alignment_Info*> (index_Nor2, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->alignChromName << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->alignChromPos << " -- " 
						<< (peAlignInfo->norAlignmentInfo_PE_2)[index_Nor2]->jumpCodeVec2Str() << endl;  
				//norAlignmentInfo_PE_1.push_back(tmpAlignmentInfo);
			}				
	}


	void pushBackPathInfo_PairedAlignments_Nor1Rcm2_unfixedTail_rcm(
		PE_Fusion_Info* peFusionInfo, PE_Read_Alignment_Info* peAlignInfo,
		Path_Info* newPathInfo, Index_Info* indexInfo, int index_Rcm2)
	{
			for(int tmpPath = 0; tmpPath < newPathInfo->finalPathVec.size(); tmpPath++)
			{
				int mapChromNameInt = (((newPathInfo->finalPathVec)[tmpPath]).first).first;
				int mapChromPosInt = (((newPathInfo->finalPathVec)[tmpPath]).first).second;
				int tmpMismatch = (newPathInfo->fixedPathMismatchVec)[tmpPath];//0;

				Alignment_Info* tmpAlignmentInfo = new Alignment_Info("+", indexInfo->chrNameStr[mapChromNameInt], 
					mapChromPosInt, (((newPathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
					tmpMismatch, indexInfo);

				(peFusionInfo->fusionAlignInfo_Rcm2_tail_rcm).push_back(
					pair<int, Alignment_Info*> (index_Rcm2, tmpAlignmentInfo));

					cout << "......tmpPath: " << tmpPath << endl << "...... " << indexInfo->chrNameStr[mapChromNameInt]
						<< " -- " << mapChromPosInt << " -- " << tmpAlignmentInfo->jumpCodeVec2Str()
						<< " ^^^^^^ " << (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromName << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->alignChromPos << " -- " 
						<< (peAlignInfo->rcmAlignmentInfo_PE_2)[index_Rcm2]->jumpCodeVec2Str() << endl;  
			}				
	}*/

	/*void detectExactFusionSite_PairedAlignments_Nor2Rcm1(
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

		if(unfixedHeadBool_Nor1)	
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

		if(unfixedTailBool_Rcm2)
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
	}*/


	/*void detectFusionSpliceSite(const string& fileName, Index_Info* indexInfo)
	{
		ifstream inputRecord_ifs(fileName.c_str());
		string line1, line2, line3, line4, line5, line6, line7,
			line8, line9, line10, line11; 
		getline(inputRecord_ifs, line1);

		while(!inputRecord_ifs.eof())
		{
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

			if((peAlignInfo->finalAlignPair_Nor1Rcm2).size() + (peAlignInfo->finalAlignPair_Nor2Rcm1).size() == 1)
			{
				if(!(peAlignInfo->allAlignmentInFinalPairCompleted()))
				{

				}	
			}
		}
	}*/
};

/*
class FusionDetection_Forw_Main
{
public:

	void detectFusion_forw_incompletePairedRead(const string& incompletePairedReadFile, Index_Info* indexInfo)
	{
		ifstream inputRecord_ifs(incompletePairedReadFile.c_str()); // unpaired alignments file
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

				
		}
		inputRecord_ifs.close();
	}
};

class FusionDetection_Forw_Info
{
public: 

	//int fusionReadType;
	bool validReadToDetectFusion;

	bool headPartToProcessForFusionDetection;
	bool TailPartToProcessForFusionDetection;


	bool Nor1Rcm2orNor2Rcm1Bool;
	bool headPartUnmapped;
	bool tailPartUnmapped;

	int unfixedHeadLength;
	int unfixedTailLength;

	int pairReadMappedPosStart;
	int pairReadMappedPosEnd;

	string chrNameStr_fusion_1stPart;
	string chrNameStr_fusion_2ndPart;

	int pairReadMappedPosStart_2ndLevelIndexNO;
	int pairReadMappedPosEnd_2ndLevelIndexNO;


	int minLength_validUnmappedSeq;

	FusionDetection_Forw_Info()
	{
		minLength_validUnmappedSeq = 50;
	}

	bool candidatePairedReadToDetectFusionBool(PE_Read_Alignment_Info* peAlignInfo)
	{
		bool validRecordToDetectFusion = false;
		if((peAlignInfo->finalAlignPair_Nor1Rcm2).size() + (peAlignInfo->finalAlignPair_Nor2Rcm1).size() == 1)
		{
			if((peAlignInfo->finalAlignPair_Nor1Rcm2).size() == 1)
			{
				Nor1Rcm2orNor2Rcm1Bool = true;
				
				headPartUnmapped = (peAlignInfo->norAlignmentInfo_PE_1)[0]->unfixedHeadExistsBool();
				tailPartUnmapped = (peAlignInfo->rcmAlignmentInfo_PE_2)[0]->unfixedTailExistsBool();
				
				unfixedHeadLength = (peAlignInfo->norAlignmentInfo_PE_1)[0]->unfixedHeadLength();
				unfixedTailLength = (peAlignInfo->rcmAlignmentInfo_PE_2)[0]->unfixedTailLength();


			}
			else
			{
				Nor1Rcm2orNor2Rcm1Bool = false;

			}
		}
		else
		{
			return false;
		}
		return validRecordToDetectFusion;
	}

}; */


/*class DetectFusionExactSite_Forw_1stPartUnfound_Info
// two candidate types to detect exact fusion site:
// 1: (Paired alignment found in the 1st part of some candidate fusion pair) Paired Alignments, fusion site in the rear end, mapping unmapped rear-end tail-seq to the 2nd part of candidate fusion pair
//seq.23/2        163     chr15   101941111       255     100M      =       101941661      0       TGAGCGGTTTCCTTTCCAAGCCTCAGTGGAGTTTGTTTTCTCCTCCAGCCCTGAGAAGGTCAAAGGCGCCTCGGCCCACACAAGGGGCCCCCTGGATGGC   *       NM:i:1  XS:A:   XF:Z:
//seq.23/1        83      chr15   101941661       255     6S22M72S  =       101941111     0GCCCCTGTCTTTCAGCAGCAGTTGGGGGCGGCCGGCCGGGGGATTCTGCTGTGTGTTCTGTGGTTAGCGTGGATAAGTGGCGCAGAGTCGCTGCTGACAG   *       NM:i:1  XS:A:   XF:Z:

// 2: (Paired alignment found in the 2st part of some candidate fusion pair) Paired Alignments, fusion site in the front end, mapping unmapped front-end head-seq to the 1st part of candidate fusion pair
//seq.23/2        163     chr15   101941111       255     80S20M      =       101941661      0       TGAGCGGTTTCCTTTCCAAGCCTCAGTGGAGTTTGTTTTCTCCTCCAGCCCTGAGAAGGTCAAAGGCGCCTCGGCCCACACAAGGGGCCCCCTGGATGGC   *       NM:i:1  XS:A:   XF:Z:
//seq.23/1        83      chr15   101941661       255     82M18S  =       101941111     0		GCCCCTGTCTTTCAGCAGCAGTTGGGGGCGGCCGGCCGGGGGATTCTGCTGTGTGTTCTGTGGTTAGCGTGGATAAGTGGCGCAGAGTCGCTGCTGACAG	GCCCCTGTCTTTCAGCAGCAGTTGGGGGCGGCCGGCCGGGGGATTCTGCTGTGTGTTCTGTGGTTAGCGTGGATAAGTGGCGCAGAGTCGCTGCTGACAG   *       NM:i:1  XS:A:   XF:Z:

{
public:
	bool validFusionReadBool;
	
	bool Pair_Nor1Rcm2_or_Nor2Rcm1_bool;

	bool unmappedSeq_End1OrEnd2;

	string otherEndChrStr;
	int otherEndMapPosInterval_start;
	int otherEndMapPosInterval_endl;
	int secondLevelIndexNO;


	int unmappedSeqStartPosInRead;
	int unmappedSeqEndPosInRead;

	vector<Alignment_Info*> otherEndAlignInfoVec;

	bool fusionSiteFound;	

	int unmappedSeqLengthMIN;


	DetectFusionExactSite_Forw_1stPartUnfound_Info()
	{
		unmappedSeqLengthMIN = 50;
	}	

	//bool checkValidFusionReadBool()
	//{}

	void generateDetectFusionExactSiteInfo(Alignment_Info* alignInfo_1, Alignment_Info* alignInfo_2)
	// (Nor_1, Rcm_2)  or  (Nor_2, Rcm_2)
	{
		if(alignInfo_1->unfixedHeadExistsBool())
			unmappedSeq_End1OrEnd2 = true;
		else
			unmappedSeq_End1OrEnd2 = false;

		if(unmappedSeq_End1OrEnd2)
		{
			otherEndChrStr = alignInfo_1->alignChromName;
			otherEndMapPosInterval_start = 
			otherEndMapPosInterval_end = 
			secondLevelIndexNO = 

			unmappedSeqStartPosInRead = 1;
			unmappedSeqEndPosInRead = alignInfo_1->unfixedHeadLength();
		}
		else
		{
			otherEndChrStr = alignInfo_1->alignChromName;
			otherEndMapPosInterval_start = 
			otherEndMapPosInterval_end = 
			secondLevelIndexNO = 

			unmappedSeqStartPosInRead = 1;
			unmappedSeqEndPosInRead = alignInfo_1->unfixedHeadLength();
		}
	}

	string validCandidateFusionRead(PE_Read_Info* peReadInfo, bool unmappedSeq_End1OrEnd2, 
		//bool unmappedSeq_NorOrRcm, 	
		int unmappedSeqStartPosInRead, int unmappedSeqEndPosInRead)
	{
		string unmappedSeq;

		unmappedSeq 
			= (peReadInfo->getReadSeq(unmappedSeq_End1OrEnd2, true)).substr(
				unmappedSeqStartPosInRead - 1, unmappedSeqEndPosInRead - unmappedSeqStartPosInRead + 1);

		return unmappedSeq;
	}

	void generateUnmappedSideAlignInfoInFusion(
		PE_Read_Info* peReadInfo, 
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		vector<unsigned int*>& secondLevelLcp,
		vector<unsigned int*>& secondLevelUp,
		vector<unsigned int*>& secondLevelDown,
		vector<unsigned int*>& secondLevelNext, 
		bool load2ndLevelIndexBool_compressedSize,
		Index_Info* indexInfo,
		bool unmappedSeq_End1OrEnd2, //bool unmappedSeq_NorOrRcm, 
		int unmappedSeqStartPosInRead, int unmappedSeqEndPosInRead,
		int secondLevelIndexNO)
	{
		string unmappedSeq = this->validCandidateFusionRead(peReadInfo, unmappedSeq_End1OrEnd2, 
			//unmappedSeq_NorOrRcm, 
			true, unmappedSeqStartPosInRead, unmappedSeqEndPosInRead);
		int unmappedSeqLength = unmappedSeq.length();

		char* unmappedSeqChar =  const_cast<char*>(unmappedSeq.c_str());
	
		Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

		if(
			((indexInfo->invalidSecondLevelIndexNOset).find(secondLevelIndexNO+1))
				!= (indexInfo->invalidSecondLevelIndexNOset).end()
			)
		{
			delete(seg2ndOriInfo);
			return;
		}

		bool unmappedSeqMapBool;

		if(!load2ndLevelIndexBool_compressedSize)
		{
			//cout << "start to mapMainSecondLevel !" << endl;
			unmappedSeqMapBool = seg2ndOriInfo->mapMainSecondLevel(
				unmappedSeqChar, 
				secondLevelSa[secondLevelIndexNO], 
				secondLevelLcp[secondLevelIndexNO], 
				secondLevelUp[secondLevelIndexNO], 
				secondLevelDown[secondLevelIndexNO], 
				secondLevelNext[secondLevelIndexNO], 
				secondLevelChrom[secondLevelIndexNO], 
				unmappedSeqLength, indexInfo->secondLevelIndexNormalSize + 1, indexInfo);
				//cout << "finish doing mapMainSecondLevel !" << endl;
				//cout << "unmapEndMapBool: " << unmapEndMapBool << endl;
		}
		else
		{
			//cout << "start to mapMainSecondLevel_compressedIndex !" << endl;
			unmappedSeqMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
				unmappedSeqChar,
				secondLevelSa[secondLevelIndexNO], 
				secondLevelLcpCompress[secondLevelIndexNO],
				secondLevelChildTab[secondLevelIndexNO],
				secondLevelChrom[secondLevelIndexNO], 
				secondLevelDetChild[secondLevelIndexNO],
				unmappedSeqLength, indexInfo);
				//cout << "finish doing mapMainSecondLevel_compressedIndex !" << endl;
				//cout << "unmapEndMapBool: " << unmapEndMapBool << endl;
		}		

		if(!unmappedSeqMapBool)
		{
			delete(seg2ndOriInfo);
			return;
		}
		
		int otherEndChrInt = indexInfo->convertStringToInt(otherEndChrStr);
		int chrPosStartIn2ndLevelIndex = indexInfo->getChrPosFromSecondLevelIndexPos(
				otherEndChrInt, secondLevelIndexNO, 1);

		Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
			otherEndMapPosInterval_start,
			otherEndMapPosInterval_endl, 
			chrPosStartIn2ndLevelIndex,
			indexInfo, otherEndChrStr);
	
		Path_Info* pathInfo = new Path_Info();
		pathInfo->getPossiPathFromSeg(segInfo);
	
		Gap_Info* gapInfo = new Gap_Info();
		gapInfo->fixGapInPath(pathInfo, segInfo, 
			indexInfo, unmappedSeq, unmappedSeqLength);

		for(int tmpPath = 0; tmpPath < pathInfo->finalPathVec.size(); tmpPath ++)
		{
			string tmpChromNameStr = otherEndChrStr;
			int tmpChromPosInt = (((pathInfo->finalPathVec)[tmpPath]).first).second;
			int tmpMismatch = (pathInfo->fixedPathMismatchVec)[tmpPath];

			string alignDirection;
			if (unmappedSeq_NorOrRcm)
			{
				alignDirection = "+";	
			}
			else
			{
				alignDirection = "-";
			}

			Alignment_Info* tmpAlignInfo = new Alignment_Info(alignDirection,
				tmpChromNameStr, tmpChromPosInt, 
				(((pathInfo->finalPathVec)[tmpPath]).second)->final_jump_code, 
				tmpMismatch, indexInfo);
			otherEndAlignInfoVec.push_back(tmpAlignInfo);
		}
					
		delete(gapInfo);
		pathInfo->memoryFree();
		delete(pathInfo);
		delete(segInfo);
		delete(seg2ndOriInfo);
	}


};*/