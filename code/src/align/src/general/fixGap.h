// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE

bool fixHead(Splice_Info* headCigarInfo, unsigned int* segmentLocInRead, char* read, char* chrom, string strand, 
	unsigned int mainFirstFragmentMapPos, int valSegStartNo, unsigned int* mapPos, const string& readString, int readLength,
	Index_Info* indexInfo)
{
	#ifdef DEBUG
	cout << "debug -- start fixHead funtion ... " << endl;
	#endif

	if(mainFirstFragmentMapPos > (indexInfo->indexSize))
	{
		////debugln("mainFirstFragmentMapPos < 0 or > (indexInfo->indexSize)");
		return false;
	}
	bool head_fixed = false;
	int headRelation = 0;
	int headLength = segmentLocInRead[valSegStartNo-1]-1;
	int extendLength = extendBack(mainFirstFragmentMapPos, segmentLocInRead[valSegStartNo-1], read, chrom, indexInfo);
	int anchorLength = ANCHOR_LENGTH;
	int headLengthAfterExtend = headLength - extendLength;
	#ifdef DEBUG
	cout << "headLengthAfterExtend = " << headLengthAfterExtend << endl;
	#endif
	bool append_fixed = false;
	string pending_seq = readString.substr(0, headLengthAfterExtend-1);

	string chrom_seq = chromString.substr(mainFirstFragmentMapPos-1, headLengthAfterExtend-1);

	size_t max_append_mismatch = (headLengthAfterExtend-1)/10;
	size_t mismatch_bits = 0;

	append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);//append first


	if(append_fixed)
	{
		//cout << "appendFixed" << endl;
		headRelation = FIX_MATCH;
		head_fixed = fixDoubleAnchorHead(headRelation, headCigarInfo, mainFirstFragmentMapPos, 
			mainFirstFragmentMapPos, headLengthAfterExtend+1, anchorLength, read, chrom, extendLength, readString, indexInfo);
		*mapPos = mainFirstFragmentMapPos;
		return head_fixed;
	}
	else 
	{

			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////   2.  directly do soft-clipping  /////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
		headRelation = FIX_SOFTCLIPPING;
		//headSoftClipping = true;
		head_fixed = fixDoubleAnchorHead(headRelation, headCigarInfo, mainFirstFragmentMapPos, 
			mainFirstFragmentMapPos, headLengthAfterExtend+1, anchorLength, read, chrom, extendLength, readString, indexInfo);
		*mapPos = mainFirstFragmentMapPos + headLengthAfterExtend;
		return head_fixed;
	}
}


bool fixTail(Splice_Info* tailCigarInfo, unsigned int* segmentLocInRead, char* read, char* chrom, string strand, 
	unsigned int mainLastFragmentMapPos, int valSegEndNo, const string& readString, int readLength, Index_Info* indexInfo)
{
	#ifdef DEBUG
	cout << "start fixTail function ... "<< endl;
	#endif
	if (mainLastFragmentMapPos > (indexInfo->indexSize))
	{
		////debugln("mainLastFragmentMapPos > (indexInfo->indexSize)");
		return false;
	}
	bool tail_fixed = false;
	int tailRelation = 0;
	int tailStartLocInRead = segmentLocInRead[valSegEndNo] - 1;
	//int readLength = READ_LENGTH;
	int tailLength = readLength - tailStartLocInRead + 1;
	int anchorLength = ANCHOR_LENGTH;
	string pending_seq = readString.substr(tailStartLocInRead, tailLength-1);
	string chrom_seq = chromString.substr
	(mainLastFragmentMapPos + tailStartLocInRead-1 , tailLength-1);
	size_t max_append_mismatch = (tailLength-1)/10;
	size_t mismatch_bits = 0;
	bool append_fixed = score_string(
		pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);

	if(append_fixed)
	{
		tailRelation = FIX_MATCH;
		tail_fixed = fixDoubleAnchorTail(tailRelation, tailCigarInfo, mainLastFragmentMapPos,
			mainLastFragmentMapPos, tailStartLocInRead-1, anchorLength, read, chrom, readString, readLength, indexInfo);
		return tail_fixed;
	}
	else
	{
		///////////////////////////////fix short tail with spliceJunctionHash////////////////////////////////
		tailRelation = FIX_SOFTCLIPPING;
		tail_fixed = fixDoubleAnchorTail(tailRelation, tailCigarInfo, mainLastFragmentMapPos,
			mainLastFragmentMapPos, tailStartLocInRead-1, anchorLength, read, chrom, readString, readLength, indexInfo);			
			
		return tail_fixed;					
	}
}

bool fixGap(unsigned int segmentNum, unsigned int* segmentLength, unsigned int* segmentLocInRead, 
	unsigned int* segmentAlignNum, 
	unsigned int mapLabel, unsigned int* segMapRangeStart, 
	unsigned int* segMapRangeEnd, unsigned int* segMapLoc, 
	char* read, char* chrom, 
	const string& readSeq,
	const string& readString, const string& alignDirection, 
	const string& readNameString, int readLength,
	vector<Alignment_Info*>& alignmentInfo,
	Index_Info* indexInfo
	)
{
	bool fix_gap=false;
	//int readLength = READ_LENGTH;
	unsigned int valSegStartNo;
	unsigned int valSegEndNo;
	unsigned int possibleMapCaseMax;
	
	
	getValSegNoShortFragIncluded(mapLabel, segMapRangeStart, 
		segMapRangeEnd, segMapLoc, &valSegStartNo, 
		&valSegEndNo, &possibleMapCaseMax);

	if(possibleMapCaseMax > POSSIBLE_MAP_CASES_MAX)
	{	
		//(*alignmentCaseNum) = 0;
		return false;
	}
	unsigned int mapCases[MAPCASES_MAX][2] = {0}; // first segment's map position and last segment's map position

	//fix all double anchor problems to get all map cases (except heads and tails) 
	unsigned int mapCasesNum = 0;
	unsigned int mapCasesNumMAX = 100;
	vector<Splice_Info*> allMapCasesCigarInfo; //(mapCasesNumMAX); // //debug: can be changed to adding elements dynamically

	for (int frag1 = 1; frag1 <= possibleMapCaseMax; frag1++)
	{	
		Jump_Code firstFragCigar(getFragmentLength(frag1, segmentNum, segmentLocInRead, segMapRangeStart, segMapRangeEnd, readLength), "M");

		Splice_Info* newSpliceInfoElement = new Splice_Info;
		newSpliceInfoElement->jump_code.push_back(firstFragCigar);

		if(*(segMapRangeEnd + (frag1-1)) == valSegEndNo)   // this fragment has covered the whole read except head and tail
		{
			mapCases[mapCasesNum][0] = segMapLoc[frag1-1];
			mapCases[mapCasesNum][1] = segMapLoc[frag1-1];

			allMapCasesCigarInfo.push_back(newSpliceInfoElement);
			
			mapCasesNum ++;

			continue;
		}
		////debugln("start to search for second fragment!!");
		int fragEndSegNo1 = *(segMapRangeEnd + (frag1-1));
		int fragMapLoc1 = segMapLoc[frag1-1]; 
		int relation = 0;
		for ( int frag2 = possibleMapCaseMax+1; frag2 <= mapLabel; frag2++)
		{
			Splice_Info* secondSpliceInfoElement = new Splice_Info;

			relation = checkRelation(*(segMapRangeStart + (frag2-1)), fragEndSegNo1, fragMapLoc1, segMapLoc[frag2-1]);		

			secondSpliceInfoElement->cleanAndCopy(*newSpliceInfoElement);

			bool fix_double_anchor = fixDoubleAnchor(relation, secondSpliceInfoElement, frag1, frag2, segMapRangeStart, segMapRangeEnd, 
				segmentLocInRead, segMapLoc, read, chrom, segmentNum, readString, readLength//, strand
				, indexInfo);
			
			if (!fix_double_anchor)
			{		
				continue;
			}
			else
			{				
				if(*(segMapRangeEnd + (frag2-1)) == valSegEndNo)
				{			
					mapCases[mapCasesNum][0] = segMapLoc[frag1-1];
					mapCases[mapCasesNum][1] = segMapLoc[frag2-1];
					allMapCasesCigarInfo.push_back(secondSpliceInfoElement);
					mapCasesNum ++;
					continue;
				}
				
				for (int frag3 = frag2 + 1; frag3 <= mapLabel; frag3++)
				{
					Splice_Info* thirdSpliceInfoElement = new Splice_Info;

					int relation2 = checkRelation(*(segMapRangeStart + (frag3-1)), *(segMapRangeEnd + (frag2-1)), 
						segMapLoc[frag2-1], segMapLoc[frag3-1]);

					thirdSpliceInfoElement->cleanAndCopy(*secondSpliceInfoElement);
					
					bool fix_double_anchor2 = fixDoubleAnchor(relation2, thirdSpliceInfoElement, frag2, frag3, segMapRangeStart, segMapRangeEnd,
					segmentLocInRead, segMapLoc, read, chrom, segmentNum, readString, readLength//, strand
					, indexInfo);

					if(!fix_double_anchor2)
					{	

						continue;
					}
					else
					{
					
						if(*(segMapRangeEnd + (frag3-1)) == valSegEndNo)
						{
							mapCases[mapCasesNum][0] = segMapLoc[frag1-1];
							mapCases[mapCasesNum][1] = segMapLoc[frag3-1];
							allMapCasesCigarInfo.push_back(thirdSpliceInfoElement);
							mapCasesNum ++;

							continue;
						}
					}
				}
			}
		}
	}

	if(allMapCasesCigarInfo.size() == 0)
	{
		fix_gap = true;
		//(*alignmentCaseNum) = 0;
		return fix_gap;
	}	
	///////After double anchor, we get N = mapCasesNum  map cases. For every mapcases we shound conduct fixHead, and fixTail../////////  
	// fix all map cases' heads

	vector<Splice_Info*> headSpliceInfo;
	vector<Splice_Info*> tailSpliceInfo;
	vector<Splice_Info*> midSpliceInfo;
	vector<Splice_Info*> wholeSpliceInfo;
	vector<unsigned int> mapPos;
	int finalMapCases = 0;
	bool headNeedFix = (valSegStartNo > 1);
	bool tailNeedFix = (valSegEndNo < segmentNum);


	if(headNeedFix && tailNeedFix)
	{	
		#ifdef DEBUG
		cout << "debug -- both head and tail need to be fixed" << endl; 
		#endif
		for(int tmpMapCasesNum = 0; tmpMapCasesNum < allMapCasesCigarInfo.size(); tmpMapCasesNum++)
		{
			Splice_Info* headCigarInfo = new Splice_Info;
			Splice_Info* tailCigarInfo = new Splice_Info;
			Splice_Info* wholeCigarInfo = new Splice_Info;
			unsigned int tmpMapPos = 0;
			unsigned int mainFirstFragmentMapPos = mapCases[tmpMapCasesNum][0];
			unsigned int mainLastFragmentMapPos = mapCases[tmpMapCasesNum][1];

			string strand = "+";

			bool fix_head = false;
			bool fix_tail = false;
			fix_head = fixHead(headCigarInfo, segmentLocInRead, read, chrom, strand, 
				mainFirstFragmentMapPos, valSegStartNo, &tmpMapPos, readString, readLength, indexInfo);

			fix_tail = fixTail(tailCigarInfo, segmentLocInRead, read, chrom, strand, 
				mainLastFragmentMapPos, valSegEndNo, readString, readLength, indexInfo);

			if (fix_head && fix_tail)
			{
				headSpliceInfo.push_back(headCigarInfo);
				tailSpliceInfo.push_back(tailCigarInfo);
				midSpliceInfo.push_back(allMapCasesCigarInfo[tmpMapCasesNum]);
				mapPos.push_back(tmpMapPos);
				wholeCigarInfo->appendJumpCode(headCigarInfo, allMapCasesCigarInfo[tmpMapCasesNum], tailCigarInfo);
				wholeSpliceInfo.push_back(wholeCigarInfo);
				finalMapCases++;
			}
		}

		for (int tmpCase = 0; tmpCase < finalMapCases; tmpCase++)
		{

			wholeSpliceInfo[tmpCase]->getFinalJumpCode();
	
			Alignment_Info* tmpAlignmentInfo 
				= new Alignment_Info(alignDirection, mapPos[tmpCase],
					wholeSpliceInfo[tmpCase]->final_jump_code, indexInfo);

			alignmentInfo.push_back(tmpAlignmentInfo);
		}

		//(*alignmentCaseNum) = finalMapCases;
		fix_gap = true;
	}
	else if(headNeedFix)
	{
		#ifdef DEBUG
		cout << endl << "debug -- only head needs to be fixed" << endl; 
		#endif
		
		for(int tmpMapCasesNum = 0; tmpMapCasesNum < allMapCasesCigarInfo.size(); tmpMapCasesNum++)
		{

			Splice_Info* headCigarInfo = new Splice_Info;
			Splice_Info* wholeCigarInfo = new Splice_Info;
			unsigned int mainFirstFragmentMapPos = mapCases[tmpMapCasesNum][0];
			unsigned int tmpMapPos = 0;
			string strand = "+";

			bool fix_head = false;
			fix_head = fixHead(headCigarInfo, segmentLocInRead, read, chrom, strand, 
				mainFirstFragmentMapPos, valSegStartNo, &tmpMapPos, readString, readLength, indexInfo);

			if (fix_head )
			{
				headSpliceInfo.push_back(headCigarInfo);
				
				midSpliceInfo.push_back(allMapCasesCigarInfo[tmpMapCasesNum]);
				mapPos.push_back(tmpMapPos);
				wholeCigarInfo->appendJumpCode(headCigarInfo, allMapCasesCigarInfo[tmpMapCasesNum]);
				wholeSpliceInfo.push_back(wholeCigarInfo);				
				finalMapCases++;											
			}
		}

		for (int tmpCase = 0; tmpCase < finalMapCases; tmpCase++)
		{
			wholeSpliceInfo[tmpCase]->getFinalJumpCode();

			Alignment_Info* tmpAlignmentInfo 
				= new Alignment_Info(alignDirection, mapPos[tmpCase],
					wholeSpliceInfo[tmpCase]->final_jump_code, indexInfo);

			alignmentInfo.push_back(tmpAlignmentInfo);
		}
		//(*alignmentCaseNum) = finalMapCases;
		fix_gap = true;
	}
	else if(tailNeedFix)
	{
		#ifdef DEBUG
		cout << "debug -- only tail needs to be fixed" << endl; 
		#endif
		for(int tmpMapCasesNum = 0; tmpMapCasesNum < allMapCasesCigarInfo.size(); tmpMapCasesNum++)
		{
			Splice_Info* tailCigarInfo = new Splice_Info;
			Splice_Info* wholeCigarInfo = new Splice_Info;
			unsigned int mainLastFragmentMapPos = mapCases[tmpMapCasesNum][1];
			unsigned int tmpMapPos = mapCases[tmpMapCasesNum][0];
			string strand = "+";
			bool fix_tail = false;

			fix_tail = fixTail(tailCigarInfo, segmentLocInRead, read, chrom, strand,
				mainLastFragmentMapPos, valSegEndNo, readString, readLength, indexInfo);
			if (fix_tail)
			{
				midSpliceInfo.push_back(allMapCasesCigarInfo[tmpMapCasesNum]);
				tailSpliceInfo.push_back(tailCigarInfo);
				mapPos.push_back(tmpMapPos);
				wholeCigarInfo->appendJumpCode(allMapCasesCigarInfo[tmpMapCasesNum], tailCigarInfo);

				wholeSpliceInfo.push_back(wholeCigarInfo);	
				finalMapCases++;
			}
		}

		for (int tmpCase = 0; tmpCase < finalMapCases; tmpCase++)
		{
			wholeSpliceInfo[tmpCase]->getFinalJumpCode();

			Alignment_Info* tmpAlignmentInfo 
				= new Alignment_Info(alignDirection, mapPos[tmpCase],
					wholeSpliceInfo[tmpCase]->final_jump_code, indexInfo);

			alignmentInfo.push_back(tmpAlignmentInfo);
		}
		//(*alignmentCaseNum) = finalMapCases;
		fix_gap = true;
	}
	else
	{
		////debugln("no head and tail ");	
		#ifdef DEBUG
		cout << "debug -- no head and tail needs to be fixed" << endl; 
		#endif	
		for(int tmpMapCasesNum = 0; tmpMapCasesNum < allMapCasesCigarInfo.size(); tmpMapCasesNum++)
		{

			allMapCasesCigarInfo[tmpMapCasesNum]->getFinalJumpCode();

			Alignment_Info* tmpAlignmentInfo 
				= new Alignment_Info(alignDirection, mapCases[tmpMapCasesNum][0],
					allMapCasesCigarInfo[tmpMapCasesNum]->final_jump_code, indexInfo);
			alignmentInfo.push_back(tmpAlignmentInfo);
		}
		//(*alignmentCaseNum) = allMapCasesCigarInfo.size();
		fix_gap = true;
	}

	// double anchor or single anchor fails, then try to append

	return fix_gap;
}

bool detAliLoc(unsigned int norSegmentNum, unsigned int* norSegmentLength, unsigned int* norSegmentLocInRead,
	unsigned int* norSegmentAlignNum, unsigned int* norSegmentAlignLoc, 
	char* read, char* chrom, 
	int* finalMapLabel, unsigned int* segMapRangeStart, 
	unsigned int* segMapRangeEnd, unsigned int* segMapLoc,
	Index_Info* indexInfo)
{
	#ifdef DEBUG
	cout << "debug -- start detAliLoc function..." << endl;
	#endif
	bool detAliLoc = false;
	unsigned int segmentNum_max = SEGMENTNUM;
	unsigned int segmentAliNum_max = CANDALILOC;
	
	unsigned int norSegMapRangeStart[SEGMENTNUM * CANDALILOC] = {0};
	unsigned int norSegMapRangeEnd[SEGMENTNUM * CANDALILOC] = {0};
	
	unsigned int norSegMapLoc[SEGMENTNUM * CANDALILOC] = {0};

	unsigned int norValSegRange[2] = {50, 0};

	map <unsigned int, unsigned int> candAliLocMap;
	map <unsigned int, unsigned int> ::iterator candAliLoc_it;
	unsigned int mapLabel = 0; // normal
	unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal

	#ifdef DEBUG
	cout << "debug -- start detAliLoc for alignment..." << endl;
	#endif
	for (unsigned int tmpSegNum = 1; tmpSegNum <= norSegmentNum; tmpSegNum++)
	{
		if ((*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
			(*(norSegmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
		{
			continue;
		}			

		for(unsigned int tmpSegAliLoc = 1; tmpSegAliLoc <= norSegmentAlignNum[tmpSegNum-1]; tmpSegAliLoc++) 
		{				
			unsigned int tmpAlignLoc = *(norSegmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
			- *(norSegmentLocInRead + tmpSegNum - 1) + 1;
			if(tmpAlignLoc > (indexInfo->indexSize))
			{
				continue;
			}
			candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

			if ((*(norSegmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
			{				
				candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
				candAliLocWeight[mapLabel] = *(norSegmentLength + tmpSegNum - 1);

				norSegMapRangeStart[mapLabel] = tmpSegNum;
				norSegMapRangeEnd[mapLabel] = tmpSegNum;
				norSegMapLoc[mapLabel] = tmpAlignLoc;
				mapLabel ++;
			}
			else if (*(norSegmentLength + tmpSegNum - 1) >= minValSegLength)// find, add corresponding weight 
			{
				candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + *(norSegmentLength + tmpSegNum - 1);
				norSegMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
			}
			else
			{
			}
		}
	}
	candAliLocMap.clear();

	if((mapLabel <= 0) || (mapLabel > (POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/)))
	{
		#ifdef DEBUG
		cout << "debug -- mapLabel too small or large : " << mapLabel << endl;
		#endif
		(*finalMapLabel) = 0;
		return false;
	}

	#ifdef DEBUG
	cout << "# of Label is : " << mapLabel << endl;
	for (unsigned int tmp_label = 0; tmp_label < mapLabel; tmp_label++)
	{
		unsigned int tmpSegMapChr, tmpSegMapChrPos;
		getChrLocation(norSegMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

		cout << "label = " << tmp_label << ", loc = " << //norSegMapLoc[tmp_label] 
		chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
			<< norSegMapRangeStart[tmp_label] << "~" << norSegMapRangeEnd[tmp_label] << endl;
	}
	#endif

	//////////////try second Level detAlignLoc

	(*finalMapLabel) = detAliLocShortSegIncluded(norSegmentNum, norSegmentLength, norSegmentLocInRead, norSegmentAlignNum, 
		norSegmentAlignLoc, mapLabel, norSegMapLoc, segMapRangeStart, segMapRangeEnd, segMapLoc, indexInfo);	

	if(((*finalMapLabel) <= 0) || ((*finalMapLabel) > POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/))
	{
		#ifdef DEBUG
		cout << "debug -- mapLabel too small or large : " << (*finalMapLabel) << endl;
		#endif
		(*finalMapLabel) = 0;
		return false;
	}	
	
	#ifdef DEBUG
	cout << "# of Label is : " << (*finalMapLabel) << endl;
	for (unsigned int tmp_label = 0; tmp_label < (*finalMapLabel); tmp_label++)
	{
		unsigned int tmpSegMapChr, tmpSegMapChrPos;
		getChrLocation(segMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

		cout << "label = " << tmp_label << ", loc = " << //segMapLoc[tmp_label] 
		chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
			<< segMapRangeStart[tmp_label] << "~" << segMapRangeEnd[tmp_label] << endl;
	}
	#endif

	detAliLoc = true;
	return detAliLoc;
}