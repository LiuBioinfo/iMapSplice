all: alignAll annotation

OPTFLAGS = -O0 -g -rdynamic -std=c++0x
OUTPUTDIR = ./

CFLAGS += $(OPTFLAGS) -fopenmp

MapSplice:
	g++  $(CFLAGS) -o $(OUTPUTDIR)MapSplice src/alignAll.cpp

MapSplice_runtime:
	g++  $(CFLAGS) -D CAL_TIME -o $(OUTPUTDIR)MapSplice_runtime src/alignAll.cpp

clean_MapSplice_runtime:
	rm MapSplice_runtime

iMapSplice:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -o $(OUTPUTDIR)iMapSplice src/alignAll.cpp

iMapSplice_runtime:
	g++ $(CFLAGS) -D CAL_TIME -D PERSONALIZED_CHR_SEQ -o $(OUTPUTDIR)iMapSplice_runtime src/alignAll.cpp

iMapSplice_twoHap:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -D TWO_HAP -o $(OUTPUTDIR)iMapSplice src/alignAll.cpp

clean_iMapSplice_runtime:
	rm iMapSplice_runtime

iMapSplice_varySNPmer:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -D VARY_SNP_MER -o $(OUTPUTDIR)iMapSplice_varySNPmer src/alignAll.cpp	

iMapSplice_twoHap:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -D TWO_HAP -o $(OUTPUTDIR)iMapSplice src/alignAll.cpp	

arrayList:
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll_arrayList src/alignAll_openMP.cpp

mpsI_arrayList:
	g++  $(CFLAGS) -D PERSONALIZED_CHR_SEQ -o $(OUTPUTDIR)mpsI_arrayList src/mpsI_arrayList.cpp	

#mps_arrayList:
#	g++  $(CFLAGS) -o $(OUTPUTDIR)mps_arrayList src/mpsI_arrayList.cpp	

alignAll_speedUp:
	g++ $(CFLAGS) -o $(OUTPUTDIR)alignAll_speedUp src/alignAll_speedUp.cpp

alignAll:
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll src/alignAll.cpp

alignAll_personalizedSeq_SNPseqIndex:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -o $(OUTPUTDIR)alignAll_personalizedSeq_SNPseqIndex src/alignAll.cpp 	

alignAll_outputUnpaired:
	g++  $(CFLAGS) -D OUTPUT_UNPAIRED_REGULAR -o $(OUTPUTDIR)alignAll_outputUnpaired src/alignAll.cpp

alignAll_phase1only:
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll_phase1only src/alignAll_phase1only.cpp

transcriptCount:
	g++	$(CFLAGS) -o $(OUTPUTDIR)transcriptCount src/alignAll_transcriptCount.cpp

geneCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)geneCount src/otherProjects/geneCount/geneCount.cpp

map2annotationOnly:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2annotationOnly src/otherProjects/map2annotationOnly/map2annotationOnly.cpp

compareIntraChromSJwithGTF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareIntraChromSJwithGTF src/MPS3SeqTools/compareIntraChromSJwithGTF.cpp

#CircRNA_begin

cirRNA:
	g++  $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)alignAll_cirRNA src/alignAll.cpp

circularRNA_post:
	g++ $(CFLAGS) -o $(OUTPUTDIR)circularRNA_post src/otherProjects/cirRNA/circularRNA_post.cpp

compareBackSpliceResultsWithValidatedList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSpliceResultsWithValidatedList src/otherProjects/cirRNA/utils/compareBackSpliceResultsWithValidatedList.cpp

simulateCirRNAalignment2IGV:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateCirRNAalignment2IGV src/otherProjects/cirRNA/utils/simulateCirRNAalignments2IGV.cpp

extractBackSpliceFromAllJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractBackSpliceFromAllJunc src/otherProjects/cirRNA/extractBackSpliceFromAllJunc.cpp

extractLariatFromBackSpliceSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLariatFromBackSpliceSJ src/otherProjects/cirRNA/utils/extractLariatFromBackSpliceSJ.cpp

extractLariatFromBackSpliceSJ_syntheticSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLariatFromBackSpliceSJ_syntheticSAM src/otherProjects/cirRNA/utils/extractLariatFromBackSpliceSJ_syntheticSAM.cpp

generateCirRNAtranscript:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCirRNAtranscript src/otherProjects/cirRNA/generateCirRNAtranscript_withBackSpliceAndNormalJunction.cpp

extractBackSpliceFromSTARfusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractBackSpliceFromSTARfusionJunc src/otherProjects/cirRNA/utils/extractBackSpliceFromSTARfusionJunc.cpp

reformatMPS2circRNA:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatMPS2circRNA src/otherProjects/cirRNA/utils/reformatMPS2circRNA.cpp

mergeCirRNAreadCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeCirRNAreadCount src/otherProjects/cirRNA/mergeCirRNAreadCount.cpp

#mergeBackSpliceJunctionWithReadCount:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount src/otherProjects/cirRNA/utils/mergeBackSpliceJunctionWithReadCount.cpp

mergeBackSpliceJunctionWithReadCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount src/otherProjects/cirRNA/utils/mergeBackSpliceJunctionWithReadCount_normalize.cpp

reformBackSpliceSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformBackSpliceSAM src/otherProjects/cirRNA/utils/reformBackSpliceSAM.cpp

compareBackSplice2gtfAnnotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSplice2gtfAnnotation src/otherProjects/cirRNA/compareBackSplice2gtfAnnotation.cpp

compareBackSplice2gtfAnnotation_geneName:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSplice2gtfAnnotation_geneName src/otherProjects/cirRNA/compareBackSplice2gtfAnnotation_geneName.cpp

compareCirRNA2gtfAnnotation_withReadCountFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareCirRNA2gtfAnnotation_withReadCountFile src/otherProjects/cirRNA/compareCirRNA2gtfAnnotation_withReadCountFile.cpp

filterBackSpliceJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterBackSpliceJunc src/otherProjects/cirRNA/backSpliceFiltering_anchorSeqSimilarity.cpp

fixHeadTail_cirRNA:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)fixHeadTail_cirRNA_withCirRNAtranscriptReadCountFile src/otherMainFiles/alignAll_fixHeadTail.cpp

normalize_foldChange_comparedTo1stSample:
	g++ $(CFLAGS) -o $(OUTPUTDIR)normalize_foldChange_comparedTo1stSample src/otherProjects/cirRNA/normalize_foldChangeCompared1stSample.cpp

countBackSpliceEncompassingReads:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countBackSpliceEncompassingReads src/otherProjects/cirRNA/countBackSpliceEncompassingReads.cpp

#CircRNA_end


initialMappingForFusionDetection:
	g++ $(CFLAGS) -D OUTPUTUNIQUEUNPAIREDBOTHENDSMAPPEDSAM -o $(OUTPUTDIR)alignAll_initialMappingForFusionDetection src/alignAll.cpp

detectFusion_post:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detectFusion_post src/otherProjects/MPS3_fusion/detectFusion_post.cpp

globalMapIncompletePairSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)globalMapIncompletePairSam src/otherProjects/MPS3_fusion/globalMapIncompletePairSam.cpp

globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion src/otherProjects/MPS3_fusion/globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion.cpp

remapAndCountEncompRead:
	g++ $(CFLAGS) -o $(OUTPUTDIR)remapAndCountEncompRead src/otherProjects/MPS3_fusion/remapAndCountEncompRead.cpp

generateFusionJunc_filter_compare2gtf:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJunc_filter_compare2gtf src/otherProjects/MPS3_fusion/generateFusionJunc_filter_compare2gtf.cpp

remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint:
	g++ $(CFLAGS) -o $(OUTPUTDIR)remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint src/otherProjects/MPS3_fusion/remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint.cpp

countUniqueUnpairedReadsEncompassingFusionBreakPoint:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countUniqueUnpairedReadsEncompassingFusionBreakPoint src/otherProjects/MPS3_fusion/countUniqueUnpairedReadsEncompassingFusionBreakPoint.cpp

analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair:
	g++ $(CFLAGS) -o $(OUTPUTDIR)analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair src/otherProjects/MPS3_fusion/analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair.cpp

separateFusionJuncStrandedOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateFusionJuncStrandedOrNot src/otherProjects/MPS3_fusion/utils/separateFusionJuncStrandedOrNot.cpp

filterFusionJunc_anchorSeqSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity src/otherProjects/MPS3_fusion/fusionJuncFiltering_anchorSeqSimilarity.cpp

filterFusionJunc_localIndexMap:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJunc_localIndexMap src/otherProjects/MPS3_fusion/fusionJuncFiltering_localIndexMap.cpp

filterNonStrandedFusionJunc_anchorSeqSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity src/otherProjects/MPS3_fusion/filterNonStrandedFusionJunc_anchorSeqSimilarity.cpp

filterOutBeersGroundTruthSJfromDetectedFusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterOutBeersGroundTruthSJfromDetectedFusionJunc src/otherProjects/MPS3_fusion/utils/filterOutBeersGroundTruthSJfromDetectedFusionJunc.cpp

reformatMPS2AnnotatedFusionResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatMPS2AnnotatedFusionResults src/otherProjects/MPS3_fusion/utils/reformatMPS2AnnotatedFusionResults.cpp

clusterAndGenerateAlterFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)clusterAndGenerateAlterFusion src/otherProjects/MPS3_fusion/clusterAndGenerateAlterFusion.cpp

generateFusionSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionSeq src/otherProjects/MPS3_fusion/utils/generateFusionSeq.cpp

fixHeadTail:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fixHeadTail src/otherMainFiles/alignAll_fixHeadTail.cpp

alignAll_PE:
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll src/alignAll_PE.cpp

multi:
	g++  $(CFLAGS) -D CHECK_MULTI -o $(OUTPUTDIR)alignAll_multi src/alignAll.cpp

debug:
	g++  $(CFLAGS) -D DEBUG -o $(OUTPUTDIR)alignAll_debug src/alignAll.cpp	
	
calTime:
	g++  $(CFLAGS) -D CAL_TIME -o $(OUTPUTDIR)alignAll_calTime src/alignAll.cpp	

#annotation:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)alignAll_annotation src/alignAll_transcript.cpp

#annotation_mapInfo:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)alignAll_annotation_mapInfo src/alignAll_transcript.cpp -D MAP_INFO

mapInfo:
	g++  $(CFLAGS) -D MAP_INFO -o $(OUTPUTDIR)alignAll_mapInfo src/alignAll.cpp

noFiltering:
	g++  $(CFLAGS) -D NO_FILTERING -o $(OUTPUTDIR)alignAll_noFiltering src/alignAll.cpp

alignment2Junc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)alignment2Junc src/MPS3SeqTools/alignmentToJunc.cpp

compareJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareJunc src/MPS3SeqTools/compareJunc.cpp

mergeSJfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeSJfile src/MPS3SeqTools/mergeMultiSJfiles2One.cpp

getAlignInferJuncHashFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile src/MPS3SeqTools/getAlignInferJuncHashFromJuncFile.cpp

getAlignInferJuncHashFromSAMfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile.cpp

getAlignInferJuncHashFromJuncFile_doAlterSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile_doAlterSpliceSite src/MPS3SeqTools/getAlignInferJuncHashFromJuncFile_doAlterSpliceSite.cpp

getAlignInferJuncHashFromSAMfile_doAlterSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile_doAlterSpliceSite.cpp

#getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifyJunc.cpp

getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity src/MPS3SeqTools/getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity.cpp

getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ src/MPS3SeqTools/getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ.cpp

getSmallExon:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSmallExon src/MPS3SeqTools/geSmallExonFromSAMfile.cpp

getSoftClipSeqSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSoftClipSeqSAM src/MPS3SeqTools/getSoftClipSeqSAM.cpp
nwdp:
	g++ $(CFLAGS) -o $(OUTPUTDIR)nwdp src/others/nw_DP_new.cpp

nwdp_oneEndOpen:
	g++ $(CFLAGS) -o $(OUTPUTDIR)nwdp_oneEndOpen src/others/nwdp_oneEndOpen.cpp

getShortAnchorSJalignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getShortAnchorSJalignment src/MPS3SeqTools/seperateTrueOrFalseShortAnchorSJalignment.cpp

refineSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM src/MPS3SeqTools/refineAlignment.cpp

refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM src/MPS3SeqTools/refineAlignment_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity.cpp

extendSJspliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extendSJspliceSite src/MPS3SeqTools/extendFromSpliceSiteInsteadOfJunction.cpp

splitSpliceSiteExtensionSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)splitSpliceSiteExtensionSJ src/MPS3SeqTools/spliceSiteExtension2SJ.cpp

compareSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSpliceSite src/MPS3SeqTools/compareSpliceSite.cpp

spliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)spliceSiteAnchorSimilarity src/MPS3SeqTools/compare_anchorSimilarity_extension_multiSpliceSite.cpp

detectIncompatibleSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detectIncompatibleSAM src/MPS3SeqTools/extractIncompatiblePairedAlignments.cpp

filterJunc_anchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterJunc_anchorSimilarity src/MPS3SeqTools/classifyJunc_getAlignInferJuncHashFromJuncFile_doAlterSpliceSite.cpp

classifySJ_withAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)classifySJ_withAnchorSimilarity src/MPS3SeqTools/classifyJunc_withAlterSpliceSiteAnchorSimilarity.cpp

getPositiveStrandFa_MPS3sam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getPositiveStrandFa_MPS3sam src/MPS3SeqTools/getPositiveStrandFastaFromMPS3sam.cpp

outputSJtypeNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)outputSJtypeNum src/MPS3SeqTools/countSJtypeNum.cpp

countSJcanonicalOrNotNumFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSJcanonicalOrNotNumFromJuncFile src/MPS3SeqTools/countSJcanonicalOrNotNumFromJuncFile.cpp

parseSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSAM src/MPS3SeqTools/parseSAM.cpp

refineSAM_SE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_SE src/MPS3SeqTools/refineSAM_SE.cpp

refineSAM_PE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_PE src/MPS3SeqTools/refineSAM_PE.cpp

softClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)softClip src/MPS3SeqTools/getSoftClippedStats.cpp

softClip_varReadLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)softClip_varReadLength src/MPS3SeqTools/getSoftClippedStats_varReadLength.cpp

insertColumn2eachLine:
	g++ $(CFLAGS) -o $(OUTPUTDIR)insertColumn2eachLine src/MPS3SeqTools/insertColumn2eachLine.cpp

SJsupportNumDiff:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJsupportNumDiff src/MPS3SeqTools/generateSJsupportNumDifference.cpp

SJsupportNumDiff_targetSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJsupportNumDiff src/MPS3SeqTools/generateSJsupportNumDifference_targetSJ.cpp

getSJsupportNum_targetSJ_differentAligner:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSJsupportNum_targetSJ_differentAligner src/MPS3SeqTools/compareJunctionSupportNumDifference.cpp

SJdiscoveryStats_supportNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJdiscoveryStats_supportNum src/MPS3SeqTools/generateSJdiscoveryStats_supportNum.cpp

countAlignmentType:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countAlignmentType src/MPS3SeqTools/countAlignmentTypeNum.cpp

checkSJanchorSimilarity_fromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJanchorSimilarity_fromJuncFile src/MPS3SeqTools/checkSJanchorSimilarity_fromJuncFile.cpp

getBothEndsReadWholeMatchFasta_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_groundTruth src/MPS3SeqTools/getBothEndsReadWholeMatchFasta_groundTruth.cpp

getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth src/MPS3SeqTools/getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth.cpp

getBothEndsReadWithoutSJfasta_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWithoutSJfasta_groundTruth src/MPS3SeqTools/getBothEndsReadWithoutSJfasta_groundTruth.cpp

convertBedtoolsResults2plotInR:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertBedtoolsResults2plotInR src/MPS3SeqTools/convertBedtoolsResults2plotInR.cpp

#sam2alignInferJuncHash_supportNum:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum src/MPS3SeqTools/sam2alignInferJuncHash_supportNum.cpp

#sam2alignInferJuncHash_supportNum_parallel:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_parallel.cpp

#sam2alignInferJuncHash_supportNum_anchorSize:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize.cpp

#sam2alignInferJuncHash_supportNum_anchorSize_parallel:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_parallel.cpp

#sam2alignInferJuncHash_supportNum_anchorSize_XM:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile.cpp

compareSam2alignInferJuncHash:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSam2alignInferJuncHash src/MPS3SeqTools/compareSam2alignInferJuncHash.cpp

countSplicedAlignmentNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSplicedAlignmentNum src/MPS3SeqTools/countSplicedAlignmentNum.cpp

separateForRevSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateForRevSam src/otherProjects/separateForRevSam.cpp

extractPeAlignmentWithInnerSoftClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractPeAlignmentWithInnerSoftClip src/otherProjects/refineInnerParts_PEreads/extractPeAlignmentWithInnerSoftClip.cpp

refineInnerSoftClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineInnerSoftClip src/otherProjects/refineInnerParts_PEreads/refineInnerSoftClip.cpp

exonSize:
	g++ $(CFLAGS) -o $(OUTPUTDIR)exonSize src/others/exonSize.cpp

compareFusionRegion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionRegion src/otherProjects/MPS3_fusion/others/compareFusionRegion.cpp

refineSAM_insertAnnotatedSJ2matchThroughEnd:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd src/MPS3SeqTools/refineSAM_insertAnnotatedSJ2matchThroughEnd.cpp

refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly:
	g++ $(CFLAGS) -D COMPLETEONLY_WHEN_INSERT_ANNOTATEDSJ -o $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly src/MPS3SeqTools/refineSAM_insertAnnotatedSJ2matchThroughEnd.cpp

countReadsStartEndMappingPosWithinOverlapRegionList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadsStartEndMappingPosWithinOverlapRegionList src/MPS3SeqTools/countReadsStartEndMappingPosWithinOverlapRegionList.cpp

compareSJfromMultiSamples:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSJfromMultiSamples src/MPS3SeqTools/compareSJfromMultiSamples.cpp

generateUnannotatedSJdistributionAlongRefGenome:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateUnannotatedSJdistributionAlongRefGenome src/MPS3SeqTools/generateUnannotatedSJdistributionAlongRefGenome.cpp

reformJuncFastaFromFusionCancer_oneLineSeqPerFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformJuncFastaFromFusionCancer_oneLineSeqPerFusion src/otherProjects/MPS3_fusion/simulate_FusionCancer/reformJuncFastaFromFusionCancer_oneLineSeqPerFusion.cpp

extractUniqueBreakPointFusionFromReformedJuncFastaFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractUniqueBreakPointFusionFromReformedJuncFastaFile src/otherProjects/MPS3_fusion/simulate_FusionCancer/extractUniqueBreakPointFusionFromReformedJuncFastaFile.cpp

simulateFusionReadsFromUniqueBreakPointFusionFasta:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateFusionReadsFromUniqueBreakPointFusionFasta src/otherProjects/MPS3_fusion/simulate_FusionCancer/simulateFusionReadsFromUniqueBreakPointFusionFasta.cpp

analysisFusionBreakPointResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)analysisFusionBreakPointResults src/otherProjects/MPS3_fusion/simulate_FusionCancer/analysisFusionBreakPointResults.cpp

generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile src/otherProjects/MPS3_fusion/utils/generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile.cpp

generateFusionJuncInfoFromAdjustedBreakPointFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile src/otherProjects/MPS3_fusion/utils/generateFusionJuncInfoFromAdjustedBreakPointFile.cpp

compareFusionJuncResultsWithGTFannotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation src/otherProjects/MPS3_fusion/utils/compareFusionJuncResultsWithGTFannotation.cpp

filterFusionJuncInfoWithSupNumThres:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJuncInfoWithSupNumThres src/otherProjects/MPS3_fusion/utils/filterFusionJuncInfoWithSupNumThres.cpp

extractFusionGeneNameFromFusionInfoFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFusionGeneNameFromFusionInfoFile src/otherProjects/MPS3_fusion/utils/extractFusionGeneNameFromFusionInfoFile.cpp

compareTwoGenePairList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareTwoGenePairList src/otherProjects/MPS3_fusion/fusionReviewNAR/compareTwoGenePairList.cpp

compareFusionGeneInfoWithGenePairList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionGeneInfoWithGenePairList src/otherProjects/MPS3_fusion/fusionReviewNAR/compareFusionGeneInfoWithGenePairList.cpp

extractGeneNamePairFromFusionGroundTruthList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractGeneNamePairFromFusionGroundTruthList src/otherProjects/MPS3_fusion/fusionReviewNAR/extractGeneNamePairFromFusionGroundTruthList.cpp

compareStrandedFusionJuncResultsWithGTFannotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareStrandedFusionJuncResultsWithGTFannotation src/otherProjects/MPS3_fusion/utils/compareStrandedFusionJuncResultsWithGTFannotation.cpp

compareFusionDetectionResultsWithFusionMapGroundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionDetectionResultsWithFusionMapGroundTruth src/otherProjects/MPS3_fusion/utils/compareFusionDetectionResultsWithFusionMapGroundTruth.cpp

countMultiMapSegMaxLength_fusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countMultiMapSegMaxLength_fusionJunc src/MPS3SeqTools/countMultiMapSegMaxLength_fusionJunc.cpp

countMultiMapSegMaxLength_regularJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countMultiMapSegMaxLength_regularJunc src/MPS3SeqTools/countMultiMapSegMaxLength_regularJunc.cpp

compareJunc2groundTruth_outputTotalInfo:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareJunc2groundTruth_outputTotalInfo src/MPS3SeqTools/compareJunc2groundTruth_outputTotalInfo.cpp

extractSplicedAlignmentsFromSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSplicedAlignmentsFromSAM src/MPS3SeqTools/extractSplicedAlignmentsFromSAM.cpp

modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount src/MPS3SeqTools/modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount.cpp

separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount src/MPS3SeqTools/separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount.cpp

#BEGIN personal genome alignment
map2syntheticSNPtransSeq_includeGlobalAlignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2syntheticSNPtransSeq_includeGlobalAlignment src/otherProjects/incorporateGenomicVariants/map2syntheticSNPtransSeq_includeGlobalAlignment.cpp

map2syntheticSNPtransSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2syntheticSNPtransSeq src/otherProjects/incorporateGenomicVariants/map2syntheticSNPtransSeq.cpp

convertBeersGeneInfo2standardGAF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertBeersGeneInfo2standardGAF src/otherProjects/incorporateGenomicVariants/utils/convertBeersGeneInfo2standardGAF.cpp

extractSJfromGAF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSJfromGAF src/otherProjects/incorporateGenomicVariants/utils/extractSJfromGAF.cpp

reformatLuanSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatLuanSNP src/otherProjects/incorporateGenomicVariants/utils/reformatLuanSNP.cpp

generateSNPhash_fromBeersSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSNPhash_fromBeersSNPfile src/otherProjects/incorporateGenomicVariants/utils/generateSNPhash_fromBeersSNPfile.cpp

simulateGenomReads_SNPpos:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateGenomReads_SNPpos src/otherProjects/incorporateGenomicVariants/simulateGenomReads_SNPpos.cpp

simulateTransReads_SNPpos:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateTransReads_SNPpos src/otherProjects/incorporateGenomicVariants/simulateTransReads_SNPpos.cpp

generateSyntheticSNPtransSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSyntheticSNPtransSeq src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPtransSeq.cpp

generateSyntheticSNPgenomSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSyntheticSNPgenomSeq src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPgenomSeq.cpp

generateDecoyAlignHash_fromSimulatedTransReadAlignResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateDecoyAlignHash_fromSimulatedTransReadAlignResults src/otherProjects/incorporateGenomicVariants/utils/generateDecoyAlignHash_fromSimulatedTransReadAlignResults.cpp

comparePSsj2SDsj:
	g++ $(CFLAGS) -o $(OUTPUTDIR)comparePSsj2SDsj src/otherProjects/incorporateGenomicVariants/utils/comparePSsj2SDsj.cpp

checkSJflankStringChanges:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringChanges src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringChanges.cpp

checkSJflankStringChanges_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringChanges_withFormattedSNPfile src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringChanges_withFormattedSNPfile.cpp

extractSamWithSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithSJ src/otherProjects/incorporateGenomicVariants/utils/extractSamWithSJ.cpp

categorizeSJwithSiteAnnotatedOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)categorizeSJwithSiteAnnotatedOrNot src/otherProjects/incorporateGenomicVariants/utils/categorizeSJwithSiteAnnotatedOrNot.cpp

generateCandiPsSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCandiPsSJ src/otherProjects/incorporateGenomicVariants/utils/generateCandiPsSJ.cpp

getBlackListSamFromSimulatedReadSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBlackListSamFromSimulatedReadSam src/otherProjects/incorporateGenomicVariants/getBlackListSamFromSimulatedReadSam.cpp

getClosestSNPbesideSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getClosestSNPbesideSpliceSite src/otherProjects/incorporateGenomicVariants/utils/getClosestSNPbesideSpliceSite.cpp

checkSJexistencInMultiJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJexistencInMultiJuncFile src/otherProjects/incorporateGenomicVariants/utils/checkSJexistencInMultiJuncFile.cpp

compareTwoAlignInferJuncHash:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareTwoAlignInferJuncHash src/otherProjects/personalGenome/utils/compareTwoAlignInferJuncHash.cpp

compareAlignInferJuncHash2annotatedSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSJ src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSJ.cpp

compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset src/otherProjects/incorporateGenomicVariants/utils/compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset.cpp

compareAlignInferJuncHash2annotatedSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSpliceSite.cpp

compareAlignInferJuncHash2annotatedSpliceSite_checkSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite_checkSNP src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSpliceSite_checkSNP.cpp

compareGeneCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareGeneCount src/otherProjects/incorporateGenomicVariants/utils/compareGeneCount.cpp

checkSNPallelExpression:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSNPallelExpression src/otherProjects/incorporateGenomicVariants/utils/checkSNPallelExpression.cpp;

compareSnpAse_refBase2MPS3sd_altBase2MPS3ps:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSnpAse_refBase2MPS3sd_altBase2MPS3ps src/otherProjects/incorporateGenomicVariants/utils/compareSnpAse_refBase2MPS3sd_altBase2MPS3ps.cpp

checkSNPbesideSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSNPbesideSJ src/otherProjects/incorporateGenomicVariants/utils/checkSNPbesideSJ.cpp

checkPsSJfrequencyInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkPsSJfrequencyInPopulation src/otherProjects/incorporateGenomicVariants/utils/checkPsSJfrequencyInPopulation.cpp

checkTotalSJfrequencyInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkTotalSJfrequencyInPopulation src/otherProjects/incorporateGenomicVariants/utils/checkTotalSJfrequencyInPopulation.cpp

parseYiPaperPsSJresults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseYiPaperPsSJresults src/otherProjects/incorporateGenomicVariants/yiPaper/parseYiPaperPsSJresults.cpp

parseLuanAlignmentImprovementResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseLuanAlignmentImprovementResults src/otherProjects/incorporateGenomicVariants/utils/parseLuanAlignmentImprovementResults.cpp

labelSJexistenceAgainstAnotherSJlist:
	g++ $(CFLAGS) -o $(OUTPUTDIR)labelSJexistenceAgainstAnotherSJlist src/otherProjects/incorporateGenomicVariants/utils/labelSJexistenceAgainstAnotherSJlist.cpp

checkSJflankStringFrequency:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringFrequency src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringFrequency.cpp

checkSJcoverageInSampleWithSNPorNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJcoverageInSampleWithSNPorNot src/otherProjects/incorporateGenomicVariants/utils/checkSJcoverageInSampleWithSNPorNot.cpp

#END personal genome alignment

extractLowSupJuncFromJuncList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLowSupJuncFromJuncList src/MPS3SeqTools/extractLowSupJuncFromJuncList.cpp

extractSamWithJuncInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithJuncInList src/MPS3SeqTools/extractSamWithJuncInList.cpp

extractSamWithReadNameInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithReadNameInList src/otherProjects/incorporateGenomicVariants/utils/extractSamWithReadNameInList.cpp

extractFaWithReadNameInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFaWithReadNameInList src/otherProjects/incorporateGenomicVariants/utils/extractFaWithReadNameInList.cpp

countSpliceSiteNumFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSpliceSiteNumFromJuncFile src/MPS3SeqTools/countSpliceSiteNumFromJuncFile.cpp

countReadOverlapWithASM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadOverlapWithASM src/otherProjects/manohar/countReadOverlapWithASM.cpp

assignRedCloverGeneName2ASM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)assignRedCloverGeneName2ASM src/otherProjects/manohar/assignRedCloverGeneName2ASM.cpp

#BEGIN MPS3_paper
count_alignmentTypeNum_softClipBaseNum_fixedLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_fixedLength.cpp

count_alignmentTypeNum_softClipBaseNum_varyLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_varyLength.cpp

count_alignmentTypeNum_softClipBaseNum_fixedLength_includeUnpaired:
	g++ $(CFLAGS) -D INCLUDE_UNPAIRED -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength_includeUnpaired src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_fixedLength.cpp

count_alignmentTypeNum_softClipBaseNum_varyLength_includeUnpaired:
	g++ $(CFLAGS) -D INCLUDE_UNPAIRED -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength_includeUnpaired src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_varyLength.cpp

compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres src/otherProjects/MPS3_paper/compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres.cpp

generateSJsupNumCorrelation_sharedByAll:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSJsupNumCorrelation_sharedByAll src/otherProjects/MPS3_paper/generateSJsupNumCorrelation_sharedByAll.cpp

convertUnmappedBasePercAfterTrimming2beforeTrimming:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertUnmappedBasePercAfterTrimming2beforeTrimming src/otherProjects/MPS3_paper/utils/convertUnmappedBasePercAfterTrimming2beforeTrimming.cpp

count_mismatchDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)count_mismatchDistribution src/otherProjects/MPS3_paper/count_mismatchDistribution.cpp

produceMappedReadSeqFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceMappedReadSeqFromSam src/otherProjects/MPS3_paper/utils/produceMappedReadSeqFromSam.cpp
#END MPS3_paper

separateMergedChromFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateMergedChromFile src/otherProjects/manohar/separateMergedChromFile.cpp

#START ALU
separateAluFamilyRegion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateAluFamilyRegion src/otherProjects/Alu/separateAluFamilyRegion.cpp

generateAluSeqFasta:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateAluSeqFasta src/otherProjects/Alu/generateAluSeqFasta.cpp

sumUpAluReadCountForEachFamily_sailfish:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sumUpAluReadCountForEachFamily_sailfish src/otherProjects/Alu/sumUpAluReadCountForEachFamily_sailfish.cpp

segmentMap2singleAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)segmentMap2singleAluFamily src/otherProjects/Alu/segmentMap2singleAluFamily.cpp

segmentMap2multiAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)segmentMap2multiAluFamily src/otherProjects/Alu/segmentMap2multiAluFamily.cpp

mergeSegMapResultsAndAssignAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeSegMapResultsAndAssignAluFamily src/otherProjects/Alu/mergeSegMapResultsAndAssignAluFamily.cpp

#END ALU

parseEquineNCBI102AnnToGenerateSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSJ src/otherProjects/equineAnn/parseEquineNCBI102AnnToGenerateSJ.cpp

addGeneSymbolAccession2SJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)addGeneSymbolAccession2SJ src/otherProjects/equineAnn/addGeneSymbolAccession2SJ.cpp

getSpliceSiteFromSJ_withGeneSymbolIdAndExonID:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSpliceSiteFromSJ_withGeneSymbolIdAndExonID src/otherProjects/equineAnn/getSpliceSiteFromSJ_withGeneSymbolIdAndExonID.cpp

#BEGIN extract SJ from gtf
extractSJfromGTF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSJfromGTF src/otherProjects/extractSJfromGTF/extractSJfromGTF.cpp
#END extract SJ from gtf

reformFeatureCount2deseq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformFeatureCount2deseq src/otherProjects/featureCounts/reformFeatureCount2deseq.cpp

reformSailfish2tpm:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformSailfish2tpm src/otherProjects/sailfish/reformSailfish2tpm.cpp

clean: clean_annotation
	rm $(OUTPUTDIR)alignAll

clean_alignAll_speedUp:
	rm $(OUTPUTDIR)alignAll_speedUp

clean_alignAll:
	rm $(OUTPUTDIR)alignAll

clean_alignAll_phase1only:
	rm $(OUTPUTDIR)alignAll_phase1only:

clean_alignAll_initialMappingForFusionDetection:
	rm $(OUTPUTDIR)alignAll_initialMappingForFusionDetection

clean_detectFusion_post:
	rm $(OUTPUTDIR)detectFusion_post

clean_cirRNA:
	rm $(OUTPUTDIR)alignAll_cirRNA

clean_circularRNA_post:
	rm $(OUTPUTDIR)circularRNA_post

clean_map2annotationOnly:
	rm $(OUTPUTDIR)map2annotationOnly

clean_compareBackSpliceResultsWithValidatedList:
	rm $(OUTPUTDIR)compareBackSpliceResultsWithValidatedList

clean_fixHeadTail:
	rm $(OUTPUTDIR)fixHeadTail

clean_simulateCirRNAalignment2IGV:
	rm $(OUTPUTDIR)simulateCirRNAalignment2IGV

clean_extractBackSpliceFromAllJunc:
	rm $(OUTPUTDIR)extractBackSpliceFromAllJunc

clean_extractLariatFromBackSpliceSJ:
	rm $(OUTPUTDIR)extractLariatFromBackSpliceSJ

clean_extractLariatFromBackSpliceSJ_syntheticSAM:
	rm $(OUTPUTDIR)extractLariatFromBackSpliceSJ_syntheticSAM

clean_generateCirRNAtranscript:
	rm $(OUTPUTDIR)generateCirRNAtranscript

clean_extractBackSpliceFromSTARfusionJunc:
	rm $(OUTPUTDIR)extractBackSpliceFromSTARfusionJunc

clean_reformatMPS2circRNA:
	rm $(OUTPUTDIR)reformatMPS2circRNA

clean_mergeCirRNAreadCount:
	rm $(OUTPUTDIR)mergeCirRNAreadCount

clean_mergeBackSpliceJunctionWithReadCount:
	rm $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount

clean_reformBackSpliceSAM:
	rm $(OUTPUTDIR)reformBackSpliceSAM

clean_compareBackSplice2gtfAnnotation:
	rm $(OUTPUTDIR)compareBackSplice2gtfAnnotation

clean_compareBackSplice2gtfAnnotation_geneName:
	rm $(OUTPUTDIR)compareBackSplice2gtfAnnotation_geneName

clean_compareCirRNA2gtfAnnotation_withReadCountFile:
	rm $(OUTPUTDIR)compareCirRNA2gtfAnnotation_withReadCountFile

clean_fixHeadTail_cirRNA:
	rm $(OUTPUTDIR)fixHeadTail_cirRNA

clean_multi:
	rm $(OUTPUTDIR)alignAll_multi

clean_debug:
	rm $(OUTPUTDIR)alignAll_debug

clean_calTime:
	rm $(OUTPUTDIR)alignAll_calTime

clean_annotation:
	rm $(OUTPUTDIR)alignAll_annotation

clean_mapInfo:
	rm $(OUTPUTDIR)alignAll_mapInfo

clean_annotation_mapInfo:
	rm $(OUTPUTDIR)alignAll_annotation_mapInfo

clean_arrayList:
	rm $(OUTPUTDIR)alignAll_arrayList

clean_mergeSJfile:
	rm $(OUTPUTDIR)mergeSJfile

clean_getAlignInferJuncHashFromJuncFile:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile

clean_getAlignInferJuncHashFromSAMfile:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile

clean_getAlignInferJuncHashFromJuncFile_doAlterSpliceSite:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile_doAlterSpliceSite

clean_getAlignInferJuncHashFromSAMfile_doAlterSpliceSite:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite

clean_getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ

clean_getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity

clean_getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ

clean_getSmallExon:
	rm $(OUTPUTDIR)getSmallExon

clean_getSoftClipSeqSAM:
	rm $(OUTPUTDIR)getSoftClipSeqSAM

clean_nwdp:
	rm $(OUTPUTDIR)nwdp

clean_nwdp_oneEndOpen:
	rm $(OUTPUTDIR)nwdp_oneEndOpen

clean_getShortAnchorSJalignment:
	rm $(OUTPUTDIR)getShortAnchorSJalignment

clean_refineSAM:
	rm $(OUTPUTDIR)refineSAM

clean_refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity

clean_extendSJspliceSite:
	rm $(OUTPUTDIR)extendSJspliceSite

clean_splitSpliceSiteExtensionSJ:
	rm $(OUTPUTDIR)splitSpliceSiteExtensionSJ

clean_compareSpliceSite:
	rm $(OUTPUTDIR)compareSpliceSite

clean_spliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)spliceSiteAnchorSimilarity

clean_detectIncompatibleSAM:
	rm $(OUTPUTDIR)detectIncompatibleSAM

clean_filterJunc_anchorSimilarity:
	rm $(OUTPUTDIR)filterJunc_anchorSimilarity

clean_getPositiveStrandFa_MPS3sam:
	rm $(OUTPUTDIR)getPositiveStrandFa_MPS3sam

clean_classifySJ_withAnchorSimilarity:
	rm $(OUTPUTDIR)classifySJ_withAnchorSimilarity

clean_outputSJtypeNum:
	rm $(OUTPUTDIR)outputSJtypeNum

clean_parseSAM:
	rm $(OUTPUTDIR)parseSAM

clean_refineSAM_SE:
	rm $(OUTPUTDIR)refineSAM_SE

clean_refineSAM_PE:
	rm $(OUTPUTDIR)refineSAM_PE

clean_softClip:
	rm $(OUTPUTDIR)softClip

clean_softClip_varReadLength:
	rm $(OUTPUTDIR)softClip_varReadLength

clean_insertColumn2eachLine:
	rm $(OUTPUTDIR)insertColumn2eachLine

clean_SJsupportNumDiff:
	rm $(OUTPUTDIR)SJsupportNumDiff

clean_SJsupportNumDiff_targetSJ:
	rm $(OUTPUTDIR)SJsupportNumDiff_targetSJ

clean_getSJsupportNum_targetSJ_differentAligner:
	rm $(OUTPUTDIR)getSJsupportNum_targetSJ_differentAligner

clean_SJdiscoveryStats_supportNum:
	rm $(OUTPUTDIR)SJdiscoveryStats_supportNum

clean_countAlignmentType:
	rm $(OUTPUTDIR)countAlignmentType

clean_checkSJanchorSimilarity_fromJuncFile:
	rm $(OUTPUTDIR)checkSJanchorSimilarity_fromJuncFile

clean_getBothEndsReadWholeMatchFasta_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_groundTruth

clean_getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth

clean_getBothEndsReadWithoutSJfasta_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWithoutSJfasta_groundTruth

clean_convertBedtoolsResults2plotInR:
	rm $(OUTPUTDIR)convertBedtoolsResults2plotInR

clean_sam2alignInferJuncHash_supportNum:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum

clean_sam2alignInferJuncHash_supportNum_parallel:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_parallel

clean_sam2alignInferJuncHash_supportNum_anchorSize:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize

clean_sam2alignInferJuncHash_supportNum_anchorSize_parallel:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_parallel

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile

clean_compareSam2alignInferJuncHash:
	rm $(OUTPUTDIR)compareSam2alignInferJuncHash

clean_countSplicedAlignmentNum:
	rm $(OUTPUTDIR)countSplicedAlignmentNum

clean_separateForRevSam:
	rm $(OUTPUTDIR)separateForRevSam

clean_extractPeAlignmentWithInnerSoftClip:
	rm $(OUTPUTDIR)extractPeAlignmentWithInnerSoftClip

clean_transcriptCount:
	rm $(OUTPUTDIR)transcriptCount

clean_refineInnerSoftClip:
	rm $(OUTPUTDIR)refineInnerSoftClip

clean_exonSize:
	rm $(OUTPUTDIR)exonSize

clean_compareFusionRegion:
	rm $(OUTPUTDIR)compareFusionRegion

clean_normalize_foldChange_comparedTo1stSample:
	rm $(OUTPUTDIR)normalize_foldChange_comparedTo1stSample

clean_refineSAM_insertAnnotatedSJ2matchThroughEnd:
	rm $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd

clean_refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly:	
	rm $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly

clean_countReadsStartEndMappingPosWithinOverlapRegionList:
	rm $(OUTPUTDIR)countReadsStartEndMappingPosWithinOverlapRegionList

clean_compareSJfromMultiSamples:
	rm $(OUTPUTDIR)compareSJfromMultiSamples

clean_generateUnannotatedSJdistributionAlongRefGenome:
	rm $(OUTPUTDIR)generateUnannotatedSJdistributionAlongRefGenome

clean_globalMapIncompletePairSam:
	rm $(OUTPUTDIR)globalMapIncompletePairSam

clean_globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion:
	rm $(OUTPUTDIR)globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion

clean_remapAndCountEncompRead:
	rm $(OUTPUTDIR)remapAndCountEncompRead

clean_generateFusionJunc_filter_compare2gtf:
	rm $(OUTPUTDIR)generateFusionJunc_filter_compare2gtf

clean_remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint:
	rm $(OUTPUTDIR)remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint

clean_countUniqueUnpairedReadsEncompassingFusionBreakPoint:
	rm $(OUTPUTDIR)countUniqueUnpairedReadsEncompassingFusionBreakPoint

clean_analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair:
	rm $(OUTPUTDIR)analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair

clean_reformJuncFastaFromFusionCancer_oneLineSeqPerFusion:
	rm $(OUTPUTDIR)reformJuncFastaFromFusionCancer_oneLineSeqPerFusion

clean_extractUniqueBreakPointFusionFromReformedJuncFastaFile:
	rm $(OUTPUTDIR)extractUniqueBreakPointFusionFromReformedJuncFastaFile

clean_simulateFusionReadsFromUniqueBreakPointFusionFasta:
	rm $(OUTPUTDIR)simulateFusionReadsFromUniqueBreakPointFusionFasta

clean_analysisFusionBreakPointResults:
	rm $(OUTPUTDIR)analysisFusionBreakPointResults

clean_generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile:
	rm $(OUTPUTDIR)generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile

clean_generateFusionJuncInfoFromAdjustedBreakPointFile:
	rm $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile

clean_compareFusionJuncResultsWithGTFannotation:
	rm $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation

clean_filterFusionJuncInfoWithSupNumThres:
	rm $(OUTPUTDIR)filterFusionJuncInfoWithSupNumThres

clean_extractFusionGeneNameFromFusionInfoFile:
	rm $(OUTPUTDIR)extractFusionGeneNameFromFusionInfoFile

clean_compareTwoGenePairList:
	rm $(OUTPUTDIR)compareTwoGenePairList

clean_compareFusionGeneInfoWithGenePairList:
	rm $(OUTPUTDIR)compareFusionGeneInfoWithGenePairList

clean_compareStrandedFusionJuncResultsWithGTFannotation:
	rm $(OUTPUTDIR)compareStrandedFusionJuncResultsWithGTFannotation

clean_clusterFusionJunc:
	rm $(OUTPUTDIR)clusterFusionJunc

clean_separateFusionJuncStrandedOrNot:
	rm $(OUTPUTDIR)separateFusionJuncStrandedOrNot

clean_filterFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity

clean_filterFusionJunc_localIndexMap:
	rm $(OUTPUTDIR)filterFusionJunc_localIndexMap

clean_filterNonStrandedFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity

clean_filterOutBeersGroundTruthSJfromDetectedFusionJunc:
	rm $(OUTPUTDIR)filterOutBeersGroundTruthSJfromDetectedFusionJunc

clean_filterBackSpliceJunc:
	rm $(OUTPUTDIR)filterBackSpliceJunc

clean_clusterAndGenerateAlterFusion:
	rm $(OUTPUTDIR)clusterAndGenerateAlterFusion

clean_compareIntraChromSJwithGTF:
	rm $(OUTPUTDIR)compareIntraChromSJwithGTF

clean_reformatMPS2AnnotatedFusionResults:
	rm $(OUTPUTDIR)reformatMPS2AnnotatedFusionResults

clean_generateFusionSeq:
	rm $(OUTPUTDIR)generateFusionSeq

clean_countMultiMapSegMaxLength_fusionJunc:
	rm $(OUTPUTDIR)countMultiMapSegMaxLength_fusionJunc

clean_countMultiMapSegMaxLength_regularJunc:
	rm $(OUTPUTDIR)countMultiMapSegMaxLength_regularJunc

clean_compareJunc2groundTruth_outputTotalInfo:
	rm $(OUTPUTDIR)compareJunc2groundTruth_outputTotalInfo

clean_geneCount:
	rm $(OUTPUTDIR)geneCount

clean_compareFusionDetectionResultsWithFusionMapGroundTruth:
	rm $(OUTPUTDIR)compareFusionDetectionResultsWithFusionMapGroundTruth

clean_extractSplicedAlignmentsFromSAM:
	rm $(OUTPUTDIR)extractSplicedAlignmentsFromSAM

clean_cluster_filterAlterFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)cluster_filterAlterFusionJunc_anchorSeqSimilarity

clean_modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount:
	rm $(OUTPUTDIR)modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount

clean_reformatBeersSNP:
	rm $(OUTPUTDIR)reformatBeersSNP

clean_generateSNPhash_fromBeersSNPfile:
	rm $(OUTPUTDIR)generateSNPhash_fromBeersSNPfile

clean_simulateGenomReads_SNPpos:
	rm $(OUTPUTDIR)simulateGenomReads_SNPpos

clean_simulateTransReads_SNPpos:
	rm $(OUTPUTDIR)simulateTransReads_SNPpos

clean_generateDecoyAlignHash_fromSimulatedTransReadAlignResults:
	rm $(OUTPUTDIR)generateDecoyAlignHash_fromSimulatedTransReadAlignResults

clean_comparePSsj2SDsj:
	rm $(OUTPUTDIR)comparePSsj2SDsj

clean_checkSJflankStringChanges:
	rm $(OUTPUTDIR)checkSJflankStringChanges

clean_extractSamWithSJ:
	rm $(OUTPUTDIR)extractSamWithSJ

clean_getBlackListSamFromSimulatedReadSam:
	rm $(OUTPUTDIR)getBlackListSamFromSimulatedReadSam

clean_getClosestSNPbesideSpliceSite:
	rm $(OUTPUTDIR)getClosestSNPbesideSpliceSite

clean_extractSamWithJuncInList:
	rm $(OUTPUTDIR)extractSamWithJuncInList

clean_countSJcanonicalOrNotNumFromJuncFile:
	rm $(OUTPUTDIR)countSJcanonicalOrNotNumFromJuncFile

clean_countReadOverlapWithASM:
	rm $(OUTPUTDIR)countReadOverlapWithASM

clean_assignRedCloverGeneName2ASM:
	rm $(OUTPUTDIR)assignRedCloverGeneName2ASM

clean_count_alignmentTypeNum_softClipBaseNum_fixedLength:
	rm $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength

clean_count_alignmentTypeNum_softClipBaseNum_varyLength:
	rm $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength

clean_compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres:
	rm $(OUTPUTDIR)compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres

clean_generateSJsupNumCorrelation_sharedByAll:
	rm $(OUTPUTDIR)generateSJsupNumCorrelation_sharedByAll

clean_separateMergedChromFile:
	rm $(OUTPUTDIR)separateMergedChromFile

clean_separateAluFamilyRegion:
	rm $(OUTPUTDIR)separateAluFamilyRegion

clean_generateAluSeqFasta:
	rm $(OUTPUTDIR)generateAluSeqFasta

clean_sumUpAluReadCountForEachFamily_sailfish:
	rm $(OUTPUTDIR)sumUpAluReadCountForEachFamily_sailfish

clean_segmentMap2singleAluFamily:
	rm $(OUTPUTDIR)segmentMap2singleAluFamily

clean_parseEquineNCBI102AnnToGenerateSJ:
	rm $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSJ

clean_parseEquineNCBI102AnnToGenerateSpliceSite:
	rm $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSpliceSite

clean_extractSJfromGTF:
	rm $(OUTPUTDIR)extractSJfromGTF

clean_reformFeatureCount2deseq:
	rm $(OUTPUTDIR)reformFeatureCount2deseq

clean_reformSailfish2tpm:
	rm $(OUTPUTDIR)reformSailfish2tpm

clean_generateCandiPsSJ:
	rm $(OUTPUTDIR)generateCandiPsSJ

clean_addGeneSymbolAccession2SJ:
	rm $(OUTPUTDIR)addGeneSymbolAccession2SJ

clean_getSpliceSiteFromSJ_withGeneSymbolIdAndExonID:
	rm $(OUTPUTDIR)getSpliceSiteFromSJ_withGeneSymbolIdAndExonID

clean_extractSpecificAluFamilySeq:
	rm $(OUTPUTDIR)extractSpecificAluFamilySeq

clean_segmentMap2multiAluFamily:
	rm $(OUTPUTDIR)segmentMap2multiAluFamily

clean_mergeSegMapResultsAndAssignAluFamily:
	rm $(OUTPUTDIR)mergeSegMapResultsAndAssignAluFamily

clean_alignAll_outputUnpaired:
	rm $(OUTPUTDIR)alignAll_outputUnpaired

clean_convertUnmappedBasePercAfterTrimming2beforeTrimming:
	rm $(OUTPUTDIR)convertUnmappedBasePercAfterTrimming2beforeTrimming

clean_count_mismatchDistribution:
	rm $(OUTPUTDIR)count_mismatchDistribution

clean_produceMappedReadSeqFromSam:
	rm $(OUTPUTDIR)produceMappedReadSeqFromSam

clean_compareTwoAlignInferJuncHash:
	rm $(OUTPUTDIR)compareTwoAlignInferJuncHash

clean_compareAlignInferJuncHash2annotatedSJ:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSJ

clean_compareAlignInferJuncHash2annotatedSpliceSite:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite

clean_compareAlignInferJuncHash2annotatedSpliceSite_checkSNP:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite_checkSNP

clean_generateSyntheticSNPtransSeq:
	rm $(OUTPUTDIR)generateSyntheticSNPtransSeq

clean_alignAll_personalizedSeq_SNPseqIndex:
	rm $(OUTPUTDIR)alignAll_personalizedSeq_SNPseqIndex

clean_convertBeersGeneInfo2standardGAF:
	rm $(OUTPUTDIR)convertBeersGeneInfo2standardGAF

clean_map2syntheticSNPtransSeq:
	rm $(OUTPUTDIR)map2syntheticSNPtransSeq

clean_map2syntheticSNPtransSeq_includeGlobalAlignment:
	rm $(OUTPUTDIR)map2syntheticSNPtransSeq_includeGlobalAlignment

clean_extractSJfromGAF:
	rm $(OUTPUTDIR)extractSJfromGAF

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ

clean_checkSNPallelExpression:
	rm $(OUTPUTDIR)checkSNPallelExpression

clean_compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset

clean_compareSnpAse_refBase2MPS3sd_altBase2MPS3ps:
	rm $(OUTPUTDIR)compareSnpAse_refBase2MPS3sd_altBase2MPS3ps

clean_checkSNPbesideSJ:
	rm $(OUTPUTDIR)checkSNPbesideSJ

clean_checkPsSJfrequencyInPopulation:
	rm $(OUTPUTDIR)checkPsSJfrequencyInPopulation

clean_checkTotalSJfrequencyInPopulation:
	rm $(OUTPUTDIR)checkTotalSJfrequencyInPopulation

clean_checkSJflankStringChanges_withFormattedSNPfile:
	rm $(OUTPUTDIR)checkSJflankStringChanges_withFormattedSNPfile

clean_parseYiPaperPsSJresults:
	rm $(OUTPUTDIR)parseYiPaperPsSJresults

clean_parseLuanAlignmentImprovementResults:
	rm $(OUTPUTDIR)parseLuanAlignmentImprovementResults

clean_labelSJexistenceAgainstAnotherSJlist:
	rm $(OUTPUTDIR)labelSJexistenceAgainstAnotherSJlist

clean_checkSJcoverageInSampleWithSNPorNot:
	rm $(OUTPUTDIR)checkSJcoverageInSampleWithSNPorNot

clean_mpsI_arrayList:
	rm $(OUTPUTDIR)mpsI_arrayList

clean_MapSplice:
	rm $(OUTPUTDIR)MapSplice

clean_iMapSplice:
	rm $(OUTPUTDIR)iMapSplice