#!/bin/sh

# IMPORTANT: always leave two empty spaces at the end of each md file,
# otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
# pandoc *.md
# -o pdf/practical.pdf
# --toc
# --variable title:"CNV calling partical - EBI cancer genomics workshops 2014"
# --variable date:"2014/06/26"
# --variable author:"Mathieu Bourgey, Ph.D"
# --variable links-as-notes
# --variable linkcolor:black
# --variable urlcolor:black
# --variable geometry:margin=3cm

pandoc README.md analysis/ASCAT.md analysis/BinCount.md analysis/DNAcopy.md solutions/1.dataDiff.md solutions/2.NgsAnalysisSummary.md solutions/3.BinCount.md solutions/4.10DnaCopy_calls.md solutions/4.11DnaCopy_output.md solutions/4.12DnaCopy_exam.md solutions/4.1DnaCopy_Ropen.md solutions/4.2DnaCopy_loadDC.md solutions/4.3DnaCopy_loadData.md solutions/4.4DnaCopy_cleanData.md solutions/4.5DnaCopy_normData.md solutions/4.6DnaCopy_LRR.md solutions/4.7DnaCopy_CNA.md solutions/4.8DnaCopy_smooth.md solutions/4.9DnaCopy_segments.md solutions/4.DnaCopy_explanation.md solutions/5.cancerChallenge.md solutions/6.SnpAnalysisSummary.md solutions/7.10ascat_saveCNA.md solutions/7.11ascat_exam.md solutions/7.1ascat_Ropen.md solutions/7.2ascat_loadM.md solutions/7.3ascat_loadD.md solutions/7.4ascat_plotRaw.md solutions/7.5ascat_segments.md solutions/7.6ascat_plotSeg.md solutions/7.7ascat_ASCP.md solutions/7.8ascat_savePlo.md solutions/7.9ascat_CNA.md -o pdf/practical.pdf --toc --variable title:"CNV calling partical - EBI cancer genomics workshops 2014" --variable date:"2014/06/26" --variable author:"Mathieu Bourgey, Ph.D" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm 


