# PID of current job: 403121
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-SetOrganism(mSet, "hsa")
geneListFile<-"replace_with_your_file_name"
geneList<-readChar(geneListFile, file.info(geneListFile)$size)
mSet<-PerformGeneMapping(mSet, geneList, "hsa", "symbol");
cmpdListFile<-"replace_with_your_file_name"
cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
mSet<-PerformCmpdMapping(mSet, cmpdList, "hsa", "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "(R)-1-Methylpiperidine-2-carboxylic acid");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "(R)-1-Methylpiperidine-2-carboxylic acid", "");
mSet<-PrepareIntegData(mSet);
mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "integ", "query");
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA, NA, NA )
mSet<-CreateIntegMatchingTable(mSet);
mSet<-PlotKEGGPath(mSet, "Lysine degradation",566, 490, "png", NULL)
mSet<-RerenderMetPAGraph(mSet, "zoom1739305829098.png",566.0, 490.0, 100.0)
mSet<-SaveTransformedData(mSet)
