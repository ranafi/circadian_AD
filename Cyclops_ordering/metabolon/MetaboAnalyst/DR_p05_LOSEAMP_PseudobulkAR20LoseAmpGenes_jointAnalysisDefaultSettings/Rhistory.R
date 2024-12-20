# PID of current job: 1090140
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-SetOrganism(mSet, "hsa")
geneListFile<-"replace_with_your_file_name"
geneList<-readChar(geneListFile, file.info(geneListFile)$size)
mSet<-PerformGeneMapping(mSet, geneList, "hsa", "symbol");
cmpdListFile<-"replace_with_your_file_name"
cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
mSet<-PerformCmpdMapping(mSet, cmpdList, "hsa", "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PrepareIntegData(mSet);
mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "integ", "query");
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA, NA, NA )
mSet<-CreateIntegMatchingTable(mSet);
mSet<-SaveTransformedData(mSet)
