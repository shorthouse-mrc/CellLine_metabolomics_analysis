# PID of current job: 424585
mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("Aniline","Pyrrole-2-carboxylic acid","1-Pyrroline-5-carboxylic acid","1-Pyrroline-2-carboxylic acid","2-Aminomuconic acid","N-Acetyl-L-glutamate 5-semialdehyde","N-Acetyl-L-aspartic acid","N-Formyl-L-glutamic acid","Phosphohydroxypyruvic acid","Argininosuccinic acid","2-Phenylethanol glucuronide","Beta-Citryl-L-glutamic acid","Hesperetin","Sphingosine 1-phosphate","CPA(18:1(11Z)/0:0)","CPA(18:1(9Z)/0:0)","Alpha-Tocopherol","DG(14:1(9Z)/15:0/0:0)","DG(15:0/14:1(9Z)/0:0)","PG(16:1(9Z)/20:4(5Z,8Z,11Z,14Z))","PG(18:2(9Z,12Z)/18:3(6Z,9Z,12Z))","PG(18:2(9Z,12Z)/18:3(9Z,12Z,15Z))","PG(18:3(6Z,9Z,12Z)/18:2(9Z,12Z))","PG(18:3(9Z,12Z,15Z)/18:2(9Z,12Z))","PG(16:0/22:5(4Z,7Z,10Z,13Z,16Z))","PG(16:0/22:5(7Z,10Z,13Z,16Z,19Z))","PG(16:1(9Z)/22:4(7Z,10Z,13Z,16Z))","PG(18:1(11Z)/20:4(5Z,8Z,11Z,14Z))","PG(18:1(9Z)/20:4(5Z,8Z,11Z,14Z))","PG(18:2(9Z,12Z)/20:3(5Z,8Z,11Z))","PG(18:2(9Z,12Z)/20:3(8Z,11Z,14Z))","PG(16:0/22:4(7Z,10Z,13Z,16Z))","PG(18:0/20:4(5Z,8Z,11Z,14Z))","PG(18:1(11Z)/20:3(5Z,8Z,11Z))","PG(18:1(11Z)/20:3(8Z,11Z,14Z))","PG(18:1(9Z)/20:3(5Z,8Z,11Z))","PG(18:1(9Z)/20:3(8Z,11Z,14Z))","PI(16:0/16:1(9Z))","PI(16:1(9Z)/16:0)","PI(16:0/16:0)","PG(18:3(6Z,9Z,12Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))","PG(18:3(9Z,12Z,15Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))","Lactosylceramide (d18:1/12:0)","PI(16:0/18:1(11Z))","PI(16:0/18:1(9Z))","PI(16:1(9Z)/18:0)","PI(18:0/16:1(9Z))","PI(18:1(11Z)/16:0)","PI(18:1(9Z)/16:0)")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA)
