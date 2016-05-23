# try to parse csv of chromatogram from gc-ms data
# want to find the difference in spectra of RJG_MS_7(WY1335 8hr DT Sup, Hour 48)
# vs RJG_MS_3 (SD+10 uM Amino Acids)

# ill plot these two chromatograms. I expect them to look reasonably
# similar. Then i calculate the difference in the two spectra to find
# what peaks have 'appeared' 

# i may need to manually retain a few peaks (26.2,30.6 for met, cys respectively)
# since those are also present in MS_3, which I'm designated as my blank for now


require(ggplot2)

setwd('X:/fast/shou_w/GC-MS/UW Chemistry/Agilent5975/20160521/CSV/')


All_Spectra<-NULL

#RJG_MS_45
{
r<-read.csv('20160521_SampleMS45.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_45',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_45=r

}

#RJG_MS_46
{
r<-read.csv('20160521_SampleMS46.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_46',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_46=r

}

#RJG_MS_47
{
r<-read.csv('20160521_SampleMS47.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_47',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_47=r

}


#RJG_MS_48
{
r<-read.csv('20160521_SampleMS48.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_48',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_48=r

}

#RJG_MS_49
{
r<-read.csv('20160521_SampleMS49.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_49',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_49=r

}

#RJG_MS_50
{
r<-read.csv('20160521_SampleMS50.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_50',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_50=r

}

#RJG_MS_51
{
r<-read.csv('20160521_SampleMS51.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_51',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_51=r

}

#RJG_MS_52
{
r<-read.csv('20160521_SampleMS52.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_52',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_52=r

}

#RJG_MS_53
{
r<-read.csv('20160521_SampleMS53.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_53',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_53=r

}

#RJG_MS_54
{
r<-read.csv('20160521_SampleMS54.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_54',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_54=r

}

#RJG_MS_55
{
r<-read.csv('20160521_SampleMS55.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_55',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_55=r

}

#RJG_MS_56
{
r<-read.csv('20160521_SampleMS56.csv',header=T)

retention<-data.frame()

for(i in 1:nrow(r))
{
  minutes<-r$RT.minutes.[i]
  print(minutes)
  t<-r[which(r$RT.minutes.==minutes),]
  t<-t[,4:ncol(t)]
  s<-sum(t[1,])
  df<-data.frame(minutes,s)
  colnames(df)<-c('Retention','Abundance')
  retention<-rbind(retention,df)
}

name<-rep('RJG_MS_56',times=nrow(retention))
retention<-cbind(retention,name)
All_Spectra$MS_56=r

}

#save.image('20160421_SpectraAll.RData')
save.image('../Analysis/20160523_SpectraAll.RData')

# maybe one thing I can try is to 'scan' for peaks at a known m/z since I know aproximately what numbers 
# to look for. For example, met has a characteristic peak at 320. Try to extract spectra with elevated levels 
# and find retention time
