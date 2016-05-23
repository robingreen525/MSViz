
require(ggplot2)

#plotting theme
{
  mytheme =   list(
    #geom_point(size=3.5),
    theme_bw(),
    theme(axis.text.x = element_text(face="bold",size=15),strip.text.x=element_text(face="bold",size=12.5),
          axis.text.y = element_text(face="bold",size=15),strip.text.y=element_text(face="bold",size=12.5),
          text=element_text(size=15,face="bold"))
    
  )
}



extractPeak<-function(obj,search_list)
{
  #obj <- the compound list that we're going to be searching for
  # search_list<- a vector that has the full GC/MS spectra names that will be searched
  
  
  for(search in search_list)
  {
    #print(search)
    
    is_57=FALSE
    is_85=FALSE
    is_159=NULL
    compound<-unique(as.character(obj$Metabolite))
    cat("Searching ", search,'for ',compound,' m/z peaks..\n')
    f<-acquireSpectra(search)
    #print(f)
    
    print(compound)
    # search for peak-57
    
    peaks<-data.frame(obj$Peak.57,obj$Peak.85,obj$Peak.159)
    #print(peaks)
    
    #get all non-zero peak-57
    {
      p57<-peaks[1]
      p57<-as.integer(p57) # m/z needs to be interger value
      #print(p57)
      strname<-paste('X',p57,sep='')
     # print(strname)
      

      peaks<-f[[strname]]
      retention<-f[["RT.minutes."]]
      spectra<-data.frame(retention,peaks)
      print(head(spectra))

      # now I want to see if there are signals in the RT spectra
      
      # look for signal greater than 100 abundance units in the RT vs. peak-57 
      query_spectra<-spectra[which(spectra$peaks >=1000),]
      if(nrow(query_spectra)>0) # if positive hit
      {
        msg<-cat("In sample",search,'there is a signal for Peak-57 for',compound,'\n')
        title=paste(compound,"Peak-57",p57,search,sep='__')
        g<-ggplot()+
          geom_point(data=spectra,aes(x=retention,y=peaks))+
          geom_line(data=spectra,aes(x=retention,y=peaks))+
          geom_text(data=query_spectra,aes(x=retention,y=peaks,label=peaks))+
          mytheme+ggtitle(title)+xlab('Rentention Time (min)')+
          ylab('Abundance')
       # print(g)
        wd<-getwd()
        path=paste(wd,'R_Plots',sep='/')
        filename<-paste(path,title,sep='/')
        #print(filename)
        spectra57<-spectra
        is_57<-TRUE
        ggsave(g,file=(paste(filename,'png',sep='.')))
        write.csv(spectra,file=(paste(filename,'csv',sep='.')))
      }
            
    }
    #get all non-zero peak-85
    {
    peaks<-data.frame(obj$Peak.57,obj$Peak.85,obj$Peak.159)    
    p85<-peaks[2]
    p85<-as.integer(p85) # m/z needs to be interger value
    #print(p85)
    strname<-paste('X',p85,sep='')
   # print(strname)
    
    
    peaks<-f[[strname]]
    retention<-f[["RT.minutes."]]
    spectra<-data.frame(retention,peaks)
    
    
    # now I want to see if there are signals in the RT spectra
    
    # look for signal greater than 100 abundance units in the RT vs. peak-57 
    query_spectra<-spectra[which(spectra$peaks >=1000),]
    if(nrow(query_spectra)>0) # if positive hit
    {
      msg<-cat("In sample",search,'there is a signal for Peak-85 for',compound,'\n')
      title=paste(compound,"Peak-85",p85,search,sep='__')
      g<-ggplot()+
        geom_point(data=spectra,aes(x=retention,y=peaks))+
        geom_line(data=spectra,aes(x=retention,y=peaks))+
        geom_text(data=query_spectra,aes(x=retention,y=peaks,label=peaks))+
        mytheme+ggtitle(title)+xlab('Rentention Time (min)')+
        ylab('Abundance')
      #print(g)
      wd<-getwd()
      path=paste(wd,'R_Plots',sep='/')
      filename<-paste(path,title,sep='/')
     # print(filename)
      spectra85<-spectra
      is_85<-TRUE
     ggsave(g,file=(paste(filename,'png',sep='.')))
    }
    
    }
    
    #get all non-zero peak-159
    {
    peaks<-data.frame(obj$Peak.57,obj$Peak.85,obj$Peak.159)    
    p159<-peaks[3]
    p159<-as.integer(p159) # m/z needs to be interger value
   # print(p159)
    strname<-paste('X',p159,sep='')
    #print(strname)
    
    
    peaks<-f[[strname]]
    retention<-f[["RT.minutes."]]
    spectra<-data.frame(retention,peaks)
    
    
    # now I want to see if there are signals in the RT spectra
    
    # look for signal greater than 100 abundance units in the RT vs. peak-57 
    query_spectra<-spectra[which(spectra$peaks >=1000),]
    if(nrow(query_spectra)>0) # if positive hit
    {
      msg<-cat("In sample",search,'there is a signal for Peak-159 for',compound,'\n')
      title=paste(compound,"Peak-159",p159,search,sep='__')
      g<-ggplot()+
        geom_point(data=spectra,aes(x=retention,y=peaks))+
        geom_line(data=spectra,aes(x=retention,y=peaks))+
        geom_text(data=query_spectra,aes(x=retention,y=peaks,label=peaks))+
        mytheme+ggtitle(title)+xlab('Rentention Time (min)')+
        ylab('Abundance')
     # print(g)
      wd<-getwd()
      path=paste(wd,'R_Plots',sep='/')
      filename<-paste(path,title,sep='/')
      #print(filename)
      spectra159<-spectra
      is_159<-TRUE
      ggsave(g,file=(paste(filename,'png',sep='.')))
    }
    
    }
    
    if( is_57 && is_85 && is_159)
    {
       print('make overlay plot')
      
       names<-rep('peak-57',times=nrow(spectra57))
       spectra57<-cbind(spectra57,names)
       print(head(spectra57))
       
       names<-rep('peak-85',times=nrow(spectra85))
       spectra85<-cbind(spectra85,names)
       print(head(spectra85))
       
       
       names<-rep('peak-159',times=nrow(spectra159))
       spectra159<-cbind(spectra159,names)
       print(head(spectra159))
       
       spectra_all<-rbind(spectra57,spectra85,spectra159)
       print(head(spectra_all))
       colnames(spectra_all)<-c('retention','peaks','type')
       
       title=paste(compound,"All",search,sep='__')
       g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
         geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
         geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
         #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=peaks,col=type))+
         mytheme+ggtitle(title)+xlab('Rentention Time (min)')+
         ylab('Abundance')+
         facet_wrap(~type,ncol=1)
       
       
       path=paste(wd,'R_Plots/compiled',sep='/')
       filename<-paste(path,title,sep='/')
       
     # print(g)
      ggsave(g,file=(paste(filename,'png',sep='.')))
      
      
      title=paste(compound,"LOG10_All",search,sep='__')
      g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
        geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
        geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
        #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=peaks,col=type))+
        mytheme+ggtitle(title)+xlab('Rentention Time (min)')+
        ylab('Abundance')+scale_y_log10()+annotation_logticks(sides='l')+
        facet_wrap(~type,ncol=1)
      
      
      path=paste(wd,'R_Plots/compiled',sep='/')
      filename<-paste(path,title,sep='/')
      
      #print(g)
      ggsave(g,file=(paste(filename,'png',sep='.')))
      
      #filter spectra
      filterMatch(spectra_all,1000,path,compound,search)
      
    }
    # want to make an overlay plot of all peak
    
    
    
    
    
    print('+++++++++++++++++++++++++++++++++++')
  }
  
}


# this function will get the selected sample of interest
# that will be searched
acquireSpectra<-function(search)
{
  spec<-All_Spectra[[search]]
  return(spec)
  
  
  
}

# want a function that looks for high signal for peak-57,peak-85,and peak-157 at same RT, filters out other peaks

filterMatch<-function(spectra_all,threshold=1000,path,compound,search,xmin=25,xmax=27)
{
  class<-c('peak-57','peak-85','peak-159')
  
  p57<-spectra_all[which(spectra_all$type==class[1]),]
  for(i in 1:nrow(p57)) # if 
  {
    peak_signal<-p57$peaks[i]
    if(peak_signal<threshold)
    {
      p57$peaks[i]<-0
     # print(threshold)
    }
    
    
  }
  
  p85<-spectra_all[which(spectra_all$type==class[2]),]
  for(i in 1:nrow(p85)) # if 
  {
    peak_signal<-p85$peaks[i]
    if(peak_signal<threshold)
    {
      p85$peaks[i]<-0
      #print(threshold)
    }
    
    
  }
  
  p159<-spectra_all[which(spectra_all$type==class[3]),]
  for(i in 1:nrow(p159)) # if 
  {
    peak_signal<-p159$peaks[i]
    if(peak_signal<threshold)
    {
      p159$peaks[i]<-0
     # print(threshold)
    }
    
    
  }
 
  spectra_all<-rbind(p57,p85,p159) 
  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=peaks,col=type))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~type,ncol=1)
 # print(g)
  
  # now only focus on RT where abundance for each m/z is greater than threshold
  
  for(i in 1:nrow(p57))
  {
    time<-p57$retention[i]
    
    a<-p57[which(p57$retention == time),]
    a<-as.numeric(a$peaks)
    
    b<-p85[which(p85$retention == time),]
    b<-as.numeric(b$peaks)
    
    c<-p159[which(p159$retention == time),]
    c<-as.numeric(c$peaks)
   
    #print(a)
    
    #if(a<threshold | b<threshold | c <threshold)
    #{
     # spectra_all[which(spectra_all$retention==time),]$peaks<-0
      
      
    #}
    
    if(!(a>=threshold && b>=threshold && c>=threshold))
    {
      spectra_all[which(spectra_all$retention==time),]$peaks<-0
      
    }
    
  }
  title=paste(compound,'filtered',search,paste('Threshold',threshold,sep ='='),sep='__')
  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=peaks,col=type))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~type,ncol=1)+
    ggtitle(title)
  
  filename<-paste(path,'filtered',title,sep='/')
  
  ggsave(g,file=(paste(filename,'png',sep='.')))
 # print(g) 
  
  
  title=paste(compound,'filtered','LOG10_ALL',paste('Threshold',threshold,sep ='='),search,sep='__')
  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=peaks,col=type))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~type,ncol=1)+
    ggtitle(title)+
    scale_y_log10()+annotation_logticks(sides='l')
  
  filename<-paste(path,'filtered',title,sep='/')
  
  ggsave(g,file=(paste(filename,'png',sep='.')))
  #print(g) 
  
  title=paste(compound,'filtered','Labeled_LOG10_ALL',paste('Threshold',threshold,sep ='='),search,sep='__')
  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_text(data=spectra_all,aes(x=retention,y=peaks,label=retention,col=type))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~type,ncol=1)+
    ggtitle(title)+
    scale_y_log10()+annotation_logticks(sides='l')
  
  filename<-paste(path,'filtered',title,sep='/')
  
  ggsave(g,file=(paste(filename,'png',sep='.')))
  #print(g) 
  
  title=paste(compound,'filtered','Labeled_LOG10_ALL_MinMax',paste('Threshold',threshold,sep ='='),search,sep='__')
  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=type))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=retention,col=type))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~type,ncol=1)+
    ggtitle(title)+
    scale_y_log10()+annotation_logticks(sides='l')+xlim(xmin,xmax)
    #print(g)
  
    filename<-paste(path,'filtered',title,sep='/')
    print(filename)
    ggsave(g,file=(paste(filename,'png',sep='.')))
#   
  
  
  
}

plotMZforSampleRT<-function(obj,samplename,retentiontime,specname,lower=0.03,upper=0.09,xlimit=200,labels=c('218','292','320'))
{

  
  min<-as.numeric(retentiontime)-lower
  max<-as.numeric(retentiontime)+upper# this will likely take some tuning depending on the peak
  cat(retentiontime,min,max,'\n')
	
  # get range of spectra
  t<-obj[which(obj$RT.minutes. >min),]
  t<-t[which(t$RT.minutes. < max),]
		
   
   plot_dat<-data.frame()
   
   print(head(t))
   for(i in 1:nrow(t))
   {
   	  time<-t[i,]$RT.minutes.
   	  temp<-t[i,4:ncol(t)]
   	  #print(temp)
   	  abundance<-sum(temp[1,])
   	  #print(abundance)
   	  df<-data.frame(time,abundance)
   	  colnames(df)<-c('time','abundance')
   	  plot_dat<-rbind(plot_dat,df)
   
   }	
#    
   g<-ggplot()+geom_point(data=plot_dat,aes(x=time,y=abundance))+
     geom_line(data=plot_dat,aes(x=time,y=abundance))+mytheme+
     xlab('Time (Minutes) in Column')+ylab('Abundance Units')
   print(g)
   ret<-g
   
   s<-g
   # plot full m/z now
   
   full_spec<-data.frame()
   temp<-t[,4:ncol(t)]
   
   val<-colSums(temp)
   name<-(colnames(temp))
   print(val)
   mz<-data.frame()
   for(i in 1:length(val))
   {
   	 # print(name[i])
   	  #masstocharge_name<-substr(as.character(name[i]),start=2,length(name[i]))
   	  #print(masstocharge_name)	
   	  #print(length(as.character(name[i])))
   	  n<-as.character(name[i])
   	 # print(n)
   	  #print(class(n))
   	  n<-(substring(n,2))
   	  n<-as.integer(n)
   	  #print(n)
   	  df<-data.frame(n,val[i])
   	  colnames(df)<-c('moverz','abundance')
   	  mz<-rbind(mz,df)
   
   }
   
   #add labels for m/z of interest
   for(i in 1:length(labels))
   {
     lab=labels[i]
     print(lab)
     
     
   }
   
   print(mz[which(mz$moverz %in% labels),])
   
   
   g<-ggplot()+#geom_point(data=mz,aes(x=moverz,y=abundance))+
   geom_segment(data=mz[which(mz$moverz >= xlimit),],aes(x=moverz,xend=moverz,y=0,yend=abundance))+mytheme+
   xlab('m/z')+ylab('Abundance Units')+geom_text(data=mz[which(mz$moverz %in% labels),],aes(x=moverz,y=1.05*abundance,label=moverz))
   s<-g
   print(g)  
   
   
   specname<-rep(specname,nrow(mz))
   #print(specname)
   
   mz<-cbind(mz,specname)
   #print(head(mz))
    
   p <- annotation_custom(grob=g,xmin = 30, xmax = 60, ymin = 60, ymax = 100)
   q<- annotation_custom(grob = s, xmin = 0.5, xmax = 0.9, ymin = 6, ymax = 10)
  
   blanktheme<- theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
      
  p = qplot(1:10, 1:10, log='x')+xlab(NULL)+ylab(NULL)+blanktheme
  s<-ggplotGrob(s)

  s<-annotation_custom(grob = s, xmin = -0.005, xmax = 1.05, ymin = 0.5, ymax = 10.5)
    
	p<-p+s
  print(p)	
	
	g<-ggplotGrob(ret)
  p<-p + annotation_custom(grob = g, xmin = 0.456, xmax = 1.04 ,ymin = 6, ymax = 10.35)
 	print(p)
  
  #ggsave(filename='temp.png',plot=p)

 	title<-paste(samplename,retentiontime,paste('',min,'-',max,'minutes',sep=''),sep='_')
 	
 	p<-p+ggtitle(title)+theme(plot.title = element_text(lineheight=2, face="bold"))
 	print(p)
 
 	path=getwd()
 	path=paste(path,'R_Plots/slices',sep='/')
 	filename=paste(path,paste(title,'pdf',sep='.'),sep='/')
 	#filename=paste(path,'temp.pdf',sep='/')
 	print(filename)
  
 	ggsave(filename=filename,plot=p)
 	
 	filename=paste(path,paste(title,'png',sep='.'),sep='/')
 	ggsave(filename=filename,plot=p)
 	
 	# temp expanding of M/z plot for power point. comment out later
 	mz<-mz[which(mz$moverz > xlimit),]
 	mz<-mz[which(mz$moverz < 400),]
 	g<-ggplot()+#geom_point(data=mz,aes(x=moverz,y=abundance))+
 	  geom_segment(data=mz,aes(x=moverz,xend=moverz,y=0,yend=abundance))+mytheme+
 	  xlab('m/z')+ylab('Abundance Units')+geom_text(data=mz[which(mz$moverz %in% labels),],aes(x=moverz,y=1.05*abundance,label=moverz))

 #	print(g)
 	ggsave("tempspectra.png")
 	
 	
 	return(mz)
   
  
   
}

compareMZ<-function(spec_1,spec_2,xlimit=0)
{
  #get relative
  
  
  
  
 	spec_2$normabundance<-spec_2$normabundance*-1
 	all_spec<-rbind(spec_1,spec_2)
 	
 	g<-ggplot()+
    geom_segment(data=all_spec[which(all_spec$moverz >= xlimit),],aes(x=moverz,xend=moverz,y=0,yend=abundance,col=specname))+mytheme+
    xlab('m/z')+ylab('abundance')+
 	 xlim(xlimit,600)
  print(g)



}

#generate plot of full spectra (abundance vs RT)
plotFullSpectra<-function(obj,samplename)
{
  print(head(obj))
  t<-obj
  plot_dat<-data.frame()
  
  for(i in 1:nrow(t))
  {
    
    time<-t[i,]$RT.minutes.
    print(time)
    temp<-t[i,4:ncol(t)]
    abundance<-sum(temp[1,])
    #print(abundance)
    df<-data.frame(time,abundance)
    colnames(df)<-c('time','abundance')
    plot_dat<-rbind(plot_dat,df)
    
  }	
  
  g<-ggplot()+geom_point(data=plot_dat,aes(x=time,y=abundance))+
    geom_line(data=plot_dat,aes(x=time,y=abundance))+mytheme+
    scale_y_log10()+annotation_logticks(sides='l')
  print(g)
  
  g<-ggplot()+geom_point(data=plot_dat,aes(x=time,y=abundance))+
    geom_line(data=plot_dat,aes(x=time,y=abundance))+mytheme
  print(g)
  
  
  g<-ggplot()+geom_point(data=plot_dat[which(plot_dat$time >15),],aes(x=time,y=abundance))+
    geom_line(data=plot_dat[which(plot_dat$time >15),],aes(x=time,y=abundance))+mytheme+
    xlab('Time (minutes) in Column') +ylab('Abundance Units')
  print(g)
  
  path<-getwd()
  path=paste(path,'R_Plots/fullspectra',sep='/')
  print(path)
  
  filename=paste(samplename,'FullSpectraPast15',sep='_')
  filename=paste(filename,'png',sep='.')
  print(filename)
  filename<-paste(path,filename,sep='/')
  ggsave(g,file=filename)
}


#background subtract m/z
BackgroundSubtractMZ<-function(spec,background)
{
  spec$abundance<-spec$abundance-background$abundance
  
  g<-ggplot()+
    geom_segment(data=spec,aes(x=moverz,xend=moverz,y=0,yend=abundance,col=specname))+mytheme+
    xlab('m/z')+ylab('abundance')
  print(g)
  
  return(spec)
}


#normalize m/z abundances by mean

NormalizeMZ<-function(spec)
{
    print(head(spec))
    mean_abundance<-mean(spec$abundance)
    print(mean_abundance)
	spec$normabundance<-spec$abundance/mean_abundance
	print(head(spec))

    g<-ggplot()+
    geom_segment(data=spec,aes(x=moverz,xend=moverz,y=0,yend=normabundance,col=specname))+mytheme+
    xlab('m/z')+ylab('relative abundance')
    print(g)
    
  return(spec)
}


SubstractSpectraforMZ<-function(spec1,spec2,mz=320)
{
	#print('subtracting spectra....')
	print(head(spec1))
	print(head(spec2))
	print('subtracting spectra....')


	# get spectra for m/z value of interest for spec1
	strname<-paste('X',mz,sep='')
	peaks<-spec1[[strname]]
	retention<-spec1[["RT.minutes."]]
    spectra1<-data.frame(retention,peaks)
	spectra1$name<-rep('Spectra1',nrow(spec1))
    print(head(spectra1))
    
    
    # get spectra for m/z value of interest for spec1
	strname<-paste('X',mz,sep='')
	peaks<-spec2[[strname]]
	retention<-spec2[["RT.minutes."]]
    spectra2<-data.frame(retention,peaks)
	spectra2$name<-rep('Spectra2',nrow(spec1))
    print(head(spectra2))
    
    spectra_all<-rbind(spectra1,spectra2)
    
    g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=retention,col=name))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~name,ncol=1)
    print(g)
   
   spectra_diff<-spectra1
   peaks<-spectra1$peaks - spectra2$peaks
   spectra_diff$peaks<-peaks
   spectra_diff$name<-rep('Difference',nrow(spec1))
   print(spectra_diff)
   
	spectra_all<-rbind(spectra1,spectra2,spectra_diff)

  g<-ggplot(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    geom_point(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    geom_line(data=spectra_all,aes(x=retention,y=peaks,col=name))+
    #geom_text(data=spectra_all,aes(x=retention,y=peaks,label=retention,col=name))+
    mytheme+xlab('Rentention Time (min)')+
    ylab('Abundance')+
    facet_wrap(~name,ncol=1)
    print(g)
   
}

IntegratePeaks<-function(obj,samplename,peaks=c(320,292,218),retentiontime=26.3,lower=0.03,upper=0.09)
  # the purpose of this function is to calculate the area under the curve of a peak
  # make up of peak-57,peak-85,and peak-159 for compounds of interest e.g met
  # to quantify
{
 # print(head(obj))
  names<-c()
  for(n in peaks)
  {
    strname<-paste('X',n,sep='')
    #print(n)
    names<-c(names,strname)
  }
  
  cols<-which(colnames(obj) %in% names)
  print(cols)
  temp<-obj[,c(1,2,cols)]
  print(head(temp))
  min<-as.numeric(retentiontime)-lower
  max<-as.numeric(retentiontime)+upper# this will likely take some tuning depending on the peak
  
  
  # get range of spectra
  t<-temp[which(temp$RT.minutes. >min),]
  t<-temp[which(t$RT.minutes. < max),]
  t<-plotMZforSubSpectrum(temp,samplename,26.3,specname='reference',upper=upper)
  print(t)
  
  a<-auc.mc(x=t$time,y=t$abundance)
  avg<-mean(a)
  sd<-sd(a)
  
  df<-data.frame(avg,sd,samplename)
  colnames(df)<-c('Mean_Area','SD_Area','Sample')
  return(df)
}

plotMZforSubSpectrum<-function(obj,samplename,retentiontime,specname,lower=0.03,upper=0.09,xlimit=200,labels=c('218','292','320'))
{
  
  
  min<-as.numeric(retentiontime)-lower
  max<-as.numeric(retentiontime)+upper# this will likely take some tuning depending on the peak
  cat(retentiontime,min,max,'\n')
  
  # get range of spectra
  t<-obj[which(obj$RT.minutes. >min),]
  t<-t[which(t$RT.minutes. < max),]
  
  
  plot_dat<-data.frame()
  
  print(head(t))
  for(i in 1:nrow(t))
  {
    time<-t[i,]$RT.minutes.
    temp<-t[i,3:ncol(t)]
    print(temp)
    abundance<-sum(temp[1,])
    #print(abundance)
    df<-data.frame(time,abundance)
    colnames(df)<-c('time','abundance')
    plot_dat<-rbind(plot_dat,df)
    
  }	
  #    
  g<-ggplot()+geom_point(data=plot_dat,aes(x=time,y=abundance))+
    geom_line(data=plot_dat,aes(x=time,y=abundance))+mytheme+
    xlab('Time (Minutes) in Column')+ylab('Abundance Units')
  print(g)
  ret<-g
  
  s<-g
  # plot full m/z now
  
  full_spec<-data.frame()
  temp<-t[,3:ncol(t)]
  
  val<-colSums(temp)
  name<-(colnames(temp))
  print(val)
  mz<-data.frame()
  for(i in 1:length(val))
  {
    # print(name[i])
    #masstocharge_name<-substr(as.character(name[i]),start=2,length(name[i]))
    #print(masstocharge_name)	
    #print(length(as.character(name[i])))
    n<-as.character(name[i])
    # print(n)
    #print(class(n))
    n<-(substring(n,2))
    n<-as.integer(n)
    #print(n)
    df<-data.frame(n,val[i])
    colnames(df)<-c('moverz','abundance')
    mz<-rbind(mz,df)
    
  }
  
  #add labels for m/z of interest
  for(i in 1:length(labels))
  {
    lab=labels[i]
    print(lab)
    
    
  }
  
  print(mz[which(mz$moverz %in% labels),])
  
  
  g<-ggplot()+#geom_point(data=mz,aes(x=moverz,y=abundance))+
    geom_segment(data=mz[which(mz$moverz >= xlimit),],aes(x=moverz,xend=moverz,y=0,yend=abundance))+mytheme+
    xlab('m/z')+ylab('Abundance Units')+geom_text(data=mz[which(mz$moverz %in% labels),],aes(x=moverz,y=1.05*abundance,label=moverz))+
    scale_y_continuous(limits=c(0,max(mz$abundance)*1.2))
  s<-g
  print(g)  
  
  
  specname<-rep(specname,nrow(mz))
  #print(specname)
  
  mz<-cbind(mz,specname)
  #print(head(mz))
  
  p <- annotation_custom(grob=g,xmin = 30, xmax = 60, ymin = 60, ymax = 100)
  q<- annotation_custom(grob = s, xmin = 0.5, xmax = 0.9, ymin = 6, ymax = 10)
  
  blanktheme<- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none",
                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())
  
  p = qplot(1:10, 1:10, log='x')+xlab(NULL)+ylab(NULL)+blanktheme
  s<-ggplotGrob(s)
  
  s<-annotation_custom(grob = s, xmin = -0.005, xmax = 1.05, ymin = 0.5, ymax = 10.5)
  
  p<-p+s
  print(p)	
  
  g<-ggplotGrob(ret)
  p<-p + annotation_custom(grob = g, xmin = 0.456, xmax = 1.04 ,ymin = 8, ymax = 10.35)
  print(p)
  
  #ggsave(filename='temp.png',plot=p)
  
  title<-paste(samplename,retentiontime,paste('',min,'-',max,'minutes_subSpectra',sep=''),sep='_')
  
  p<-p+ggtitle(title)+theme(plot.title = element_text(lineheight=2, face="bold"))
  print(p)
  
  path=getwd()
  path=paste(path,'R_Plots/slices',sep='/')
  filename=paste(path,paste(title,'pdf',sep='.'),sep='/')
  #filename=paste(path,'temp.pdf',sep='/')
  print(filename)
  
  ggsave(filename=filename,plot=p)
  
  filename=paste(path,paste(title,'png',sep='.'),sep='/')
  ggsave(filename=filename,plot=p)
  
  # temp expanding of M/z plot for power point. comment out later
  mz<-mz[which(mz$moverz > xlimit),]
  mz<-mz[which(mz$moverz < 400),]
  g<-ggplot()+#geom_point(data=mz,aes(x=moverz,y=abundance))+
    geom_segment(data=mz,aes(x=moverz,xend=moverz,y=0,yend=abundance))+mytheme+
    xlab('m/z')+ylab('Abundance Units')+geom_text(data=mz[which(mz$moverz %in% labels),],aes(x=moverz,y=1.05*abundance,label=moverz))+
    scale_y_continuous(limits=c(0,150000*2))
  
   #print(g)
  #ggsave("tempspectra.png")
  
  
  return(plot_dat)
  
  
  
}


