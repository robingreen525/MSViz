# the purpose of this script is to search GC/MS spectra
# for peaks of a specific m/z ratio. The logic is that specific compounds
# will show specific m/z peaks (such as peak-57,peak-85,and peak -159).
# See http://www.sigmaaldrich.com/technical-documents/articles/reporter-us/the-derivatization.html

# I will write a script that searches through the previously generarated
# 'All_Spectra' data object, where each entry is was a csv file of
# raw m/z ratios at a given rentention time on the column.

# I will extract the points in each sample that shows this charactersitc 
# peak and look up the rentention time. I will then investigate
# these retention time peaks against the NIST library manually.

rm(list=ls())

require(ggplot2)
require(png)
require(flux)
# directory where my data are
setwd('X:/fast/shou_w/GC-MS/UW Chemistry/Agilent5975/20160521/Analysis')


#load RData file with my compiled spectra
load('20160523_SpectraAll.RData')

#call functions file
source('extractPeaks_FUNCTIONS.R')




#read list of compounds of interest
r<-read.csv('compoundsofinterest_metONLY.csv',header=T)

#sink('20160503_SpectralAnalysisLOG.txt')
search_list<-c('MS_45','MS_46','MS_47','MS_48','MS_49','MS_50','MS_51','MS_52',
               'MS_53','MS_54','MS_55','MS_56')



for(i  in 1:nrow(r))
{
  print('***************************')
  obj<-NULL
  obj<-r[i,]

    
  #extractPeak(obj,search_list)
  print('***************************')
  

  
}


# 
# obj<-acquireSpectra('MS_45')
# ref<-plotMZforSampleRT(obj,'MS_45',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_45')
# 
# obj<-acquireSpectra('MS_46')
# ref<-plotMZforSampleRT(obj,'MS_46',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_46')
# 
# obj<-acquireSpectra('MS_47')
# ref<-plotMZforSampleRT(obj,'MS_47',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_47')
# 
# obj<-acquireSpectra('MS_48')
# ref<-plotMZforSampleRT(obj,'MS_48',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_48')
# 
# obj<-acquireSpectra('MS_49')
# ref<-plotMZforSampleRT(obj,'MS_49',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_49')
# 
# obj<-acquireSpectra('MS_50')
# ref<-plotMZforSampleRT(obj,'MS_50',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_50')
# 
# obj<-acquireSpectra('MS_51')
# ref<-plotMZforSampleRT(obj,'MS_51',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_51')
# 
# obj<-acquireSpectra('MS_52')
# ref<-plotMZforSampleRT(obj,'MS_52',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_52')
# 
# obj<-acquireSpectra('MS_53')
# ref<-plotMZforSampleRT(obj,'MS_53',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_53')
# 
# obj<-acquireSpectra('MS_54')
# ref<-plotMZforSampleRT(obj,'MS_54',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_54')
# 
# obj<-acquireSpectra('MS_55')
# ref<-plotMZforSampleRT(obj,'MS_55',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_55')
# 
# obj<-acquireSpectra('MS_56')
# ref<-plotMZforSampleRT(obj,'MS_56',26.3,specname='reference',xlimit=200)
# plotFullSpectra(obj,'MS_56')

areas<-data.frame()

obj<-acquireSpectra('MS_53')
area<-IntegratePeaks(obj,'MS_53',peaks=c(320,292,218),upper=0.065)

search_list<-c('MS_45','MS_46','MS_47','MS_48','MS_49','MS_50','MS_51','MS_52',
               'MS_53','MS_54','MS_55','MS_56')
for(sample in search_list)
{
  obj<-acquireSpectra(sample)
  area<-IntegratePeaks(obj,sample,peaks=c(320,292,218),upper=0.065)
  areas<-rbind(areas,area)
  
}

MetConc<-t(data.frame(10,10,3,3,1,1,0.3,0.3,10,10,1,1))
colnames(MetConc)<-c('MetConc')
areas<-cbind(areas,MetConc)

Media<-t(data.frame('1/10x SD','1/10x SD','1/10x SD','1/10x SD','1/10x SD','1/10x SD','1/10x SD','1/10x SD',
                    'Water','Water','Water','Water'))
colnames(Media)<-c('Media')
areas<-cbind(areas,Media)

g<-ggplot()+
  geom_pointrange(data=areas,aes(x=MetConc,y=Mean_Area,col=Media,
                                 ymax=Mean_Area+2*SD_Area,
                                 ymin=Mean_Area-2*SD_Area),size=1.5,position = position_jitter())+
  mytheme+xlab('uM Met of Stock')+ylab('Integrated Area (Trapizodial Rule)')
ggsave(g,file='20160523_Areas.pdf')
ggsave(g,file='20160523_Areas.png')
print(g)