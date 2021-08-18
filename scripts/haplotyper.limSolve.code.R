blah = function(par,tw){
	ww = exp(-unlist(tw)/(2*par^2))
	ww = ww/sum(ww)
	NumberSites = 500
	portionOfWeight = 0.50
	# the idea is to choose sigma, such that the "NumberSites" closest sites to pos
	# account for "portionOfWeight" of the weight.  In yeast I used 50:0.50
	# the trade of is bigger and bigger windows cut down resolution so this should sort of scale
	# with haplotype block size
	# maybe in flies I can use 100:0.50...
	(sum(sort(ww,decreasing=TRUE)[NumberSites:length(ww)])-portionOfWeight)^2
	}
 
runscan <- function(Cpool, Tpool, OutName, folder, SNPtable, foundernames){

	#  Cpool, Tpool = Control and Treatment pools to compare
	#  OutName = prefix for output file
	#  folder = output folder for results
	#  SNPtable = file of SNP frequencies in pool and founders.  The structure of the file matters
	#     the 1st column is a chromosome name, and the 2nd a position, there can be additional columns
	#     eventually you want paired columns of the 1st is the frequency in some sample and the 2nd the coverage
	#     I have a python script that makes these files from "vcf" files, vcf are a very standard format
	#     this file should be preprocessed with GrannySNP.R
	#  foundernames = file with the names of the founder haplotypes (and a short name) 

	poolroot=Cpool
	yy = read.table(file=foundernames,header=FALSE)
	founderroots = as.character(yy[,1])
	foundershortnames = as.character(yy[,2])
	filename=paste(folder,"/",OutName,"_hap_freq.txt",sep='')
	
	if (!file.exists(filename)){
		xx = read.table(file=SNPtable)
		NCpool <-  match(paste("freq_",Cpool,sep=''),names(xx))
		NTpool <-  match(paste("freq_",Tpool,sep=''),names(xx))
		founders <- match(paste("freq_",founderroots,sep=''),names(xx))
		Nf = length(founders)
		stepSize = 10000
		WINSIZE = 100000
		names = paste(foundershortnames, collapse=";")
		lppp = 0
		for(chr in c("chrX","chr2L","chr2R","chr3L","chr3R","chr4")){
			Maxpos = max(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			Minpos = min(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			ppp = seq(Minpos+WINSIZE,Maxpos-WINSIZE,stepSize)
			lppp = lppp + length(ppp)
			}
		ddd = data.frame(poolroot=rep("",lppp),chr=rep("",lppp),pos=rep(0,lppp),NSNPs=rep(0,lppp),foundernames=rep("",lppp), cutree=rep("",lppp),numberlevs=rep(0,lppp),Cfounderfreqs=rep("",lppp),Cadjfounderfreqs=rep("",lppp),Tfounderfreqs=rep("",lppp),Tadjfounderfreqs=rep("",lppp),diffByGroup=rep("",lppp),SS=rep(0,lppp),CC=rep(0,lppp),lpF=rep(0,lppp),RSSf=rep(0,lppp),RSSr=rep(0,lppp),df1=rep(0,lppp),df2=rep(0,lppp),stringsAsFactors=FALSE)
		i=1
		for(chr in c("chrX","chr2L","chr2R","chr3L","chr3R","chr4")){
			Maxpos = max(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			Minpos = min(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			ppp = seq(Minpos+WINSIZE,Maxpos-WINSIZE,stepSize)
			for(pos in ppp){
				#pos=9700
				#pos=50700
				# median non-recombined block size is 14kb
				#  so step size should be much smaller than that (500bp = 1/30th)
				#  I would guess optimum window size is +/- 7kb (= 14 total) ish.  
				#  but markers are limiting...
				MynumberSNP<- NA
				Groups = rep(NA,Nf)
				groups = NA
				nGroups = NA
				CFreqs = rep(NA,Nf)
				CGroupAvFreqs = rep(NA,Nf)
				TFreqs = rep(NA,Nf)
				TGroupAvFreqs = rep(NA,Nf)
				CC = NA
				lpF = NA
				RSSf = NA
				RSSr = NA
				df1 = NA
				df2 = NA

				if(sum(xx$CHROM == chr & xx$POS > pos - WINSIZE & xx$POS < pos + WINSIZE)>100){
					predictors <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=founders)
					CY <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=NCpool)
					TY <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=NTpool)
					tw = (subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=POS)-pos)^2
					MynumberSNP <- nrow(predictors)
					# pick a sigma, near the ends of chromosomes sigma has to be bigger than in the middle
					predictNotMissing = apply(predictors,1,function(x) sum(is.na(x))==0) & (!is.na(CY)) & (!is.na(TY))
					predictors = predictors[predictNotMissing,]
					CY = CY[predictNotMissing,]
					TY = TY[predictNotMissing,]
					tw = tw[predictNotMissing,]
					out <- optimize(blah,c(20000,100000),tw)
					sigma=out$minimum
					weights = exp(-tw/(2*sigma^2))

					Groups <- cutree(hclust(dist(t(predictors*unlist(weights)))),h=1)
					nGroups = nlevels(as.factor(Groups))
					groups = paste(as.numeric(Groups),collapse=";")
					d = ncol(predictors)
					A = predictors[weights>0.01,]
					E = t(matrix(rep(1,d)))
					F = 1
					G = diag(rep(1,d))
					H = matrix(rep(0.0003,d))
					if(nrow(A)>=100){
						outC = lsei(A=A,B=CY[weights>0.01],E=E,F=F,G=G,H=H,verbose=TRUE)
						outT = lsei(A=A,B=TY[weights>0.01],E=E,F=F,G=G,H=H,verbose=TRUE)
						outR = lsei(A=rbind(A,A),B=c(CY[weights>0.01],TY[weights>0.01]),E=E,F=F,G=G,H=H,verbose=TRUE)
						if(!outC$IsError & !outT$IsError & !outR$IsError){
							CFreqs <- outC$X
							Cmeanbygroup <- tapply(CFreqs,as.numeric(Groups),mean)
							CGroupAvFreqs <- Cmeanbygroup[as.numeric(Groups)]      # freqs replaced by group means
							TFreqs <- outT$X
							Tmeanbygroup <- tapply(TFreqs,as.numeric(Groups),mean)
							TGroupAvFreqs <- Tmeanbygroup[as.numeric(Groups)]      # freqs replaced by group means
							RFreqs <- outR$X
							Rmeanbygroup <- tapply(RFreqs,as.numeric(Groups),mean)
							RGroupAvFreqs <- Rmeanbygroup[as.numeric(Groups)]
							# mm = (TGroupAvFreqs + CGroupAvFreqs)/2
							RSSCFree = sum((CY[weights>0.01] - as.matrix(A) %*% CGroupAvFreqs)^2)
							# RSSCRest = sum((CY[weights>0.01] - as.matrix(A) %*% mm)^2)
							RSSTFree = sum((TY[weights>0.01] - as.matrix(A) %*% TGroupAvFreqs)^2)
							# RSSTRest = sum((TY[weights>0.01] - as.matrix(A) %*% mm)^2)
							RSSr = sum((c(CY[weights>0.01],TY[weights>0.01]) - as.matrix(rbind(A,A)) %*% RGroupAvFreqs)^2)
							RSSf = RSSCFree + RSSTFree
							df1 = ncol(A)
							df2 = 2*nrow(A)-2*ncol(A)
							CC = RSSr/RSSf
							F = ((RSSr - RSSf)*df2)/(RSSf*df1)
							lpF = -(1/log(10))*pf(F,df1,df2,log.p=TRUE,lower.tail=FALSE)
							} # model could be fitted
						} # nrow >= 100
				} # >100 raw SNPs present
				Cfreqs = paste(round(as.numeric(CFreqs),4),collapse=";")
				Cafreqs = paste(round(as.numeric(CGroupAvFreqs),4),collapse=";")
				Tfreqs = paste(round(as.numeric(TFreqs),4),collapse=";")
				Tafreqs = paste(round(as.numeric(TGroupAvFreqs),4),collapse=";")
				diffByGroup = paste(round(as.numeric(TGroupAvFreqs - CGroupAvFreqs),4),collapse=";")
				SS = sum((TGroupAvFreqs - CGroupAvFreqs)^2)
				ddd[i,] = c(poolroot, chr, pos, MynumberSNP, names, groups, nGroups, Cfreqs, Cafreqs, Tfreqs, Tafreqs, diffByGroup, SS, CC, lpF, RSSf, RSSr, df1, df2)
				i=i+1
				}  # over positions
			}  # chr
		write.table(ddd,filename,row.names=FALSE,sep="\t",quote=FALSE)
		} # file does not exist
	} # function
			
			
			
			
			
			
			
			

