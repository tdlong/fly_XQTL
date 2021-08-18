cutter = function(x,NC){
	cut=x[2*NC+1]
	strand=x[2*NC+2]
	if(strand){
		out=c(x[(1):(cut)],x[(NC+cut+1):(2*NC)])
		}else{
		out=c(x[(NC+1):(NC+cut)],x[(cut+1):(NC)])
		}
	out
	}

bang_out_a_replicate = function(TASK,myseed){
myseed=as.numeric(as.character(myseed))
TASK=as.numeric(as.character(TASK))
set.seed(myseed)
library(tidyverse)
libby = read.table("helperfile/libby.txt")
RILnames = unique(libby$RIL)
ALLPOS = unique(libby$pos)
SS = c(200,400,600,800)
SampleSize=SS[(TASK %% 4) + 1]    # Number of DSPR RILs
RRR = (TASK %/% 4) + 1
NN=20000
#  Pure Reps * SampleSize * Intensity * NNind * reps of experiment (3) * number of positions
#  colnames(MM) = c("Pure_Rep_number", "N_RILs_base", "Selection_Intensity", "Ne_Fly_DNA_pool", "Pos_On_3L", "N_exp_reps", "LOD", "nLOD")
MM = matrix(nrow=3*3*3*length(ALLPOS)*3,ncol=9)
cc=1
for(intensity in c(0.05,0.10,0.20)){
	for(NNind in c(150,300,600)){
	for(NgenRandomMate in c(4,8,16)){

# debug		
#SampleSize=200
#intensity=0.05
#NNind=150
		
		cat(SampleSize,"\t",intensity,"\t",NNind,"\n")
		g = rnorm(8) 
		RILsample = libby[libby$RIL %in% sample(RILnames,SampleSize,replace=FALSE),]
		temp = RILsample %>% pivot_wider(names_from=pos,values_from=founder)
		MMh = temp[,3:ncol(temp)]
		MMe = data.frame(matrix(rnorm(ncol(MMh)*nrow(MMh)),nrow=nrow(MMh)))
		T3D = data.frame(hap=rep(NA,length(ALLPOS)*8*2*12),trt=rep(NA,length(ALLPOS)*8*2*12), RR=rep(NA,length(ALLPOS)*8*2*12),pos=rep(NA,length(ALLPOS)*8*2*12),afreq=rep(NA,length(ALLPOS)*8*2*12))
		T3D2 = data.frame(hap=rep(NA,length(ALLPOS)*8*2*12),trt=rep(NA,length(ALLPOS)*8*2*12), RR=rep(NA,length(ALLPOS)*8*2*12),pos=rep(NA,length(ALLPOS)*8*2*12),afreq=rep(NA,length(ALLPOS)*8*2*12))
		start = 1
		NC = ncol(MMh)
		for(mrep in 1:12){
			# I now blow the population up from SampleSize to 20000
			popsize = SampleSize*2.6			
			S1 = sample(1:SampleSize,popsize,replace=TRUE)
			S2 = sample(1:SampleSize,popsize,replace=TRUE)
			S3 = sample(1:SampleSize,popsize,replace=TRUE)
			MMh_old_female = cbind(MMh[S1,],MMh[S2,])
			MMe_old_female = cbind(MMe[S1,],MMe[S2,])
			MMh_old_male = MMh[S3,]
			MMe_old_male = MMe[S3,]
			for(i in 1:NgenRandomMate){
				cut = sample(2:(NC-1),popsize,replace=TRUE)
				strand = sample(c(TRUE,FALSE),popsize,replace=TRUE)
				MMh_old_female$cut = cut
				MMe_old_female$cut = cut
				MMh_old_female$strand = strand
				MMe_old_female$strand = strand
				MMh_female = t(apply(MMh_old_female,1,function(x) cutter(as.numeric(x),NC)))
				MMe_female = t(apply(MMe_old_female,1,function(x) cutter(as.numeric(x),NC)))
				MMh_male = MMh_old_male
				MMe_male = MMe_old_male
					
				popsize = min(floor(popsize*2.6),NN)
				S1 = sample(1:nrow(MMh_female),popsize,replace=TRUE)
				S2 = sample(1:nrow(MMh_female),popsize,replace=TRUE)
				S3 = sample(1:nrow(MMh_female),popsize,replace=TRUE)
				MMh_old_female = data.frame(cbind(MMh_female[S1,],MMh_male[S2,]))
				MMe_old_female = data.frame(cbind(MMe_female[S1,],MMe_male[S2,]))
				MMh_old_male = data.frame(MMh_female[S3,])
				MMe_old_male = data.frame(MMe_female[S3,])
				}
			# position 60 is causative
			D = data.frame(BG = apply(MMe_old_female[,1:NC],1,sum) + apply(MMe_old_female[,(NC+1):(2*NC)],1,sum))
			D$BG = (sqrt(8) * (D$BG - mean(D$BG))/sd(D$BG))                 # var of 8
			D$G = g[MMh_old_female[,60]] + g[MMh_old_female[,60+NC]]
			D$G =  (sqrt(2) * (D$G - mean(D$G))/sd(D$G))                 		# var of 2
			D$P = D$G + D$BG + rnorm(popsize,0,sqrt(9*(var(D$G) + var(D$BG))))		# since heritability = 50, but other chromosomes are noise

			oo = order(D$P,decreasing=TRUE)[1:(popsize*intensity)]
			Roo = sample(1:popsize,NNind)
			Roo2 = sample(1:popsize,NNind)
			Soo = sample(oo,NNind)
			cccc = 1
			for (ppp in ALLPOS){
				f_full      = (table(factor(MMh_old_female[Soo ,cccc],levels=1:8)) + table(factor(MMh_old_female[Soo, NC+cccc],levels=1:8)))/(2 * NNind)
				f_sample    = (table(factor(MMh_old_female[Roo ,cccc],levels=1:8)) + table(factor(MMh_old_female[Roo, NC+cccc],levels=1:8)))/(2 * NNind)
				f_sample_c  = (table(factor(MMh_old_female[Roo2,cccc],levels=1:8)) + table(factor(MMh_old_female[Roo2,NC+cccc],levels=1:8)))/(2 * NNind)
				TC=cbind(1:8, rep("C",8), rep(mrep,8), rep(ppp,8), asin(sqrt(f_full)))
				TS=cbind(1:8, rep("S",8), rep(mrep,8), rep(ppp,8), asin(sqrt(f_sample)))
				TS2=cbind(1:8, rep("S2",8), rep(mrep,8), rep(ppp,8), asin(sqrt(f_sample_c)))
				T3D[start:(start+15),] = rbind(TC,TS)	
				T3D2[start:(start+15),] = rbind(TS,TS2)	
				start = start + 16
				cccc=cccc+1
				}  # positions
			} # mrep			
		T3D$hap = as.factor(T3D$hap)
		T3D$trt = as.factor(T3D$trt)
		T3D$pos = as.numeric(as.character(T3D$pos))
		T3D$RR = as.numeric(as.character(T3D$RR))
		T3D$afreq = as.numeric(as.character(T3D$afreq))
		T3D2$hap = as.factor(T3D2$hap)
		T3D2$trt = as.factor(T3D2$trt)
		T3D2$pos = as.numeric(as.character(T3D2$pos))
		T3D2$RR = as.numeric(as.character(T3D2$RR))
		T3D2$afreq = as.numeric(as.character(T3D2$afreq))
		for(ppp in ALLPOS){
			for(mrep in c(4,8,12)){
				D=T3D[T3D$pos==ppp & T3D$RR <= mrep,]
				out = anova(lm(afreq~trt+hap+hap:trt,data=D))
				LODfq = -pf(out[3,3]/out[4,3],out[3,1],out[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)	
				D2=T3D2[T3D2$pos==ppp & T3D2$RR <= mrep,]
				out2 = anova(lm(afreq~trt+hap+hap:trt,data=D2))
				LODfq2 = -pf(out2[3,3]/out2[4,3],out[3,1],out2[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)	
				MM[cc,] = c(RRR, SampleSize, intensity, NNind, NgenRandomMate, ppp, mrep, LODfq, LODfq2)
				cc = cc+1
				} # mrep				
			} # positions				
		}} # NNind, #NgenRandomMate		
	}  # intensity
MM = data.frame(MM)
# Ne_Fly_DNA_pool is an effective number that includes errors due to DNA prep sampling (like size of flies, etc.)
#      plus the uncertainty in the haplotype estimator, this is generally unknown
colnames(MM) = c("Pure_Rep_number", "N_RILs_base", "Selection_Intensity", "Ne_Fly_DNA_pool", "N_Gen_Randon_Mate", "Pos_On_3L", "N_exp_reps", "LOD", "nLOD")
fileout=paste("powerdir2/powersim.",myseed,".txt",sep="")
write.table(MM,fileout,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

