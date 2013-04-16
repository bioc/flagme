gatherInfo<-function(pD,obj,newind=NULL,method=c("apex"),findmzind=TRUE,useTIC=FALSE,top=NULL,intensity.cut=.05) {
  # pind - indices of peaks
  # rind - indices of runs
  # newind ... comes from an imputation step, if at all.
  method<-match.arg(method)
  if (is(obj,"multipleAlignment")) {
    pind<-obj@betweenAlignment@ind
    rind<-obj@betweenAlignment@runs
    groups<-obj@betweenAlignment@groups
  } else if (is(obj,"progressiveAlignment")) {
    #n<-length(obj@merges[[1]])
    #	cat(n,"\n")
    pind<-obj@merges[[1]]$ind   # for example, if obj=ma$pas[1]
    rind<-obj@merges[[1]]$runs
    groups<-rep(names(obj),length(rind))
  }
  out<-matrix(NA,nrow(pind),ncol(pind))
  for(j in 1:ncol(pind)) {
    data<-pD@peaksrt[[rind[j]]]
    keep<-which(!is.na(pind[,j]))
	out[keep,j]<-data[pind[keep,j]]
  }
  colnames(out)<-paste(groups,rind,sep=".")
  rownames(out)<-1:nrow(out)
  peaks<-vector("list",nrow(pind))
  for(i in 1:length(peaks)) {
    if (findmzind) {
	  mzind<-rep(0,length(pD@mz))
      for(j in 1:ncol(pind)) {
        data<-pD@peaksdata[[rind[j]]]
	    if(!is.na(pind[i,j]))
	      mzind<-mzind+data[,pind[i,j]]
	  }
	  #mzind<-which(mzind>0)
	  mzind<-which(mzind>intensity.cut*max(mzind))
	} else {
	  mzind<-1:length(pD@mz)
	}
    if (useTIC)
      this.d<-matrix(,nrow=1,ncol=ncol(pind))
    else
      this.d<-matrix(,nrow=length(mzind),ncol=ncol(pind))
		
	for(j in 1:ncol(pind)) {
	  r<-rind[j]
  	  overlap<-NULL
	  if(!is.na(pind[i,j]))
		if(method=="apex") {
		  if (useTIC)
	        this.d[1,j]<-sum(pD@rawdata[[r]][,pD@peaksind[[r]][pind[i,j]]]) ##
		  else
	        this.d[,j]<-pD@rawdata[[r]][mzind,pD@peaksind[[r]][pind[i,j]]]
		  } else { # method="area"
		    #nn<-length(pD@peaksind.start[[r]])
			#w<-which(round(pD@peaksind.start[[r]][-1]) < round(pD@peaksind.end[[r]][-nn]))
			#overlap<-cbind(w,w+1)
			#v<-which(pind[i,j]==overlap,arr.ind=TRUE)
			#if(nrow(v)>0) {
			  # do mixture thing. 
		      #pks<-overlap[v[1,1],]
			  #cat("overlap -- run",r,"-- peaks",pks,"--",v[1,2],"\n")
			  #if(r==12&pks[1]==4)
			  #  return(list(mzind=mzind))
			  #this.d[,j]<-partitionWithNormalMixture2(pD,r,pks[1],pks[2],mzind,useTIC=useTIC)[,v[1,2]]
		      #if(j==ncol(pind)) cat("\n")
			#} else {
			  # just simple sum of intensity
 		      #st<-pD@peaksind.start[[r]][pind[i,j]]
		      #en<-pD@peaksind.end[[r]][pind[i,j]]
		      #if (useTIC)
	          #  this.d[1,j]<-sum(pD@rawdata[[r]][,st:en])  ##
			  #else
	          #  this.d[,j]<-rowSums(pD@rawdata[[r]][mzind,st:en])
			#}
		  }
	    if(!is.null(newind))
		  if (!is.na(newind$apex[i,j])) {
  		    if(method=="apex") {
			  if (useTIC)
	            this.d[1,j]<-sum(pD@rawdata[[r]][,newind$apex[i,j]])  ##
			  else
	            this.d[,j]<-pD@rawdata[[r]][mzind,newind$apex[i,j]]
		      out[i,j]<-pD@rawrt[[r]][newind$apex[i,j]]
		    } else { # method="area"
		      #st<-newind$start[i,j]
		      #en<-newind$end[i,j]
			  #cat(r,i,j,st,en,"\n")
			  #if (useTIC)
	          #  this.d[1,j]<-sum(pD@rawdata[[r]][,st:en])  ##
			  #else
	          #  this.d[,j]<-rowSums(pD@rawdata[[r]][mzind,st:en])
			}
		  }
	  }
	  colnames(this.d)<-paste(groups,rind,sep=".")
	  peaks[[i]]<-list(rt=out[i,],mz=pD@mz[mzind],data=this.d)
	}
	if( !is.null(top) ) {
	  for(i in seq(along=peaks)) {
	    o<-sort(order(-rowSums(peaks[[i]]@data,na.rm=TRUE))[1:top])
	    peaks[[i]]@mz<-peaks[[i]]@mz[o]
	    peaks[[i]]@data<-peaks[[i]]@data[o,]
	  }
	}
	return(peaks)
}
