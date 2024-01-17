# Run on delta.
# Generate a synteny plot for a list of species from TOGA results.

#===========Usage.
#source("toga2synteny.r")

## Parameters.
#spes<-c("HLcarPer2","HLgloSor2","HLphyHas1","HLphyDis3","HLdesRot8A","HLdiaYou1","HLdipEca1","HLmicMeg1","HLmyoMyo6","HLrhiFer5")
#ref<-"HLdesRot8A"
#TOGA_dir<-"/projects/hillerlab/genome/gbdb-HL/hg38/TOGA"
#minChrLen<-5e6
#palette(palette.colors(palette="R4",alpha=0.5))

## Get synteny.
#mySynteny<-toga2synteny(spes,ref,TOGA_dir,minChrLen)

## Plot synteny.
#plot.toga2synteny(mySynteny,ref,arrange=FALSE)

## Plot synteny by rearranging the chromosomes.
#plot.toga2synteny(mySynteny,ref,arrange=TRUE)

#===========Define synteny blocks for a list of species based on orthologous gene alignments.
# By Shenglin Liu, on Aug 27, Sat, 2022.

# Caution: the algorithm overlooks fine-scale (<= thr.shft) disruption of collinearity (normally due to inversions);
#	it also overlooks strand of the genes, assuming the strandedness fits with collinearity.

#-----------Utility functions for plotting.
convergence<-function(b1,b2,chrs1)
{
	# NA values are forbiden in b1 and b2.
	sapply(chrs1,function(chr1)
	{
		index<-b1[,1]==chr1
		if(sum(index)==0)return(NA)
		w<-b1[index,4]
		p<-sapply(split(w,b2[index,1]),sum)/sum(w)
		-sum(p*log(p))
	})
}
containment<-function(b1,b2,chrs1)
{
	# NA values are forbiden in b1 and b2.
	sapply(chrs1,function(chr1)
	{
		index1<-b1[,1]==chr1
		if(sum(index1)==0)return(NA)
		chrs2<-unique(b2[index1,1])
		mean(sapply(chrs2,function(chr2)
		{
			index2<-b2[,1]==chr2
			w<-b2[index2,4]
			w.joint<-b2[index1&index2,4]
			p<-sum(w.joint)/sum(w)
			if(p<1)p<-c(p,1-p)
			-sum(p*log(p))
		}))
	})
}
convert.coord<-function(chr,pos,chrs,lengths)
{
	chr<-as.character(chr)
	chrs<-as.character(chrs)
	
	pos<-as.numeric(pos)
	lengths<-as.numeric(lengths)
	
	n.chrs<-length(chrs)
	increments<-c(0,cumsum(lengths[-n.chrs]))
	names(increments)<-chrs
	
	pos+increments[chr]
}
reverse.coord<-function(a,chrs,lengths)
{
	names(lengths)<-chrs
	index<-which(a[,1]%in%chrs)
	if(length(index)>0)
	{
		a[index,2]<-lengths[as.character(a[index,1])]-a[index,2]+1
		a[index,3]<-lengths[as.character(a[index,1])]-a[index,3]+1
	}
	a
}
move.chrs<-function(b1,b2,chrs1,lengths1,chrs2,lengths2)
{
	# NA values are forbiden in b1 and b2.
	x0<-convert.coord(b1[,1],b1[,2],chrs1,lengths1)
	loc<-sapply(chrs2,function(chr2)
	{
		index<-b2[,1]==chr2
		w<-b2[index,4]
		sum(x0[index]*w)/sum(w)
	})
	index<-order(loc,na.last=TRUE)
	list(chrs=chrs2[index],lengths=lengths2[index])
}
flip.chrs<-function(b1,b2,chrs1,lengths1,chrs2,lengths2)
{
	# NA values are forbiden in b1 and b2.
	orientation<-sapply(chrs2,function(chr2)
	{
		index<-b2[,1]==chr2
		w<-b2[index,4]
		same<-(b1[index,3]>b1[index,2])==(b2[index,3]>b2[index,2])
		sum(same*w)/sum(w)
	})
	index<-which(orientation<0.5)
	list(chrs.rev=chrs2[index],lengths.rev=lengths2[index])
}

#-----------Core functions.
# Function for getting synteny blocks.
get.blocks<-function(a.list,ref,thr.shft=100,min.nGen=5)
{
	# a.list: a named list of data.frames, each with 3 columns recording gene coordinates for a species;
	#	the data.frames MUST be aligned;
	#	the data.frames MUST NOT contain NA or missing values;
	#	for a gene in a genome, if Column 3 is larger than Column 2, it is on + strand, otherwise on - strand.
	# ref: name of the reference genome; must be contained in the names of a.list.
	# thr.shft: When the neighboring gene ranks on a query have difference larger than this, a block boundary is defined;
	#	defaults to 100.
	# min.nGen: Blocks with gene number lower than this will be discarded; this value must be lower than thr.shft;
	#	defaults to 5.
	# Returns a list of data.frames with the same length as a.list;
	#	each data.frame has 4 columns, defining the coordinates and gene counts of synteny blocks on the corresponding genome.

	spes<-names(a.list)
	n<-length(spes)
	i.ref<-which(spes==ref)

	# Order the genes according to reference.
	a1<-a.list[[i.ref]]
	index<-order(a1[,1],a1[,2])	#,a1[,3]
	a.list<-lapply(a.list,function(x)x[index,])

	# Generate blocks.
	ENs<-lapply(a.list[-i.ref],function(a2){
		ordr<-order(order(a2[,1],a2[,2]))	#,a2[,3]	# Get the gene ranks on each query.
		shft<-c(ordr[-1],-thr.shft)-ordr
		which(abs(shft)>thr.shft)
	})
	EN<-sort(unique(unlist(ENs)))
	ST<-c(1,EN[-length(EN)]+1)

	# Remove genes in small blocks and remake blocks.
	# Particularly useful for cases where one or a few genes jumped to a location far away (side effect?).
	index<-which((EN-ST+1)<min.nGen)
	if(length(index)>0)
	{
		# Remove genes in small blocks.
		index<-unlist(lapply(index,function(x)ST[x]:EN[x]))
		a.list<-lapply(a.list,function(x)x[-index,])

		# Remake blocks.
		ENs<-lapply(a.list[-i.ref],function(a2){
			ordr<-order(order(a2[,1],a2[,2]))	#,a2[,3]
			shft<-c(ordr[-1],-thr.shft)-ordr
			which(abs(shft)>thr.shft)
		})
		EN<-sort(unique(unlist(ENs)))
	}

	# Chromosome change.
	index<-lapply(a.list,function(x)which(c(x[-1,1],"")!=x[,1]))
	EN<-sort(unique(c(EN,unlist(index))))
	ST<-c(1,EN[-length(EN)]+1)

	# Remove small blocks.
	index<-which((EN-ST+1)<min.nGen)
	if(length(index)>0)
	{
		ST<-ST[-index]
		EN<-EN[-index]
	}

	# Get coordinates for the blocks.
	get.coords<-function(a2,ST,EN)
	{
		c2<-character();s2<-numeric();e2<-numeric();n2<-integer()
		for(i in 1:length(ST))
		{
			c2<-c(c2,a2[ST[i],1])
			if(a2[EN[i],2]>a2[ST[i],2])
			{
				s2<-c(s2,min(a2[ST[i],-1]))
				e2<-c(e2,max(a2[EN[i],-1]))
			}else
			{
				s2<-c(s2,max(a2[ST[i],-1]))
				e2<-c(e2,min(a2[EN[i],-1]))
			}
			n2<-c(n2,EN[i]-ST[i]+1)
		}
		data.frame(chr=c2,start=s2,end=e2,wt=n2)
	}
	lapply(a.list,function(x)get.coords(x,ST,EN))
}
# Function for plotting synteny blocks.
plot.blocks<-function(b.list,c.list,l.list,ref,chrGap=5e6,chrWid=0.2,arrange=FALSE)
{
	# b.list, c.list, l.list must be named.
	# Names of b.list must be contained by names of c.list and l.list.
	
	spes<-names(b.list)
	n<-length(spes)
	i.ref<-which(spes==ref)
	c.list<-c.list[spes]
	l.list<-l.list[spes]
	chrExp<-chrWid/2
	original.ref.chrs<-c.list[[i.ref]]
	c.col.list<-lapply(c.list,function(chrs)
	{
		cols<-rep("dimgrey",length(chrs))
		names(cols)<-chrs
		cols
	})

	# Rearrange chromosomes to maximize visual interpretability; flip if necessary.
	if(arrange)
	{
		# Reorder the chromosomes of the reference species by containment.
		b1<-b.list[[i.ref]]
		chrs1<-c.list[[i.ref]]
		index<-intersect(1:n,c(-1,1)+i.ref)
		cnt.ref<-rowMeans(sapply(b.list[index],function(b2)containment(b1,b2,chrs1)))
		index<-order(cnt.ref,na.last=TRUE)
		c.list[[i.ref]]<-c.list[[i.ref]][index]
		l.list[[i.ref]]<-l.list[[i.ref]][index]

		# Reorder and flip the query chromosomes.
		if(i.ref>1)
		for(i in (i.ref-1):1)
		{
			b1<-b.list[[i+1]]
			chrs1<-c.list[[i+1]]
			lengths1<-l.list[[i+1]]
			b2<-b.list[[i]]
			chrs2<-c.list[[i]]
			lengths2<-l.list[[i]]
			# Reorder the chromosomes of the query species.
			tmp<-move.chrs(b1,b2,chrs1,lengths1,chrs2,lengths2)
			c.list[[i]]<-tmp$chrs
			l.list[[i]]<-tmp$lengths
			# Get chromosomes to flip in the query species and flip them.
			tmp<-flip.chrs(b1,b2,chrs1,lengths1,chrs2,lengths2)
			if(length(tmp$chrs.rev)>0)
			{
				b.list[[i]]<-reverse.coord(b2,tmp$chrs.rev,tmp$lengths.rev)
				c.col.list[[i]][tmp$chrs.rev]<-"grey"
			}
		}
		if(i.ref<n)
		for(i in (i.ref+1):n)
		{
			b1<-b.list[[i-1]]
			chrs1<-c.list[[i-1]]
			lengths1<-l.list[[i-1]]
			b2<-b.list[[i]]
			chrs2<-c.list[[i]]
			lengths2<-l.list[[i]]
			# Reorder the chromosomes of the query species.
			tmp<-move.chrs(b1,b2,chrs1,lengths1,chrs2,lengths2)
			c.list[[i]]<-tmp$chrs
			l.list[[i]]<-tmp$lengths
			# Get chromosomes to flip in the query species and flip them.
			tmp<-flip.chrs(b1,b2,chrs1,lengths1,chrs2,lengths2)
			if(length(tmp$chrs.rev)>0)
			{
				b.list[[i]]<-reverse.coord(b2,tmp$chrs.rev,tmp$lengths.rev)
				c.col.list[[i]][tmp$chrs.rev]<-"grey"
			}
		}
	}

	# Make gaps between chromosomes.
	c.list<-lapply(c.list,function(chrs)as.vector(rbind(chrs,NA)))
	l.list<-lapply(l.list,function(lengths)as.vector(rbind(lengths,chrGap)))
	
	# Synteny polygons.
	x1<-x2<-numeric()
	for(spe in spes)
	{
		blks<-b.list[[spe]]
		chrs<-c.list[[spe]]
		lengths<-l.list[[spe]]
		xx1<-with(blks,convert.coord(chr,start,chrs,lengths))
		xx2<-with(blks,convert.coord(chr,end,chrs,lengths))
		x1<-rbind(x1,xx1,xx1)
		x2<-rbind(x2,xx2,xx2)
	}
	x<-as.numeric(rbind(x1,x2[nrow(x2):1,],NA))
	yy<-as.numeric(rbind(1:n-chrExp,1:n+chrExp))
	y<-c(yy,rev(yy),NA)
	xy<-cbind(x,y)

	# Plot.
	plot(1,type="n",xlim=c(1,max(sapply(l.list,sum))),ylim=c(n+chrExp,1-chrExp),xlab=NA,ylab=NA,axes=F,xaxs="i",yaxs="i")
	polygon(xy,col=factor(b.list[[i.ref]]$chr,original.ref.chrs),border=NA)	#,na.omit(c.list[[i.ref]])

	# Add chromosomes.
	for(i in 1:n)
	{
		chrs<-c.list[[i]]
		lengths<-l.list[[i]]
		rect(1,i-chrExp,sum(lengths),i+chrExp,border=NA,col="white")
		index<-!is.na(chrs)
		x2<-convert.coord(chrs[index],lengths[index],chrs,lengths)
		x1<-x2-lengths[index]+1
		y1<-i-chrExp*0.7
		y2<-i+chrExp*0.7
		cols<-c.col.list[[i]][chrs[index]]
		rect(x1,y1,x2,y2,col=cols,border=NA)
	}
	axis(2,at=1:n,labels=spes,tick=FALSE,las=1)
}

#===========Get orthologous gene alignments for a list of species from TOGA results (delta).
genes_one2one<-function(spes,TOGA_dir="/projects/hillerlab/genome/gbdb-HL/hg38/TOGA")
{
	f_toga_isoforms<-sprintf("%s/toga.isoforms.tsv",TOGA_dir)
	toga_isoforms<-read.table(f_toga_isoforms,header=T,sep="\t",stringsAsFactors=F)
	genes<-sort(unique(toga_isoforms[,1]))

	alignments<-list()
	for(spe in spes)
	{
		f_orthology<-sprintf("%s/vs_%s/orthology_classification.tsv",TOGA_dir,spe)
		orthology<-read.table(f_orthology,header=T,sep="\t",stringsAsFactors=F)
		index<-orthology[,5]=="one2one"
		orthology<-unique(orthology[index,c(1,3)])
		ref2que<-orthology[,2]
		names(ref2que)<-orthology[,1]

		f_gene_spans<-sprintf("%s/vs_%s/query_gene_spans.bed",TOGA_dir,spe)
		gene_spans<-read.table(f_gene_spans,sep="\t",stringsAsFactors=F)
		que2crd<-gene_spans[,1:3]
		index<-gene_spans[,6]=="-"
		tmp<-que2crd[index,2]
		que2crd[index,2]<-que2crd[index,3]
		que2crd[index,3]<-tmp
		rownames(que2crd)<-gene_spans[,4]

		index<-ref2que[genes]
		alignments<-c(alignments,list(que2crd[index,]))
	}
	names(alignments)<-spes
	list(genes=genes,alignments=alignments)
}

#===========Get chromosomes and their lengths for a list of species (delta).
# minChrLen: scaffolds shorter than this will be excluded.
get.chromosomes<-function(spes,minChrLen=5e6)
{
	files<-sprintf("/projects/hillerlab/genome/gbdb-HL/%s/chrom.sizes",spes)
	names(files)<-spes
	lapply(files,function(file)
	{
		a<-read.table(file)
		index<-a[,2]>=minChrLen
		a[index,]
	})
}

#===========Wrapper (delta).
toga2synteny<-function(spes,ref,TOGA_dir,minChrLen)
{
	tmp<-genes_one2one(spes,TOGA_dir)
	a.list<-tmp$alignments

	tmp<-get.chromosomes(spes,minChrLen)
	c.list<-lapply(tmp,function(x)x[,1])
	l.list<-lapply(tmp,function(x)x[,2])

	# Remove genes containing NA or on chromosomes < minChrLen.
	index<-integer()
	for(spe in spes)
	{
		a<-a.list[[spe]]
		chrs<-c.list[[spe]]
		index<-c(index,which(rowSums(is.na(a))>0))
		index<-c(index,which(!(a[,1]%in%chrs)))
	}
	index<-sort(unique(index))
	if(length(index)>0)a.list<-lapply(a.list,function(x)x[-index,])

	b.list<-get.blocks(a.list,ref)
	echo<-sprintf("Generated %d synteny blocks based on %d genes.\n",nrow(b.list[[ref]]),sum(b.list[[ref]][,4]))
	cat(echo)
	list(b.list=b.list,c.list=c.list,l.list=l.list)
}
plot.toga2synteny<-function(from_toga2synteny,ref,arrange=FALSE)
{
	with(from_toga2synteny,
		plot.blocks(b.list,c.list,l.list,ref,chrGap=5e6,chrWid=0.2,arrange))
}
