assignSeq2REsite <-
function(input.RD, REmap.RD, cut.offset=1, seq.length=36, allowed.offset=5, min.FragmentLength=60, max.FragmentLength=300,partitionMultipleRE=c("unique","average","estimate", "best", "random"))
{
		cat(date(), "Validating input ...\n");
		partitionMultipleRE = match.arg(partitionMultipleRE)
		if (missing(input.RD)) {
        stop("Missing required argument input.RD!")
    }
    if (class(input.RD) != "GRanges") {
        stop("No valid input.RD passed in. It needs to be GRanges object")
    }
    if (missing(REmap.RD)) {
        stop("Missing requirement argument REmap.RD!")
		}
    if (class(REmap.RD) != "GRanges") {
        stop("No valid REmap.RD passed in. It needs to be GRanges")
    }		
	allChr.map = unique(as.character(seqnames(REmap.RD)))
	allChr.input = unique(as.character(seqnames(input.RD)))
	cat(date(), "Prepare map data ...\n");
	REmap.RD = REmap.RD[as.character(seqnames(REmap.RD)) %in% allChr.input,]
	r1 = do.call(rbind,lapply(allChr.input,function(chr)
	{
               	   if (chr %in% allChr.map) {
			cat(date(), "Align to chromosome ", chr, " ...\n");
			thisInput = input.RD[as.character(seqnames(input.RD)) == chr,]
			thisMap = REmap.RD[as.character(seqnames(REmap.RD))== chr,]
			Input.IR = IRanges(start=start(thisInput),end=end(thisInput),
			  names = names(thisInput))
			Map.IR = IRanges(start=start(thisMap),end=end(thisMap),
			   names = names(thisMap))
			this.strand = as.character(strand(thisInput))
			rm(thisInput)
			rm(thisMap)
 	                nearestRE = Map.IR[nearest(Input.IR, Map.IR)]	
			Distance = start(Input.IR) - start(nearestRE);                  
			data.frame(SEQid= names(Input.IR), REid = names(nearestRE),Chr = c(chr), strand=this.strand,SEQstart=start(Input.IR), SEQend=end(Input.IR),REstart=start(nearestRE), REend=end(nearestRE), Distance=Distance, Weight=c(1))
                }
     }))
	cat(date(), "Finished 1st round of aligning! Start the 2nd round of aligning ...\n");
	r1.keep = r1[(r1$Distance <= (cut.offset + allowed.offset) & as.character(r1$strand) %in% c("1", "+") & r1$Distance >=0) | (r1$Distance >= (cut.offset-seq.length - allowed.offset) & as.character(r1$strand) %in% c("-1", "-") & r1$Distance <=0),]
	newMap = unique(data.frame(REid=r1.keep$REid,Chr=r1.keep$Chr, REstart=r1.keep$REstart, REend=r1.keep$REend))
	newInput = r1[!r1$SEQid %in% r1.keep$SEQid,]
	rm(r1)
	allChr.map = unique(as.character(newMap$Chr))
	allChr.input = unique(as.character(newInput$Chr))
	r2 = do.call(rbind,lapply(allChr.input,function(chr)
		{
               if (chr %in% allChr.map) {
				 cat(date(), "Align to chromosome ", chr, " ... \n");
				thisInput = newInput[as.character(newInput$Chr) == chr,]
				thisMap = newMap[as.character(newMap$Chr)== chr,]
				Input.IR = IRanges(start=thisInput$SEQstart,end=thisInput$SEQend, names=thisInput$SEQid)
				id.strand = cbind(as.character(thisInput$SEQid),as.character(thisInput$strand))
				colnames(id.strand) = c("SEQid", "strand")
				Map.IR = IRanges(start=thisMap$REstart,end=thisMap$REend,
				names=thisMap$REid)
            			matches = findOverlaps(Input.IR, Map.IR, maxgap = max.FragmentLength, select ="all")
                		matchmatrix = as.matrix(matches)
                		qname = Input.IR[matchmatrix[, 1]]
                		tname = Map.IR[matchmatrix[, 2]]
				Distance = start(qname) - start(tname)
				SEQid = cbind(names(qname))
				colnames(SEQid) = "SEQid"
				this.strand = merge(id.strand, SEQid, by="SEQid")
                		data.frame(SEQid = names(qname), REid = names(tname), Chr = c(chr), strand = this.strand[,2], SEQstart=start(qname), SEQend=end(qname),REstart=start(tname), REend=end(tname), Distance=Distance, Weight= c(1))         
                }
     }))
	cat(date(), "Start filtering ... \n");
	r.passed = r2[(r2$Distance <= 0 & abs(r2$Distance) <= max.FragmentLength  & abs(r2$Distance) >= min.FragmentLength) |(r2$Distance >= 0 & abs(r2$Distance) <= (max.FragmentLength -seq.length) & abs(r2$Distance) >= (min.FragmentLength -seq.length)),]
	r.notpassed = r2[!r2$SEQid %in% r.passed$SEQid,]
	
	seqID.count = as.data.frame(table(r.passed$SEQid))
	colnames(seqID.count) = c("SEQid", "count")

	seqID.mRE = as.character(seqID.count[seqID.count$count>1,]$SEQid)
	r.passed.mRE = r.passed[as.character(r.passed$SEQid) %in% seqID.mRE,]
	r.passed.uRE = r.passed[!as.character(r.passed$SEQid) %in% seqID.mRE,]
	rm(seqID.mRE)
	#save(r.passed,file="rpassed.rda")
	#save(r1.keep,file="r1keep.rda")
	#save(r.notpassed,file="notpassed.rda")
	r.unique = rbind(r1.keep,r.passed.uRE)
	if(partitionMultipleRE=="unique")
	{
		list(passed.filter=r.unique,notpassed.filter=rbind(r.notpassed,r.passed.mRE)) 	}
	else 
	{
		cat(date(),"Partitioning reads over RE sites within ", max.FragmentLength, "...\n");
		if (partitionMultipleRE == "average" || partitionMultipleRE == "random")
		{
			nMapped = rowsum(r.passed.mRE$Weight, group=r.passed.mRE$SEQid)
			nMapped = cbind(rownames(nMapped),nMapped[,1])
			colnames(nMapped) = c("SEQid", "nMapped")
			r.passed.mRE.1 = merge(r.passed.mRE,nMapped, by="SEQid")
			if (partitionMultipleRE == "average")
			{
				r.passed.mRE.1$Weight = 1/as.numeric(as.character(r.passed.mRE.1$nMapped))
				list(passed.filter=rbind(r.unique,r.passed.mRE.1[,1:10]), notpassed.filter=r.notpassed,mREwithDetail=r.passed.mRE.1)
			}
			else
			{
				r.passed.mRE.1$Weight = 0
				i = 1
				while (i < dim(r.passed.mRE.1)[1])
				{
					j = sample(seq(1,as.numeric(as.character(r.passed.mRE.1$nMapped[i]))), 1)
					r.passed.mRE.1$Weight[i + j -1] =1
					i = i + as.numeric(as.character(r.passed.mRE.1$nMapped[i]))				
				}
				r.passed.mRE.1 = r.passed.mRE.1[r.passed.mRE.1$Weight ==1,]
				list(passed.filter=rbind(r.unique,r.passed.mRE.1[,1:10]), notpassed.filter=r.notpassed,mREwithDetail=r.passed.mRE.1)
			}
		}
		else
		{
  			cat(date(), "get count for each RE ...\n")
			reID.count = as.data.frame(table(r1.keep$REid))
			colnames(reID.count) = c("REid", "count")
			r.passed.mRE.1 = merge(as.data.frame(r.passed.mRE), as.data.frame(reID.count),by="REid")
			temp = rowsum(r.passed.mRE.1$count, group=r.passed.mRE.1$SEQid)
			temp = cbind(rownames(temp), temp[,1])
			colnames(temp) = c("SEQid", "total.count")
			r.passed.mRE = merge(r.passed.mRE.1, temp, by="SEQid")
  			r.passed.mRE$Weight = r.passed.mRE$count/as.numeric(as.character(r.passed.mRE$total.count))
			if (partitionMultipleRE == "estimate")
			{
				list(passed.filter=rbind(r.unique,r.passed.mRE[,1:10]), notpassed.filter=r.notpassed,mREwithDetail=r.passed.mRE)
			}
			else
			{
				nMapped = as.data.frame(table(r.passed.mRE$SEQid))
				colnames(nMapped) = c("SEQid", "nMapped")
				r.passed.mRE = merge(r.passed.mRE,nMapped, by="SEQid")
				r.passed.mRE$Weight = 0
				i = 1
				while (i < dim(r.passed.mRE)[1])
				{
					r.passed.mRE$Weight[r.passed.mRE$count == max(r.passed.mRE$count[i:(i + as.numeric(as.character(r.passed.mRE$nMapped[i]))-1)]) & r.passed.mRE$SEQid == r.passed.mRE$SEQid[i]] =1
					i = i + as.numeric(as.character(r.passed.mRE$nMapped[i]))				
				}
				r.passed.mRE = r.passed.mRE[r.passed.mRE$Weight ==1,]
				list(passed.filter=rbind(r.unique,r.passed.mRE[,1:10]), notpassed.filter=r.notpassed,mREwithDetail=r.passed.mRE)	
			}
		}
	}
}

