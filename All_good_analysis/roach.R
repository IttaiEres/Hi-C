#!/software/R-3.2-el6-x86_64/bin/R

#R file for creating tabulation DFs on midway. Verified as working on midway's version of R (R-3.2-el6-x86_64)

###Grab necessary packages
library(plyr)
library(tidyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(plotly)
library(dplyr)
library(Hmisc)

###Define all necessary functions
#Simple function that finds the sets of bins in each block of reads that will yield the greatest sum of reads, regardless of read inclusion.
optimize.bin.sums <- function(coverage_df){
  
  #define some necessary variables based on the coverage df; initiliaze some variables
  bin.size <- coverage_df$end[1]-coverage_df$start[1]
  window.slide <- coverage_df$start[2]-coverage_df$start[1]
  index.slide <- (bin.size/window.slide) #obtain the value for how many indices to slide forward for the sliding window!
  loc.index <- 1 #initialize an index for moving through the coverage_df
  coverage_df$bin <- 0 #intiailize another column in the coverage_df storing whether a given range ends up being a bin (1 for yes)
  
  while(loc.index<=nrow(coverage_df)){#iterate through all the positions in the coverage df
    
    if(coverage_df$count[loc.index]==0){ #move on to the next index if the given row has no reads in that bin--it's irrelevant.
      loc.index <- loc.index+1 #increment the index
      next}# move on to iterate this while loop again
    
    #If we see a switch from 0 to >0, the window has slid to include a new read in it. First, define the full range of the reads we can find here:
    end.index <- loc.index + match(0, coverage_df$count[loc.index:nrow(coverage_df)])-2 #pulls out just b/f the first instance of 0 in the count column--i.e. end.index is the last index at which a non-zero value is found
    tmp.indices <- seq(loc.index, end.index, index.slide) #grab first set of possible bins from this read region
    tmp.max <- sum(coverage_df$count[tmp.indices]) #store the sum of the read counts from this set as the temporary maximum
    tmp.max.indices <- tmp.indices #initialize the tmp.max.indices variable with these first indices, for proper bin assignation after the loop
    
    #Include something here to deal with the case when the initial tmp.indices is only a single element? I.e. where a region doesn't have options?
    #Check on the other regions where this is happening--the issue might not be that, but, rather, that there are ties!
    
    while(head(tmp.indices, 1)<=(loc.index+index.slide-1)){#keep going until we've considered all the possible bin combinations in this block of reads--check by stopping if the first value in tmp.indices is greater than index.slide above the loc.index
      tmp.indices <- tmp.indices + 1 #increment to look at the next set of bins
      tmp.indices <- tmp.indices[which(tmp.indices<=end.index)] #get rid of indices in tmp.indices that go past the end.index (i.e. move into other read blocks that are separated by a 0!)
      
      if(!length(tmp.indices)){break} #If tmp.indices turns into an empty vector (cases where there are not multiple options for bin combos), break this loop and move on (not gonna get higher tmp.max).
      
      if(sum(coverage_df$count[tmp.indices])>tmp.max){#If we have found a new way to maximize total read count in a bin set,
        tmp.max.indices <- tmp.indices #store these given indices in the tmp.max.indices variable, for bin assignation after the loop
        tmp.max <- sum(coverage_df$count[tmp.indices]) #make this the new tmp.max
      }
    }
    coverage_df$bin[tmp.max.indices] <- 1 #assign a status of being a chosen bin to the indices that gave us the max sum!
    loc.index <- tail(tmp.max.indices, 1)+index.slide #Move the loc index to the next possible bin after the chosen set (tmp.max.indices)
  }
  return(coverage_df[which(coverage_df$bin==1),-5]) #return a df of all the bins and their counts. minus the "bin" column
}

#A function that tabulates the frequency of contact between different bins on a position dataframe.
basic_tabulation <- function(coverage_df, pos_df, distance1, distance2){
  colnames(pos_df) <- c("pos1", "pos2") #ensure the column names of the pos_df work with the next lines
  
  pos1_new <- ifelse(pos_df$pos1<pos_df$pos2, pos_df$pos1, pos_df$pos2) #reorders columns so pos1 < pos2. important for not repeating bin pairs
  pos2_new <- ifelse(pos_df$pos1<pos_df$pos2, pos_df$pos2, pos_df$pos1) #If true, returns second arg; if false, returns third.
  pos_df$pos1 <- pos1_new #assignation of the new columns
  pos_df$pos2 <- pos2_new #assignation of the new columns
  
  pos_df$diff <- abs(pos_df$pos1-pos_df$pos2) #make a third column--the difference between the first two
  pos_df <- pos_df[pos_df$diff<(distance2)&pos_df$diff>=(distance1),1:2] #pull out pairs in the given distance range
  
  if(!is.null(nrow(pos_df))){#If we've found pairs at the given distance, tabulate contact counts!
    
    bins.df <- optimize.bin.sums(coverage_df)
    bins.vec <- unique(sort(c(bins.df$start, bins.df$end))) #Pulls out the vector of bins for the given chromosome at given bin size.
    
    pos_df <- pos_df %>% mutate(pos1_bin = cut2(pos1, bins.vec, onlycuts=TRUE, oneval=FALSE, digits=9)) %>% mutate(pos2_bin = cut2(pos2, bins.vec, onlycuts=TRUE, oneval=FALSE, digits=9)) #make two new columns in the pos_df, with bin assignments for pos1 and pos2
    pos_df %>% group_by(pos1_bin, pos2_bin) %>% summarise(count=n()) #tabulate the bin contacts that are repeated to find total bin-bin contacts
  }
  else(return(0))#print("No contacts at that distance on this chromosome!"))
}

#Creates a massive long-form dataframe containing information on count, mean, locus 1 and locus 2, and chromosome for each distance bin--called on individual chromosomes at a time.
full.chrom.tabulator <- function(coverage_df, pos_df, start, end, by, chr){
  distances <- seq(start, end, by) #Get the set of distances being looked at
  tabulation.list <- vector("list", length(distances)) #Initialize a blank list of length distances to store dfs of contact count dist'ns for different distance bins in.
  for(distance in 1:length(distances)){ #Iterate through all the user-input distances
    if(is.null(nrow(basic_tabulation(coverage_df, pos_df, distances[distance], distances[distance]+by)))){next} #Check that there are contacts at that distance on this chromosome. If not, move on.
    tmp.table <- basic_tabulation(coverage_df, pos_df, distances[distance], distances[distance]+by) #If contacts found, tabulate!
    tmp.table <- tmp.table[as.character(tmp.table$pos1_bin)!=as.character(tmp.table$pos2_bin),] #Remove any contacts that are between loci within the same bin. Meaningless! Or will this ruin the background model estimation? We shall see! Had to add as.character() to each of the values to resolve issue with level sets of factors being different.
    tabulation.list[[distance]] <- data.frame(chr=rep(chr, length(tmp.table$pos1_bin)),  dist.bin=rep(distances[distance], length(tmp.table$pos1_bin)), loc1_bin=tmp.table$pos1_bin, loc2_bin=tmp.table$pos2_bin, count=tmp.table$count) #Put a data frame with all the info about contacts at this distance into the list
  }
  tabulation.df <- do.call(rbind, tabulation.list) #Combine all the lists into a single long-form dataframe.
  return(tabulation.df)
}

#Function to obtain the values that represent the contact count # threshold for the top 'percents' of each distribution.
sig.thresh.finder <- function(full_tab_df, percents){
  
  bins <- unique(full_tab_df$dist.bin) #Grab the bins
  percents <- percents*.01 #Re-scale the percents to work with quantile in the loop below.
  
  values.df <- data.frame(matrix(nrow=length(bins), ncol=length(percents))) #Initialize the df to be returned
  colnames(values.df) <- percents #Initialize the df's colnames for iterating through
  rownames(values.df) <- bins #Initialize the df's rownames for iterating through
  
  for(bin in bins){
    for(percent in percents){
      values.df[as.character(bin), as.character(percent)] <- quantile(full_tab_df$count[full_tab_df$dist.bin==bin], 1-percent)[[1]]#To pull out the value above which are the top percent % of loci
    }
  }
  return(values.df)
}

#Function that takes a homer file, a specified chromosome, and a full_tab_df for a given locus size to produce a matrix of significant hit overlaps b/t my data and homer data at a variety of bin distances and percentages.
signif.comparison <- function(chr, full_info, locus_size, percents=c(20, 15, 10, 5, 2.5, 1, 0.5, 0.1, 0.01)){
  #Subset down the full data frame info to the chromosome of interest.
  chr_df <- full_info[full_info$chr==chr,]
  colnames(chr_df) <- c("chr", "dist.bin", "loc1_bin", "loc2_bin", "count") #Important for proper functioning of call to merge() later on.
  
  #Obtain the bins for assigning bins to homer paired loci.
  characters <- (unique(nchar(as.character(chr_df$loc1_bin)))-3)/2
  chr.bins <- unique(sort(c(as.numeric(substr(chr_df$loc1_bin, 2, characters+1)), as.numeric(substr(chr_df$loc1_bin, characters+3, 2*characters+2)), as.numeric(substr(chr_df$loc2_bin, 2, characters+1)), as.numeric(substr(chr_df$loc2_bin, characters+3, 2*characters+2)))))
  
  #Get a data frame of the values above which we would consider hits "significant" in the different distance bins.
  sig.values <- sig.thresh.finder(chr_df, percents)
  
  #Read in homer file, change column names to work with next lines
  homer.sig <- fread(paste("/project/gilad/ittai/roach/homer.signif/sig_", locus_size, "kb/", chr, "_", locus_size, "kb_siglocs.txt", sep=""), header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  colnames(homer.sig) <- c("start1", "end1", "start2", "end2")
  
  ##If the homer.sig df is empty, rename the columns appropriately and then skip all the next reformatting steps:
  if(nrow(homer.sig)==0){
    homer.sig <- homer.sig[,1:3]
    colnames(homer.sig) <- c("distance", "loc1_bin", "loc2_bin")
  }
  
  else{#If the homer.sig df is not empty, do all this necessary reformatting!
    
    ##Reformat the homer file to be easier to work with for matching hereafter.
    start1_new <- ifelse(homer.sig$start1<homer.sig$start2, homer.sig$start1, homer.sig$start2) #reorders columns so pos1 < pos2. important for not repeating bin pairs
    start2_new <- ifelse(homer.sig$start1<homer.sig$start2, homer.sig$start2, homer.sig$start1) #If true, returns second arg; if false, returns third.
    end1_new <-ifelse(homer.sig$start1<homer.sig$start2, homer.sig$end1, homer.sig$end2)
    end2_new <- ifelse(homer.sig$start1<homer.sig$start2, homer.sig$end2, homer.sig$end1)
    
    #Assign the new columns.
    homer.sig$start1 <- start1_new
    homer.sig$start2 <- start2_new
    homer.sig$end1 <- end1_new
    homer.sig$end2 <- end2_new
    
    #Get individual locations as midpoints of loci identified as significant by homer, so I can assign them to my set of bins for the given chr.
    homer.sig$loc1 <- rowMeans(homer.sig[,1:2]) #Get midpoint of start1 and end1, for loc1
    homer.sig$loc2 <- rowMeans(homer.sig[,3:4]) #Same for loc2 (start2 and end2)
    
    #Get another column that is distance between the two loci, for subsetting down later on.
    homer.sig$distance <- abs(homer.sig$loc1-homer.sig$loc2)
    
    
    #If the homer sig df isn't empty, assign the individual locations of homer-significant loci to my set of bins for the given chr.
    homer.sig <- homer.sig %>% mutate(loc1_bin=cut2(loc1, chr.bins, onlycuts=TRUE, oneval=FALSE, digits=9)) %>% mutate(loc2_bin=cut2(loc2, chr.bins, onlycuts=TRUE, oneval=FALSE, digits=9))
    
    #Subset down to just the pairs of significant hits bins from the homer data.
    homer.sig <- homer.sig[,7:9] #Gets rid of columns 1-6, column 7 with distance between loci and columns 8 and 9 with position bin assignations are all that matter.
    
    #Remove any contacts that are between loci within the same bin. Meaningless!I think. We shall see!
    homer.sig <- homer.sig[as.character(homer.sig$loc1_bin)!=as.character(homer.sig$loc2_bin),] 
  }
  
  #Grab the distance between bins quickly.
  slide <- as.numeric(unique(chr_df$dist.bin)[2]-unique(chr_df$dist.bin)[1])
  
  #Initialize a list to add the data frames to. Will convert to long-form df at the end!
  contacts.list <- vector("list", dim(sig.values)[1]*dim(sig.values)[2])
  
  #Get the names for the contacts.list to be compiled at the end. Critical for inserting dfs into it in the loop below.
  names(contacts.list) <- do.call("paste", c(melt(t(sig.values))[,1:2], sep=", "))
  
  #Iterate across the significance threshold data frame, creating individual dfs for each distance bin-percent combination, then adding a column indicating whether there's overlap with homer hits and finally merging them all together.
  for(bin in rownames(sig.values)){
    for(percent in colnames(sig.values)){
      tmp.my.sigs <- chr_df[(chr_df$dist.bin==bin)&((chr_df$count)>=(sig.values[bin, percent])),] #Pull out df of my sig hits for this bin and percentage.
      tmp.homer.sigs <- homer.sig[(homer.sig$distance)>=(as.numeric(bin))&(homer.sig$distance)<(as.numeric(bin)+slide), 2:3] #Same for homer's significant hits.
      
      #Add in a few key columns to tmp.my.sigs to make it easier to get the final df after
      tmp.my.sigs$percentile <- (1-as.numeric(percent))*100 #Indicate which percentile category these hits fall under
      tmp.my.sigs$homer.hit <- 0 #Indicate status WRT whether we can find this pair of contacts in homer significant hits. 0=no, 1=yes!
      tmp.my.sigs$my.hits.tot <- nrow(tmp.my.sigs) #How many total significant hits do I find?
      tmp.my.sigs$homer.hits.tot <- nrow(tmp.homer.sigs) #How many total significant hits does homer find?
      tmp.my.sigs$overlap <- nrow(merge(tmp.my.sigs[,3:4], tmp.homer.sigs)) #What's the overlap of those hits?
      
      tmp.my.sigs.vec <- do.call("paste", tmp.my.sigs[,3:4]) #Create a vector of all the rows for comparing
      tmp.homer.sigs.vec <- do.call("paste", tmp.homer.sigs) #Create a vector of all the rows for comparing
      
      matching.indices <- which(tmp.my.sigs.vec %in% tmp.homer.sigs.vec) #Extract row indices from my tmp.sigs that match homer sigs!
      
      tmp.my.sigs$homer.hit[matching.indices] <- 1 #For the rows that are matches with homer hits, indicate this!
      
      contacts.list[[paste(percent, bin, sep=", ")]] <- tmp.my.sigs
    }
  }
  
  final_df <- rbind.fill(contacts.list) #Turn the list into one long-form DF.
  rownames(final_df) <- NULL #Those are useless.
  return(final_df)
}

do.it.all <- function(locus.bin.size, start=0, end=250000, by=5000, variable.list=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")){
  coverage <- fread(paste("/project/gilad/ittai/roach/cis_only_cov_pos/coverage_cis_only/coverage_", locus.bin.size, "kb_5kb.bed", sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE, data.table=FALSE) #Pull out coverage for the given size
  colnames(coverage) <- c("chr", "start", "end", "count") #Rename the columns of the coverage df to work with other functions
  
  #Initialize named list for storing information on different chromosomes
  full.tab.list <- vector("list", length(variable.list))
  names(full.tab.list) <- variable.list
  
  #Begin iterating through all chromosomes!
  for(chromosome in variable.list){
    
    #Grab the coverage df and midpoint positions df for this chromosome
    tmp.coverage <- coverage[coverage$chr==chromosome,]
    tmp.pos.df <- fread(paste("/project/gilad/ittai/roach/cis_only_cov_pos/paired_read_midpoint_positions/midpoint.pairs.", substr(chromosome, 4, 5), sep=""), header=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
    
    #Get all the tabulation dfs for all chromosomes at this bin size. Form: tab_df.chrX.Xkb
    full.tab.list[[chromosome]] <- full.chrom.tabulator(tmp.coverage, tmp.pos.df, start, end, by, chromosome)
    print(paste("Tab list for chr ", chromosome, " completed with ", (nrow(full.tab.list[[chromosome]])), " rows.", sep=""))
  }
  
  #Once all the tabulation list has been created, assign it into one long-form df!
  full_tab_df <- rbind.fill(full.tab.list)
  
  print(paste("Full tab df completed for ", locus.bin.size, "kb, with dimensions ", dim(full_tab_df)))
  
  #Write this long-form df out into a separate file!
  write.table(full_tab_df, file=paste("/project/gilad/ittai/roach/full_data_frames/full_tab_df.", locus.bin.size, "kb", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
  
  #Initialize named list for storing homer overlap information on different chromosomes
  full.homer.overlap.list <- vector("list", length(variable.list))
  names(full.homer.overlap.list) <- variable.list
  
  #Iterate through dem chromosomes to do dat homer signif comparison!
  for(chromosome in variable.list){
    full.homer.overlap.list[[chromosome]] <- signif.comparison(chromosome, full_tab_df, locus.bin.size, percents=c(20, 10, 5, 1, 0.5, 0.1))
  }
  
  #Once all the homer overlap list has been created, assign it into one long-form df!
  full_homer_overlap_df <- rbind.fill(full.homer.overlap.list)
  
  #Write this long-form df out into a separate file!
  write.table(full_homer_overlap_df, file=paste("/project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.", locus.bin.size, "kb", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

##Actually do it!
for(locus.bin.size in c(seq(10, 100, 10), c(125, 150, 200))){
  do.it.all(locus.bin.size = locus.bin.size)
}