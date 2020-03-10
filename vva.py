#!/usr/bin/env python
import os
import sys
import gzip
import time
import string

try:
	import tabix
except ImportError:
	print "Please install pytabix module [ pip install pytabix ]"
try:
	from optparse import OptionParser
except ImportError:
	print "Please Install optparse module [ pip install optparse ]"
try:
	import multiprocessing
except ImportError:
	print "Please install multprocessing module [ pip install multiprocessing ]"

parser = OptionParser()
parser.add_option("-a", "--case",	dest="casevcf",		help="input case VCF",					default=False)
parser.add_option("-b", "--control",	dest="ctrlvcf",		help="input control VCF",				default=False)
parser.add_option("-g", "--gene",	dest="genelist",	help="target gene list",				default=False)
parser.add_option("-c", "--consevation",dest="conservation",	help="mark conserved variant, default is \"FALSE\"",	default=False)
parser.add_option("-o", "--output",	dest="output",		help="out directory name, default is \"OUTPUT\"",	default="OUTPUT")
parser.add_option("-n", "--process",	dest="n_processes",	help="number of processes to use",			default=4)
(options, args) = parser.parse_args()

##### ##### ##### Messages ##### ##### #####
def printMsg():
        print """
######################################################
##  ##   ##  ##   ##   #####   #                     #
##  ##   ##  ##   ##  ##   ##  #  Variant            #
##   ## ##    ## ##   #######  #  Visualization and  #
##   ## ##    ## ##   ##   ##  #  Annotation         #
##    ###      ###    ##   ##  #                     #
######################################################

Usage: ./vva -a <caseVCF> -b <controlVCF> [-options]

Options:   -g/--gene           <string>      target gene list or gene symbols with comma seperator
           -c/--conservation   <string>      show conservation score (GERP++), defualt if FALSE
           -o/--output         <string>      output directory name, default is OUTPUT
           -n/--process        <number>      number of processes to use, default is 4
"""

##### ##### ##### R Code ##### ##### #####
R_code = """
library(Cairo)

FET <- function(dat){
	x <- matrix(as.numeric(dat), byrow=T, nrow=2)

	obj <- try(fisher.test(x), silent=TRUE)
	if(is(obj, "try-error")) {
		value <- NA
	} else {
		value <- obj$p.value
	}
	return(value)
}

rsidUpdate <- function(dat){
	x <- matrix(dat)
	re <- x[2,]
	if(is.na(re)){
		re <- x[1,]
	}	
	return(re)
}

getPos <- function(dat){
        pos <- as.numeric(strsplit(dat, ":")[[1]][2])
        return(pos)
}

all_variants <- read.table(paste(output, "/genes/", gene_symbol, ".txt", sep=""), header=TRUE)

if(dim(all_variants)[1] > 0){
	png_file_name = paste(output, "/", gene_symbol, ".png", sep="")

	png_width = 600 + (dim(all_variants)[1]*20)
	png_height = 900

	Cairo(width = png_width, height = png_height, file = png_file_name, type = "png", bg = "white", canvas = "white", units = "px", dpi = 200)
	
	par(xpd=NA)
	par(oma=c(2,0,0,2))

	colors <- c(rep("grey80",dim(all_variants)[1]))

	names <- apply(all_variants[,c("Pos","rsID")], 1, rsidUpdate)
	all_variants$Pos2 <- apply(matrix(all_variants$Pos), 1, getPos)
	
	fisher_table <- cbind(case_num - all_variants$CASE, all_variants$CASE, ctrl_num - all_variants$CTRL, all_variants$CTRL)
	if(dim(fisher_table)[1] > 1){
                all_variants$p <- apply(as.matrix(fisher_table[, c(1:4)]), 1, FET)
        } else {
                all_variants$p <- apply(as.matrix(fisher_table), 1, FET)
        }

	case_variants <- c(rep(" ",dim(all_variants)[1]))
	case_variants[which(all_variants$CASE + all_variants$CTRL > 0)] <- "*"

	highlights <- c(rep(" ",dim(all_variants)[1]))
	highlights[which(all_variants$p < 0.05)] <- "*"

	reported_snvs <- length(all_variants$CASE)
	observed_snvs <- length(which(all_variants$CASE + all_variants$CTRL > 0))

	if(dim(all_variants)[1] < 2){
                figure_main = paste(gene_symbol,"(",gene_strand,")", sep="")
        } else {
                figure_main = paste(gene_symbol,"(",gene_strand,") ",gene_location, sep="")
        }

        b <- barplot(as.vector(all_variants$SIFT), main=figure_main, ann=FALSE, yaxt="n", names.arg=paste(names,all_variants$Base,sep=" "), cex.names=0.5, las=3, ylim=c(0,1), col=colors , border=NA)

	if(isTRUE(gerp_score)){
		conservation <- c(rep(" ",dim(all_variants)[1]))
		conservation[which(all_variants$GERP > 2)] <- "+"
		text(b, 1, label=conservation, pos = 3, cex = 0.7, col = "black")
	}

	text(b, -0.07, labels=case_variants, col="black")
	text(b, -0.07, labels=highlights, col="red")

	ticks <- c()

	if(dim(b)[1]>2){

        	for(j in 2:dim(b)[1]-1){
                	alpha <- floor((all_variants$Pos2[j+1]-all_variants$Pos2[j])/10)
	                if(alpha > 0){
        	                by <- (b[j+1,1]-b[j,1])/alpha
                	        interval <- seq(from=(b[j,1])+(by/2.0), to=b[j+1,1], by=by)
                        	ticks <- c(ticks, interval)
	                }
        	}
        	axis(side = 1, at = ticks, labels=FALSE)
	}

	if(length(b) > 1){
		lines(b, all_variants$AFR/661, type="l", col="red", ylim=c(0,1), lwd=0.8)

		lines(b, all_variants$AMR/347, type="l", col="blue", ylim=c(0,1), lwd=0.8)

		lines(b, all_variants$EUR/503, type="l", col="green", ylim=c(0,1), lwd=0.8)

		lines(b, all_variants$SAS/489, type="l", col="gold", ylim=c(0,1), lwd=0.8)					

		lines(b, all_variants$EAS/504, type="l", col="orange", ylim=c(0,1), lwd=0.8)					

		lines(b, all_variants$CTRL/ctrl_num, type="l", col="black", ylim=c(0,1), lwd=0.8)

		lines(b, all_variants$CASE/case_num, type="l", col="black", ylim=c(0,1), lty=2)

	} else {
		points(b, all_variants$AFR/661, col="red", ylim=c(0,1), pch=20)

		points(b, all_variants$AMR/347, col="blue", ylim=c(0,1), pch=20)

		points(b, all_variants$EUR/503, col="green", ylim=c(0,1), pch=20)

		points(b, all_variants$SAS/489, col="gold", ylim=c(0,1), pch=20)                                     

		points(b, all_variants$EAS/504, col="orange", ylim=c(0,1), pch=20)                                     

		points(b, all_variants$CTRL/ctrl_num, col="black", ylim=c(0,1), pch=20)

		points(b, all_variants$CASE/case_num, col="black", ylim=c(0,1), pch=17, cex=1)

	}

	axis(2, at=c(seq(0,1,0.2)), label=c(seq(0,100,20)), las=2)
	mtext("Allele carrier frequency (%)", side = 2, line=3)

	axis(4, at=c(seq(0,1,0.2)), label=c("0","0.2","0.4","0.6","0.8","1"), las=2, col="black",col.axis="black")
	mtext("SIFT Score", side = 4, line=3, col="black")

	mtext(paste(observed_snvs," / ", reported_snvs," SNVs reported in VCFs",sep=""), side = 1, line=6, las=1, cex=0.8)

	if(length(b) > 1){
		legend("topleft", legend=c(paste("CASE,N=",case_num,sep=""),paste("CTRL,N=",ctrl_num,sep=""),"AFR,N=661","AMR,N=347","EAS,N=504","SAS,N=489","EUR,N=503"), lwd=c(1,1,1,1,1,1,1), lty=c(2,1,1,1,1,1,1), col=c("black","black","red","blue","orange","gold","green"), bty="n",cex=0.5)
	} else { 
		legend("topleft", legend=c(paste("CASE,N=",case_num,sep=""),paste("CTRL,N=",ctrl_num,sep=""),"AFR,N=661","AMR,N=347","EAS,N=504","SAS,N=489","EUR,N=503"), pch=c(17,20,20,20,20,20,20), col=c("black","black","red","blue","orange","gold","green"), bty="n",cex=0.5)
	}

	dev.off()
}
quit("no")
"""

##### ##### ##### Mean ##### ##### #####
def mean(alist):
        num = 0
        for n in alist:
                num += len(n)
        return float(num / len(alist))

##### ##### ##### Check input options ##### ##### #####
def checkArguments(options):
	global gcasevcf, gctrlvcf
	if not options.casevcf or not os.path.exists(options.casevcf):
 		print("ERROR::Invalid case VCF file %s"%(options.casevcf))
                return False
	else:
		gcasevcf = options.casevcf
        if not options.ctrlvcf or not os.path.exists(options.ctrlvcf):
                print("ERROR::Invalid ctrl VCF file %s"%(options.ctrlvcf))
                return False
	else:
		gctrlvcf = options.ctrlvcf

        try:
                int(options.n_processes)
        except Exception, e:
                print("ERROR::Invalid process defined %s"%(options.n_processes))
                return False

        return True

##### ##### ##### VCF validity ##### ##### #####
def checkVCF(filename):
	sample_count = False
	infile = gzip.open(filename, 'rb')
	try:
		for line in infile.readlines():
			if line.startswith("#CHROM"):
				items = line.strip().split("\t")
				sample_count = len(items) - 9
				break
		infile.close()
	except:
		infile.close()

	return sample_count

##### ##### ##### Is TXT or gene symbol(s) ##### ##### #####
def isTxt(filename):
	if os.path.isfile(filename):
		s = open(filename).read(512)
		text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
		_null_trans = string.maketrans("", "")

		if not s: # Empty files are considered text
        	        return True

		if "\0" in s: # Files with null bytes are likely binary
			return False

		t = s.translate(_null_trans, text_characters)

		if float(len(t))/float(len(s)) > 0.30:
			return False
	else:
		return False

	return True

##### ##### ##### Query VCF and retrun genotype summary ##### ##### #####
def queryVCF(filename, query):
	return_dic = {}
	tb = tabix.open(filename)
	try:
		output_call = tb.querys(query)
		sample_genotype = []
		for items in output_call:
			genotypes = []
                        chr = items[0]
                        pos = items[1]
                        ref = items[3]
                        alt = items[4].split(",")

			if len(ref) > 1 or mean(alt) != 1:
				continue

                        length = len(alt)
                        cnt = len(items[9:])

                        for j in range(9, cnt + 9):
				buff = items[j].split(":")[0]
				if buff == "." or buff == "./.":
					buff = "0/0"
                                genotypes.append(buff)

                        if length > 1:
                                for index in range(0, length):
					key = chr + ":" + str(pos) + ref + ">" + alt[index]
                                        gts = list(genotypes)
                                        for jndex in range(0, cnt):
                                                if genotypes[jndex] != ".":
                                                        buff = genotypes[jndex].split("/")
                                                        gt = []
                                                        for b in buff:
                                                                if int(b) == 0:
                                                                        gt.append('0')
                                                                elif int(b) != int(index+1):
                                                                        gt.append('0')
                                                                elif int(b) == int(index+1):
                                                                        gt.append('1')
                                                        gts[jndex] = '/'.join(gt)
                                        sample_genotypes = list(gts)
					
					wt = 0
					het = 0
					hom = 0

					for gt in sample_genotypes:
						alleles = gt.split("/")
						if alleles[0] == "0" and alleles[1] == "0":
							wt += 1
						elif alleles[0] != alleles[1]:
							het += 1
						else:
							hom += 1
					return_dic[key] = str(wt) + ":" + str(het) + ":" + str(hom)

                        else:
				key = chr + ":" + str(pos) + ref + ">" + alt[0]
                                sample_genotypes = list(genotypes)

				wt = 0
				het = 0
				hom = 0
				for gt in sample_genotypes:
					alleles = gt.split("/")
					if alleles[0] == "0" and alleles[1] == "0":
						wt += 1
					elif alleles[0] != alleles[1]:
						het += 1
					else:
						hom += 1
				return_dic[key] = str(wt) + ":" + str(het) + ":" + str(hom)

	except tabix.TabixError:
		print "ERROR::Query failed - " + filename
		pass
	
	return return_dic
				
##### ##### ##### Get T1GP Data ##### ##### #####
def queryT1gp(gene, keys):
	return_dic = {}
	infile = gzip.open("./resource/" + gene + ".vcf.gz", 'rb')
	for line in infile.readlines():
		if not line.startswith("#"):
			items = line.strip().split("\t")
			chr, pos, rsid, ref, alt = items[:5]
			key = chr + ":" + str(pos) + ref + ">" + alt
			filt = items[6]
			info = items[7].split(";")

			if key in keys or filt == "T1GP":
				return_dic[key] = info
	infile.close()
	return return_dic


##### ##### ##### VVA Start ##### ##### #####
def VVAStart(gene, query, strand, queue):
	print "VVA::" + gene 

	outfile = open(options.output + "/genes/" + gene + ".txt", 'w')
	outfile.write("Pos\tBase\trsID\tSIFT\tGERP2\tCASE\tCTRL\tAFR\tAMR\tEAS\tEUR\tSAS\n")
	
	cases = queryVCF(gcasevcf, query)
	ctrls = queryVCF(gctrlvcf, query)
	keys = list(set(list(cases.keys()) + list(ctrls.keys())))
	t1gp = queryT1gp(gene, keys)
	keys = list(set(list(cases.keys()) + list(ctrls.keys()) + list(t1gp.keys())))
	for key in sorted(keys):
		case_carrier = 0
		ctrl_carrier = 0
		t1gp_info = {}

		if key in cases:
			case_gt = cases[key].split(":")
			case_carrier = int(case_gt[1]) + int(case_gt[2])			
		if key in ctrls:
			ctrl_gt = ctrls[key].split(":")
			ctrl_carrier = int(ctrl_gt[1]) + int(ctrl_gt[2])			
		if key in t1gp:
			t1gp_buff = t1gp[key]
			for item in t1gp_buff:
				key2 = item.split("=")[0]
				value = item.split("=")[1]
				if key2 in ["rsID", "GERP++"]:
					if value == ".":
						value = "NA"
				if key2 in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
					if not value.startswith("NA"):
						value_tmp = value.split(":")
						value = int(value_tmp[1]) + int(value_tmp[2])
					else:
						value = 0
				t1gp_info[key2] = value
		else:
			t1gp_info['rsID'] = "NA"
			t1gp_info['SIFT'] = "NA"
			t1gp_info['GERP++'] = "NA"
			t1gp_info['AFR'] = 0
			t1gp_info['AMR'] = 0
			t1gp_info['EAS'] = 0
			t1gp_info['EUR'] = 0
			t1gp_info['SAS'] = 0

		if t1gp_info['SIFT'] != "NA":
			outfile.write(key[:-3] + "\t" + key[-3:] + "\t" + t1gp_info['rsID'] + "\t" + t1gp_info['SIFT'] + "\t" + t1gp_info['GERP++'] + "\t" + str(case_carrier) + "\t" + str(ctrl_carrier) + "\t" +  str(t1gp_info['AFR']) + "\t" + str(t1gp_info['AMR']) + "\t" + str(t1gp_info['EAS']) + "\t" + str(t1gp_info['EUR']) + "\t" + str(t1gp_info['SAS']) + "\n")

	outfile.close()

	outfile = open(options.output + "/" + gene + ".R", 'w')
	outfile.write("gene_symbol <- \"" + gene + "\"\n")
	n_strand = "+"
	if strand < 0:
		n_strand = "-"
	outfile.write("gene_strand <- \"" + n_strand + "\"\n")
	outfile.write("gene_location <- \"" + query + "\"\n")
	outfile.write("case_num <- " + str(gcasevcf_cnt) + "\n")
	outfile.write("ctrl_num <- " + str(gctrlvcf_cnt) + "\n")
	conserv_t_f = "FALSE"
	if options.conservation:
		conserv_t_f = "TRUE"
	outfile.write("gerp_score <- " + conserv_t_f + "\n")
	outfile.write("output <- \"" + options.output + "\"\n")
	outfile.write(R_code)
	outfile.close()

	os.system("R CMD BATCH " + options.output + "/" + gene + ".R")
	os.system("mv " + options.output + "/" + gene + ".R " + options.output + "/genes/")
	os.system("mv " + gene + ".Rout " + options.output + "/genes/")

	return

##### ##### ##### Main ##### ##### #####
def main():
	if not checkArguments(options):
		printMsg()
                sys.exit(2)

	##### ##### ##### Check VCFs ##### ##### #####
	global gcasevcf_cnt, gctrlvcf_cnt
	gcasevcf_cnt = checkVCF(gcasevcf)
	gctrlvcf_cnt = checkVCF(gctrlvcf)

	if not gcasevcf_cnt or not gctrlvcf_cnt:
		print("ERROR::Invalid VCF(s) files - require: VCF header, compressed by bgzip")
		sys.exit(2)

	##### ##### ##### Output directiry setting ##### ##### #####
	working_dir = os.getcwd()
	new_wd = working_dir + "/" + options.output
	if not os.path.exists(new_wd + "/genes"):
		os.makedirs(new_wd + "/genes")
		os.chmod(new_wd, 0777)
		os.chmod(new_wd + "/genes", 0777)

	##### ##### ##### Get target gene list ##### ##### #####
	genes_dic = {}
	genes_tmp = []
	infile = open("./resource/HUGO_protein_coding_genes.txt")
	for line in infile.readlines():
		items = line.strip().split("\t")
		symbol = items[0]
		position = items[1]
		strand = items[2]

		genes_tmp.append(symbol)
		genes_dic[symbol] = position + "$$" + str(strand)
	infile.close()

	genes_input = []
	if options.genelist:
		istext_flag = isTxt(options.genelist)

	if options.genelist and istext_flag:
		infile = open(options.genelist)
		for line in infile.readlines():
			genes_input.append(line.strip())
	elif options.genelist and not istext_flag:
		genes_input = options.genelist.split(",")
	else:
		genes_input = genes_tmp

	not_hugo = set(genes_input).difference(set(genes_tmp))
	if len(not_hugo) > 0:
		print("WARNING::Invalid gene symbol(s) - " + " ".join(sorted(not_hugo)))

	genes = set(genes_input).intersection(set(genes_tmp)) 

	##### ##### ##### Multiprocessing ##### ##### #####
	manager = multiprocessing.Manager()
	queue = manager.Queue()
	pool = multiprocessing.Pool(int(options.n_processes))

	jobs = []

	for gene in genes:
		query = genes_dic[gene].split("$$")[0]
		strand = genes_dic[gene].split("$$")[1]

		job = pool.apply_async(VVAStart, (gene, query, strand, queue))
		jobs.append(job)

	for job in jobs:
		job.get()

	queue.put('STOP')
	pool.close()
	pool.join()

if __name__ == "__main__":
#       main(sys.argv[1:])
        main()

## VVA v1.0 (@'3(o.o);
## EOF 2017.06.26 16:28
