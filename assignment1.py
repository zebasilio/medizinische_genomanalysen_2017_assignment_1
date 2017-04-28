#! /usr/bin/env python3

'''
To run the script, in addition to the imported modules, please install samtools and bedtools: sudo apt-get install bedtools; sudo apt-get install samtools
Before running the script, change the path of the bam file in the bam variable (line 27)
to run the script, in the command line:
$chmod +x assignment1.py
$./assignment1.py > assignment1.txt 

Running this script as mentioned above will produce the following files:
mygene.txt
gene_proper_paired_reads.txt
flagstat.txt
assignment1.txt (where the assigment summary will be written)  

  

'''
import mysql.connector
import pysam
import pybedtools
import subprocess


__author__ = 'José Basílio'

bam = '/home/jose/MedGenAn/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam'


class Assignment1:
    def __init__(self):
        #Gene of interest
        self.gene = "AMPD3"

	#give path of the bam file
        self.bam = bam

	#convert bam to sam:
        self.sam = pysam.AlignmentFile(self.bam, 'rb')    
    
    def fetch_gene_coordinates(self, genome_reference, file_name):
        print("\n+++++++++++++++++++\nConnecting to UCSC and Fetching Data:")
        
        # Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        # Get cursor
        cursor = cnx.cursor()
        
        # Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        # Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        # Execute query
        cursor.execute(query)
        
        # Write to file
        
        self.gene_info = []
        with open(file_name, "w") as fh:
            for row in cursor:
                if row[0] == "AMPD3":
                    for attribute in row:
                        fh.write(str(attribute) + "\n")
                        self.gene_info.append(attribute)
 
    	##Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("\n+++++++++++++++++++\nDone Fetching Data:")
        
    	# get header of whole file
                
    def get_sam_header(self):
        print ("\n+++++++++++++++++++\nSam Header:")
        cmd = ["samtools_0.1.18 view -H {}".format(self.bam)]
        subprocess.call(cmd, shell=True)
        print()


 # http://pysam.readthedocs.io/en/latest/api.html
    # For the function, .bam file will be indexed (samtools index HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam ) and the proper paired reads of the gene will be written into a new  file "gene_proper_paired_reads.txt"
    def get_properly_paired_reads_of_gene(self):
        print ("\n+++++++++++++++++++\nproperly paired reads of gene:")
        # call samtools via command line; indexes the bam file
        cmd = ["samtools_0.1.18 index {}".format(self.bam)]
        subprocess.call(cmd, shell=True)
        samfile = pysam.AlignmentFile(self.bam, "rb")
        properpairedreads = pysam.AlignmentFile("gene_proper_paired_reads.txt", "w", template=samfile)
        for read in samfile.fetch("11", self.gene_info[3], self.gene_info[4]):
            if read.is_proper_pair:
             		properpairedreads.write(read)
        samfile.close()
        #counts the number of properly_paired_reads_of_gene and prints it
        file = open("gene_proper_paired_reads.txt", "r")
        read_count = 0
        reads = []
        for line in file:
            if line.startswith("SRR"):
               read_count += 1
               reads.append(line)
        print(read_count)
        return read_count
        
        

    def get_properly_paired_reads(self):
        print ("\n+++++++++++++++++++\nProperly paired reads:")

# call samtools via command line and save it in an extra file (flagsat.txt)
        cmd = ["samtools_0.1.18 flagstat {} > flagstat.txt".format(self.bam)]
        subprocess.call(cmd, shell=True)
        file = open("flagstat.txt", "r")
        rows = []
        for row in file:
            rows.append(row)
        target_row = rows[6]
        properly_paired_reads = target_row.split(" ")[0]
    
        print(properly_paired_reads)
        return properly_paired_reads

#Specifies and returns the number of reads with indels, CIGAR field either as I (Insert) or D (Delition)
	
    def get_gene_reads_with_indels(self):
        print ("\n+++++++++++++++++++\nGene Reads with Indels:")
        read_count = 0
        reads = []
        for read in self.sam.fetch("11", self.gene_info[3], self.gene_info[4]):
            columns = str(read).split("\t")
            if "I" in str(columns[5]) or "D" in str(columns[5]):
                read_count += 1
                reads.append(read)
        print(read_count)
        return reads
        
        
        
    def calculate_total_average_coverage(self):
        print("\n+++++++++++++++++++\nTotal average coverage:")
        bed = pybedtools.BedTool(self.bam)
        coverage = bed.genome_coverage(bg=True)
        line_count = 0 #number of regions
        total_coverage = 0 #added coverage across all regions
        for line in coverage:
            line_count += 1
            total_coverage += int(line[3])

        total_average_coverage = total_coverage / line_count
        print("%.2f" % round(total_average_coverage,2))
        return total_average_coverage



    def calculate_gene_average_coverage(self):
        print("\n+++++++++++++++++++\nGene Average Coverage:")

        bed = pybedtools.BedTool(self.bam)
        coverage = bed.genome_coverage(bg=True)

        line_count = 0 #number of regions
        total_coverage = 0 #added coverage across all regions
        for line in coverage:
            if int(line[1]) >= int(self.gene_info[3]) and int(line[2]) <= int(self.gene_info[4]):
                line_count += 1
                total_coverage += int(line[3])
                
        gene_average_coverage = total_coverage / line_count
        print("%.2f" % round(gene_average_coverage,2))
        return gene_average_coverage
        
        
    def get_number_mapped_reads(self):
        print("\n+++++++++++++++++++\nNumber of Mapped Reads:")
        file = open("flagstat.txt", "r")
        rows = []
        for row in file:
            rows.append(row)
        target_row = rows[2]
        number_mapped_reads = target_row.split (" ")[0]

        print(number_mapped_reads)
        return number_mapped_reads
	
        
    def get_gene_symbol(self):
        print("\n+++++++++++++++++++\nGene Symbol:")
	
        print(self.gene_info[0])
        return self.gene_info[0]
        
    def get_region_of_gene(self):
        print("\n+++++++++++++++++++\nRegion of Gene:")
        
        print("Chromosome: {}\nStart: {}\nEnd: {}\n".format(self.gene_info[2], self.gene_info[3], self.gene_info[4]))
        return [self.gene_info[2], self.gene_info[3], self.gene_info[3], self.gene_info[4]]
        
    def get_number_of_exons(self):
        print("\n+++++++++++++++++++\nNumber of Exons:")

        print(self.gene_info[6])
        return self.gene_info[6]

	
    
    def print_summary(self):
        self.fetch_gene_coordinates("hg19", "mygene.txt")
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_properly_paired_reads()	
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()

if __name__ == '__main__':
    print ("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
