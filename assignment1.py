import mysql.connector

__author__ = 'XXX'


class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.gene = ""

    
    def fetch_gene_coordinates(self, genome_reference, file_name):
        
        print "Connecting to UCSC to fetch data"
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
    
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print "Done fetching data"
                
    def get_sam_header(self):
        print "todo"
        
    def get_properly_paired_reads_of_gene(self):
        print "todo"
        
    def get_gene_reads_with_indels(self):
        print "todo"
        
    def calculate_total_average_coverage(self):
        print "todo"
        
    def calculate_gene_average_coverage(self):
        print "todo"
        
    def get_number_mapped_reads(self):
        print "todo"
        
    def get_gene_symbol(self):
        print "todo"
        
    def get_region_of_gene(self):
        print "ads"
        
    def get_number_of_exons(self):
        print "ads"
    
    def print_summary(self):
        print "Print all results here"
    
        
if __name__ == '__main__':
    print "Assignment 1"
    assignment1 = Assignment1()
    assignment1.print_summary()
    
    assignment1.fetch_gene_coordinates("hg19", "MYFILE.TXT") ## TODO change filename
    
    

