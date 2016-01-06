__author__ = 'dsvet'


"""
Script converts BED region to all possible substitutions

"""

NUCLEOTIDES = ['A', 'C', 'T', 'G']

import bx.seq.twobit
import sys
import optparse

###----Read 2bit file---###

fmt = optparse.IndentedHelpFormatter( indent_increment=2, max_help_position=9, width=79, short_first=1 )
usage = "Usage: %prog [options]"
parser = optparse.OptionParser( usage = usage, version = "%prog v1.0", formatter = fmt )
parser = optparse.OptionParser()
parser.add_option( "-f", "--fasta_2bit", action = "store", type = "string", dest = "fasta_2bit", help = '2bit file with fasta' )
parser.add_option( "-i", "--input_file", action = "store", type = "string", dest = "input_file", help = 'input bed file' )
(options, args) = parser.parse_args()


PATH_TO_FASTA_FILE = options.fasta_2bit
PATH_TO_BED_FIEL = options.input_file

def select_fasta_for_region(chrID_fasta_dict, region_bed_format_list):
    chr_name = region_bed_format_list[0]
    start = int(region_bed_format_list[1])
    end = int(region_bed_format_list[2])
    fasta_sequence=chrID_fasta_dict[chr_name].get(start,end).upper()
    return fasta_sequence


chrID_fasta_dict = bx.seq.twobit.TwoBitFile( open( PATH_TO_FASTA_FILE, 'rb' ) )

with open( PATH_TO_BED_FIEL, 'r' ) as my_file:
    for line in my_file:
        line = line.split()
        if len(line)==3:
            ###---iterate over nucleotides in bed file---###
            region_start = int(line[1])
            region_end = int(line[2])
            for nucl_position in range( region_start,region_end ):
                position_to_get = [ line[0], nucl_position, nucl_position+1 ]
                ref_fasta = select_fasta_for_region(chrID_fasta_dict, position_to_get)
                ###---generate substitutions---###
                subst_nucletids = [ nuc for nuc in NUCLEOTIDES if nuc not in ref_fasta ]
                for nucleotide in subst_nucletids:
                    string_to_print = "\t".join( [ line[0], str(nucl_position), str( nucl_position+1 ), ref_fasta, nucleotide  ] )
                    sys.stdout.write(string_to_print + "\n")
        if len(line) > 3:
            ###---iterate over nucleotides in bed file---###
            region_start = int(line[1])
            region_end = int(line[2])
            for nucl_position in range( region_start,region_end ):
                position_to_get = [ line[0], nucl_position, nucl_position+1 ]
                ref_fasta = select_fasta_for_region(chrID_fasta_dict, position_to_get)
                ###---generate substitutions---###
                subst_nucletids = [ nuc for nuc in NUCLEOTIDES if nuc not in ref_fasta ]
                for nucleotide in subst_nucletids:
                    string_to_print = "\t".join( [ line[0], str(nucl_position), str( nucl_position+1 ), ref_fasta, nucleotide, str(line[3])  ] )
                    sys.stdout.write(string_to_print + "\n")
















