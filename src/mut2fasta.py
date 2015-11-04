import bx.seq.twobit
from collections import OrderedDict
import optparse
import sys



LEFT=400
RIGHT=400

def select_fasta_for_region_SNV(chrID_fasta_dict, region_bed_format_list):
    chr_name = region_bed_format_list[0]
    start = int(region_bed_format_list[1])
    end = int(region_bed_format_list[2])
    fasta_sequence=chrID_fasta_dict[chr_name].get(start,end).upper()
    return fasta_sequence


def read_file_with_mutations(path_to_file):
    chr_start_end_ref_mut_dict=OrderedDict()
    with open(path_to_file,'r') as my_file:
        for line in my_file:
            line = line.split()
            region_ID="_".join([line[0],line[1],line[2],line[3],line[4]])
            chr_start_end_ref_mut_dict[region_ID]  = region_ID

    return chr_start_end_ref_mut_dict


def create_sequence_for_region(twobit_file_object,mutations_dict):
    for elem in mutations_dict:
        data = elem.split("_")
        ###----Get infor is it insertion, deletion of  SNV-----###
        if len(data[3])==1 and len(data[4])==1:
            ##---SNV---###
            chr = data[0]
            ref_nuc = data[3]

            ##---Take only first variant (NOTE these calls are usually not very confident)---###
            mut_nuc = data[4]
            mut_nuc = mut_nuc.split(",")[0]

            start_left = int(data[1])- LEFT
            end_right = int(data[1]) + 1 + RIGHT
            left_seq_to_take = select_fasta_for_region_SNV(chrID_fasta_dict, [chr, start_left, int(data[1]) ])
            right_seq_to_take =  select_fasta_for_region_SNV(chrID_fasta_dict, [chr, int(data[1])+1, end_right ])

            ###--merge--
            mut_sequence = left_seq_to_take + mut_nuc + right_seq_to_take
            ref_sequence = left_seq_to_take + ref_nuc + right_seq_to_take

            # print "SNVs"
            # print elem
            # print "mut_sequence",mut_sequence
            # print "ref_sequence", ref_sequence

            file_to_save_ref.write(">" +elem + "\n" )
            file_to_save_ref.write(ref_sequence + "\n" )


            file_to_save_mut.write(">" +elem + "\n" )
            file_to_save_mut.write(mut_sequence + "\n" )

        if len(data[3])==1 and len(data[4])>1:
            ###----Make reference and mutant fasta for insertion----###
            chr = data[0]
            ref_nuc = data[3]
            mut_nuc = data[4]
            start_left = int(data[1])- LEFT
            end_right = int(data[1]) + 1 + RIGHT
            left_seq_to_take = select_fasta_for_region_SNV(chrID_fasta_dict, [chr, start_left, int(data[1]) ])
            right_seq_to_take =  select_fasta_for_region_SNV(chrID_fasta_dict, [chr, int(data[1])+1, end_right ])

            ###--merge--
            mut_sequence = left_seq_to_take + mut_nuc + right_seq_to_take
            ref_sequence = left_seq_to_take + ref_nuc + right_seq_to_take

            # print "INS"
            # print elem
            # print "mut_sequence",mut_sequence
            # print "ref_sequence", ref_sequence


            file_to_save_ref.write(">" +elem + "\n" )
            file_to_save_ref.write(ref_sequence + "\n" )


            file_to_save_mut.write(">" +elem + "\n" )
            file_to_save_mut.write(mut_sequence + "\n" )

        if len(data[3])>1 and len(data[4])==1:
            ###----Make reference and mutant fasta for insertion----###
            chr = data[0]
            mut_nuc = data[3]
            ref_nuc = data[4]
            start_left = int(data[1])- LEFT
            end_right_ref = int(data[1]) + 1 + RIGHT
            end_right_mut  = int(data[1]) + 1 + RIGHT + len(mut_nuc) - 1
            left_seq_to_take = select_fasta_for_region_SNV(chrID_fasta_dict, [chr, start_left, int(data[1]) ])
            right_seq_to_take_mut =  select_fasta_for_region_SNV(chrID_fasta_dict, [chr, int(data[1]) +len(mut_nuc) , end_right_mut ])
            right_seq_to_take_ref = select_fasta_for_region_SNV(chrID_fasta_dict, [chr, int(data[1])+1, end_right_ref ])
            ###--merge--
            mut_sequence = left_seq_to_take + ref_nuc + right_seq_to_take_mut
            ref_sequence = left_seq_to_take + ref_nuc + right_seq_to_take_ref




            # print "DEL"
            # print elem
            # print "mut_sequence",mut_sequence
            # print "ref_sequence", ref_sequence


            file_to_save_ref.write(">" +elem + "\n" )
            file_to_save_ref.write(ref_sequence + "\n" )


            file_to_save_mut.write(">" +elem + "\n" )
            file_to_save_mut.write(mut_sequence + "\n" )





if __name__ == '__main__':
    fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
    usage = "Usage: %prog [options]"

    parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)

    parser = optparse.OptionParser()
    parser.add_option("-r", "--ref_save", action = "store", type = "string", dest = "ref_save", help = 'File with feature names')
    parser.add_option("-m", "--mut_save", action = "store", type = "string", dest = "mut_save", help = 'file with mutations (insertions)')
    parser.add_option("-f", "--fasta_path"   , action = "store", type = "string", dest = "fasta_path", help = 'Path to save results')
    parser.add_option("-v", "--mut_path", action = "store", type = "string", dest = "mut_path", help = 'Path to .2bit file')

    (options, args) = parser.parse_args()



    PATH_TO_FASTA_FILE=options.fasta_path
    PATH_TO_MUTS=options.mut_path
    PAHT_TO_SAVE_REF=options.ref_save
    PAHT_TO_SAVE_MUT=options.mut_save


    file_to_save_ref=open(PAHT_TO_SAVE_REF,'w')
    file_to_save_mut=open(PAHT_TO_SAVE_MUT,'w')

    # Check if we have an expression matrix filea FASTA or twobit file is given as input.
    if ( (options.ref_save is None) or (options.mut_save is None) or (options.fasta_path is None)  \
                 or (options.mut_path is None) ):
        parser.print_help()
        print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n'
        sys.exit(1)


    chrID_fasta_dict = bx.seq.twobit.TwoBitFile( open(PATH_TO_FASTA_FILE, 'rb') )
    chr_start_end_ref_mut_dict = read_file_with_mutations(PATH_TO_MUTS)
    create_sequence_for_region(chrID_fasta_dict,chr_start_end_ref_mut_dict)

    file_to_save_ref.close
    file_to_save_mut.close

