#!/usr/bin/env python
#
# cbust.py
#
# Copyright (C) 2012 - Gert Hulselmans
#
# Based on cbust.sh written by Bram Van de Sande.



import os
import sys
import subprocess
import optparse
import tempfile
from signal import signal, SIGPIPE, SIG_DFL


def which(program):
    """
    Check if a certain program exists and that it is executable.
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def prg_in_path(program):
    """
    Check if a certain program exists and that it is executable.
    If the program doesn't exist, display an error message and exit this script.
    """

    if which(program) is None:
        print >> sys.stderr, '\nProgram "{0:s}" is could not be found in $PATH.\n'.format(program)
        sys.exit(1)

    return None


def run_cmd(cmd):
    """
    Run the program with the provided arguments.

      - When the program succeeds, return stdout and stderr.
      - When the program fails, print stdout and stderr messages.
    """

    try:
        pid = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        stdoutdata, stderrdata = pid.communicate(None)
    except OSError, msg:
        print >> sys.stderr, "\nERROR during execution of '" + ' '.join(cmd) + "': " + str(msg)
        sys.exit(1)

    if (pid.returncode != 0):
        print >> sys.stderr, "\nERROR during execution of '" + ' '.join(cmd) + "':"
        print >> sys.stderr, "\nStandard output:\n----------------"
        print >> sys.stderr, stdoutdata
        print >> sys.stderr, "\nStandard error:\n---------------"
        print >> sys.stderr, stderrdata
        sys.exit(1)

    return [stdoutdata, stderrdata]



def fix_genomic_coordinates_cbust(twobit, bed_start, bed_end, bed_strand, feature_start, feature_end, feature_strand,
                                  feature_type):
    """
    Fix start and end coordinates so that they refer to the position on the genome, rather than the position in the
    fasta sequences given to Cluster-Buster after intersection of a whole genome fasta file with a BED file. Also fix
    the strand information.

    Caveats:
    --------

    1. The locations of the regions are given in 0-based half-open intervals (end position is not included in the
       interval). In contrast, the CRM locations provided by ClusterBuster are provided as 1-based closed intervals.
    2. The strand of a CRM or motif is derived based on the following rules:
          + CRM/motif & + region => + CRM/motif
          - CRM/motif & + region => - CRM/motif
          + CRM/motif & - region => - CRM/motif
          - CRM/motif & - region => + CRM/motif
    3. For sequences located on the negative strand, twoBitToFa (from the UCSC Kent Toolkit) returns the reverse
       complemented nucleotide sequence. The CRM locations provided by ClusterBuster on these negative strands need to
       be corrected. When we start from fasta files and use fastaFromBed (BEDtools), we don't have this problem.
    """

    if twobit is not None:

        if bed_strand == '+':
            fixed_feature_start  = int(bed_start) + (int(feature_start) - 1)
            fixed_feature_end    = int(bed_start) + int(feature_end)
            fixed_feature_strand = str(feature_strand)
        elif bed_strand == '-':
            fixed_feature_start  = int(bed_end) - int(feature_end)
            fixed_feature_end    = int(bed_end) - (int(feature_start) - 1)

            if feature_strand == '+':
                fixed_feature_strand = '-'
            elif feature_strand == '-':
                fixed_feature_strand = '+'
    else:
        fixed_feature_start  = int(bed_start) + (int(feature_start) - 1)
        fixed_feature_end    = int(bed_start) + int(feature_end)

        if feature_type == 'CRM':
            fixed_feature_strand = str(bed_strand)
        else:
            fixed_feature_strand = feature_strand

    return [fixed_feature_start, fixed_feature_end, fixed_feature_strand]


def cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, BEDinput, twobit, TFname, normalize_scores):
    """
    Parse Cluster-Buster output and write BED and/or GFF output to a file or stdout.
    If normalize_scores is True, also normalize the CRM/motif scores for usage with the UCSC browser.
    """

    name = 'cbust'

    # Dictionaries to store the almost final BED and GFF results.
    # Only the CRM/motif score needs to be modified inf normalize_scores is True.
    parsed_output_dict_bed = dict()
    parsed_output_dict_gff = dict()

    parsed_output_dict_bed_key = 0
    parsed_output_dict_gff_key = 0

    # Set minimum and maximum CRM score to values that will be overruled by a real CRM score.
    min_CRM_score = float(999999999999)
    max_CRM_score = float(-999999999999)

    # Set minimum and maximum motif score to values that will be overruled by a real motif score.
    min_motif_score = float(999999999999)
    max_motif_score = float(-999999999999)


    for line in cbust_output:
        line = line.rstrip()

        if line.startswith('CLUSTER '):
            feature_cluster = line.split(" ")[1]
        elif line.startswith('>'):
            feature_seqname = line.split(" ")[0][1:]
        elif line.startswith('Location:'):
            (location, feature_start, to, feature_end) = line.split(" ")
        elif line.startswith('Score:'):
            feature_score = line.split(" ")[1]

            # Now we have all info to make a CRM line.

            # See if we found a new minimum or maximum CRM score.
            min_CRM_score = min(min_CRM_score, float(feature_score))
            max_CRM_score = max(max_CRM_score, float(feature_score))


            if BEDinput is None:
                if fh_bed_output is not None:
                    parsed_output_dict_bed[parsed_output_dict_bed_key] = (
                        'CRM',
                        '\t'.join([str(feature_seqname), str(feature_start), str(feature_end),
                         'CRM[' + TFname + '@' + str(feature_seqname) + ']']),
                        feature_score,
                        '+'
                        )

                    parsed_output_dict_bed_key += 1

                if fh_gff_output is not None:
                    parsed_output_dict_gff[parsed_output_dict_gff_key] = (
                        'CRM',
                        '\t'.join([str(feature_seqname), 'Cbust', 'CRM', str(feature_start), str(feature_end)]),
                        feature_score,
                        '\t'.join(['+', '.', 'cluster "' + str(name) + '-cluster-' + str(feature_cluster) + '"'])
                        )

                    parsed_output_dict_gff_key += 1

            else:
                (bed_chrom, bed_start, bed_end, bed_strand) = bed_input_dict[str(feature_seqname)]


                # Fix start and end coordinates so that they refer to the position on the genome, rather than the
                # position in the fasta sequences given to Cluster-Buster after intersection of a whole genome
                # fasta file with a BED file. Also fix the strand information.
                (fixed_feature_start, fixed_feature_end, fixed_feature_strand) = fix_genomic_coordinates_cbust(
                    twobit=twobit, bed_start=int(bed_start), bed_end=int(bed_end), bed_strand=str(bed_strand),
                    feature_start=int(feature_start), feature_end=int(feature_end), feature_strand='+',
                    feature_type='CRM')


                if fh_bed_output is not None:
                    parsed_output_dict_bed[parsed_output_dict_bed_key] = (
                        'CRM',
                        '\t'.join([str(bed_chrom), str(fixed_feature_start), str(fixed_feature_end),
                                   'CRM[' + TFname + '@' + str(feature_seqname) + ']']),
                        feature_score,
                        str(fixed_feature_strand)
                        )

                    parsed_output_dict_bed_key += 1

                if fh_gff_output is not None:
                    parsed_output_dict_gff[parsed_output_dict_gff_key] = (
                        'CRM',
                        '\t'.join([str(bed_chrom), 'Cbust', 'CRM', str(fixed_feature_start), str(fixed_feature_end)]),
                        feature_score,
                        '\t'.join([fixed_feature_strand, '.', 'seqname "' + str(feature_seqname) + '"; cluster "' +
                                   str(name) + '-cluster-' + feature_cluster + '"'])
                        )

                    parsed_output_dict_gff_key += 1


        else:
            cbust_columns = line.split('\t')

            if len(cbust_columns) == 6:
                # This is a motif line.

                (feature_id, feature_start, feature_end, feature_strand, feature_score, feature_site) = cbust_columns

                # See if we found a new minimum or maximum motif score.
                min_motif_score = min(min_motif_score, float(feature_score))
                max_motif_score = max(max_motif_score, float(feature_score))


                if BEDinput is None:
                    if fh_bed_output is not None:
                        parsed_output_dict_bed[parsed_output_dict_bed_key] = (
                            'motif',
                            '\t'.join([str(feature_seqname), str(feature_start), str(feature_end),
                                       'motif[' + TFname + '@' + str(feature_seqname) + ']']),
                            feature_score,
                            str(feature_strand)
                            )

                        parsed_output_dict_bed_key += 1

                    if fh_gff_output is not None:
                        parsed_output_dict_gff[parsed_output_dict_gff_key] = (
                            'motif',
                            '\t'.join([str(feature_seqname), 'Cbust', 'motif', str(feature_start), str(feature_end)]),
                            feature_score,
                            '\t'.join([feature_strand, '.',
                                       'id "' + feature_id + '"; site "' + feature_site + '"; cluster "'
                                       + str(name) + '-cluster-' + feature_cluster + '"'])
                            )

                        parsed_output_dict_gff_key += 1

                else:
                    # Fix start and end coordinates so that they refer to the position on the genome, rather than the
                    # position in the fasta sequences given to Cluster-Buster after intersection of a whole genome
                    # fasta file with a BED file. Also fix the strand information.
                    (fixed_feature_start, fixed_feature_end, fixed_feature_strand) = fix_genomic_coordinates_cbust(
                        twobit=twobit, bed_start=int(bed_start), bed_end=int(bed_end), bed_strand=str(bed_strand),
                        feature_start=int(feature_start), feature_end=int(feature_end), feature_strand=feature_strand,
                        feature_type='motif')


                    if fh_bed_output is not None:
                        parsed_output_dict_bed[parsed_output_dict_bed_key] = (
                            'motif',
                            '\t'.join([str(bed_chrom), str(fixed_feature_start), str(fixed_feature_end),
                                       'motif[' + TFname + '@' + str(feature_seqname) + ']']),
                            feature_score,
                            str(fixed_feature_strand))

                        parsed_output_dict_bed_key += 1


                    if fh_gff_output is not None:
                        parsed_output_dict_gff[parsed_output_dict_gff_key] = (
                            'motif',
                            '\t'.join([str(bed_chrom), 'Cbust', 'motif', str(fixed_feature_start), str(fixed_feature_end)]),
                            feature_score,
                            '\t'.join([fixed_feature_strand, '.',
                                       'seqname "' + str(feature_seqname) + '"; id "' + feature_id + '"; site "'
                                       + feature_site + '"; cluster "' + str(name) + '-cluster-' + feature_cluster + '"'])
                            )

                        parsed_output_dict_gff_key += 1



    if fh_bed_output is not None:
        # Write BED file.

        i = 0
        
        while i < parsed_output_dict_bed_key:
            (feature, bed_part1, feature_score, bed_part2) = parsed_output_dict_bed[i]
            
            if normalize_scores is True:
                # If normalize_scores is True, also normalize the CRM/motif scores for usage with the UCSC browser.
                if feature == 'CRM':
                    print >> fh_bed_output, '\t'.join([bed_part1,
                                                       str((float(feature_score)-min_CRM_score)/(max_CRM_score-min_CRM_score) * 1000.0),
                                                       bed_part2])
                else:
                    print >> fh_bed_output, '\t'.join([bed_part1,
                                                       str((float(feature_score)-min_motif_score)/(max_motif_score-min_motif_score) * 1000.0),
                                                       bed_part2])
            else:
                # Write original CRM/motif scores in the BED file.
                print >> fh_bed_output, '\t'.join([bed_part1, str(feature_score), bed_part2])
            
            i += 1

    if fh_gff_output is not None:
        # Write GFF file.

        i = 0

        while i < parsed_output_dict_gff_key:
            (feature, gff_part1, feature_score, gff_part2) = parsed_output_dict_gff[i]

            if normalize_scores is True:
                # If normalize_scores is True, also normalize the CRM/motif scores for usage with the UCSC browser.
                if feature == 'CRM':
                    print >> fh_gff_output, '\t'.join([gff_part1,
                                                       str((float(feature_score)-min_CRM_score)/(max_CRM_score-min_CRM_score) * 1000.0),
                                                       gff_part2])
                else:
                    print >> fh_gff_output, '\t'.join([gff_part1,
                                                       str((float(feature_score)-min_motif_score)/(max_motif_score-min_motif_score) * 1000.0),
                                                       gff_part2])
            else:
                # Write original CRM/motif scores in the GFF file.
                print >> fh_gff_output, '\t'.join([gff_part1, str(feature_score), gff_part2])
            
            i += 1
    
            
    return None


def main():
    """
    cbust.py main function.
    """

    # Ignore SIG_PIPE and don't throw exceptions on it:
    #   http://docs.python.org/library/signal.html
    signal(SIGPIPE, SIG_DFL)

    # Build command line options.
    fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1 )

    usage = "Usage: %prog [options]"

    parser = optparse.OptionParser(usage=usage, version="%prog v1.0", formatter=fmt)

    parser.add_option("-f", "--fasta",
        action="store", type="string",
        dest="fastaFile", help="input FASTA file")
    parser.add_option("-2", "--2bit",
        action="store", type="string",
        dest="twobitFile", help="input 2bit file")
    parser.add_option("-b", "--bed",
        action="store", type="string",
        dest="BEDinput", help="input BED file")
    parser.add_option("-C", "--cb",
        action="store", type="string",
        dest="cbFile", help="input Cluster-Buster matrix motif file")
    parser.add_option("-B", "--bedoutput",
        action="store", type="string",
        dest="BEDoutput", help="output BED file")
    parser.add_option("-G", "--gffoutput",
        action="store", type="string",
        dest="GFFoutput", help="output GFF file")
    parser.add_option("-t", "--tfname",
        action="store", type="string",
        dest="TFname", help="transcription factor name (default = Cluster-Buster matrix motif filename)")
    parser.add_option("-n", "--normalize-score",
        action="store_true", default=False,
        dest="normalize_scores", help="Normalize scores between range of 0 to 1000 (useful for viewing in UCSC genome browser)")
    parser.add_option("-c", "--clusterscorethreshold",
        action="store", type="int", default=5,
        dest="clusterscorethreshold", help="Cluster score threshold (default = 5)")
    parser.add_option("-m", "--motifscorethreshold",
        action="store", type="int", default=6,
        dest="motifscorethreshold", help="Motif score threshold (default = 6)")
    parser.add_option("-g", "--gap",
        action="store", type="int", default=35,
        dest="gap", help="Expected gap (bp) between neighboring motifs in a cluster (default = 35)")
    parser.add_option("-r", "--range",
        action="store", type="int", default=100,
        dest="range", help="Range in bp for counting local nucleotide abundances (default = 100)")
    parser.add_option("-l", "--maskrepeat",
        action="store_true", default=False,
        dest="maskrepeat", help="Mask lowercase letters in the sequences")
    parser.add_option("-p", "--pseudocount",
        action="store", type="float", default=0.375,
        dest="pseudocount", help="Pseudocount (default = 0.375)")


    (options, args) = parser.parse_args()


    # Check if a FASTA or twobit file is given as input.
    if ( (options.fastaFile is None) and (options.twobitFile is None) ):
        parser.print_help()
        print >> sys.stderr, '\nNo FASTA file (--fasta, -f) or 2bit file (--2bit, -2) specified.\n'
        sys.exit(1)
    elif (options.twobitFile is not None):
        if not os.path.isfile(options.twobitFile):
            print >> sys.stderr, '\n2bit file "{0:s}" could not be found.\n'.format(options.twobitFile)
            sys.exit(1)
    elif (options.fastaFile is not None):
        if not os.path.isfile(options.fastaFile):
            print >> sys.stderr, '\nFASTA file "{0:s}" could not be found.\n'.format(options.fastaFile)
            sys.exit(1)


    # Check if a Cluster-Buster matrix motif file is given.
    if (options.cbFile is None):
        parser.print_help()
        print >> sys.stderr, "\nNo Cluster-Buster matrix motif file specified (--cb, -i).\n"
        sys.exit(1)
    elif (not os.path.isfile(options.cbFile)):
        print >> sys.stderr, '\nCluster Buster matrixfile "{0:s}" could not be found.\n'.format(options.cbFile)
        sys.exit(1)
    elif options.TFname is not None:
        # Set TFname to specified value.
        TFname = options.TFname
    else:
        # If no TFname is specified, use the Cluster-Buster matrix motif filename.
        TFname = os.path.basename(options.cbFile)

        # If the extention is ".cb", remove it.
        if TFname.endswith('.cb'):
            TFname = TFname[:-3]


    # Check that at least one output format (BED and/or GFF) is chosen.
    if options.BEDoutput is None and options.GFFoutput is None:
        parser.print_help()
        print >> sys.stderr, '\nChoose at least one output format (BED and/or GFF).\n'
        sys.exit(1)


    # Look for Cluster Buster.
    prg_in_path('cbust')

    bed_input_dict = dict()

    # Check if there is a BED file specified.
    if (options.BEDinput is not None):
        if not os.path.isfile(options.BEDinput):
            print >> sys.stderr, '\nBED file "{0:s}" could not be found.\n'.format(options.BEDinput)
            sys.exit(1)

        with open(options.BEDinput, 'r') as fh_bed_input:
            for bed_line in fh_bed_input:
                if bed_line.startswith('#'): continue

                bed_columns = bed_line.rstrip().split('\t')

                nbr_bed_columns = len(bed_columns)

                if nbr_bed_columns < 4:
                    print >> sys.stderr, '\nBED file "{0:s}" has less than 4 columns per line.\n'.format(
                        options.BEDinput)
                    sys.exit(1)
                elif nbr_bed_columns <= 5:
                    (bed_chrom, bed_start, bed_end, bed_name) = bed_columns[0:4]
                    bed_strand = '+'
                else:
                    (bed_chrom, bed_start, bed_end, bed_name, bed_score, bed_strand) = bed_columns[0:6]


                if bed_input_dict.get(bed_name) is not None:
                    print >> sys.stderr, '\nThe name (4th) column of the BED file must have unique values.\n'
                    sys.exit(1)

                bed_input_dict[bed_name] = (bed_chrom, bed_start, bed_end, bed_strand)

        fh_cbust_fasta = tempfile.NamedTemporaryFile(prefix='cbust-fasta-', dir='/tmp')
        cbust_fasta = fh_cbust_fasta.name

        if (options.twobitFile is not None):
            # Look for twoBitToFa (Kent tools).
            prg_in_path('twoBitToFa')

            # Create temporary BED file for usage with twoBitToFa:
            #   - needs at least 6 columns, so twoBitToFa knows the strand.
            #   - the score needs to be a integer, but is not really used so we can use '0'.
            fh_twobitofa_bed = tempfile.NamedTemporaryFile(prefix='twoBitToFa-bed-', dir='/tmp')
            twobitofa_bed = fh_twobitofa_bed.name

            for bed_name, bed_input_values in bed_input_dict.iteritems():
                (bed_chrom, bed_start, bed_end, bed_strand) = bed_input_values
                print >> fh_twobitofa_bed, '{0:s}\t{1:s}\t{2:s}\t{3:s}\t0\t{4:s}\n'.format(bed_chrom, bed_start, bed_end, bed_name, bed_strand)

            # Flush the BED file for twoBitToFa, so twoBitToFa is able to read the whole file.
            fh_twobitofa_bed.flush()

            twoBitToFa_cmd = ['twoBitToFa', '-bed=' + twobitofa_bed, options.twobitFile, cbust_fasta]

            run_cmd(twoBitToFa_cmd)
        elif (options.fastaFile is not None):
            # Look for fastaFromBed (BEDtools).
            prg_in_path('fastaFromBed')

            fastaFromBed_cmd = ['fastaFromBed', '-name', '-fi', options.fastaFile, '-bed', options.BEDinput, '-fo', cbust_fasta]

            run_cmd(fastaFromBed_cmd)
    elif (options.twobitFile is not None):
        # Convert whole 2bit file to FASTA.

        # Look for twoBitToFa (Kent tools).
        prg_in_path('twoBitToFa')

        fh_cbust_fasta = tempfile.NamedTemporaryFile(prefix='cbust-fasta-', dir='/tmp')
        cbust_fasta = fh_cbust_fasta.name

        twoBitToFa_cmd = ['twoBitToFa', options.twobitFile, cbust_fasta]

        run_cmd(twoBitToFa_cmd)
    else:
        # Use provided FASTA file as it is.
        cbust_fasta = options.fastaFile



    # Build Cluster-Buster commandline and filter out "None" values (introduced when the "maskrepeat" option was not set).
    cbust_cmd = filter(None,
                       [ 'cbust',
                         '-f0',
                         '-c{0:f}'.format(options.clusterscorethreshold),
                         '-m{0:f}'.format(options.motifscorethreshold),
                         '-g{0:d}'.format(options.gap),
                         '-r{0:d}'.format(options.range),
                         '-p{0:f}'.format(options.pseudocount),
                         '-l' if options.maskrepeat is True else None,
                         options.cbFile,
                         cbust_fasta
                       ])

    # Run Cluster-Buster.
    (cbust_stdoutdata, cbust_stderrdata) = run_cmd(cbust_cmd)


    # Make an array (each element is one one) of the Cluster-Buster output.
    cbust_output = cbust_stdoutdata.split('\n')


    # Parse Cluster-Buster output and write BED and/or GFF output to a file or stdout.

    fh_bed_output = None
    fh_gff_output = None

    if options.BEDoutput is not None:
        # Output the BED file content.
        if options.BEDoutput == '-' or options.BEDoutput == 'stdout':
            # BED file output to stdout.
            fh_bed_output = sys.stdout

            if options.GFFoutput is not None:
                # Output the GFF file content.
                if options.GFFoutput == '-' or options.GFFoutput == 'stdout':
                    # GFF file output to stdout.
                    fh_gff_output = sys.stdout
                    cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
                else:
                    # GFF file output to file.
                    with open(options.GFFoutput, 'w') as fh_gff_output:
                        cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
            else:
                # Don't output the GFF file content.
                cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
        else:
            # BED file output to file.
            if options.GFFoutput is not None:
                if options.GFFoutput == '-' or options.GFFoutput == 'stdout':
                    # GFF file output to stdout.
                    fh_gff_output = sys.stdout

                    with open(options.BEDoutput, 'w') as fh_bed_output:
                        cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
                else:
                    # GFF file output to file.
                    with open(options.BEDoutput, 'w') as fh_bed_output, open(options.GFFoutput, 'w') as fh_gff_output:
                        cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
            else:
                # Don't output the GFF file content.
                with open(options.BEDoutput, 'w') as fh_bed_output:
                    cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
    else:
        # Don't output the BED file content.
        if options.GFFoutput is not None:
            # Output the GFF file content.
            if options.GFFoutput == '-' or options.GFFoutput == 'stdout':
                # GFF file output to stdout.
                fh_gff_output = sys.stdout
                cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)
            else:
                # GFF file output to file.
                with open(options.GFFoutput, 'w') as fh_gff_output:
                    cbust_parse_output(fh_bed_output, fh_gff_output, cbust_output, bed_input_dict, options.BEDinput, options.twobitFile, TFname, options.normalize_scores)



    sys.exit(0)


if __name__ == "__main__":
    main()
