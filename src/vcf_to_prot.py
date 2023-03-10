"""
    This script generates a protein fasta file out of a variant call format (vcf) file.
    Author: Georges BEDRAN, gbadran_90@live.com

"""
import argparse
import tempfile
import warnings

from pyteomics import parser

import numpy as np
import pandas as pd
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=FutureWarning)

CONV = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': '', '.': ''}

AA = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V', ''
]


def splitDataFrameList(df, target_column, separator):
    """
        splits a DataFrame by a separator on a specific column and adds the split
        to new duplicated rows.

        Arguments:
        ----------
            - df: (pd.DataFrame) dataframe to split,
            - target_column: (str) the column containing the values to split
            - separator: (Str) the symbol used to perform the split

        Return:
        -------
            - pd.DataFrame with each entry for the target column separated, with
              each element moved into a new row. The values in the other columns
              are duplicated across the newly divided rows.
    """

    def splitListToRows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(
        splitListToRows, axis=1, args=(new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


def read_gtf(gtf_file, filter_feature="CDS"):
    """
        Reads GTF and returns a double dict with the following structure:
        double_dict[transcript_id][exon_number]=[start, end, frame].

        Arguments:
        ----------
            - gtf_file: (File handle) to a gtf file
            - filter_feature: (str) a feature to filter the gft on (3rd field in
                            a GTF file)

        Return:
        ------
            - dict
    """
    desc = {}
    for line in gtf_file:
        if line.startswith('#'):
            continue
        line = line.strip('\n')
        splitted = line.split('\t')
        if filter_feature == splitted[2]:
            desc_dict = {}
            start = splitted[3]
            end = splitted[4]
            frame = splitted[7]
            desc_splitted = splitted[8].replace('"', '').replace(
                '; ', ';').split(';')
            for element in desc_splitted:
                element = element.split(' ')
                try:
                    desc_dict[element[0]] = element[1]
                except IndexError:
                    continue
            if desc_dict['transcript_id'] not in desc:
                desc[desc_dict['transcript_id']] = {}
            desc[desc_dict['transcript_id']][int(desc_dict['exon_number'])] = [
                int(start), int(end), int(frame)
            ]
    return desc


def prepare_vcf(lvcf_path):
    """
        Prepares a vcf for pybedtools
    """
    vcf_df = pd.read_csv(
        lvcf_path,
        usecols=[0, 1, 2, 3, 4],
        comment='#',
        sep='\t',
        header=None,
        low_memory=False)
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    vcf_df.CHROM = vcf_df['CHROM'].str.replace('chr', '')
    vcf_df = splitDataFrameList(vcf_df, 'ALT', ',')
    vcf_df.drop_duplicates(inplace=True)
    vcf_df.reset_index(inplace=True, drop=True)
    cols = ['CHROM', 'POS', 'POS', 'ID', 'REF', 'ALT']
    vcf = pybedtools.BedTool.from_dataframe(vcf_df[cols])
    return vcf


def prepare_gtf(lgtf_path):
    """
        Prepares a gtf for pybedtools
    """
    gtf = pybedtools.BedTool(lgtf_path)  # GTf file

    gtf_df = gtf.to_dataframe(comment='#', low_memory=False)
    # gtf_df.start = gtf_df.start.map(int) - 1
    gtf_df = gtf_df[gtf_df.feature == 'CDS'][[
        'seqname', 'start', 'end', 'strand', 'frame', 'attributes'
    ]]
    gtf = pybedtools.BedTool.from_dataframe(gtf_df)
    return gtf


def get_cds_pos(exons, exon_number, row):
    """
        Calculates a CDS position given a dict of exons, the exon number
        within the transcript, the transcript strand (+ or -) and a position
        within the exon number.

        Arguments:
        ----------
            - exons: (dict) dict[int]= [start, end, frame]
            - exon_number: (int) the number of the exon where the mutation takes place
            - row (pd.DataFrame row)

        Returns:
        --------
            - cds_pos (int)
    """
    cds_pos = 0
    if min(exons.keys()) < exon_number:
        for index in range(min(exons.keys()), exon_number):
            cds_pos = cds_pos + exons[index][1] - exons[index][0] + 1
    if row.strand == '+':
        # genomic pos  - start pos + 1 + frame
        cds_pos = cds_pos + int(row.pos) - exons[exon_number][0]
    else:
        # end pos - genomic pos - frame
        cds_pos = cds_pos + exons[exon_number][1] - int(row.pos)
    return cds_pos


def introduce_mut2cds(row, cds_pos, str_seq, reversed_=False):
    """
        Takes a string (CDS) and introduces a mutation given a row
        (pd.DataFrame row) containing a strand (+ or -) and alt (alternative seq)
        and a cds_pos (int): the position in the string to introduce the mutation.

        Arguments:
        ---------
            - row (pd.DataFrame row)
            - cds_pos (int)
            - str_seq (str)
            - reversed_: (bool) True if row.ref has been converted to it's
                       reverse complement equivalent. False if not.

        Return:
        -------
            - (str) the mutated sequence
    """
    if not reversed_ and row.strand == "-":
        row.ref = ''.join([CONV[x] for x in row.ref[::-1]])
        row.alt = ''.join([CONV[x] for x in row.alt[::-1]])

    # in case deletion
    if len(row.ref) - len(row.alt) > 0:
        if row.strand == '+':
            str_seq = str_seq[:cds_pos] + row.alt + str_seq[cds_pos +
                                                            len(row.ref):]
        else:
            str_seq = str_seq[:cds_pos - len(row.ref) +
                              1] + row.alt + str_seq[cds_pos + 1:]
    else:
        # in case insertion and substitution
        str_seq = str_seq[:cds_pos] + row.alt + str_seq[cds_pos + 1:]
    return str_seq


def get_aminoacid_ref_pos_alt(row, str_seq, str_seq_mut, cds_pos, to_stop):
    """
        Returns the reference amino acid sequence on the variation position,
        the variation position and the alternative amino acid sequence after
        introducing the mutation (protein level).

        Arguments:
        ----------
            - row (pd.DataFrame row)
            - str_seq (str) normal CDS sequence
            - str_seq_mut (str) mutated CDS sequence
            - cds_pos (int) position of the mutation start within the cds seq

        Return:
        ------
            - (aa_ref, aa_pos, aa_alt): tuple (str, int, str)
    """
    # if 6%3==0 ==> first base of the codon is 6-2
    first_base_codon = {0: 2, 2: 1, 1: 0}
    # consider if an indel without fs is intra-codon or inter-codon
    # used to consider reporting an extra amino acid in the header
    last_base_codon = {0: 3, 2: 3, 1: 0}

    aa_frame = (cds_pos + 1) % 3
    first_base_pos = cds_pos - first_base_codon[aa_frame]
    codon_ref = str_seq[first_base_pos:first_base_pos + 3]
    aa_ref = str(Seq(codon_ref).translate(to_stop=to_stop))
    aa_pos = (first_base_pos + 3) / 3
    # in case substitution or indel without frameshit
    if abs(len(row.ref) - len(row.alt)) % 3 == 0:
        codon_alt = str_seq_mut[first_base_pos:first_base_pos + 3]
        aa_alt = str(Seq(codon_alt).translate())
        if len(row.ref) - len(row.alt) > 0:
            if row.strand == '+':
                aa_frame = (cds_pos + 2) % 3
                first_base_pos = cds_pos + 1 - first_base_codon[aa_frame]
                aa_pos = (first_base_pos + 3) / 3
                codon_ref = str_seq[first_base_pos:first_base_pos +
                                    len(row.ref) - 1 +
                                    last_base_codon[aa_frame]]
                aa_ref = str(Seq(codon_ref).translate())
                aa_alt = 'del'
            else:
                aa_frame = (cds_pos + 1 - len(row.ref) + 1) % 3
                first_base_pos = cds_pos - (
                    len(row.ref) - 1) - first_base_codon[aa_frame]
                aa_pos = (first_base_pos + 3) / 3
                codon_ref = str_seq[first_base_pos:first_base_pos +
                                    len(row.ref) - 1 +
                                    last_base_codon[aa_frame]]
                aa_ref = str(Seq(codon_ref).translate())
                aa_alt = 'del'
        elif len(row.ref) - len(row.alt) < 0:
            aa_frame = (cds_pos + 2) % 3
            first_base_pos = cds_pos + 1 - first_base_codon[aa_frame]
            codon_alt = str_seq_mut[first_base_pos:first_base_pos +
                                    len(row.alt) - 1 +
                                    last_base_codon[aa_frame]]
            aa_alt = str(Seq(codon_alt).translate()) + '-ins'
    else:
        # in case ins with frameshift
        if len(row.alt) > len(row.ref):
            aa_alt = 'ins-fs'
        # in case del with frameshift
        else:
            aa_alt = 'del-fs'
    return (aa_ref, aa_pos, aa_alt)


def list_to_dict(list_):
    # list_ = list of string, pair length
    dict_ = {}
    for index in range(0, len(list_), 2):
        try:
            dict_[list_[index]] = list_[index + 1].replace('"', '')
        except IndexError:
            pass
    return (dict_)


def genomic_to_cds_coords(lintersect_df,
                          transcript_dict,
                          CDS_DB_file,
                          loutput_file,
                          include_id=False,
                          to_stop=True):
    """
        Converts genomic coordinates to CDS coordinates.

        Arguments:
        ---------
            - lintersect_df: (pd.DataFrame) with the columns 'chr', 'pos',
              'id', 'ref', 'alt', 'strand', 'frame', 'ann'. The column ann is
              the "attribute"
                             field of the GTF.
            - transcript_dict: (dict) double_dict[transcript_id][exon_number]=[start, end, frame]
                                output of read_gtf
            - CDS_DB_file: (str) path to CDS fasta
            - loutput_file = (str) path to the output_fasta
            - include_id : (bool) whether to include the vcf id value in the
                           header or not

        Return:
        -------
            - void
    """
    # removing the transcript version from the cds file
    warnings = {}
    output = open(loutput_file, 'w')
    tmp_cds_db = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(CDS_DB_file) as fh:
        for line in fh:
            if line.startswith('>'):
                if '.' in line.split(' ', 1)[0]:
                    line = line.split('.', 1)[0] + '\n'
                tmp_cds_db.write(line)
            else:
                tmp_cds_db.write(line)
    tmp_cds_db.close()
    CDS_DB = SeqIO.index(tmp_cds_db.name, 'fasta')
    lintersect_df.ann = lintersect_df.ann.map(lambda x: x[:-1])
    lintersect_df.ann = lintersect_df.ann.str.split('; | ').map(
        lambda x: list_to_dict(x))
    lintersect_df.reset_index(inplace=True, drop=True)
    tot = lintersect_df.shape[0]
    for index_, row in tqdm(
            lintersect_df.iterrows(), desc='Processing mutation', total=tot):
        # assign some variable for readability
        protein_id = row.ann.get('protein_id', '')
        transcript_name = row.ann.get('transcript_name', 'no-trans-name')
        gene_id = row.ann.get('gene_id', 'no-gene-id')
        ref = row.ref
        alt = row.alt
        transcript_id = row.ann.get('transcript_id', 'No-Trans-id')
        transcript_version = row.ann.get('transcript_version', '')
        try:
            exon_number = int(row.ann['exon_number'])
        except KeyError:
            print('KeyError for row number {}'.format(index_))
            continue
        exons = transcript_dict[transcript_id]
        # convert genomic position to CDS position
        cds_pos = get_cds_pos(exons, exon_number, row)
        # if - strand get the reverse complement for ref and alt
        if row.strand == "-":
            row.ref = ''.join([CONV[x] for x in row.ref[::-1]])
            row.alt = ''.join([CONV[x] for x in row.alt[::-1]])
        if row.alt == '.':
            row.alt = ''
        try:
            # the transcript versions between the GTF and CDS fasta aren't all
            # matching ==> ignore the transcript version.
            # seq=CDS_DB[transcript_id +'.' + transcript_version].seq
            seq = CDS_DB[transcript_id].seq
        except KeyError:
            print(transcript_id)
            continue
        str_seq = str(seq)

        # accounting for N letters at the beginning of the CDS
        # to correct CDS pos
        N = 0
        for base in str_seq:
            if base == 'N':
                N += 1
            else:
                break
        cds_pos += N

        # checking if everything is alright by extracting the CDS sequence and
        # comparing it to the ref sequence in the VCF.
        # if they match ==> the cds_pos calculation is correct.
        # note: be careful in case of mutations spanning exon-intron junctions
        if str_seq[cds_pos:cds_pos +
                   len(row.ref)] != row.ref and row.strand == '+':
            warning = "Mutation spanning exon junction ignored "
            warning += f"{row.chr} {row.pos} {row.ref} {row.alt} strand: {row.strand}"
            warnings[warning] = 1
            continue
        if str_seq[cds_pos - len(row.ref) + 1:cds_pos +
                   1] != row.ref and row.strand == '-':
            warning = "Mutation spanning exon junction ignored "
            warning += f"{row.chr} {row.pos} {row.ref} {row.alt} strand: {row.strand}"
            warnings[warning] = 1
            continue
        # introducing the mutation to the CDS sequence
        str_seq_mut = introduce_mut2cds(row, cds_pos, str_seq, reversed_=True)
        # generating the ref, position and alt amino acid for the fasta header
        aa_ref, aa_pos, aa_alt = get_aminoacid_ref_pos_alt(
            row, str_seq, str_seq_mut, cds_pos, to_stop=to_stop)

        seq = Seq(str_seq_mut)
        # ignoring synonymous mutations and translate to protein
        if aa_ref == aa_alt:
            continue
        try:
            if (len(row.ref) - len(row.alt)) % 3 == 0:
                protein_seq = str(seq.translate())
            else:
                protein_seq = str(seq.translate(to_stop=to_stop))
        except:
            warning = "mutation %s %s %s %s translation failed "
            warning += f"{row.chr} {row.pos} {row.ref} {row.alt} strand: {row.strand}"
            warnings[warning] = 1
            continue
        if include_id is True and row.id != '.':
            header = str(row.id) + '_' + protein_id
        else:
            header = protein_id
        header = header + '_c.' + str(cds_pos) + ':' + row.ref + '>' + row.alt
        header = header + '_p.' + aa_ref + str(int(aa_pos)) + aa_alt
        header = header + ' ' + transcript_id + '.' + transcript_version
        header = header + '|' + transcript_name + '|' + gene_id
        header = header + '|' + row.strand + '|'
        header = header + 'chr' + str(row.chr) + ':g.' + str(row.pos)
        header = header + ':' + ref + '>' + alt
        entry = ">%s\n%s\n" % (header, protein_seq)
        output.write(entry)
    output.close()
    for k in warnings.keys():
        print(k)


def stop_codon_permutation(protein_fasta,
                           output_fasta,
                           aa_list,
                           window_size=21,
                           cleavage=True,
                           **kwards):
    database = SeqIO.index(protein_fasta, 'fasta')
    output = open(output_fasta, "w")
    for key, value in database.items():
        seq = str(value.seq).strip('*')
        position = None
        if '*' in seq:
            positions = np.where([x == '*' for x in seq])[0]
            try:
                if positions[1] - positions[0] < (window_size - 1) / 2:
                    seq = seq[:positions[1]]
            except IndexError:
                pass
            position = positions[0]
            if position - (window_size - 1) / 2 + 1 < 0:
                start = 0
            else:
                start = int(position - (window_size - 1) / 2)
            if len(seq) - position + 1 < (window_size - 1) / 2:
                end = len(seq)
            else:
                end = int(position + (window_size - 1) / 2 + 1)
            new_seqs = []
            cleaved_towrite = []
            for aa in aa_list:
                new_seq = seq[start:position] + aa + seq[position + 1:end]
                new_pos = position - start
                header = '>' + f'{aa}_{position+1}:{start+1}-{end}_' + key
                new_seqs.append(header)
                new_seqs.append(new_seq)
                if cleavage:
                    cleaved_seqs = unspecific_cleavage(new_seq, **kwards)
                    for c_seq in cleaved_seqs:
                        c_start = new_seq.index(c_seq)
                        c_end = c_start + len(c_seq)
                        if c_start <= new_pos < c_end:
                            if aa == '':
                                cheader = '>' + f'skip_{position+1}:{c_start + 1}-{c_end}_' + key
                            else:
                                cheader = '>' + f'{aa}_{position+1}:{c_start + 1}-{c_end}_' + key
                            cleaved_towrite.append(cheader)
                            cleaved_towrite.append(c_seq)
            # output.write('\n'.join(new_seqs) + '\n')
            if cleavage:
                output.write('\n'.join(cleaved_towrite) + '\n')
    output.close()


def unspecific_cleavage(sequence,
                        num_missed_cleavage=11,
                        peptide_min_length=8,
                        rule="([^\\*])"):
    """
        Cleaves protein sequences into a list of unique peptides.
        for more about pyteomics.parser.cleave and cleavage rules:
        <https://pythonhosted.org/pyteomics/api/parser.html>
    """
    parser.expasy_rules['new_enzyme'] = rule
    peptide_set = parser.cleave(
        sequence=sequence,
        rule=parser.expasy_rules['new_enzyme'],
        missed_cleavages=num_missed_cleavage,
        min_length=peptide_min_length)
    return list(peptide_set)


def main(arguments):
    gtf_path = arguments.input_gtf[0]  # GTf file
    vcf_path = arguments.input_vcf[0]  # VCF file
    CDS_DB_file = arguments.input_cds[0]  # CDS fasta
    output_file = arguments.output[0]  # output fasta
    transcript_dict = read_gtf(open(gtf_path))
    py_vcf = prepare_vcf(vcf_path)
    py_gtf = prepare_gtf(gtf_path)
    intersect = py_vcf.intersect(py_gtf, wb=True)
    # intersect_df = intersect.to_dataframe(header=-1, low_memory=False)
    intersect_df = intersect.to_dataframe(
        header=None, disable_auto_names=True, low_memory=False)

    intersect_df = intersect_df[[0, 1, 3, 4, 5, 9, 10, 11]]
    intersect_df.columns = [
        'chr', 'pos', 'id', 'ref', 'alt', 'strand', 'frame', 'ann'
    ]
    genomic_to_cds_coords(
        intersect_df,
        transcript_dict,
        CDS_DB_file,
        output_file,
        include_id=True)


def main_nmd(arguments):
    gtf_path = arguments.input_gtf[0]  # GTf file
    vcf_path = arguments.input_vcf[0]  # VCF file
    CDS_DB_file = arguments.input_cds[0]  # CDS fasta
    output_file = arguments.output[0]  # output fasta
    transcript_dict = read_gtf(open(gtf_path))
    py_vcf = prepare_vcf(vcf_path)
    py_gtf = prepare_gtf(gtf_path)
    intersect = py_vcf.intersect(py_gtf, wb=True)
    # intersect_df = intersect.to_dataframe(header=-1, low_memory=False)
    intersect_df = intersect.to_dataframe(
        header=None, disable_auto_names=True, low_memory=False)

    intersect_df = intersect_df[[0, 1, 3, 4, 5, 9, 10, 11]]
    intersect_df.columns = [
        'chr', 'pos', 'id', 'ref', 'alt', 'strand', 'frame', 'ann'
    ]
    genomic_to_cds_coords(
        intersect_df,
        transcript_dict,
        CDS_DB_file,
        output_file,
        include_id=True,
        to_stop=False)
    perm_output = '/'.join(
        output_file.split('/')[:-1]) + 'perumutated_cleaved_'
    perm_output += output_file.split('/')[-1]
    stop_codon_permutation(output_file, perm_output, AA)


def args():
    parser = argparse.ArgumentParser(
        description=
        'generates a protein fasta file out of a variant call format (vcf) file'
    )
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-input_gtf',
        nargs=1,
        type=str,
        help='Path to input GTF fasta',
        required=True)
    requiredNamed.add_argument(
        '-input_cds',
        nargs=1,
        type=str,
        help='Path to input cds fasta',
        required=True)

    requiredNamed.add_argument(
        '-input_vcf',
        nargs=1,
        type=str,
        help='Path to input vcf variant file',
        required=True)

    requiredNamed.add_argument(
        '-permutation',
        nargs=1,
        type=int,
        choices=[0, 1],
        help='generate aa permutations for novel stop codons',
        required=True)

    requiredNamed.add_argument(
        '-output',
        nargs=1,
        type=str,
        help='Path to output fasta',
        required=True)
    return parser


if __name__ == '__main__':
    arguments = args().parse_args()
    permutation = arguments.permutation[0]

    if permutation:
        main_nmd(arguments)
    else:
        main(arguments)
