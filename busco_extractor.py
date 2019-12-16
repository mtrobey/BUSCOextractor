import glob
import os
import sys
import subprocess
from typing import Dict
from pathlib import Path
from collections import defaultdict
import argparse

from Bio import SeqIO

def _parse_arguments() -> argparse.Namespace:

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_dir', type=str, default='./input',
        help="""Directory containing protein fasta files""")
    parser.add_argument('-hmm', '--hmm_dir', type=str, default='./hmms',
        help="""Directory containing BUSCO HMMs""")
    parser.add_argument('-o', '--output_dir', type=str, default='./output',
        help="""Directory that will contain output""")
    parser.add_argument('-c', '--cores', type=int, default=os.cpu_count(),
        help="""Number of CPUs for protein alignment""")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
        help="""Displays debug messages""")
    parser.add_argument('-e', '--e_value_min', type=float, default=1E-10,
        help="""Min E value for hmmsearch results""")
    parser.add_argument('-l', '--seq_length_min', type=float, default=100,
        help="""Min protein sequence length for hmmsearch results""")
    parser.add_argument('-g', '--gap_threshold', type=float, default=0.9,
        help="""Gap threshold for trimal""")
    args = parser.parse_args()

    return args


def hmmsearch(faa: str, hmm_dir: str, output_dir: str, verbose: bool, cores: int, e_value_min: float, seq_length_min: int, use_cache: bool = True) -> None:
    """Extracts BUSCO genes using hmmsearch

        args:
            faa (str): path to genome protein file
            hmm_dir (str): path to HMM directory
            output_dir (str): output directory
            verbose (bool): prints debugging messages
            cores (int): number of CPU cores
            e_value_min (float): min E value for hmmsearch result
            seq_length_min (int): min protein length for hmmsearch result
        
        returns:
            None

    """
    missed = []
    for hmm in glob.glob(os.path.join(hmm_dir, '*.hmm')):
        faa_name = os.path.splitext(os.path.basename(faa))[0]
        
        hmm_name = os.path.splitext(os.path.basename(hmm))[0]
        domtable_name = f'{os.path.join(output_dir, hmm_name, faa_name)}.domtable' 
        busco_faa_name = os.path.join(output_dir, hmm_name, f'{faa_name}.faa')
        if not os.path.exists(os.path.join(output_dir, hmm_name)):
            os.mkdir(os.path.join(output_dir, hmm_name))
        hmmsearch_cmd = f'hmmsearch --domtblout {domtable_name} --noali --cpu {cores} {hmm} {faa}'
        if use_cache and os.path.exists(domtable_name) and os.path.exists(busco_faa_name):
            continue
        if verbose:
            print(hmmsearch_cmd)
        subprocess.check_output(hmmsearch_cmd, shell=True)
        proteins = []
        with open(domtable_name, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    row = line.split(' ')
                    row = [i for i in row if i]
                    try:
                        start = int(row[19])
                        end = int(row[20])
                        e_value = float(row[6])
                        if e_value < e_value_min and ((end-start) > seq_length_min):
                            proteins.append({
                                'name': row[0],
                                'e-value': e_value,
                                'start': start,
                                'end': end,
                            })
                    except ValueError:
                        sys.exit(f'Error parsing hmmsearch domtable row {row}') 
        if len(proteins) == 0:
            missed.append(os.path.basename(hmm))
            protein_seq = '-'
        else:
            proteins = sorted(proteins, key=lambda x: x['e-value'])
            protein_seq = get_protein(proteins[0]['name'], faa, proteins[0]['start'], proteins[0]['end'])
            if not protein_seq:
                print(f'Warning: error retrieving protein sequence {faa_name} {hmm}.')
                continue
        with open(busco_faa_name, 'w') as f:
            fasta_header = '>{}'.format(faa_name)
            f.write('\n'.join([fasta_header, str(protein_seq)]))
    
    with open(os.path.join(output_dir, faa_name+'_missed.txt'), 'w') as f:
        f.write('\n'.join(missed))


def get_protein(protein_name: str, faa: str, start: int, end: int):
    """Retrieves a protein sequence by fasta header
    
        args:
        protein_name (str): name of protein fasta header
        faa (str): path to protein fasta file
        start (int): start coordinate of protein to extract
        end (int): end coordinate of protein to extract
        
        returns:
            Bio.Seq protein object
            
    """
    if not os.path.exists(faa):
        sys.exit(f'Unable to find protein fasta file {faa} after hmmsearch')
    with open(faa, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.name == protein_name and record.seq[start:end]:
                return record.seq[start:end]
        sys.exit(f'Unable to find protein sequence {protein_name} in record')


def align_busco_genes(output_dir: str, verbose: bool, cores: int, genome_number: int, gap_threshold: float) -> None:
    """Aligns BUSCO protein fasta files
    
        args:
            output_dir (str): output directory
            verbose (bool): for displaying debug messages
            cores (int): CPU cores
            genome_number (int): Number of genomes. Used to filter out BUSCO genes not in all genomes.
            gap_threshold (float): gap threshold for trimal. This is the required fraction of non-gapped characters in each column.

        returns:
            None
    """

    aligned_hmms = []
    for subdir, dirs, files in os.walk(output_dir):
        for dir in dirs:
            merged_name = os.path.join(subdir, dir, f'{dir}_merged.faa')
            faa_files = glob.glob(os.path.join(subdir, dir, '*.faa'))
            if len(faa_files) != genome_number:
                continue
            aligned_hmms.append(dir)
            files = []
            for faa in faa_files:
                with open(faa, 'r') as f:
                    files.append(f.read())
            files = '\n'.join(files)
            with open(merged_name, 'w') as f:
                f.write(files)
            muscle_cmd = f'mafft --auto {merged_name} {merged_name}_mafft.faa'
            if verbose:
                print(muscle_cmd)
            subprocess.check_output(muscle_cmd, shell=True)

            trim_alignments(f'{merged_name}_mafft.faa', gap_threshold)

    if len(aligned_hmms) > 0:
        with open(os.path.join(output_dir, 'aligned_hmms.txt'), 'w') as f:
            f.write('\n'.join(aligned_hmms))
    else:
        sys.exit('No BUSCO genes were found. Try adjusting the hmmsearch parameters or remove genomes with missed BUSCOs')

def trim_alignments(fasta: str, gap_threshold: float) -> None:
    """Removes gappy regions from alignments.

        args:
            fasta (str): path to alignment, in fasta format
            gap_threshold (float): gap threshold for trimal. This is the required fraction of non-gapped characters in each column.
        
        returns:
            Nones
    """

    cmd = f'trimal -in {fasta} -out {fasta}.trimal -gt {gap_threshold}'
    subprocess.check_output(cmd, shell=True)

def concat_alignments(output_dir: str) -> None:
    muscle_files = list(Path(output_dir).glob('**/*.trimal'))
    if len(muscle_files) == 0:
        sys.exit(f'No aligned files were found in {output_dir}')
    seqs = defaultdict(list)
    for f in muscle_files:
        with open(f) as f:
            for record in SeqIO.parse(f, 'fasta'):
                seqs[record.name].append(str(record.seq))
    with open(os.path.join(output_dir, 'merged_buscos.faa'), 'w') as f:
        for key, val in seqs.items():
            joined_seqs = ''.join(val)
            f.write('\n'.join(['>'+key, joined_seqs])+'\n')

def main() -> None:
    args = _parse_arguments()
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    genomes = list(
            glob.glob(os.path.join(args.input_dir, '*.faa'))
        ) + list (
            glob.glob(os.path.join(args.input_dir, '*.fasta'))
        )

    if len(genomes) == 0:
        sys.exit('No genomes in input directory. Genomes must be in .faa or .fasta format')
    
    print(f'Running BUSCO extraction on {len(genomes)} genomes')


    print('Running hmmsearch')

    for genome in genomes:
        if args.verbose:
            print(genome)
        
        hmmsearch(
            genome,
            args.hmm_dir,
            args.output_dir,
            args.verbose,
            args.cores,
            args.e_value_min,
            args.seq_length_min
            )
        
    print('Aligning BUSCO protein sequences')

    align_busco_genes(args.output_dir, args.verbose, args.cores, len(genomes), args.gap_threshold)
    
    print('Concatenating alignments')

    concat_alignments(args.output_dir)

    print('Finished!')


if __name__ == '__main__':
    main()