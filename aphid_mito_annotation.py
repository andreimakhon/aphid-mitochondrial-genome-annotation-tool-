from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Data import CodonTable
from io import StringIO
from os import listdir
import os
import shutil
import subprocess
import pandas as pd
import sys, getopt


def find_rna_position(align_seq, seq_len):
    '''
    function find find position of RNA in genome
    if 75% of same nucleotides alignmented its start or stop position

    Parameters
    ----------
    align_seq : Alignmented sequences 

    Returns
    -------
    start : start rna 
    end : end rna

    '''
    seqs_rev = [seq.reverse_complement() for seq in align_seq]
    len_seqs = len(align_seq)
    start = 0
    end = 0
    #find start position
    for i,s in enumerate(zip(*align_seq)):
        num_allign = len_seqs - s.count('-') - 1    #number of alignmented nucleotides in same position
        if num_allign/len_seqs >= 0.75 and not start:
            start = i + 1
    for i,s in enumerate(zip(*seqs_rev)):
        num_allign = len_seqs - s.count('-') - 1
        if num_allign/len_seqs >= 0.75 and not end:
            end = seq_len - i  + 1
    return start,end


def rna_main(main_rec,path_to_rna_db):
    all_rna = []
    for f in path_to_rna_db:
        db = open(f'dbr/{f}')
        to_allignm = [db_rec for db_rec in SeqIO.parse(db,'fasta')]
        to_allignm.append(main_rec)
        muscle_cline = MuscleCommandline()
        handle = StringIO()
        SeqIO.write(to_allignm, handle, "fasta")
        data = handle.getvalue()
        stdout, stderr = muscle_cline(stdin=data)
        allignmented_seqs = [align.seq for align in 
                             SeqIO.parse(StringIO(stdout),'fasta')]
        start_rna,end_rna = find_rna_position(allignmented_seqs,len(main_rec.seq))
        seq_rec = SeqRecord(
            main_rec.seq[start_rna -1 : end_rna],
            id=f[:-4],name='',
            annotations={'location':f'{start_rna}-{end_rna}',
                         'strand' : '-1', 
                         'size' : len(main_rec.seq[start_rna -1 : end_rna]),})
        all_rna.append(seq_rec)
    return all_rna
    

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    if start >= 0:
                        answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    return answer    
    
    
def find_probable_cds(record,table,min_pro_len):
    prob_cds = []
    orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
    for i,(start, end, strand, pro) in enumerate(orf_list):
        if strand > 0:
            s = record.seq[start:end].translate(table=5)
            seq_rec = SeqRecord(s,id=f'prot_{str(i)}',name='',
                            description=f'{start}-{end}')
        else:
            s = record.seq[start:end].reverse_complement().translate(table=5)
            seq_rec = SeqRecord(s,id=f'prot_{str(i)}',name='',
                            description=f'c({start}-{end})')
        prob_cds.append(seq_rec)
    return prob_cds
    
    
def main_cds(record,prob_cds):
    os.mkdir('temp')
    proteins = []
    genome_seq = record.seq
    for count,seq_rec in enumerate(prob_cds):
        coord = seq_rec.description
        SeqIO.write(seq_rec,
                f'temp/a.fasta', "fasta")
        blastp_cline = NcbiblastpCommandline(
            query=f'temp/a.fasta',
            db="db/all_genes_mito", 
            evalue=3e-20, outfmt=5, 
            out=f'temp/{count}.xml')
        stdout, stderr = blastp_cline()
        blast_res = open(f'temp/{count}.xml')
        blast_record = NCBIXML.read(blast_res)
        query_f,query_t = [],[]
        name = ''
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if not name:
                    name = alignment.title.split()[1].split('_')[0].lower()
                query_f = hsp.query_start
                query_t = hsp.query_end
                if 'c' in coord:
                    strand = -1
                    stop = int(coord[:-1].split('-')[1])
                    positions = [stop-(query_t*3)+1,stop-((query_f*3)-3)]
                    prob_start = genome_seq[positions[1]-3:positions[1]].reverse_complement()
                    prob_end = genome_seq[positions[0]-1:positions[0]+2].reverse_complement()
                    while prob_start not in start_codons:
                        positions[1] -= 3
                        prob_start = genome_seq[positions[1]-3:positions[1]].reverse_complement()
                    while prob_end not in stop_codons and name not in exeptions:
                        positions[0] -= 3
                        prob_end = genome_seq[positions[0]-1:positions[0]+2].reverse_complement()    
                else:
                    strand = 1
                    start = int(coord.split('-')[0])
                    positions = [start+1+((query_f-1)*3),start+(query_t*3)]
                    prob_start = genome_seq[positions[0]-1:positions[0]+2]
                    prob_end = genome_seq[positions[1]-3:positions[1]]
                    while prob_start not in start_codons:
                        positions[0] += 3
                        prob_start = genome_seq[positions[0]-1:positions[0]+2]
                    while prob_end not in stop_codons and name not in exeptions:
                        positions[1] += 3
                        prob_end = genome_seq[positions[1]-3:positions[1]]
                proteins.append(SeqRecord(
                genome_seq[positions[0]-1:positions[1]],
                id=name,name='',
                annotations = {'location' : f'{positions[0]}-{positions[1]}',
                'strand' : strand, 
                'size': len(genome_seq[positions[0]-1:positions[1]]),
                'start_stop_codons': f'{prob_start}/{prob_end}'}))
                break
            break  
    shutil.rmtree('temp')
    return proteins


def trna_main(all_genes):
    
    def merge(l):
       i = 0
       len_coords = len(coords)
       while i < len_coords:
           try:
               if coords[i][1] > coords[i+1][0]:
                   coords[i:i+2] = [[coords[i][0],coords[i+1][1]]]
                   len_coords -= 1
               else:
                   i+=1
           except:
               break
    all_genes = sorted(all_genes,key= lambda pos: 
                           int(pos.annotations['location'].split('-')[0]))
    trna = []                      
    coords = [[int(gene.annotations['location'].split('-')[0]), 
               int(gene.annotations['location'].split('-')[1])] for gene in all_genes]
    merge(coords)
    coords_range = [list(range(i[0], i[1]+1))for i in coords]
    result_arwen = subprocess.run(['arwen/./arwen',
                                   main_file_path,
                                   '-gcinvert', '-w','-rp'], 
                                  stdout=subprocess.PIPE).stdout.decode('utf-8')
    low_score_trna = [trna for trna in list(StringIO(result_arwen))[2:] if '*' in trna]
    hight_score_trna = [trna for trna in list(StringIO(result_arwen))[2:] if '*' not in trna]
    
    for trna_list in (hight_score_trna,low_score_trna):
        for record_trna in trna_list:
            pos = record_trna.split()[-3]
            end = int(pos.split(',')[1][:-1])
            if 'c' in pos:
                start = int(pos.split(',')[0][2:])
                strand = -1
            else:
                start = int(pos.split(',')[0][1:])
                strand = 1
            len_list = len(coords_range)
            counter = 0
            integrate = 0
            while counter < len_list:
                if start in coords_range[counter] and end in coords_range[counter]:
                    break
                elif start not in coords_range[counter] and end not in coords_range[counter]:
                    counter += 1
                    continue
                else:
                    if start in coords_range[counter]:
                        integrate += coords_range[counter][-1] - start
                    try:
                        if end in coords_range[counter+1]:
                            integrate += end - coords_range[counter+1][0]
                    except:
                        counter += len_list
                    counter += len_list
            else:
                if integrate/(end-start+1) < 0.25:
                    if '*' in record_trna.split()[1]:
                        name = record_trna.split()[1][:-1]
                    else:
                        name = record_trna.split()[1]
                    trna.append(SeqRecord(record.seq[start-1:end], id=name,name='',
                    annotations = {'location' : f'{start}-{end}',
                      'strand' : strand, 
                      'size': (end-start+1),
                      'anticodon' : record_trna.split()[-1][1:-1].upper(),
                      'anticodon position': int(record_trna.split()[-2])+start}))
        coords.extend([int(gene.annotations['location'].split('-')[0]), 
               int(gene.annotations['location'].split('-')[1])] for gene in trna)
        coords.sort(key= lambda pos: pos[0])
        merge(coords)
        coords_range = [list(range(i[0], i[1]+1))for i in coords]
    return trna


def non_coding_regions(left_gene,right_gene,non_coding_name):
    left_pos,left_index = 0,0
    right_pos,right_index = 0,0
    for index, gene in enumerate(all_genes):
        if left_gene in gene.id:
            left_index = index
            left_pos = int(gene.annotations['location'].split('-')[1]) + 1
        if right_gene in gene.id:
            right_index = index
            right_pos = int(gene.annotations['location'].split('-')[0]) -1
    if left_pos > right_pos:
        right_pos = len(record.seq)
        right_index = len(all_genes)
    del all_genes[left_index+1:right_index]
    all_genes.insert(left_index+1,SeqRecord(
        record.seq[left_pos-1:right_pos],
        id=non_coding_name,name='',annotations={
            'location':f'{left_pos}-{right_pos}',
            'strand':'',
            'size':right_pos-left_pos-1}))    

main_file_path = ''
ofile = 'Annotation Table.xlsx'
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")
for o, a in myopts:
    if o == '-i':
        main_file_path=a
    elif o == '-o':
        ofile=a
        
all_genes = []
pd_dataframe = {'Gene':[],'Strand':[],'Position':[],'Size(bp)':[],
                'Anticodon':[],'Anticodon Positions':[],'Start/Stop Codons':[],
                'Intergenic Nucleotides':[] }
mito_table = CodonTable.unambiguous_dna_by_id[5]
start_codons = mito_table.start_codons
stop_codons = mito_table.stop_codons
exeptions = ['nad4']
rna_db_folder = listdir('dbr/')    #path to rna database
with open(main_file_path) as mainf:
    for record in SeqIO.parse(mainf,'fasta'):
        for rna in rna_main(record, rna_db_folder):
            all_genes.append(rna)
        table = 5
        min_pro_len = 45
        prob_cds = find_probable_cds(record, table, min_pro_len)
        for cds in main_cds(record,prob_cds):
            all_genes.append(cds)
        for trna in trna_main(all_genes):
            all_genes.append(trna)
        all_genes = sorted(all_genes,key= lambda pos: 
                           int(pos.annotations['location'].split('-')[0]))
        non_coding_regions('mtRNA-Glu','mtRNA-Phe', 'repeat region')
        non_coding_regions('rrnS','mtRNA-Ile', 'control region')
        for index, gene in enumerate(all_genes):
            pd_dataframe['Gene'].append(gene.id)
            if gene.annotations['strand'] == -1:
                pd_dataframe['Strand'].append('N')
            else: pd_dataframe['Strand'].append('J')
            pd_dataframe['Position'].append(gene.annotations['location'])
            pd_dataframe['Size(bp)'].append(gene.annotations['size'])
            if 'anticodon' in gene.annotations:
                pd_dataframe['Anticodon'].append(gene.annotations['anticodon'])
            else: pd_dataframe['Anticodon'].append('')
            if 'anticodon position' in gene.annotations:
                pd_dataframe['Anticodon Positions'].append(
                    gene.annotations['anticodon position'])
            else: pd_dataframe['Anticodon Positions'].append('')
            if 'start_stop_codons' in gene.annotations:
                pd_dataframe['Start/Stop Codons'].append(
                    gene.annotations['start_stop_codons'])
            else: pd_dataframe['Start/Stop Codons'].append('')
            try:
                intg_reg = (int(all_genes[index+1].annotations['location'].split('-')[0])
                -1-int(gene.annotations['location'].split('-')[1]))
                if intg_reg:
                    pd_dataframe['Intergenic Nucleotides'].append(intg_reg)
                else: 
                    pd_dataframe['Intergenic Nucleotides'].append('-')
            except:
                intg_reg = len(record.seq)-int(gene.annotations['location'].split('-')[1])
                if intg_reg:
                    pd_dataframe['Intergenic Nucleotides'].append(intg_reg)
                else: 
                    pd_dataframe['Intergenic Nucleotides'].append('-')
        df1 = pd.DataFrame(pd_dataframe, index=pd_dataframe['Gene'])
        df1.to_excel(ofile,index=False)  
print('Job done!')
            
            
            
        
            
        
        
        