#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO, pairwise2

REFS = {'A': 'NC_038235.1', 'B': 'NC_001781.1'}

def fetch_genbank_record(acc, email, api_key=None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(handle, 'genbank')
    handle.close()
    return rec

def get_gene_coords(rec, gene_name):
    for feat in rec.features:
        if feat.type == 'CDS':
            if 'gene' in feat.qualifiers and gene_name.lower() in [x.lower() for x in feat.qualifiers['gene']]:
                loc = feat.location
                start = int(loc.start)
                end = int(loc.end)
                strand = int(loc.strand or 1)
                return start, end, strand
    raise ValueError(f'Gene {gene_name} not found in {rec.id}')

def map_slice(ref_seq, q_seq, start, end):
    aln = pairwise2.align.globalms(ref_seq.upper(), q_seq.upper(), 2, -1, -5, -1, one_alignment_only=True)[0]
    ref_aln, q_aln = aln.seqA, aln.seqB
    ref_i = 0
    q_i = 0
    q_start = None
    q_end = None
    for ra, qa in zip(ref_aln, q_aln):
        if ra != '-':
            if ref_i == start and q_start is None:
                q_start = q_i
            ref_i += 1
            if ref_i == end:
                q_end = q_i
                break
        if qa != '-':
            q_i += 1
    if q_start is None:
        q_start = 0
    if q_end is None:
        q_end = q_i
    return q_start, q_end

def aa_diffs(ref_aa, q_aa):
    diffs = []
    L = min(len(ref_aa), len(q_aa))
    for i in range(L):
        r = ref_aa[i]
        q = q_aa[i]
        if r != q:
            diffs.append((i+1, r, q))
    return diffs

def main():
    ap = argparse.ArgumentParser(description='Report AA mutations in G and F genes relative to RSV references')
    ap.add_argument('--email', required=False, help='Entrez email, or set env NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set env NCBI_API_KEY')
    ap.add_argument('--subtype', choices=['A','B'], required=True)
    ap.add_argument('--genomes', required=True, help='FASTA with one or more complete RSV genomes')
    ap.add_argument('--out_tsv', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL to comply with NCBI usage policy')
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    ref_acc = REFS[args.subtype]
    ref_gb = fetch_genbank_record(ref_acc, email, api_key)
    ref_seq = ref_gb.seq
    g_start, g_end, g_strand = get_gene_coords(ref_gb, 'G')
    f_start, f_end, f_strand = get_gene_coords(ref_gb, 'F')
    ref_g_nt = ref_seq[g_start:g_end] if g_strand == 1 else ref_seq[g_start:g_end].reverse_complement()
    ref_f_nt = ref_seq[f_start:f_end] if f_strand == 1 else ref_seq[f_start:f_end].reverse_complement()
    ref_g_aa = ref_g_nt.translate()
    ref_f_aa = ref_f_nt.translate()

    recs = list(SeqIO.parse(args.genomes, 'fasta'))
    with open(args.out_tsv, 'w') as out:
        out.write('sample	gene	pos	ref_aa	alt_aa
')
        for rec in recs:
            gs, ge = map_slice(ref_seq, rec.seq, g_start, g_end)
            fs, fe = map_slice(ref_seq, rec.seq, f_start, f_end)
            g_nt = rec.seq[gs:ge]
            f_nt = rec.seq[fs:fe]
            if g_strand == -1:
                g_nt = g_nt.reverse_complement()
            if f_strand == -1:
                f_nt = f_nt.reverse_complement()
            g_aa = g_nt.translate()
            f_aa = f_nt.translate()
            for gene, refaa, qaa in [('G', ref_g_aa, g_aa), ('F', ref_f_aa, f_aa)]:
                for pos, r, q in aa_diffs(refaa, qaa):
                    out.write(f'{rec.id}	{gene}	{pos}	{r}	{q}
')

if __name__ == '__main__':
    main()
