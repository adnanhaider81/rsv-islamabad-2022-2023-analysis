#!/usr/bin/env python3
import argparse, subprocess, shlex, sys
from pathlib import Path

def run(cmd):
    print('+', cmd, file=sys.stderr)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print(p.stderr, file=sys.stderr)
        sys.exit(p.returncode)
    return p.stdout

def main():
    ap = argparse.ArgumentParser(description='Find closest nt matches for contigs using BLAST remote')
    ap.add_argument('--in', dest='contigs', required=True)
    ap.add_argument('--max_hits', type=int, default=50)
    ap.add_argument('--out_tsv', required=True)
    args = ap.parse_args()
    fmt = '6 qseqid sacc pident length bitscore staxids stitle'
    cmd = f"blastn -remote -db nt -query {shlex.quote(args.contigs)} -outfmt '{fmt}' -max_target_seqs {args.max_hits}"
    out = run(cmd)
    Path(args.out_tsv).write_text(out, encoding='utf-8')
    print(f'Wrote {args.out_tsv}')

if __name__ == '__main__':
    main()
