import numpy as np
import pandas as pd
import altair as alt
import duckdb
from pathlib import Path
import gzip
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor
import re
import os
import time
import argparse
from collections import defaultdict


def pigz_file(path, threads=1):
    subprocess.run(['pigz', '-p', str(threads), '-f', path], check=True)
    return path + '.gz'


def pigz_all(files, threads):
    # largest first so the big job gets full parallelism immediately
    files = sorted(files, key=lambda p: os.path.getsize(p), reverse=True)
    for p in files:
        subprocess.run(['pigz', '-p', str(threads), '-f', p], check=True)

def parse_bc_from_header(header):
    """Extract barcode from a FASTQ header containing BC=..."""
    match = re.search(r'BC=([ACGTN]*)', header)
    if match:
        return match.group(1)
    return None


def read_fastq_records(filepath):
    """Generator yielding (header, seq, qual) tuples."""
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().rstrip('\n')
            if not header:
                break
            seq = f.readline().rstrip('\n')
            f.readline()  # '+' line
            qual = f.readline().rstrip('\n')
            yield header, seq, qual


def demux_paired_fastq(fwd_fq, rev_fq, bc_to_ref, output_dir, sample, out_rpt, threads, wt_only):
    """
    Demultiplex paired fastq files by barcode -> refname and output one pair of files per refname. 
    """
    # remove any prior output for this run
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    os.makedirs(output_dir, exist_ok=True)

    # want to buffer output files to increase efficency

    buffers_fwd = defaultdict(list)
    buffers_rev = defaultdict(list)

    unmatched_fwd = []
    unmatched_rev = []

    FLUSH_THRESHOLD = 10000 # how many reads to accumulate in buffers before writing to output file

    stats = defaultdict(int)
    total_reads = 0
    unmatched_count = 0

    def format_fq_lines(header, seq, qual):
        return f'{header}\n{seq}\n+\n{qual}\n'

    def flush_refname(refname, sample):
        fwd_out = os.path.join(output_dir, f'{refname}_{sample}.fwd.fq')  
        rev_out = os.path.join(output_dir, f'{refname}_{sample}.rev.fq')   

        with open(fwd_out, 'a') as outf:
            outf.writelines(buffers_fwd[refname])
        with open(rev_out, 'a') as outf:
            outf.writelines(buffers_rev[refname])
        
        buffers_fwd[refname] = []
        buffers_rev[refname] = []
    

    fwd_reader = read_fastq_records(fwd_fq)
    rev_reader = read_fastq_records(rev_fq)

    for (fh, fs, fq), (rh, rs, rq) in zip(fwd_reader, rev_reader):
        total_reads += 1

        bc_fwd = parse_bc_from_header(fh)
        bc_rev = parse_bc_from_header(rh)

        if bc_fwd != bc_rev:
            raise ValueError(f"BC mismatch at read {total_reads}: fwd={bc_fwd} rev={bc_rev}")

        if bc_fwd and bc_fwd in bc_to_ref:
            refname = bc_to_ref[bc_fwd]
            stats[refname] += 1

            buffers_fwd[refname].append(format_fq_lines(fh, fs, fq))
            buffers_rev[refname].append(format_fq_lines(rh, rs, rq))

            if len(buffers_fwd[refname]) >= FLUSH_THRESHOLD:
                flush_refname(refname, sample)
        else:
            unmatched_count += 1
            unmatched_fwd.append(format_fq_lines(fh, fs, fq))
            unmatched_rev.append(format_fq_lines(rh, rs, rq))

            if len(unmatched_fwd) >= FLUSH_THRESHOLD and not wt_only:
                with open(os.path.join(output_dir, f'{sample}_unmatched.fwd.fq'), 'a') as f:
                    f.writelines(unmatched_fwd)
                with open(os.path.join(output_dir, f'{sample}_unmatched.rev.fq'), 'a') as f:
                    f.writelines(unmatched_rev)
                unmatched_fwd, unmatched_rev = [], []

        if total_reads % 10000000 == 0:
            print(f" {sample} {total_reads:,} pairs | Matched: {total_reads - unmatched_count:,} | Unmatched: {unmatched_count:,}")

    # Final flush and write table
    demux_report = {
            'libname': [],
            'refname': [],
            'fastq_fwd': [],
            'fastq_rev': [],
            'nreads': [],
            'star_ref': []
        }
    for refname in list(buffers_fwd.keys()):
        if buffers_fwd[refname]:
            flush_refname(refname, sample)
        
        demux_report['libname'].append(sample)
        demux_report['refname'].append(refname)
        demux_report['fastq_fwd'].append(os.path.join(output_dir, f'{refname}_{sample}.fwd.fq'))
        demux_report['fastq_rev'].append(os.path.join(output_dir, f'{refname}_{sample}.rev.fq'))
        demux_report['nreads'].append(stats[refname])
        demux_report['star_ref'].append(refname)

    if unmatched_fwd and not wt_only:
        with open(os.path.join(output_dir, f'{sample}_unmatched.fwd.fq'), 'a') as f:
            f.writelines(unmatched_fwd)
        with open(os.path.join(output_dir, f'{sample}_unmatched.rev.fq'), 'a') as f:
            f.writelines(unmatched_rev)

    # gzip files
    print(f"\nDone: {sample} {total_reads:,} pairs | Matched: {total_reads - unmatched_count:,} | Unmatched: {unmatched_count:,}")
    print(f"Refnames with reads: {len(stats)}")
    
    start_time = time.perf_counter()
    files_to_gzip = []
    unmatched_to_gzip = []
    for refname in demux_report['refname']:
        # only the ones that exist
        for p in (os.path.join(output_dir, f'{refname}_{sample}.fwd.fq'),
                  os.path.join(output_dir, f'{refname}_{sample}.rev.fq')):
            if os.path.exists(p):
                files_to_gzip.append(p)

    # unmatched files (may or may not exist)
    for p in (os.path.join(output_dir, f'{sample}_unmatched.fwd.fq'),
              os.path.join(output_dir, f'{sample}_unmatched.rev.fq')):
        if os.path.exists(p):
            unmatched_to_gzip.append(p)

    pigz_all(unmatched_to_gzip, threads)

    files_to_gzip = sorted(files_to_gzip, key=os.path.getsize, reverse=True)
    with ThreadPoolExecutor(max_workers=threads) as ex:
        list(ex.map(pigz_file, files_to_gzip))

    end_time = time.perf_counter()

    print(f"{sample} gzip took {end_time - start_time} seconds to complete.")
    # update report paths to point at the .gz files
    demux_report['fastq_fwd'] = [p + '.gz' for p in demux_report['fastq_fwd']]
    demux_report['fastq_rev'] = [p + '.gz' for p in demux_report['fastq_rev']]

    pd.DataFrame(demux_report).to_csv(out_rpt, sep='\t', index=False)
   

    return stats


# ============ RUN ============
def main():
    # Step 1: Load full lookup into memory

    parser = argparse.ArgumentParser(description='load in paired fastq, a pairing duckdb, output fastqs separated by ref based on pairing bcs and table of samples')

    parser.add_argument('--fq_fwd',  dest='fq_fwd' )
    parser.add_argument('--fq_rev',  dest='fq_rev' )
    parser.add_argument('--pairing_db',  dest='pairing_db' )
    parser.add_argument('--pairing_tbl_name', dest='pairing_tbl_name', default='haps')
    parser.add_argument('--fq_outdir', dest='fq_outdir', help='directory to put demuxed fastq files')
    parser.add_argument('--libname', dest='libname')
    parser.add_argument('--out_rpt',  dest='out_rpt' )
    parser.add_argument('--threads', dest='threads', type=int, default=1)
    parser.add_argument('--select_wt_only', dest='select_wt_only', action='store_true', help='if included then will only output reads that are WT based on bc')

    args = parser.parse_args()

    
    print("Loading barcode lookup...")
    con = duckdb.connect(args.pairing_db, read_only=True)
    
    query = f'SELECT readgroupid, refname FROM {args.pairing_tbl_name}'
    if args.select_wt_only:
        query += ' WHERE n_variants_passing = 0'
        print("Filtering for WT only (n_variants_passing == 0)")

    bc_to_ref = dict(
        con.execute(query).fetchall()
    )

    con.close()
    print(f"Loaded {len(bc_to_ref):,} barcodes")

    # Step 2: Demux
    stats = demux_paired_fastq(
        fwd_fq=args.fq_fwd,
        rev_fq=args.fq_rev,
        bc_to_ref=bc_to_ref,
        output_dir=args.fq_outdir,
        sample = args.libname, 
        out_rpt = args.out_rpt,
        threads=args.threads,
        wt_only=args.select_wt_only
    )


if __name__ == "__main__":
   main()