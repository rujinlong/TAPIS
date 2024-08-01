#!/usr/bin/env python
#
#    Align PacBio reads to reference genome and fix indels/mismatches
#
#    Copyright (C) 2022  Jianqiang Sun
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import sys
import os
import subprocess
import argparse



def filename(fpath):
    return os.path.splitext(os.path.basename(fpath))[0]



def gamp_align(fasta_fname, fasta_fpath, index_dpath, index_name, max_intron, procs, tmp_dpath, epoch, overwrite, verbose):
    align_result_fpath = os.path.join(tmp_dpath,'%s_r%d.sam' % (fasta_fname, epoch + 1))
    
    if overwrite or (not os.path.exists(align_result_fpath)):

        _ = '%s %s | gmap -D %s -d %s --no-chimeras --cross-species --expand-offsets 1 -B 5 -K %s -f samse -n 1 -t %d > %s 2> %s'
        cmd = _ % ('zcat' if fasta_fpath.endswith('gz') else 'cat',
                   fasta_fpath,
                   index_dpath, index_name,
                   max_intron, procs, 
                   align_result_fpath,
                   os.path.splitext(align_result_fpath)[0] + '.log')
    
        if verbose: sys.stderr.write('Executing: %s\n' % cmd)
        status = subprocess.call(cmd, shell=True)
    
        assert status == 0, 'The process of GMAP failed to run.'
    
    else:
        if verbose: sys.stderr.write('Skip GMAP alignment. Results are existed.\n')



def sam2bam(sam_fpath, procs, overwrite, verbose):
    
    bam_fpath = os.path.splitext(sam_fpath)[0] + '.bam'
    
    if overwrite or (not os.path.exists(bam_fpath)):
        
        # SAM to BAM
        cmd = 'samtools view -@ %d -Sbh %s > %s ' % (procs, sam_fpath, bam_fpath)
        if verbose: sys.stderr.write('Executing: %s\n' % cmd)
        status = subprocess.call(cmd, shell=True)
        assert status == 0, 'samtools view failed.'

        # sorting
        sorted_bam_fpath = os.path.splitext(bam_fpath)[0] + '.sorted.bam'
        cmd = 'samtools sort -@ %d -m %dG -o %s %s' % (procs, 2, sorted_bam_fpath, bam_fpath)
        if verbose: sys.stderr.write('Executing: %s\n' % cmd)
        status = subprocess.call(cmd, shell=True)
        assert status == 0, 'samtools sort failed.'

        cmd = 'mv %s %s' % (sorted_bam_fpath, bam_fpath)
        status = subprocess.call(cmd, shell=True)

        # indexing
        cmd = 'samtools index %s' % (bam_fpath)
        if verbose: sys.stderr.write('Executing: %s\n' % cmd)
        status = subprocess.call(cmd, shell=True)
        assert status == 0, 'samtools index failed.'
        
        # remove SAM file
        if os.path.exists(sam_fpath):
            os.remove(sam_fpath)
    
    else:
        if verbose: sys.stderr.write('Skip SAM/BAM conversion. Results are existed.\n')


 

def clean_alignments(tmp_dpath, fasta_fname, reference, iterations, epoch, edit_dist_ratio, overwrite, verbose):
    
    output_files = os.path.join(tmp_dpath, '%s_%s_r%d.%s' % (fasta_fname, '%s', epoch + 1, '%s'))
    
    if overwrite or (not os.path.exists(os.path.join(output_files % ('fixed', 'fa')))):
    
        _ = 'cleanAlignments.py -e %f -t %d -f %s -j %s -s %s -u %s -r %s %s %s %s'
        cmd = _ % (edit_dist_ratio,
                   0 if epoch + 1 < iterations else 40, 
                   output_files % ('fixed', 'fa'),
                   output_files % ('junctions', 'fa'),
                   output_files % ('fixed', 'bam'),
                   output_files % ('unaligned', 'fa'),
                   output_files % ('filtered', 'fa'),
                   reference,
                   os.path.join(tmp_dpath, '%s_r%d.bam' % (fasta_fname, epoch + 1)),
                   '-a 10' if epoch == 0 else '')
        if verbose: cmd = cmd + ' -v'
    
        if verbose: sys.stderr.write('Executing: %s\n' % cmd)
        status = subprocess.call(cmd, shell=True)
    
        assert status == 0, 'The process of cleanAlignments.py failed to run.'
    
    else:
        if verbose: sys.stderr.write('Skip cleanAlignments.py. Results are existed.\n')
 


def main(fasta_fpath, output_dpath,
         iterations, index_dpath, index_name, reference,
         max_intron, edit_dist_ratio, procs, overwrite=False, verbose=True):
    
    tmp_dpath = os.path.join(output_dpath, 'tmp')
    if not os.path.exists(output_dpath): os.makedirs(output_dpath)
    if not os.path.exists(tmp_dpath): os.makedirs(tmp_dpath)

    fasta_fname = filename(fasta_fpath)

    # iterative alignment correction
    for epoch in xrange(iterations):
        if os.path.getsize(fasta_fpath) == 0: break
        if verbose: sys.stderr.write('Starting iteration: %d\n' % (epoch+1))
        
        gamp_align(fasta_fname, fasta_fpath, index_dpath, index_name, max_intron, procs, tmp_dpath, epoch, overwrite, verbose)
        sam2bam(os.path.join(tmp_dpath,'%s_r%d.sam' % (fasta_fname,epoch+1)), procs, overwrite, verbose)
        clean_alignments(tmp_dpath, fasta_fname, reference, iterations, epoch, edit_dist_ratio, overwrite, verbose)
        
        fasta_fpath = os.path.join(tmp_dpath, '%s_fixed_r%d.fa' % (fasta_fname, epoch + 1))
    
    
    
    if verbose: sys.stderr.write('Merging filtered alignments...\n')
    cmd = 'samtools merge %s %s/%s_fixed_r*.bam' % (os.path.join(output_dpath, 'aligned.bam'), tmp_dpath, fasta_fname)
    if verbose: sys.stderr.write('Executing: %s\n' % (cmd))
    status = subprocess.call(cmd, shell=True)
    
    assert status == 0, 'The final step failed to run.'

    if verbose: sys.stderr.write('All processes have been completed successfully.\n')




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Iteratively fix aligned reads using reference genome')
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action='store_true', default=True,
                        help='Verbose mode')
    parser.add_argument('-w', '--overwrite', dest='overwrite',
                        action='store_true', default=False,
                        help='Overwrite the existed files.')

    parser.add_argument('-i', '--iterations', dest="iterations",
                        action='store', type=int, default=3,
                        help='Number of aligment iterations, default=3')
    parser.add_argument('-e', '--edr', dest='edr',
                        action='store', type=float, default=0.10,
                        help='Edit distance ratio, default=10')
    parser.add_argument('-o', '--outdir', dest="outdir",
                        action='store', type=str, default='filtered',
                        help='Output directory, default=./filtered')
    parser.add_argument('-p', '--procs', dest="procs",
                        action='store', type=int, default=1,
                        help='Number of processors, default=1')
    parser.add_argument('-K', '--maxIntron', dest="maxIntron",
                        action='store', type=int, default=8000,
                        help='maximum intron length for gmap, default=8000')

    parser.add_argument("indexesDir",
                        action='store', type=str,
                        help='directory to gmap indexes')
    parser.add_argument("indexName",
                        action='store', type=str,
                        help='name of gmap index')

    parser.add_argument('reference', action='store', 
                        type=str, help='Reference sequence')
    parser.add_argument('fasta', action='store', 
                        type=str, help='Reads to align')

    args = parser.parse_args()
    
    main(args.fasta, args.outdir,
         args.iterations, args.indexesDir, args.indexName, args.reference,
         args.maxIntron, args.edr, args.procs, args.overwrite, args.verbose)



