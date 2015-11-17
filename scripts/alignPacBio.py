#! /usr/bin/env python

#
#    Align PacBio reads to reference genome and fix indels/mismatches
#
#    Copyright (C) 2015  Michael Hamilton
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

import sys, os, subprocess
from argparse import ArgumentParser, ArgumentTypeError
parser = ArgumentParser(description='Iteratively fix aligned reads using reference genome')
parser.add_argument('-v', '--verbose', dest="verbose",
                    action='store_true', default=False,
                    help='Verbose mode')

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

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
             
base = os.path.basename(args.fasta).split('.')[0]
gmapIn = args.fasta


# fix to make generic
GMAP    = 'gmap -D %s -d %s --no-chimeras --cross-species --expand-offsets 1 -B 5 -K %s -f samse -n 1 -t %d %s > %s 2> %s'
CONVERT = 'convertSam.py %s'
CLEAN   = 'cleanAlignments.py -e %f -t %d -f %s -j %s -s %s -u %s -r %s %s %s'
 
for epoch in xrange(args.iterations):
    if os.path.getsize(gmapIn) == 0:
        break
    if args.verbose:
        sys.stderr.write('Starting iteration: %d\n' % (epoch+1))
    cmd = GMAP % (args.indexesDir, args.indexName, args.maxIntron, args.procs, gmapIn, 
                  os.path.join(args.outdir,'%s_r%d.sam' % (base,epoch+1)), 
                  os.path.join(args.outdir,'%s_r%d.log' % (base,epoch+1)))
    if args.verbose:
        sys.stderr.write('Executing: %s\n' % cmd)
    status = subprocess.call(cmd,shell=1)
    if args.gmap or epoch == 0:
        cmd = CONVERT % (os.path.join(args.outdir,'%s_r%d.sam' % (base,epoch+1)) )
    else:
        cmd = CONVERT % (os.path.join(args.outdir,'%s_r%dAligned.out.sam' % (base,epoch+1)) )
    if args.verbose:
        sys.stderr.write('Executing: %s\n' % cmd)
        cmd += ' -v'

    status = subprocess.call(cmd,shell=1)
    if args.gmap or epoch == 0:
        cmd = CLEAN % ( args.edr,
                        0 if epoch+1 < args.iterations else 40, 
                        os.path.join(args.outdir,'%s_fixed_r%d.fa' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_junctions_r%d.fa' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_fixed_r%d.sam' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_unaligned_r%d.fa' % (base, epoch+1)),
                        os.path.join(args.outdir,'%s_filtered_r%d.fa' % (base, epoch+1)),
                        args.reference,
                        os.path.join(args.outdir,'%s_r%d.bam' % (base, epoch+1)), 
                    )
    else:
        cmd = CLEAN % ( args.edr,
                        0 if epoch+1 < args.iterations else 40, 
                        os.path.join(args.outdir,'%s_fixed_r%d.fa' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_junctions_r%d.fa' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_fixed_r%d.sam' % (base, epoch+1)), 
                        os.path.join(args.outdir,'%s_unaligned_r%d.fa' % (base, epoch+1)),
                        os.path.join(args.outdir,'%s_filtered_r%d.fa' % (base, epoch+1)),
                        args.reference,
                        os.path.join(args.outdir,'%s_r%dAligned.bam' % (base, epoch+1)), 
                    )

    if epoch == 0:
        cmd += ' -a 10'
    if args.verbose:
        sys.stderr.write('Executing: %s\n' % cmd)
        cmd += ' -v'
    status = subprocess.call(cmd,shell=1)
    gmapIn = os.path.join(args.outdir,'%s_fixed_r%d.fa' % (base,epoch+1))

