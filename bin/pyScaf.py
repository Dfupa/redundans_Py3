#!/usr/bin/env python
desc="""Perform scaffolding of contigs using information from (in this order):
- paired-end (PE) and/or mate-pair (MP) libraries (!!!NOT IMPLEMENTED YET!!!)
- long reads
- synteny to reference genome
More info at: https://github.com/lpryszcz/pyScaf
"""
epilog="""Author:
l.p.pryszcz@gmail.com
Warsaw, 12/03/2016
"""

import math, os, sys
import subprocess, resource, subprocess
from datetime import datetime
from FastaIndex import FastaIndex

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values. 
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @return - the percentile of the values
    From http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def median(N, key=lambda x:x):
    """median is 50th percentile."""
    return percentile(N, 0.5, key)

def mean(data):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    #if n < 1:
    #    raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) 

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    #if n < 2:
    #    raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5

class Graph(object):
    """Graph class to represent scaffolds. It shouldn't be invoked directly,
    rather use its children: PairedGraph or SyntenyGraph.
    import pyScaf as ps
    
    #####
    # Scaffolding using PE/MP libraries
    reload(ps); s = ps.ReadGraph(fasta, mapq=10, limit=19571, log=sys.stderr);
    s.add_library(fastq, name="lib1", isize=600, stdev=100, orientation="FR"); print s
    s.add_library(fastq, name="lib2", isize=5000, stdev=1000, orientation="FR"); print s
    s.save(out=open(fasta+".scaffolds.fa", "w"))
    
    ######
    # Reference-based scaffolding
    reload(ps); s = ps.SyntenyGraph('test/run1/contigs.reduced.fa', 'test/ref.fa')
    s.save(out=open(fasta+".scaffolds.ref.fa", "w"))
    """
    def __init__(self, mingap=15, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "Graph"
        self.log = log
        self.printlimit = printlimit
        self.ilinks  = 0

    def logger(self, mssg="", decorate=1):
        """Logging function."""
        head = "\n%s"%("#"*50,)
        timestamp = "\n[%s]"% datetime.ctime(datetime.now())
        memusage  = "[%5i Mb] "%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, )
        if self.log:
            if decorate:
                self.log.write("".join((head, timestamp, memusage, mssg)))
            else:
                self.log.write(mssg)
        
    def shorter(self, v, i=4, sep="_"):
        """Return shortened contig name"""
        return sep.join(v.split(sep)[:i])
        
    def __str__(self):
        """Produce string representation of the graph"""
        out = '%s\n%s contigs: %s\n\n%s links:\n' % (self.name, len(self.contigs), list(self.contigs.keys())[:self.printlimit], self.ilinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                '''if not self.links[ref1][end1]:
                    continue
                ref2, end2, links, gap = self.links[ref1][end1]'''
                for (ref2, end2), (links, gap) in list(self.links[ref1][end1].items()):
                    # skip if v2 longer than v1
                    if self.contigs[ref1] < self.contigs[ref2]:
                        continue
                    out += ' %s (%s) - %s (%s) with %s links; %s bp gap\n' % (self.shorter(ref1), end1, self.shorter(ref2), end2, links, gap)
                    # print up to printlimit
                    i += 1
                    if i > self.printlimit:
                        break
        return out
        
    def _init_storage(self, genome):
        """Load sequences from genome, their sizes and init links"""
        # load fasta into index
        self.sequences = FastaIndex(genome)
        self.seq = self.sequences
        # prepare storage
        self.contigs = {c: self.seq.id2stats[c][0] for c in self.seq}
        self.links   = {c: [{}, {}] for c in self.contigs}
        self.ilinks  = 0
        
    def _add_line(self, ref1, ref2, end1, end2, links, gap):
        """Add connection between contigs. """
        # store connection details
        self.links[ref1][end1][(ref2, end2)] = (links, gap)
        # update connection counter 
        self.ilinks += 1
        
    def _get_seqrecord(self, name, scaffold, orientations, gaps):
        """"Return name & seq for given scaffold"""
        # add empty gap at the end
        gaps.append(0)
        seqs = []
        for c, reverse, gap in zip(scaffold, orientations, gaps):
            seq = self.seq.get_sequence(c, reverse)
            # adjust gap size
            if gap and gap < self.mingap:
                strip = int(gap - self.mingap)
                seq = seq[:strip]
                gap = self.mingap
            seqs.append(str(seq)+"N"*gap)
        return name, "".join(seqs)

    def _format_fasta(self, name, seq, linebases=60):
        """Return formatted FastA entry"""
        seq = '\n'.join(seq[i:i+linebases] for i in range(0, len(seq), linebases))
        return ">%s\n%s\n"%(name, seq)
        
    def save(self, out, format='fasta'):
        """Resolve & report scaffolds"""
        # generate scaffolds
        self._get_scaffolds()
        # report
        self.logger("Reporting scaffolds...\n")
        # log scaffold structure
        log = open(out.name+".tsv", "w")
        logline = "%s\t%s\t%s\t%s\t%s\t%s\n"
        log.write("# name\tsize\tno. of contigs\tordered contigs\tcontig orientations (0-forward; 1-reverse)\tgap sizes (negative gap size = adjacent contigs are overlapping)\n")
        totsize = 0
        added = set()
        for i, (scaffold, orientations, gaps) in enumerate(self.scaffolds, 1):
            # scaffold00001
            name = str("scaffold%5i"%i).replace(' ','0')
            # save scaffold
            name, seq = self._get_seqrecord(name, scaffold, orientations, gaps)
            out.write(self._format_fasta(name, seq))
            # report info
            log.write(logline%(name, len(seq), len(scaffold), " ".join(scaffold),
                               " ".join(map(str, (int(x) for x in orientations))), # x may be bool!
                               " ".join(map(str, (x for x in gaps)))))
            totsize += len(seq)
            # store info about added
            for c in scaffold:
                added.add(c)
        # add contigs that were not scaffolded
        for contig in [x for x in self.seq if x not in added]:
            name = "unscaffolded.%s"%contig
            seq = self.seq.get_sequence(contig)
            out.write(self._format_fasta(name, seq))
        # close output & loge
        out.close()
        log.close()
        self.logger(" %s bp in %s scaffolds. Details in %s\n"%(totsize, len(self.scaffolds), log.name), 0)
        self.logger("Scaffolds saved to: %s\n"%out.name, 0)

    def save_dotplot(self, query, readstdin=False): # open('/dev/null','w')): #
        """Produce query to reference dotplot"""
        ### TBD: generate dotplot only if less than several hundreds of seqs in query
        outfn = "%s.%s"%(query, self.dotplot)
        self.logger("Saving dotplots to: %s\n"%outfn) 
        args = ["last-dotplot", "-", outfn]
        if readstdin:
            proc = subprocess.Popen(args, stdin=subprocess.PIPE, stderr=self.log)
        else:
            proc = subprocess.Popen(args, stdin=self._lastal(query), stderr=self.log)
        return proc

    def _lastal(self, queries=[]):
        """Start LAST in local mode and with FastQ input (-Q 1)."""
        # build db
        if not os.path.isfile(self.ref+".suf"):
            os.system("lastdb %s %s" % (self.ref, self.ref))
        # decide on input
        if not queries:
            queries = self.fastq
        # convert filename to list of filenames
        if type(queries) is str:
            queries = [queries, ]
        # run LAST aligner, split and maf-convert in pipe
        seqformati = 1
        args0 = ["cat", ] + queries
        if queries[0].endswith('.gz'):
            args0[0] = "zcat"
            seqformati += 1
        # deduce sequence format
        seqformat = queries[0].split(".")[-seqformati].lower()
        if seqformat in ("fasta", "fa"):
            seqformatcode = "0" # FASTA
        elif seqformat in ("fastq", "fq"):
            seqformatcode = "1" # FastQ
        else:
            self.logger("[WARNING] Unsupported sequence format `%s` in %s\n"%(seqformat, queries[0]), 0)
            sys.exit(1)
        # combine processes
        proc0 = subprocess.Popen(args0, stdout=subprocess.PIPE, stderr=sys.stderr)
        args1 = ["lastal", "-Q", seqformatcode, "-P", str(self.threads), self.ref, "-"]
        proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stdin=proc0.stdout, stderr=sys.stderr)
        args2 = ["last-split", "-"]
        proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stdin=proc1.stdout, stderr=sys.stderr)
        args3 = ["maf-convert", "tab", "-"]
        proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stdin=proc2.stdout, stderr=sys.stderr)
        #print " ".join(args1), queries
        return proc3.stdout

    def _lastal_global(self, query=''):
        """Start LAST in overlap mode. Slightly faster than local mode,
        but much less sensitive."""
        # build db
        if not os.path.isfile(self.ref+".suf"):
            os.system("lastdb %s %s" % (self.ref, self.ref))
        # select query
        if not query:
            query = self.genome
        # run LAST
        args = ["lastal", "-T", "1", "-f", "TAB", "-P", str(self.threads), self.ref, query]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=sys.stderr)
        return proc.stdout

    ###
    # SCAFFOLDING PART
    def _populate_scaffold(self, links, pend, sid, scaffold, orientations, gaps, porientation):
        """Add links to scaffold representation. Experimental!
        !!!Plenty of work needed here!!!
        """
        # self.links[ref1][end1][(ref2, end2)] = (links, gap)
        #ref, end, links, gap = links
        ## there may be many connections now, but only one is processed so far!!!
        #print links
        links_sorted = sorted(iter(links.items()), key=lambda x: x[1][0], reverse=1)
        if len(links_sorted)>1:
            self.logger(" multi connections: %s %s\n"%(scaffold, links), 0)
        (ref, end), (links, gap) = links_sorted[0]
        #for (ref, end), (links, gap) in links.iteritems(): break
        # skip if already added
        if ref in self.contig2scaffold:
            return scaffold, orientations, gaps, porientation
        # get orientation - get forward/reverse-complement signal by XOR
        ## if previous is 1 (end) & current is (0) start & orientation is 0 (forward) --> keep forward orientation
        orientation = (pend != end) != porientation
        # store at the end if previous contig was F and pend 1
        if porientation != pend:
            scaffold.append(ref)
            orientations.append(not orientation)
            gaps.append(gap)
        else:
            scaffold.insert(0, ref)
            orientations.insert(0, not orientation)
            gaps.insert(0, gap)
        # update contigs2scaffold info
        self.contig2scaffold[ref] = sid
        # populate further connections from another end
        links = self.links[ref][abs(end-1)]
        # skip if not links
        if not links:
            return scaffold, orientations, gaps, orientation
        # populate further connections from another end
        return self._populate_scaffold(links, end, sid, scaffold, orientations, gaps, orientation)
        
class ReadGraph(Graph):
    """Graph class to represent scaffolds derived from NGS libraries.
    
    It stores connections from paired alignments in self.readlinks. 
    """
    def __init__(self, genome, mapq=10, load=0.2, threads=4, dotplot="png", 
                 frac=1.5, mingap=15, ratio=0.75, minlinks=5,
                 log=sys.stderr, printlimit=10):
        self.name = "ReadGraph"
        self.log = log
        self.printlimit = printlimit
        self.threads = threads
        # load fasta into index
        self.ref = genome
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.mapq  = mapq
        self.limit = load * sum(self.contigs.values())
        # scaffolding parameters
        self.ratio    = ratio
        self.minlinks = minlinks
        self.mingap   = mingap
        self.frac     = frac
        
    def _set_isize(self, isize=0, stdev=0, orientation=0):
        """Estimate isize stats if necessary"""
        # insert size
        self.isize = isize
        self.stdev = stdev
        self.orientation = orientation
        
    # read-specific part
    def show(self):
        """Produce string representation of the graph"""
        # header
        out = '%s\n%s contigs: %s\n\n%s links' % (self.name, len(self.contigs), list(self.contigs.keys())[:self.printlimit], self.ireadlinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                for ref2 in sorted(self.readlinks[ref1][end1], key=lambda x: self.contigs[x], reverse=1):
                    for end2 in range(2):
                        # skip if v2 longer than v1
                        if self.contigs[ref1] < self.contigs[ref2]:
                            continue
                        # get positions of links
                        positions = self.readlinks[ref1][end1][ref2][end2]
                        if not positions:
                            continue
                        # print up to printlimit
                        i += 1
                        if i > self.printlimit:
                            break
                        out += ' %s (%s) - %s (%s) with %s links: %s\n' % (self.shorter(ref1), end1, self.shorter(ref2), end2, len(positions), str(positions[:self.printlimit]))
        return out

    def _present(self, c):
        """Return True if v in present in the graph"""
        if c in self.readlinks:
            return True
        self.logger("[WARNING] %s not in contigs!\n"%c, 0)
            
    def _get_distance_FR(self, ref, pos, flag):
        """Return True if distance from the v end is smaller that frac * isize"""
        # FR - /1 and /2 doesn't matter, only F / R
        # reverse -> beginning of the contig
        if flag&16:
            return pos, 0
        # forward -> end of the contig
        else:
            return self.contigs[ref]-pos, 1

    def _get_distance_RF(self, ref, pos, flag):
        """Return distance and end from which the read is originating"""
        # RF - /1 and /2 doesn't matter, only F / R
        # forward -> beginning of the contig
        if not flag&16:
            return pos, 0
        # reverse -> end of the contig
        else:
            return self.contigs[ref]-pos, 1
                
    def add_readline(self, ref1, ref2, pos1, pos2, flag1, flag2):
        """Add a line from v1 to v2"""
        if not self._present(ref1) or not self._present(ref2):
            return
        # check if distance is correct
        ## note, distance is corrected by pair orientation
        d1, end1 = self.get_distance(ref1, pos1, flag1)
        d2, end2 = self.get_distance(ref2, pos2, flag2)
        if d1 + d2 > self.maxdist:
            return
        # add connection
        if ref2 not in self.readlinks[ref1][end1]:
            self.readlinks[ref1][end1][ref2] = [[],[]]
        if ref1 not in self.readlinks[ref2][end2]:
            self.readlinks[ref2][end2][ref1] = [[],[]]
        # store connection details
        self.readlinks[ref1][end1][ref2][end2].append(float('%s.%s'%(d1, d2)))
        self.readlinks[ref2][end2][ref1][end1].append(float('%s.%s'%(d2, d1)))
        # update connection counter 
        self.ireadlinks += 1

    def sam_parser(self, handle):
        """Return tuple representing paired alignment from SAM format."""
        q1 = q2 = ""
        for l in handle:
            l = l.strip()
            if not l or l.startswith('@'):
                continue
            sam = l.split('\t')
            #first in pair
            if int(sam[1]) & 64:
                #skip multiple matches
                if sam[0] == q1:
                    continue
                q1, flag1, ref1, start1, mapq1, len1 = sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])
            else:
                #skip multiple matches
                if sam[0] == q2:
                    continue
                q2, flag2, ref2, start2, mapq2, len2 = sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])
            #report
            if q1 == q2:
                yield q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2
        
    def _bwamem(self, fastqs, bwalog=open("/dev/null", "w")):
        """Return pipe to stdout of BWA MEM subprocess"""
        # build db if index missing .sa .bwt .ann .amb .pac
        if not os.path.isfile(self.ref+".sa") or not os.path.isfile(self.ref+".bwt"):
            os.system("bwa index %s > %s 2>&1" % (self.ref, bwalog.name))
        # run BWA MEM
        args = ['bwa', 'mem', '-S', '-t', str(self.threads), self.ref] + fastqs
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=bwalog)
        return proc.stdout
            
    def load_from_SAM(self, handle):
        """Populate graph with readlinks from SAM file.
        Note, SAM file has to be ordered by read name, F read always before R.
        Multiple alignments are allowed, but only the first (assuming best)
        alignment is taken into accound. 
        """
        # parse SAM
        i = j = k = pq1 = 0
        for q1, flag1, r1, s1, mapq1, len1, q2, flag2, r2, s2, mapq2, len2 in self.sam_parser(handle):
            i   += 1
            # skip
            if self.limit and i > self.limit:
                break
            if not self.rlen:
                self.rlen = int(len1+len2)/2
            # skip alignments with low quality
            if self.mapq:
                if mapq1 < self.mapq or mapq2 < self.mapq:
                    continue  
            j   += 1
            #skip self matches
            if r1==r2:
                #isizes.append(abs(s2-s1))
                continue
            k += 1
            #print output
            self.add_readline(r1, r2, s1, s2, flag1, flag2)
        #return isizes
        
    def _get_major_link(self, links, c1, end1):
        """Return major readlink if any of the links full fill quality criteria.
        So far this is too simplistic!
        It works for PE, but for MP you need to allow for multi-joins, 
        as many contigs will be shorter than i.e. 5kb insert... 
        """
        if not links:
            return []
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c])]
        # sort starting from contig with most links
        best  = sorted(links, key=lambda x: x[-1], reverse=1)[0]
        # skip if not enough links or many links to more than one contigs
        sumi = sum(len(pos) for c, e, pos in links)#; print sumi, best
        if len(best[-1]) < self.minlinks or len(best[-1]) < self.ratio*sumi:
            return []
        return [best,]

    def _cluster_links(self, links, maxdist=200):
        """Merge readlinks from the same regions. Experimental!
        !!!Plenty of work needed here!!!
        """
        alllinks = {}
        for (c, e, pos) in sorted(links, key=lambda x: len(x[-1]), reverse=1):
            # get median pos - maybe use mean & stdev?
            mpos = median([int(str(p).split('.')[0]) for p in pos]) 
            if not alllinks:
                alllinks[mpos] = [(c, e, pos)]
                continue
            # sort by distance #filter(lambda x: abs(x-mpos)<maxdist, alllinks)
            closest = sorted(list(alllinks.keys()), key=lambda x: abs(x-mpos))[0]
            if abs(closest-mpos)<maxdist:
                alllinks[closest].append((c, e, pos))
            else:
                alllinks[mpos] = [(c, e, pos)]
        
        for mpos, links in list(alllinks.items()):
            yield links
        
    def _select_links(self, links): 
        """Return best readlinks between set of contigs"""
        if not links:
            return []
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c]) if pos]
        # split links from different regions
        bests = []
        alllinks = self._cluster_links(links)
        for links in alllinks:
            best = links[0] #sorted(links, key=lambda x: len(x[-1]), reverse=1)[0]
            # skip if not enough links or many links to more than one contigs
            sumi = sum(len(pos) for c, e, pos in links)
            if len(best[-1]) < self.minlinks or len(best[-1]) < self.ratio*sumi:
                continue
            ### also check if the same connection in reciprocal bests - no need, will do it during scaffolding
            bests.append(best)
        # make sure best links are not overlapping too much
        
        return bests

    def _calculat_gap_size(self, positions):
        """Return estimated size of the gap"""
        # unload positions
        dists = [self.isize - sum(map(int, str(pos).split('.'))) for pos in positions]
        return median(dists)-self.rlen
        
    def _filter_links(self):
        """Filter readlinks by removing contigs with too many connections ie. repeats"""
        # first identify contigs with too many matches
        
    def get_links(self):
        """Generate connections contig connections for given library."""
        # memory profilling
        #from pympler import asizeof; print "Links: %s bases\n"%asizeof.asizeof(self.links)
        # filter links
        self._filter_links()
        # process starting from the longest 
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                # best connections
                for ref2, end2, pos in self._select_links(self.readlinks[ref1][end1]): # self._get_major_link(self.readlinks[ref1][end1], ref1, end1):
                    # get distance
                    gap = self._calculat_gap_size(pos)
                    # add connection
                    self._add_line(ref1, ref2, end1, end2, len(pos), gap)

    def add_library(self, fastqs, isize=0, stdev=0, orientation=0):
        """Add connections from sequencing library to the graph"""
        self.logger("Adding library...\n")
        # for read info
        self.readlinks  = {c: [{},{}] for c in self.contigs}
        self.ireadlinks = 0
        # set isize
        self._set_isize(isize, stdev, orientation)
        # set distance function
        if   orientation in (1, "FR"):
            self.get_distance = self._get_distance_FR
        elif orientation in (2, "RF"):
            self.get_distance = self._get_distance_RF
        else:
            self.logger("[WARNING] Provided orientation (%s) is not supported!\n"%orientation, 0)
        # set parameters
        self.maxdist = self.isize + self.frac * self.stdev
        self.rlen = None
        # run alignments; parse alignments & get links
        #self.logger(" Library stats...\n", 0)
        handle = self._bwamem(fastqs)
        self.load_from_SAM(handle)
        self.get_links()
        
    # scaffolding part
    def _simplify(self):
        """Simplify scaffold graph by resolving circles"""
        # remove non-reciprocal connections
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                if not self.links[ref1][end1]:
                    continue
                for (ref2, end2), (links, gap) in list(self.links[ref1][end1].items()): 
                    if not (ref1, end1) in self.links[ref2][end2]:
                        # remove link and update counter
                        #self.logger("Removing connection %s -> %s\n"%(str((ref1, end1)), str((ref2, end2))), 0)
                        self.links[ref1][end1].pop((ref2, end2))
                        self.ilinks -= 1
                        
        ## make sure there are no circles

    def _get_scaffolds(self):
        """Resolve & report scaffolds"""
        # simplify graph
        self._simplify()
        
        # build scaffolds
        self.scaffolds = []
        self.contig2scaffold = {}
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            # skip if already added
            if ref1 in self.contig2scaffold:
                continue
            # get scaffold id
            sid = len(self.scaffolds)
            self.contig2scaffold[ref1] = sid
            # store scaffold and its orientation (0-forward, 1-reverse-complement) and gaps
            ## consider blist instead of list!
            scaffold, orientations, gaps = [ref1], [0], []
            # populate scaffold with connections from both ends
            for end1 in range(2):
                links = self.links[ref1][end1]
                if not links:
                    continue
                scaffold, orientations, gaps, porientation = self._populate_scaffold(links, end1, sid, scaffold, orientations, gaps, 0)
            # store
            self.scaffolds.append((scaffold, orientations, gaps))                    

class LongReadGraph(Graph):
    """Graph class to represent scaffolds derived from long read information"""
    def __init__(self, genome, fastq, identity=0.51, overlap=0.66, norearrangements=0, 
                 threads=4, dotplot="png", mingap=15, maxgap=0, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "ReferenceGraph"
        self.log = log
        self.printlimit = printlimit
        # vars
        self.genome = genome
        self.ref = self.genome
        self.fastq = fastq
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.identity = identity
        self.overlap  = overlap
        self.threads  = threads
        self.dotplot  = dotplot
        # scaffolding options
        self.mingap  = mingap
        self._set_maxgap(maxgap)
        self.maxoverhang = 0.1
        # store long links
        self.longlinks = {c: [{}, {}] for c in self.contigs}
        self.ilonglinks = 0
        
    def _set_maxgap(self, maxgap=0, frac=0.01, min_maxgap=10000):
        """Set maxgap to 0.01 of assembly size, 0.01 of assembly size"""
        # set to 0.01 of assembly size
        if not maxgap:
            maxgap = int(round(frac * sum(self.contigs.values())))
        # check if big enough
        if maxgap < min_maxgap:
            maxgap = min_maxgap
        # set variable
        self.maxgap = maxgap
        self.logger(" maxgap cut-off of %s bp\n"%self.maxgap, 0)
        
    def _get_hits(self):
        """Resolve & report scaffolds"""
        # maybe instead of last-split, get longest, non-overlapping matches here
        q2hits, q2size = {}, {}
        for l in self._lastal(): # open('contigs.reduced.fa.tab7'): #
            # strip leading spaces
            l = l.lstrip()        
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            if q not in q2hits:
                q2hits[q] = []
                q2size[q]  = qsize
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                # reverse -> adjust start
                strand = 1 
                qstart = qsize - qstart - qalg
            q2hits[q].append((qstart, qalg, strand, t, tstart, talg))

        return q2hits, q2size
        
    def _get_best_global_match(self, hits):
        """Return best, longest match for given q-t pair"""
        newHits = [[hits[0]]]
        for hit in hits[1:]:
            # break synteny if too large gap
            #if hit[2]=='scaffold22|size195699': print hit, hit[0]-newHits[-1][-1][0]
            if hit[0]-newHits[-1][-1][0] > self.maxgap:
                newHits.append([])
            newHits[-1].append(hit)
        # sort by the longest consecutive alg
        newHits = sorted(newHits, key=lambda x: sum(y[1] for y in x), reverse=1)
        return newHits[0]
        
    def _add_longread_line(self, ref1, ref2, end1, end2, gap):
        """Add connection between contigs from long reads."""
        if (ref2, end2) not in self.longlinks[ref1][end1]:
            self.longlinks[ref1][end1][(ref2, end2)] = []
            self.longlinks[ref2][end2][(ref1, end1)] = []
        # store connection details
        self.longlinks[ref1][end1][(ref2, end2)].append(gap)
        self.longlinks[ref2][end2][(ref1, end1)].append(gap)
        # update connection counter 
        self.ilonglinks += 1

    def _contained_hits(self, hits):
        """Return True if cointaned hits ie A, B, A"""
        added = set([hits[0][3]])
        for i, (qstart, qalg, strand, t, tstart, talg) in enumerate(hits[1:]):
            if t!=hits[i][3] and t in added:
                return True
            added.add(t)
        
    def _hits2longlinks(self, q2hits, q2size, score=0):
        """Filter alignments and populate links.
        Skip: 
        - long reads aligning to only one contig
        - check read overlap
        - mixed alignments ie c1, c2, c1, c2
        - clearly wrong alignment ie c1s - c2s
        """
        for q in list(q2hits.keys()):
            qsize, hits = q2size[q], q2hits[q]
            
            # check if more than 2 contigs aligned
            targets = set(t for qstart, qalg, strand, t, tstart, talg in hits)
            if len(targets)<2:
                continue
                
            # check if enough overlap
            aligned = sum(qalg for qstart, qalg, strand, t, tstart, talg in hits)
            if aligned < self.overlap*qsize:
                continue
            
            hits.sort()

            # check if contained hit (ie A, B, A)
            if self._contained_hits(hits):
                self.logger("contained hits" + "\n".join("\t".join(map(str, (qsize, qstart, qalg, strand, "_".join(t.split("_")[:2]), tstart, talg, self.contigs[t]))) for qstart, qalg, strand, t, tstart, talg in hits) + "\n", 0)
                continue
            
            # combine hits for the same pair
            t2hits = {}
            for qstart, qalg, strand, t, tstart, talg in hits:
                if t not in t2hits:
                    t2hits[t] = []
                t2hits[t].append((tstart, talg, q, qstart, qalg, strand, score))

            uhits = []
            for t in t2hits:
                # get best target hit
                hits = self._get_best_global_match(sorted(t2hits[t]))
                
                # and report rearrangements
                
                # get strand correctly - by majority voting
                talg   = sum(x[1] for x in hits)
                strand = int(round(1.0 * sum(x[1]*x[5] for x in hits) / talg))
                #if strand>1: self.logger("[ERROR] Wrong strand from majority voting: %s %s\n"%(strand, str(hits)), 0)

                if strand:
                    qstart = hits[-1][3]
                    qend   = hits[0][3] + hits[0][4]
                else:
                    qstart = hits[0][3]
                    qend   = hits[-1][3] + hits[-1][4]
                # get global start & end
                tstart = hits[0][0]
                tend   = hits[-1][0] + hits[-1][1]
                uhits.append((q, qstart, qend, strand, t, tstart, tend))

            uhits = sorted(uhits, key=lambda x: x[1])
            self.logger("\t".join(map(str, uhits[0]))+"\n", 0)
            for i, (q, qstart, qend, strand, t, tstart, tend) in enumerate(uhits[1:]):
                dist = qstart - uhits[i][2]
                c1, c2 = t, uhits[i][4]
                # get contig orientation
                end1, end2 = 0, 1
                pos1 = tstart
                if strand:
                    end1 = 1
                    pos1 = self.contigs[c1] - tend
                pos2 = self.contigs[c2] - uhits[i][6]
                if uhits[i][3]: 
                    end2 = 0
                    pos2 = uhits[i][5]
                # calculate gap
                gap = dist - pos1 + pos2
                overhang = gap - dist
                self.logger("\t".join(map(str, (q, qstart, qend, strand, t, tstart, tend, gap)))+"\n", 0)
                self.logger(" %s:%s %s -> %s:%s %s  %s  %s bp\n"%(c1, pos1, end1, c2, pos2, end2, dist, gap), 0)
                # skip if too big overhang on contig edge
                if gap > self.maxgap or overhang > self.maxoverhang*(self.contigs[c1]+self.contigs[c2]):
                    self.logger(" too big contig overhang (%s) or gap (%s)!\n\n"%(overhang, gap), 0)
                    continue
                self._add_longread_line(c1, c2, end1, end2, gap)
            self.logger("\n", 0)
                
        # get links
        self._get_links()
        self.logger("%s %s\n"%(self.ilonglinks, self.ilinks), 0)
                
    def _get_links(self):
        """Combine longlinks into links"""
        for ref1 in self.longlinks:
            for end1, data in enumerate(self.longlinks[ref1]):
                for (ref2, end2), gaps in data.items():
                    if ref1 > ref2:
                        continue
                    links = len(gaps)
                    gap = int(round(mean(gaps)))
                    gapstd = pstdev(gaps)
                    if gap and gapstd>100 and gapstd / gap > 0.25:
                        self.logger("highly variable gap size at %s %s -> %s %s: %s +- %.f %s\n"%(ref1, end1, ref2, end2, gap, gapstd, str(gaps)), 0)
                    self._add_line(ref1, ref2, end1, end2, links, gap)
                    self._add_line(ref2, ref1, end2, end1, links, gap)
                
    def _get_scaffolds(self):
        """Resolve & report scaffolds"""
        self.logger("Aligning long reads on contigs...\n")
        # get best ref-match to each contig
        q2hits, q2size = self._get_hits()

        # get simplified global alignments
        #self.longlinks = {c: [{}, {}] for c in self.contigs}
        self._hits2longlinks(q2hits, q2size)
        
        # build scaffolds
        self.scaffolds = []
        self.contig2scaffold = {}
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            # skip if already added
            if ref1 in self.contig2scaffold:
                continue
            # get scaffold id
            sid = len(self.scaffolds)
            self.contig2scaffold[ref1] = sid
            # store scaffold and its orientation (0-forward, 1-reverse-complement) and gaps
            ## consider blist instead of list!
            scaffold, orientations, gaps = [ref1], [0], []
            # populate scaffold with connections from both ends
            for end1 in range(2):
                links = self.links[ref1][end1]
                if not links:
                    continue
                scaffold, orientations, gaps, porientation = self._populate_scaffold(links, end1, sid, scaffold, orientations, gaps, 0)
            # store
            self.scaffolds.append((scaffold, orientations, gaps))                    

            
class SyntenyGraph(Graph):
    """Graph class to represent scaffolds derived from synteny information"""
    def __init__(self, genome, reference, identity=0.51, overlap=0.66, norearrangements=0, 
                 threads=4, mingap=15, maxgap=0, nofilter_overlaps=0, 
                 dotplot="png", printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "ReferenceGraph"
        self.log = log
        self.printlimit = printlimit
        # vars
        self.genome = genome
        # don't load reference genome - maybe we can avoid that
        self.reference = reference
        self.ref = self.reference
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.identity = identity
        self.overlap  = overlap
        self.threads  = threads
        self.dotplot  = dotplot
        # 0-local alignment; 
        if norearrangements:
            # 1-global/overlap - simpler and faster
            self._get_hits = self._get_hits_global
            self._lastal   = self._lastal_global # needed by save_dotplot
        # scaffolding options
        self.mingap  = mingap
        self._set_maxgap(maxgap)
        self.nofilter_overlaps = nofilter_overlaps
        
    def _set_maxgap(self, maxgap=0, frac=0.01, min_maxgap=10000):
        """Set maxgap to 0.01 of assembly size, 0.01 of assembly size"""
        # set to 0.01 of assembly size
        if not maxgap:
            maxgap = int(round(frac * sum(self.contigs.values())))
        # check if big enough
        if maxgap < min_maxgap:
            maxgap = min_maxgap
        # set variable
        self.maxgap = maxgap
        self.logger(" maxgap cut-off of %s bp\n"%self.maxgap, 0)

    def _clean_overlaps_on_reference(self, _t2hits):
        """Remove hits that overlap on reference too much"""
        t2hits = {}
        for t, hits in _t2hits.items():
            t2hits[t] = []
            # remove hits overlapping too much # sort by r pos
            for tstart, tend, q, qstart, qend, strand in sorted(hits):
                # overlap with previous above threshold
                if t2hits[t] and t2hits[t][-1][1]-tstart > self.overlap*self.contigs[q]:
                    phit = t2hits[t][-1]
                    # do nothing if previous hit is better
                    if tend - tstart <= phit[1]-phit[0]:
                        continue
                    # remove previous match
                    t2hits[t].pop(-1)
                # add match only if first,
                # no overlap with previous or better than previous
                t2hits[t].append((tstart, tend, q, qstart, qend, strand))
        return t2hits

    def _get_hits_global(self):
        """Resolve & report scaffolds"""
        dotplot = None
        if self.dotplot:
            dotplot = self.save_dotplot(self.genome, readstdin=True)
        ## consider splitting into two functions
        ## to facilitate more input formats
        t2hits = {}
        t2size = {}
        q2hits = {}
        for l in self._lastal_global():
            if dotplot:
                try:
                    dotplot.stdin.write(l)
                except:
                    self.logger("[WARNING] dotplot generation failed!\n", 0)
                    dotplot = None
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            #get score, identity & overlap
            identity = 1.0 * score / qalg
            overlap  = 1.0 * qalg / qsize
            #filter by identity and overlap
            if identity < self.identity or overlap < self.overlap:
                continue
            # keep only best match to ref for each q
            if q in q2hits:
                pt, pscore = q2hits[q]
                # skip if worse score
                if score < q2hits[q][1]:
                    continue
                # update if better score
                q2hits[q] = (t, score)
                # remove previous q hit
                # normally should be at the end of matches to pt
                t2hits[pt].pop(-1)
            else:
                q2hits[q] = (t, score)
            # store
            if t not in t2hits:
                t2hits[t] = []
                t2size[t] = tsize            
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                # reverse
                strand = 1 
                qstart = qsize - qstart - qalg
            qend, tend = qstart + qalg, tstart + talg
            t2hits[t].append((tstart, tend, q, qstart, qend, strand))

        # remove q that overlap too much on t
        t2hits = self._clean_overlaps_on_reference(t2hits)
        return t2hits, t2size

    def _get_best_global_match(self, hits):
        """Return best, longest match for given q-t pair"""
        newHits = [[hits[0]]]
        for hit in hits[1:]:
            # break synteny if too large gap
            #if hit[2]=='scaffold22|size195699': print hit, hit[0]-newHits[-1][-1][0]
            if hit[0]-newHits[-1][-1][0] > self.maxgap:
                newHits.append([])
            newHits[-1].append(hit)
        # sort by the longest consecutive alg
        newHits = sorted(newHits, key=lambda x: sum(y[1] for y in x), reverse=1)
        return newHits[0]
        
    def _calculate_global_algs(self, t2hits):
        """Return simplified, global alignments"""
        t2hits2 = {}
        for t in t2hits:
            if t not in t2hits2:
                t2hits2[t] = []
            for q in t2hits[t]:
                # sort by r pos
                hits = self._get_best_global_match(sorted(t2hits[t][q]))
                #get score, identity & overlap
                score = sum(x[-1] for x in hits)
                qalg  = sum(x[4] for x in hits)
                identity = 1.0 * score / qalg # this local identity in hits ignoring gaps
                overlap  = 1.0 * qalg / self.contigs[q]
                #filter by identity and overlap
                if identity < self.identity or overlap < self.overlap:
                    continue
                # get global start & end
                ## this needs work and bulletproofing!!!
                tstart = hits[0][0]
                tend   = hits[-1][0] + hits[-1][1]
                qstart = hits[0][3]
                qend   = hits[-1][3] + hits[-1][4]
                # and report rearrangements
                
                # get strand correctly - by majority voting
                strand = int(round(1.0 * sum(x[4]*x[5] for x in hits) / qalg))
                t2hits2[t].append((tstart, tend, q, qstart, qend, strand))
                #print t, tstart, tend, q, qstart, qend, strand, len(t2hits[t][q]), identity, overlap
                
        return t2hits2
        
    def _get_hits(self):
        """Resolve & report scaffolds"""
        dotplot = None
        if self.dotplot:
            dotplot = self.save_dotplot(self.genome, readstdin=True)
        ## consider splitting into two functions
        ## to facilitate more input formats
        t2hits, t2size = {}, {}
        q2hits = {}
        for l in self._lastal(self.genome):
            # strip leading spaces
            l = l.lstrip()
            if dotplot:
                try:
                    dotplot.stdin.write(l)
                except:
                    self.logger("[WARNING] dotplot generation failed!\n", 0)
                    dotplot = None
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            # prepare storage
            if t not in t2hits:
                t2hits[t] = {q: []}
                t2size[t] = tsize
            elif q not in t2hits[t]:
                t2hits[t][q] = []
            if q not in q2hits:
                q2hits[q] = []
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                # reverse -> adjust start
                strand = 1 
                qstart = qsize - qstart - qalg
            t2hits[t][q].append((tstart, talg, q, qstart, qalg, strand, score))
            q2hits[q].append((qstart, qalg, strand, t, tstart, talg))

        # get simplified global alignments
        t2hits = self._calculate_global_algs(t2hits)
        #for t, hits in t2hits.iteritems(): print t, hits

        # clean overlaps on reference
        if not self.nofilter_overlaps:
            t2hits = self._clean_overlaps_on_reference(t2hits)

        #print "after clean-up"
        #for t, hits in t2hits.iteritems(): print t, hits
        return t2hits, t2size

    def _estimate_gap(self, data, pdata):
        """Return estimated gap size"""
        # current start - previous end
        # this needs to be corrected by qstart and qend !!
        gap = data[0] - pdata[1]
        return gap
        
    def _get_scaffolds(self):
        """Resolve & report scaffolds"""
        self.logger("Aligning contigs on reference...\n")
        # get best ref-match to each contig
        t2hits, t2size = self._get_hits()

        # store scaffold structure
        ## [(contigs, orientations, gaps), ]
        self.scaffolds = []
        added = set()
        for t in sorted(t2size, key=lambda x: t2size[x], reverse=1):
            # skip one-to-one matches
            if len(t2hits[t])<2:
                continue
            # add empty scaffold
            scaffold, orientations, gaps = [], [], []
            #print t, t2size[t]
            for data in t2hits[t]:
                gap = 0
                tstart, tend, q, qstart, qend, strand = data
                # calculate gap only if scaffold has at least one element
                if scaffold:
                    gap = self._estimate_gap(data, pdata)
                    if gap > self.maxgap:
                        # add to scaffolds
                        if len(scaffold)>1:
                            self.scaffolds.append([scaffold, orientations, gaps])
                            added.update(scaffold)
                        # reset storage
                        scaffold, orientations, gaps = [], [], []
                    else:    
                        gaps.append(gap)
                # store contig & orientation
                scaffold.append(q)
                orientations.append(strand)
                # keep track of previous data
                pdata = data
                
            # add to scaffolds
            if len(scaffold)>1:
                self.scaffolds.append([scaffold, orientations, gaps])
                added.update(scaffold)
                
        # add missing
        for c in [x for x in self.contigs if x not in added]:
            self.scaffolds.append([(c,),(0,),[]])

def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return "".join(p.stdout.readlines())

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in dependencies.items():
        out = _check_executable(cmd)
        if "not found" in out:
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            out = subprocess.getoutput("%s --version"%cmd)
            curver = out.split()[-1]
            if not curver.isdigit():
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            elif int(curver)<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/pyScaf#dependencies !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)
            
def main():
    import argparse
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-f", "--fasta", required=1, 
                        help="assembly FASTA file")
    parser.add_argument("-o", "--output", default="scaffolds.fa", type=argparse.FileType('w'), 
                        help="output stream [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, 
                        help="max no. of threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), 
                        help="output log to [stderr]")
    parser.add_argument("--dotplot", default="png", choices=["", "png", "gif", "pdf"],
                        help="generate dotplot as [%(default)s]")
    
    refo = parser.add_argument_group('Reference-based scaffolding options')
    refo.add_argument("-r", "--ref", "--reference", default='', 
                      help="reference FastA file")
    refo.add_argument("--identity", default=0.33, type=float,
                      help="min. identity [%(default)s]")
    refo.add_argument("--overlap", default=0.66, type=float,
                      help="min. overlap  [%(default)s]")
    refo.add_argument("-g", "--maxgap", default=0, type=int,
                      help="max. distance between adjacent contigs [0.01 * assembly_size]")
    refo.add_argument("--norearrangements", action='store_true', 
                      help="high identity mode (rearrangements not allowed)")
    refo.add_argument("-nfo", "--nofilter_overlaps", action='store_true', 
                      help="don't filter overlaps on reference ie. polyploid or transcriptome scaffolding")
    
    scaf = parser.add_argument_group('long read-based scaffolding options (EXPERIMENTAL!)')
    scaf.add_argument("-n", "--longreads", nargs="+",
                      help="FastQ/FastA file(s) with PacBio/ONT reads")
    
    scaf = parser.add_argument_group('NGS-based scaffolding options (!NOT IMPLEMENTED!)')
    scaf.add_argument("-i", "--fastq", nargs="+",
                      help="FASTQ PE/MP files")
    scaf.add_argument("-j", "--joins",  default=5, type=int, 
                      help="min pairs to join contigs [%(default)s]")
    scaf.add_argument("-a", "--linkratio", default=0.7, type=float,
                       help="max link ratio between two best contig pairs [%(default)s]")    
    scaf.add_argument("-l", "--load",  default=0.2, type=float, 
                      help="align subset of reads [%(default)s]")
    scaf.add_argument("-q", "--mapq",    default=10, type=int, 
                      help="min mapping quality [%(default)s]")
    
    # standard
    #parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.12a')
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    #if o.verbose:
    o.log.write("Options: %s\n"%str(o))

    # check if all executables exists & in correct versions
    dependencies = {'lastal': 700, 'lastdb': 700, }
    _check_dependencies(dependencies)
    
    # check logic
    if not o.ref and not o.fastq and not o.longreads:
        sys.stderr.write("Provide FastQ files, reference genome and/or long reads!\n")
        sys.exit(1)
        
    # check if input files exists
    fnames = [o.fasta, o.ref]
    if o.fastq:
        fnames += o.fastq
    if o.longreads:
        fnames += o.longreads
    for fn in fnames: 
        if fn and not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)
            
    fasta = o.fasta
    out = o.output
    log = o.log

    # perform scaffolding using long reads
    if o.longreads:
        s = LongReadGraph(fasta, o.longreads, identity=o.identity, overlap=o.overlap, \
                          maxgap=o.maxgap, threads=o.threads, dotplot=o.dotplot, 
                          norearrangements=o.norearrangements, log=log)
        # save output
        s.save(out)
        if o.dotplot:
            s.save_dotplot(out.name).wait()
        #fasta = fasta+".scaffolds.fa"
            
    # perform PE/MP scaffolding if NGS provided
    elif o.fastq:
        # NOT IMPLEMENTED YET
        sys.stderr.write("NGS-based scaffolding is not implemented yet! Stay tuned :)\n"); sys.exit(1)
        # get library statistics
        # init
        s = ReadGraph(fasta, mapq=o.mapq, load=o.load, threads=o.threads, dotplot=o.dotplot,
                      log=log)
        # add library
        s.add_library(o.fastq, isize=600, stdev=100, orientation="FR"); s.show()
        # add another library

        # save output
        s.save(out)
        if o.dotplot:
            s.save_dotplot(out.name).wait()
        # update fasta at the end
        #fasta = fasta+".scaffolds.fa"
    
    # perform referece-based scaffolding only if ref provided
    elif o.ref:
        # init
        s = SyntenyGraph(fasta, o.ref, identity=o.identity, overlap=o.overlap, \
                         maxgap=o.maxgap, threads=o.threads, dotplot=o.dotplot, \
                         nofilter_overlaps=o.nofilter_overlaps, 
                         norearrangements=o.norearrangements, log=log)
        # save output
        s.save(out)
        if o.dotplot:
            s.save_dotplot(out.name).wait()
        s.logger("Done!\n")

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    
