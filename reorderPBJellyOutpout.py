#!/usr/bin/env python
from __future__ import division
import sys, os, re ,json, argparse, tempfile, subprocess
import math


class DNASequence(object):
    """
    An object to hold a DNA sequence with some methods to manipulate it
    """
    
    alphabet = set('ACTGNactgn')
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',\
                  'a': 't','c' : 'g', 'g': 'c','t': 'a','n':'n'}
    
    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self.sequence = sequence
        
        for letter in self.sequence:
            if letter not in self.alphabet:
                raise IncorrectSequenceLetter(letter,\
                                              self.__class__.__name__)
    def __len__(self):
        return len(self.sequence)
    
    def __cmp__(self, other):
        """
        return longer sequence of two
        """
        if len(self)==len(other):   return 0
        elif len(self)<len(other):  return -1
        elif len(self)>len(other):  return 1

    def get_identifier(self):
        return self.identifier

    def get_sequence(self):
        return self.sequence

    def get_fasta_string(self):
        """
        return a nicely formated fasta string
        """
        splitted_seq = []
        for x in xrange(0,len(self.sequence),80):
            splitted_seq.append(self.sequence[x:x+80])
        return ">%s\n%s\n" %(self.identifier,"\n".join(splitted_seq))
    
    def get_complement(self):
        """
        return a DNASequence instance of the complement of the current sequence
        """
        complementSeq = "".join([ self.complement[letter] for letter in self.sequence ])
        return DNASequence(identifier = self.identifier+"_complement",\
                           sequence = complementSeq[::-1] )


class IncorrectSequenceLetter(ValueError):

    def __init__(self, letter, classname):
        self.message = "The sequence item %s is not found in the alphabet of class %s\n" %(letter, classname)
    
class Psl(object):
    ''' Class to represent PSL output from BLAT.
    Provides basic methods such as score() and calcPercentIdentity().
    Argument: String in PSL format
    #note by LK: This class was not written by me, but I can't find
    #the source.
    '''
    
    def __init__(self, s):
        
        # split and tokenize input
        fields = s.strip().split()
        num_fields = len(fields)
        matches, mismatches, repmatches, ncount, qnuminsert, qbaseinsert, \
            tnuminsert, tbaseinsert, strand, qname, qsize, qstart, qend, \
            tname, tsize, tstart, tend, blockcount, blocksizes, qstarts, \
            tstarts = fields[0:21]
        
        self.matches = int(matches)
        self.mismatches = int(mismatches)
        self.repmatches = int(repmatches)
        self.ncount = int(ncount)
        self.qnuminsert = int(qnuminsert)
        self.qbaseinsert = int(qbaseinsert)
        self.tnuminsert = int(tnuminsert)
        self.tbaseinsert = int(tbaseinsert)
        self.strand = strand
        self.qname = qname
        self.qsize = int(qsize)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.tname = tname
        self.tsize = int(tsize)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.blockcount = int(blockcount)
        self.blocksizes = [int(x) for x in blocksizes.split(',')[0:-1]]
        self.qstarts = [int(x) for x in qstarts.split(',')[0:-1]]
        self.tstarts = [int(x) for x in tstarts.strip().split(',')[0:-1]]
        
    ## Private methods
        
    def __lenmul(self):
        ''' Determine length multiplier.
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        if self.__isProtein:
            return 3
        else:
            return 1

    def __isProtein(self):
        ''' Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        lastblock = self.blockcount - 1
        return (self.strand[1:1] == '+' and \
                    self.tend == (self.tstarts[lastblock] + (3 * self.blocksizes[lastblock]))) or \
                    ((self.strand[1:1] == '-') and \
                         (self.tstart == (self.tsize - (self.tstarts[lastblock] + 3*self.blocksizes[lastblock]))))
    
    def __calcMilliBad(self, ismrna):
        ''' Return number of non-identical matches. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        qalisize = self.__lenmul() * self.qspan()
        alisize = min(qalisize, self.tspan())
        millibad = 0
        if alisize <= 0: return 0
        sizediff = alisize - self.tspan()
        if sizediff < 0:
            if ismrna:
                sizediff = 0
            else:
                sizediff = -sizediff
        insertfactor = self.qnuminsert
        if not ismrna: insertfactor += self.tnuminsert
        total = self.__lenmul() *\
            (self.matches + self.repmatches + self.mismatches)
        if total != 0:
            millibad = (1000 * (self.mismatches * self.__lenmul() + insertfactor + \
                                    round(3*math.log(1 + sizediff)))) / total
        return millibad

    # Public methods
    
    def qspan(self):
        ''' Span of alignment on query sequence '''
        return self.qend - self.qstart
    
    def tspan(self):
        ''' Span of alignment on target sequence '''
        return self.tend - self.tstart
    
    def score(self):
        ''' Score as calculated by web-BLAT. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        return self.matches + (self.repmatches / 2) - self.mismatches - \
            self.qnuminsert - self.tnuminsert
    
    def calcPercentIdentity(self):
        ''' Percent identity as calculated by web-BLAT. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        return 100.0 - self.__calcMilliBad(True) * 0.1
    
def loadScaffolOrder(fasta):
    """
    Parses a fasta to return a list of its headers in the order they appear
    in the assembly
    """
    
    with open(fasta) as f:
        return [line.strip(">").rstrip() for line in f if line.startswith(">")]
    
def loadFasta(fasta):
    """
    A generator function to parse a fasta and return it as
    DNASequence instances
    """
    
    with open(fasta) as f:
        sequence = ""
        for line in f:
            if line[0] == ">":
                if len(sequence) > 0:
                    yield DNASequence(identifier, sequence)
                identifier = line.strip(">").rstrip()
                sequence = ""
            else:
                sequence += line.rstrip()

        if len(sequence) > 0:
            yield DNASequence(identifier, sequence)
            
def getOrientationBlat(query, reference, prePath):
    """
    Given two sequences, run blat in order to determine
    their relative orientaiton
    """
    output = os.path.join(prePath,"out")
    blatCall = "blat {} {} {} -out=psl -noHead".format(query,reference,output)
    
    n = subprocess.Popen([blatCall], shell=True)
    nReturnCode = n.wait()
    
    pslOut = output +".bestHit"
    pslTrash = output+".pslTrash"
    
    pslRepsCall = "pslReps {} {} {} -singleHit -nohead".format(output, pslOut, pslTrash)
    f = subprocess.Popen([pslRepsCall], shell=True, stdout = subprocess.PIPE, \
                        stderr = subprocess.PIPE)
    
    fReturnCode = f.wait()
    
    #check if blat executed correctly
    if nReturnCode != 0:
        sys.stderr.write("blat crashed with exit code {}\n".format(nReturnCode))
        return None
        #delete the temporary files:
        os.remove(output)
        os.remove(pslOut)
        os.remove(pslTrash)
    #check is pslReps executed correctly
    elif fReturnCode != 0:
        sys.stderr.write("psl crashed with exit code {}\n".format(fReturnCode))
        os.remove(output)
        os.remove(pslOut)
        os.remove(pslTrash)
    else:
        
        #assign scores based on percentage identity * average alignment length
        #to alignments on both strands and return higher scoring strand
        
        scoreS = {"+": 0.0, "-": 0.0}
        with open(pslOut) as o:
            for line in o:
                p = Psl(line)
                percID = p.calcPercentIdentity()
                alnLength = (p.qspan() + p.tspan())/2
                score = alnLength * percID
                scoreS[p.strand] += score
                
        #delete the temporary file:
        os.remove(output)
        os.remove(pslOut)
        os.remove(pslTrash)
        
        return max(scoreS, key=scoreS.get)


def getOrientationNucmer(query, reference, prePath):
    """
    Given two sequences, this function runs nucmer in the background to
    determine their relative orientation
    """
    
    nucmerOutput =  os.path.join(prePath,"out")
    nucmerCall = "nucmer -p {} {} {}".format(nucmerOutput,query,reference)
    showCoordsCall = "show-coords -HlTd {}".format(nucmerOutput+".delta")    
    n = subprocess.Popen([nucmerCall], shell=True, stdout = subprocess.PIPE, \
                        stderr = subprocess.PIPE)
    nReturn = n.wait()
    
    s = subprocess.Popen([showCoordsCall], shell=True, stdout = subprocess.PIPE, \
                        stderr = subprocess.PIPE)
    sReturn = s.wait()
    
    #check if one of them crashed 
    if nReturn != 0 or sReturn != 0:
        return None
    
    scoreS = {"+": 0.0, "-": 0.0}
    
    #parse the alignments
    for line in s.stdout:
        line = line.rstrip().split('\t')
        alnLenUpg, alnLenRef, identity, orientationUpg,orientationRef= \
        line[4],line[5],line[6],line[9], line[10]
        #score the alignment
        score = float(((int(alnLenUpg)+int(alnLenUpg))/2)*float(identity))
        orientation = '+' if (orientationUpg==orientationRef) else '-'
        scoreS[orientation] += score
    
    #delete the temporary file:
    os.remove(nucmerOutput+".delta")
    
    return max(scoreS, key=scoreS.get)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json',  help='Path to the liftOver.json generated by PBJelly',
                    required=True)
    parser.add_argument('-gapInfo', help='Path to the *.gapInfo.bed file created by the PBJelly setup stage',
                    required=True)
    parser.add_argument('-originalRef', help='Path to the originl reference (renamed to reference.fasta.original by PBJelly)',
                    required=True)
    parser.add_argument('-jellyOutRef', help='Path to the output assembly of PBJelly (may be before or after correcting with quiver)',
                    required=True)
    parser.add_argument('-blat', help='Run blat to do the alignments. Nucmer will be used by default.',
                    required=False, action='store_true')
    args = parser.parse_args()    
    
    OUTFILEREF = "upgradeReference.ordered.oriented.fa"
    OUTFILELOOKUP = "lookupTable.txt"
    
    if args.blat:
        aligner = "blat"
    else:
        aligner = "nucmer"
    
    #create 2 lookup tables for the headers:
    #one for ref\d(7) -> Contig\d
    #one for originalHeader = ref\d(7)
    #create dict for final names (ref concatenations instea of "ContigD")
    lookUpRC = {}
    lookUpOR = {}
    lookUpCR = {}
    lookUpRO = {}
    finalNames = {}
    
    refd7Pattern = re.compile("ref\d{7}")
    
    
    liftOver = json.load(open(args.json))
    for key in liftOver:
        refs = set()
        for refId, strand, size in liftOver[key]:
            refs.update(refd7Pattern.findall(refId))
        newName = []
        for ref in refs:
            newName.append(ref)
            lookUpRC[ref] = key
        finalNames[key] = "_".join(newName)
        lookUpCR[key] = refs
            
    #Create originalHeader -> ref\d(7) lookup table
    with open(args.gapInfo) as g:
        for line in g:
            entry = line.split()
            entry = entry[0].split("|")
            refId = entry[-1]
            originalHeader = "|".join(entry[:-1])
            lookUpOR[originalHeader] = refId
            lookUpRO[refId] = originalHeader
            
    #load in the 2 assemblies:
    originalReference = {}
    jellyOutputReference = {}
    
    for scaffold in loadFasta(args.originalRef):
        originalReference[scaffold.get_identifier()] = scaffold
    
    for scaffold in loadFasta(args.jellyOutRef):
        #in case assembly is quiver corrected, get rid of "|quiver" in the name
        k = scaffold.get_identifier().split("|")
        if len(k) > 1:
            lookUpID = "|".join(k[:-1])
        else:
            lookUpID = k[0]
        jellyOutputReference[lookUpID] = scaffold    
    
    #get set of missing scaffolds, in case they were dropped by quiver
    missingScaffoldsC = set()
    missingScaffoldsR = set()
    
    for scaffoldC in lookUpCR:
        if scaffoldC not in jellyOutputReference:
            missingScaffoldsC.add(scaffoldC)
            for scaffoldR in lookUpCR[scaffoldC]:
                missingScaffoldsR.add(scaffoldR)
    
    #get the order of the original assembly:
    scaffoldOrder = loadScaffolOrder(args.originalRef)
    
    #go through each scaffold in the scaffold order list, write the original and
    #the upgrade to a tmp file, align them to get the orientation

    cwd = os.getcwd()
    tmpdir = os.path.join(cwd,"OOtmp")
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    else:
        sys.stderr.write("Directory {} already exists. Rename or delete it.\nRefusing to touch it...Exiting\n".format(tmpdir))
        sys.exit(1)
    
    outfileRef =  os.path.join(cwd,OUTFILEREF)
    outfileLT  =  os.path.join(cwd,OUTFILELOOKUP)
    if os.path.exists(outfileRef):
        sys.stderr.write("File {} already exists. Rename or delete it.\nRefusing to touch it...Exiting\n".format(outfileRef))
        sys.exit(1)
    if os.path.exists(outfileLT):
        sys.stderr.write("File {} already exists. Rename or delete it.\nRefusing to touch it...Exiting\n".format(outfileLT))
        sys.exit(1)    
    
    fR = open(outfileRef,'w')
    
    #keep track of scaffoldC that we already checked
    scaffoldCChecked = set()
    
    for scaffoldO in scaffoldOrder:
        #various scaffoldO may point to the same scaffoldC
        scaffoldR = lookUpOR[scaffoldO]
        if lookUpRC[scaffoldR] in scaffoldCChecked:
            continue
        
        #check if scaff was dropped
        if not scaffoldR in missingScaffoldsR:
            #check if output mapps to more then one input scaffold and select longest
            scaffoldfRecordsR = set()
            for i in lookUpCR[lookUpRC[scaffoldR]]:
                scaffoldfRecordsR.add(originalReference[lookUpRO[i]])
                
            scaffoldRecordR = max(scaffoldfRecordsR)
            scaffoldRecordC = jellyOutputReference[lookUpRC[scaffoldR]]
            
            #f1 = tempfile.NamedTemporaryFile()
            #f2 = tempfile.NamedTemporaryFile()
            tempseqA = os.path.join(tmpdir+"/sequenceA.fa")
            tempseqB = os.path.join(tmpdir+"/sequenceB.fa")
            
            if os.path.exists(tempseqA) or os.path.exists(tempseqB):
                sys.stderr.write("File {} or {}already exists. Rename or delete it.\nRefusing to touch it...Exiting\n".format(tempseqA,tempseqB))
                sys.exit(1) 
            
            f1 = open(tempseqA, 'w')
            f2 = open(tempseqB, 'w')
            f1.write(scaffoldRecordR.get_fasta_string())
            f2.write(scaffoldRecordC.get_fasta_string())
            f1.flush()
            f2.flush()
            f1.close()
            f2.close()
            
            if aligner == 'nucmer':
                orientation = getOrientationNucmer(f1.name,f2.name,tmpdir)
            elif aligner == 'blat':
                orientation = getOrientationBlat(f1.name,f2.name,tmpdir)
                
            if orientation == "-":
                revComp = scaffoldRecordC.get_complement()
                fR.write(revComp.get_fasta_string())
                
            elif orientation == "+":
                fR.write(scaffoldRecordC.get_fasta_string())
                
            elif orientation is None:
                sys.stderr.write("Nucmer crashed trying to align {} and {}, asuming + orientation\n".\
                                 format(scaffoldO,lookUpRC[scaffoldR]))
                fR.write(scaffoldRecordC.get_fasta_string())

            #keep track of scaffolds that we allready output    
            scaffoldCChecked.add(lookUpRC[scaffoldR])
            #remove the tempfiles
            os.remove(tempseqA)
            os.remove(tempseqB)
            
        else: #in case a scaffold was dropped by quiver 
            #sys.stdout.write(originalReference[lookup])
            sys.stderr.write("Scaffold {} ({}) is absent in the jelly output. Reinserting original scaffold\n".\
                             format(lookUpRC[scaffoldR],scaffoldO))
            fR.write(originalReference[scaffoldO].get_fasta_string())
        
        
            
    fR.close()
    os.rmdir(tmpdir)
    sys.exit(0)
