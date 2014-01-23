#! /usr/bin/python
import os, sys

#gffTidy takes a GFF formatted file as input and produces a tidied up
#gff file as output, kepping only those features listed below.

usage = [ "python gffTidy.py [options] <input> <output>" ]
example = [ "python gffTidy.py -nc ref_Amel_4.5_top_level.gff3 apis.new.gff3" ]

options = {'-nc':['Include non-coding features ',False], '-features':['Include features from text file', False], '-quiet':['Turn off GFF summary statistics', False]}

features = ['CDS', 'exon', 'gene', 'misc_RNA', 'intron', "5'UTR", "3'UTR", 'mRNA', 'region', 'transcript', 'primary_transcript']

attributes = ['Dbxref', 'gene', 'genome', 'mol_type', 'strain', 'product', 'ID', 'Parent']

noncoding = ['ncRNA', 'miRNA', 'tRNA', 'rRNA', 'lincRNA']

ids = {'gene':0, 'rna':0, 'id':0, 'cds':0}

seq_region = {0:0}

def main():
    if len(sys.argv) < 3:
        print "Please specify input and output gff3 files."
        printUsage()
        sys.exit()
    for x in range(0, len(sys.argv)):
        if sys.argv[x] == '-nc': options['-nc'][1] = True
        if sys.argv[x] == '-features': 
            options['-features'][1] = True
            featureFile = sys.argv[x+1]
        if sys.argv[x] == '-quiet': options['-quiet'][1] = True
        if sys.argv[x] == '-help' or sys.argv[x] == '-h':
            printUsage()
            sys.exit()
    inputFile = sys.argv[-2]
    outputFile = sys.argv[-1]
    if not inputFile.startswith('-') and not outputFile.startswith('-'):
        if not gffCheck(inputFile):
            print "Input file was not valid GFF!"
            sys.exit()
    try:
        f = open(inputFile)
    except IOError as e:
        print "Can't open file: " + inputFile
        printUsage()
        sys.exit()

    try:
        o = open(outputFile, 'w')
    except IOError as e:
        print "Can't open file: " + outputFile
        print Usage()
        sys.exit()

    if options['-quiet'][1] == False: 
        gffStat(inputFile)
    if options['-features'][1] == True:
        try:
            f1 = open(featureFile)
        except IOError as e:
            print "Can't open feature file: " + featureFile
        for line in f1:
            feat = line.strip()
            for x in features: features.remove(x)
            if feat not in features: features.append(feat)

    if options['-nc'][1] == True:
        for x in noncoding:
            if x not in features: features.append(x)

    start, geneFeatures, region = False, [], False
    for line in f:
        if not line.startswith('#'):
            splitter = line.strip().split('\t')
            if splitter[2] == 'region': region = True
            #if splitter and len(splitter) == 9 and splitter[2] in features:
            if splitter:
                if splitter[2] == 'gene' or region:
                    if start and geneFeatures: 
                        writeGene(o, geneFeatures)
                        geneFeatures = [] 
                    else: start = True
            geneFeatures.append(splitter)
            region = False
        else:
            if line.startswith('##FASTA') or line.startswith('>'): break
            elif 'sequence-region' in line:
                 o.write(line)
                 region = True
            elif not line.startswith('\n'): 
                region = False
                o.write(line)
    if geneFeatures:
        writeGene(o, geneFeatures)
    f.close()
    o.close()
    gffStat(outputFile)

def writeGene(o, geneFeatures):
    o1 = open('misfits.txt', 'a')
    coding,block, geneSplit = True, '', []
    rnas = ['mRNA', 'ncRNA', 'transcript', 'misc_RNA', 'tRNA']
    for x in geneFeatures:
        newL = ''
        splitee = x[8].split(';')
        if not options['-nc'][1]:
            for nc in noncoding:
                if nc in x[8] or nc in x[1] or nc in x[2]: coding = False
        if x[2] in features and coding:
            for i in x[:8]:
                newL += i + '\t'
            for each in splitee:
                for attribute in attributes:
                    if attribute in each and attribute not in newL:
                        newL += each + ';'
            if newL[-1] == ';': newL = newL[:-1]
            if newL[-1] == '\t': newL = newL[:-1]
            if newL[-1] != '\n': newL += '\n'
            block += newL
    o.write(block)

def gffStat(fileHandle):
    f = open(fileHandle)
    genes, sources, seqs, types, geneCount, seqids = {}, {}, {}, {}, 0, 0
    for line in f:
        if not line.startswith('#') and not line.startswith('\n'):
            x = line.split()
            
            seqid= x[0]
            source = x[1]
            type = x[2]
            if seqid not in seqs:
                seqs[seqid] = 0
            seqs[seqid] += 1
            if source not in sources:
                sources[source] = 0
            sources[source] += 1
            if type not in types:
                types[type]= 0
            types[type] += 1
            if type == 'gene':
                geneCount += 1
                ID = x[8].split(';')[0].split('ID=')[1]
                genes[ID] = {'rna':0, 'exon':0}
            if 'RNA' in type:
                genes[ID]['rna'] += 1
            elif 'exon' in type:
                genes[ID]['exon'] += 1

    sumRNA, sumExons = 0,0
    for x,y in genes.iteritems():
        sumRNA += y['rna']
        sumExons += y['exon']
    print "GFF Report: " + fileHandle
    print "    %d Sequences" % len(seqs)
    print "    %d Genes" % geneCount
    print "    %.4f Average mRNAs per Gene" % ( float(sumRNA) / float(geneCount) )
    print "    %.4f Average exons per Gene" % ( float(sumExons) / float(geneCount) )
    srcList, typeList = '\t', '\t'
    for x,y in sources.iteritems(): srcList += x +'('+ str(y)+')' + ','
    print "    sources:\n",srcList[:-1]
    for x,y in types.iteritems(): typeList += x + '(' + str(y) + ')' + ','
    print "    types:\n", typeList[:-1]
    print '\n'


def printUsage():
    print "Usage: ", usage[0]
    print "[options]"
    for x,y in options.iteritems():
         if not x == '-features': print '  '+x+'\n    '+ y[0] + '\tdefault=' + str(y[1])
         else: print '  '+x+' <text file>\n    '+y[0] + '\tdefault=' + str(y[1])

def gffCheck(inputFile):
    try:
        f = open(inputFile)
    except IOError as e:
        print "Can't open file: " + inputFile
        printUsage()
        sys.exit()
    gffHeader = False
    nineFields = False
    num = 0
    for line in f:
        num += 1
        if line.startswith('#'):
            if line.startswith('#gff-version 3') or line.startswith('##gff-version 3') or 'gff-version' in line:
                gffHeader = True
        elif len(line.split()) == 9:
            nineFields = True
            break
        if num > 100: break
    if gffHeader or nineFields: return True
    else: 
        if not gffHeader: print "No #gff-version 3 line found."
        return False 
     


main()   
