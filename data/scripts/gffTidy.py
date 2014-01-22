#! /usr/bin/python
import os, sys

#gffTidy takes a GFF formatted file as input and produces a tidied up
#gff file as output, kepping only those features listed below.

usage = [ "python gffTidy.py [options] <input> <output>" ]
example = ["python gffTidy.py -nc ref_Amel_4.5_top_level.gff3 apis.new.gff3"]

options = {'-nc':['Include non-coding features ',False], '-features':['Include features from text file', False], '-quiet':['Turn off GFF summary statistics', False]}

features = ['CDS', 'exon', 'gene', 'misc_RNA', 'intron', "5'UTR", "3'UTR", 'mRNA', 'region']

attributes = ['Dbxref', 'gene', 'genome', 'mol_type', 'strain', 'product']

noncoding = ['ncRNA', 'microRNA', 'tRNA', 'rRNA']

ids = {'gene':0, 'rna':0, 'id':0}

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
            if feat not in features: features.append(feat)
    if options['-nc'][1] == True:
        for x in noncoding:
            if x not in features: features.append(x)

    start, geneFeatures=False, []
    for line in f:
        if not line.startswith('#'):
            splitter = line.split()
            if splitter and len(splitter) == 9 and splitter[2] in features:
                if splitter[2] == 'gene':
                    if start: 
                        writeGene(o, geneFeatures)
                        geneFeatures = [] 
                    else: start = True
            geneFeatures.append(splitter)
        else:
            if line.startswith('##FASTA') or line.startswith('>'): break
            if not line.startswith('\n'): o.write(line)
    writeGene(o, geneFeatures)
    f.close()
    o.close()
    gffStat(outputFile)

def writeGene(o, geneFeatures):
    never = False
    mRNA = False
    region = False
    coding = True
    newL = ''
    for x in geneFeatures:
        if not options['-nc'][1]:
            for nc in noncoding:
                if nc in x[8] or nc in x[1] or nc in x[2]: coding = False
        if coding and x[2] in features:
            splitee = x[8].split(';')
            if x[2] == 'gene':
                geneSplit = x
                start, stop = x[3], x[4]
                ids['gene'] += 1
                geneID = 'gene' + str(ids['gene'])
                mRNA = True
                attrib = 'ID='+geneID + ';'
            elif x[2] == 'mRNA' or x[2] == 'ncRNA':
                if mRNA: mRNA = False
                ids['rna'] += 1
                attrib = 'ID=rna'+str(ids['rna']) + ';Parent=gene'+str(ids['gene'])+';'
            elif x[2] == 'CDS' or x[2] == 'exon':
                ids['id'] += 1
                attrib = 'ID=id' + str(ids['id']) + ';Parent=rna'+str(ids['rna']) + ';'
            elif x[2] in features:
                ids['id'] += 1
                attrib = 'ID=id' + str(ids['id']) + ';'
            for y in x[:8]: newL += y + '\t'
            for column in splitee:
                for attribute in attributes:
                    if column.startswith(attribute): attrib += column + ';'
            newL += attrib + '\n'
    if coding:
        o.write(newL)
        
def gffStat(fileHandle):
    f = open(fileHandle)
    genes, sources, seqs, types, geneCount, seqids = {}, {}, {}, {}, 0, 0
    for line in f:
        if not line.startswith('#'):
            x = line.split()
            seqid= x[0]
            source = x[1]
            type = x[2]
            if seqid not in seqs:
                seqs[seqid] = 1
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
    print "    %d Scaffolds" % len(seqs)
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
