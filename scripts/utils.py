# from https://www.biostars.org/p/344555/

import pysam

# Read in a BAM file
def readBAM():
    bamfile = pysam.AlignmentFile("../data/mapping_BK_SRR_sorted.bam","rb")
    for line in bamfile:
        print(line)

# Read in a SAM file 
def readSAM():
    samfile = pysam.AlignmentFile("outputFile.sam","r")
    for line in samfile:
        print(line)

# Read in a CRAM file. 
def readCRAM(): 
    CRAMfile = pysam.AlignmentFile("outputfile.cram", "rc")


def readVCF(filename):    
    data = []
    header = []
    with open(filename, "r") as f:
        for line in f:
            # get header
            if line.startswith('##'):
                pass
            elif line.startswith('#'):
                line = line[1:-1]
                header = line.split('\t')
            # get data
            else:
                raw = line[0:-1].split('\t')
                dic = {header[i]: raw[i] for i in range(len(raw))}
                data.append(dic)
    return data

def readMpileup(filename):      
    data = []
    header = []
    with open(filename, "r") as f:
        for line in f:
            # get header
            if line.startswith('##'):
                pass
            elif line.startswith('#'):
                line = line[1:-1]
                header = line.split('\t')
            # get data
            else:
                raw = line[0:-1].split('\t')
                dic = {header[i]: raw[i] for i in range(len(raw))}
                data.append(dic)
    return data


def filterVCFByPos(data, rangePos):
    dataPos = []
    for d in data:
        if rangePos[0] < float(d['POS']) < rangePos[1]:
            dataPos.append(d)
            print(d)
    return dataPos

def vcfAnalysis():
    data = readVCF('../data/BK_SRR_sorted_v3.vcf')

    # /product="18S ribosomal RNA", positions 4008..5877
    print('\n\nproduct="18S ribosomal RNA", positions 4008..5877:\n')
    dataPos = filterVCFByPos(data, [4008, 5877])

    # /product="5.8S ribosomal RNA", positions 6878..7034
    print('\n\nproduct="5.8S ribosomal RNA", positions 6878..7034:\n')
    dataPos = filterVCFByPos(data, [6878, 7034])

    # /product="28S ribosomal RNA", positions 8123..12852
    print('\n\nproduct="28S ribosomal RNA", positions 8123..12852:\n')
    dataPos = filterVCFByPos(data, [8123, 12852])

def runCommand(command):
    import subprocess
    proc = subprocess.Popen(command.split(' '), 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT)
    stdout,stderr = proc.communicate()
    print(stdout)
    print(stderr)


def mpileupAnalysis():
    inputFile = '../data/BK_SRRsorted_v2.mpileup'
    quality_thresholds = [20, 25, 30]
    poly_threshold = 1.0
    threads = 4


    for q in quality_thresholds:
        # /product="18S ribosomal RNA", positions 4008..5877
        print('\n\nproduct="18S ribosomal RNA", positions 4008..5877:\n')
        command = './mpileup-parser-v3.py -f %s -t %s -q %s -s 4008 -e 5877 -p %f' % (inputFile, threads, q, poly_threshold)
        runCommand(command)

        # /product="5.8S ribosomal RNA", positions 6878..7034
        print('\n\nproduct="5.8S ribosomal RNA", positions 6878..7034:\n')
        command = './mpileup-parser-v3.py -f %s -t %s -q %s -s 6878 -e 7034 -p %f' % (inputFile, threads, q, poly_threshold)
        runCommand(command)

        # /product="28S ribosomal RNA", positions 8123..12852
        print('\n\nproduct="28S ribosomal RNA", positions 8123..12852:\n')
        command = './mpileup-parser-v3.py -f %s -t %s -q %s -s 8123 -e 12852 -p %f' % (inputFile, threads, q, poly_threshold)
        runCommand(command)



mpileupAnalysis()