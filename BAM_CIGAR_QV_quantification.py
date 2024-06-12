#!/usr/bin/env python
'''

Purporse: 

Given a target VCF and BAM file, calculate the QV of each read and mutation info

NOTE: require python>3.7, bcftools in root path

'''


import pysam, re, argparse, logging
import numpy as np
from pysam import VariantFile

logging.basicConfig(
    level=logging.DEBUG,
    format="[%(levelname)s][%(asctime)s]%(message)s",
    handlers=[
        logging.FileHandler('QV_quantification.py.log', mode='w'),
        logging.StreamHandler()
    ]
)

def VCF(vcf):
    vcf_in=VariantFile(vcf)

    for rec in vcf_in:
        print(rec.pos)
        
def var_count(cigar_tuples):
    snp_count=0
    indel_count=0

    for c in cigar_tuples:
        operation=int(c[0])
        length=int(c[1])
        if operation == 8:
            snp_count=snp_count+length
        if operation in [1,2]:
            indel_count=indel_count+length

    return((snp_count, indel_count))

def getOverlap(a,b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def is_overlap(start, end):
    '''
       gene coordinates: 
       288-1829
       2270-5173
    '''
    reg1=getOverlap([start, end],[288,1829]) * 100 / (1829-288)
    reg2=getOverlap([start, end],[2270,5173]) * 100 / (5173-2270)
    if (reg1 > 95) | (reg2 > 95):
        return True
    return False

def write_table(d, header=None, outfn='test.txt'):
    '''
        header: list of columns, e.g.['A','B','C']
    '''
    with open(outfn,"w") as fh:
        if header:
            fh.write('\t'.join(header)+'\n')
        for c1 in d:
            for c2 in d[c1]:
                val=d[c1][c2]
                res=[c1,c2,val]
                fh.write('\t'.join(map(str,res))+'\n')

def QV(bam,pos,prefix=None):
    samfile=pysam.AlignmentFile(bam, 'rb') 

    SNP_QV_res=dict()
    Indel_QV_res=dict()
    SNP_pass_res=dict()
    Indel_pass_res=dict()
    QV_pass_res=dict()


    for read in samfile:
         QV=round(np.mean(read.query_qualities),0)
         SNPn, indel_n=var_count(read.cigartuples)
         overlapped=is_overlap(read.reference_start, read.reference_end)
         #print("Read: {}".format(read.query_name))
         #print("QV: {}".format(QV))
         pass_n=None
         if read.has_tag('np'):
             pass_n=read.get_tag('np')
         #print("is_overlap: {}".format(overlapped))
         #print("np: {}".format(pass_n))

         ### output
         if overlapped:
             if(QV not in SNP_QV_res): 
                 SNP_QV_res[QV]=dict()
                 Indel_QV_res[QV]=dict()


             SNP_QV_res[QV][SNPn]= 1 if SNPn not in SNP_QV_res[QV] else SNP_QV_res[QV][SNPn] + 1 
             Indel_QV_res[QV][indel_n]= 1 if indel_n not in Indel_QV_res[QV] else Indel_QV_res[QV][indel_n] + 1

             if pass_n not in SNP_pass_res:
                 SNP_pass_res[pass_n]=dict()
                 Indel_pass_res[pass_n]=dict()
                 QV_pass_res[pass_n]=dict()

             SNP_pass_res[pass_n][SNPn]= 1 if SNPn not in SNP_pass_res[pass_n] else SNP_pass_res[pass_n][SNPn] + 1
             Indel_pass_res[pass_n][indel_n]= 1 if indel_n not in Indel_pass_res[pass_n] else Indel_pass_res[pass_n][indel_n] + 1
             QV_pass_res[pass_n][QV]= 1 if QV not in QV_pass_res[pass_n] else QV_pass_res[pass_n][QV] + 1

             
             ### write to file
             write_table(SNP_QV_res, header=['QV','SNPn','ReadN'], outfn='.'.join([prefix,'SNP_QV','txt']))
             write_table(Indel_QV_res, header=['QV','INDELn','ReadN'], outfn='.'.join([prefix,'INDEL_QV','txt']))
             write_table(SNP_pass_res, header=['PASSn','SNPn','ReadN'], outfn='.'.join([prefix,'PASS_SNP','txt']))
             write_table(Indel_pass_res, header=['PASSn','INDELn','ReadN'], outfn='.'.join([prefix,'PASS_INDEL','txt']))
             write_table(QV_pass_res, header=['PASSn','QV','ReadN'], outfn='.'.join([prefix,'PASS_QV','txt']))

def main():

    parser = argparse.ArgumentParser(description='Pacbio BAM')
    parser.add_argument('-b', '--bam', help='BAM', required=True)
    parser.add_argument('-pre', '--prefix', help='prefix',required=True, default='test' )

    #parser.add_argument('-v' , '--vcf', help='VCF', required=True)

    args=parser.parse_args()
    print("[BAM] {}".format(args.bam))
    #print("[VCF] {}".format(args.vcf))
    #pos=VCF(args.vcf)
    pos=None
    QV(args.bam, pos, args.prefix)



if __name__ == "__main__":
    main()




