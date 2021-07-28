import pysam
import numpy as np
import cv2
from PIL import Image
from multiprocessing import Process, Queue
import random
import glob

class Query:
    def __init__(self,name,cigar,query,qual,mapq,position,indent,ref):
        self.name=name
        self.ref=ref
        self.cigar=cigar
        self.query=query
        self.qual=qual
        self.mapq=mapq
        self.position=position
        self.indent=indent

def making_query_dict(path,pos_chr,pos_start,insertion_length=0):
    term=32
    samfile = pysam.AlignmentFile(path, "rb" )
    pileupdata=samfile.pileup(pos_chr, pos_start, pos_start+1)
    insertion_dict={}
    query_dict={}
    veryfirst=0
    start_index=0
    for pileupcolumn in samfile.pileup(pos_chr, pos_start-term, pos_start+1+term):
        if(pileupcolumn.pos>=pos_start-term and pileupcolumn.pos<pos_start+1+term):
            for pileupread in pileupcolumn.pileups:
                if(pileupread.alignment.query_name not in query_dict.keys()):
                    isinsertion=0
                    cigar=pileupread.alignment.cigartuples
                    cigar_info=''
                    cigar_pos=0
                    query_sequence=pileupread.alignment.query_sequence
                    position=pileupread.alignment.get_reference_positions()
                    if(veryfirst==0):veryfirst=position[0]
                    first_pos=0
                    softclip=0
                    for i,e in enumerate(cigar):
                        if(e[0]==0):
                            if(first_pos==0):first_pos=position[0]
                            cigar_info=cigar_info+'M'*e[1]
                            cigar_pos=cigar_pos+e[1]
                        elif(e[0]==4):
                            softclip=e[1]
                            if(first_pos==0):
                                first_pos=position[0]-e[1]
                                for sftidx in range(0,e[1]):position=[first_pos+sftidx]+position
                            elif(cigar_pos>=len(position)):
                                position=position+[position[-1]+1]
                            else:
                                for sftidx in range(0,e[1]):position=position[0:cigar_pos]+[position[cigar_pos]+1]+position[cigar_pos:len(position)]
                            cigar_info=cigar_info+'S'*e[1]
                            cigar_pos=cigar_pos+e[1]

                        elif(e[0]==5):
                            cigar_info=cigar_info+'H'*e[1]
                            cigar_pos=cigar_pos+e[1]
                        elif(e[0]==1):
                            cigar_info=cigar_info+'I'*e[1]
                            if(position[0]+cigar_pos - veryfirst not in insertion_dict):insertion_dict[pileupread.alignment.query_name]=position[0]+cigar_pos-veryfirst
                            cigar_pos=cigar_pos+1
                            isinsertion=1
                        elif(e[0]==2):
                            cigar_info=cigar_info+'D'*e[1]
                            query_sequence=query_sequence[0:cigar_pos]+'D'*e[1]+query_sequence[cigar_pos:len(query_sequence)]
                            position=position[0:cigar_pos]+[0]*e[1]+position[cigar_pos:len(query_sequence)]
                            cigar_pos=cigar_pos+e[1]
                    indent=0
                    query=''
                    start_softclip,end_softclip=0,0
                    start_deletion,end_deletion=0,0
                    ref=pileupread.alignment.get_reference_sequence()
                    qual=str(pileupread.alignment.qual)
                    for q in range(0,len(query_sequence)):
                        nt=query_sequence[q]
                        if(cigar_info[q]=='S'):
                            query+='S'
                            if(q==start_softclip):start_softclip+=1
                            else:end_softclip+=1
                        elif(cigar_info[q]=='I'):continue
                        elif(cigar_info[q]=='D'):
                            query+=nt
                            qual=qual[0:q]+' '+qual[q:len(qual)]
                        else:query+=nt
                    ref=' '*(start_softclip)+ref+' '*(end_softclip)
                    if(first_pos-veryfirst>=0):
                        query=' '*(first_pos-veryfirst)+query
                        ref=' '*(first_pos-veryfirst)+ref
                        qual=' '*(first_pos-veryfirst)+qual
                        for i in range(first_pos-veryfirst,0):position=[position[0]-1]+position
                    else:
                        query=query[veryfirst-first_pos:len(query)]
                        ref=ref[veryfirst-first_pos:len(ref)]
                        qual=qual[veryfirst-first_pos:len(qual)]
                        position=position[veryfirst-first_pos:len(position)]
                    if(start_index==0):
                        if(pos_start not in position):start_index=position.index(pos_start-term)
                        else:start_index=position.index(pos_start)-term
                    qual=qual[start_index:start_index+term*2]
                    query=query[start_index:start_index+term*2]
                    ref=ref[start_index:start_index+term*2]
                    # print(query)
                    if(isinsertion==1):
                        if(insertion_dict[pileupread.alignment.query_name]<start_index or insertion_dict[pileupread.alignment.query_name]>start_index+30):del(insertion_dict[pileupread.alignment.query_name])
                        else:insertion_dict[pileupread.alignment.query_name]=insertion_dict[pileupread.alignment.query_name]-start_index
                    query_dict[pileupread.alignment.query_name]=Query(pileupread.alignment.query_name,cigar_info,query,qual,pileupread.alignment.mapq,position,indent,ref)
    return([query_dict,insertion_dict,veryfirst])


def making_bam_array(query_dict,insertion_dict):
    bam_array = np.zeros([64,64,3])

    scoreDict = {'A':110, 'T': 140, 'G': 170, 'C': 200, 'D' : 230, ' ': 0, 'S':20, 'N':100, '*' : 100, "I" : 80 , "D" : 10}

    readct = 0
    for read in query_dict.keys():
        for breakpoint in list(set(insertion_dict.values())):
            if(breakpoint<len(query_dict[read].query) and query_dict[read].query[breakpoint]!=' '):
                if(read in insertion_dict.keys()):
                    if(insertion_dict[read]==breakpoint):
                        query_dict[read].query=query_dict[read].query[0:breakpoint]+'I'+query_dict[read].query[breakpoint:len(query_dict[read].query)]
                        query_dict[read].ref=query_dict[read].ref[0:breakpoint]+'*'+query_dict[read].ref[breakpoint:len(query_dict[read].ref)]
                        query_dict[read].qual=query_dict[read].qual[0:breakpoint]+' '+query_dict[read].qual[breakpoint:len(query_dict[read].qual)]
                    else:
                        query_dict[read].query=query_dict[read].query[0:breakpoint]+'*'+query_dict[read].query[breakpoint:len(query_dict[read].query)]
                        query_dict[read].ref=query_dict[read].ref[0:breakpoint]+'*'+query_dict[read].ref[breakpoint:len(query_dict[read].ref)]
                        query_dict[read].qual=query_dict[read].qual[0:breakpoint]+' '+query_dict[read].qual[breakpoint:len(query_dict[read].qual)]
                else:
                    query_dict[read].query=query_dict[read].query[0:breakpoint]+'*'+query_dict[read].query[breakpoint:len(query_dict[read].query)]
                    query_dict[read].ref=query_dict[read].ref[0:breakpoint]+'*'+query_dict[read].ref[breakpoint:len(query_dict[read].ref)]
                    query_dict[read].qual=query_dict[read].qual[0:breakpoint]+' '+query_dict[read].qual[breakpoint:len(query_dict[read].qual)]
            else:
                query_dict[read].query=query_dict[read].query[0:breakpoint]+' '+query_dict[read].query[breakpoint:len(query_dict[read].query)]
                query_dict[read].ref=query_dict[read].ref[0:breakpoint]+' '+query_dict[read].ref[breakpoint:len(query_dict[read].ref)]
                query_dict[read].qual=query_dict[read].qual[0:breakpoint]+' '+query_dict[read].qual[breakpoint:len(query_dict[read].qual)]

    for read in query_dict.keys():
        for pos in range(0, len(query_dict[read].query) ):
            bam_array[readct, pos,0] = scoreDict[query_dict[read].query[pos].upper()]
            bam_array[readct, pos,1] = scoreDict[query_dict[read].ref[pos].upper()]
        readct += 1
    return bam_array




def work(id, sliced_var_list, result):
    for term in sliced_var_list :
        try:
            random_var = random.sample(list(range(0,64)),term[0])
            if(term[0] == 2): output_path = "/home/sylee/GOTCHA/training/" + "snv"
            elif(term[0] == 20): output_path = "/home/sylee/GOTCHA/training/" + "del"
            elif(term[0] == 50): output_path = "/home/sylee/GOTCHA/training/" + "ins"
            for r in random_var:
                pos_start = int(term[-1][1]) + r
                print(term[4])
                mother= making_query_dict(term[3],term[-1][0],pos_start)
                father= making_query_dict(term[2],term[-1][0],pos_start)
                proband= making_query_dict(term[1],term[-1][0],pos_start)

                family_dict={**mother[1],**father[1],**proband[1]}

                print(mother[2],father[2],proband[2])
                mother_array=making_bam_array(mother[0],family_dict)
                proband_array=making_bam_array(proband[0],family_dict)
                father_array=making_bam_array(father[0],family_dict)
                
                cv2.imwrite(output_path + "/father/%s_s.%d_vp.%s"%(term[-1][0],pos_start,term[-1][1]) + ".png",father_array)
                cv2.imwrite(output_path + "/mother/%s_s.%d_vp.%s"%(term[-1][0],pos_start,term[-1][1]) + ".png",mother_array)
                cv2.imwrite(output_path + "/proband/%s_s.%d_vp.%s"%(term[-1][0],pos_start,term[-1][1]) + ".png",proband_array)
        except ValueError:
            continue
        except IndexError:
            continue
    return id


class inputData:
    def __init__(self, proband, father, mother, bam_path):
        self.proband = proband
        self.father = father
        self.mother = mother
        self.bam_path = bam_path
        self.snv_ti = ''
        self.ins = ''
        self.snv_tv = ''
        self.deletion = ''
        

    def show_info(self):
        print(self.proband, self.father, self.mother, self.bam_path)

    def path_match(self, variant_file_list):
        for file in variant_file_list:
            if "Ti" in file : self.snv_ti = open(file,'r').readlines()
            elif "Tv" in file : self.snv_tv = open(file,'r').readlines()
            elif "DEL" in file : self.deletion = open(file,'r').readlines()
            elif "INS" in file : self.ins = open(file,'r').readlines()


if __name__ == "__main__":
    print("###########################################################################")
    table_list = ["/EXTDATA/mhchoi/Denovo_mutation_project/Control_trio_vs_ABS_Analysis/Directly_blood_draw/2020/2nd" ,
    "/EXTDATA/mhchoi/Denovo_mutation_project/Control_trio_vs_ABS_Analysis/Directly_blood_draw/2021/1st",
    "/EXTDATA/mhchoi/Denovo_mutation_project/Control_trio_vs_ABS_Analysis/RareDisease_dataSet",
    "/EXTDATA/mhchoi/Denovo_mutation_project/Atomic_Bomb_Survivors_Analysis/2021"]

    family_dnmfile_dict = dict()

    var_list = []

    for path in table_list:
        for specified_path in glob.glob(path + '/*'):
            bamlist = glob.glob(specified_path + "/**/*.final.bam", recursive = True)
            for more_specified_path in glob.glob(specified_path +'/*.T*'):
                

                if 'UNIST' in more_specified_path : 
                    tag = '_'.join(more_specified_path.split('/')[-3:])
                    more_specified_path = more_specified_path + '/1.Processing'
                    variant_file_list = glob.glob(more_specified_path + '/dSNV_dINDEL_final_table/*TRUE.table' )
                
                elif '2021' in more_specified_path:
                    tag = '_'.join(more_specified_path.split('/')[-4:])
                    variant_file_list = glob.glob(more_specified_path + "/dSNV_dINDEL_final_table/*TRUE.table")
                
                else : 
                    tag = '_'.join(more_specified_path.split('/')[-4:])
                    variant_file_list = glob.glob(more_specified_path + '/Reproducibility_test/dSNV_dINDEL_final_table/*TRUE.table' )
                
                #print(variant_file_list)

                ped = open(glob.glob(more_specified_path + '/*.ped')[0], 'r').readlines()[-1].split('\t')
                
                mom_path = ""
                dad_path = ""
                proband_path = ""
                if(len(bamlist) != 0) :
                    for path in bamlist:
                        if ped[1] in path: proband_path = path
                        elif ped[2] in path : dad_path = path
                        elif ped[3] in path : mom_path = path
                    
                    for var in variant_file_list:
                        file = open(var,'r').readlines()
                        if ".dSNV" in var :
                            for line in file:
                                var_list.append([2, proband_path, dad_path, mom_path, line.split('\t')[0:2]])
                        elif ".dDEL" in var :
                            for line in file:
                                var_list.append([20, proband_path, dad_path, mom_path, line.split('\t')[0:2]])
                        elif ".dINS" in var :
                            for line in file:
                                var_list.append([50, proband_path, dad_path, mom_path, line.split('\t')[0:2]])

    print(len([x for x in var_list if x[0] == 2]))
    print(len([x for x in var_list if x[0] == 20]))
    print(len([x for x in var_list if x[0] == 50]))

    
    result=Queue()
    procs=[]
    START,END=0,len(var_list)
    interval = END // 10
    for term in range(0,10):
        print(interval * term)
        procs.append(Process(target=work, args=(term, var_list[interval * term : interval * (term + 1)], result)))
        procs[-1].start()
    for q in range(0,10):
        procs[q].join()
