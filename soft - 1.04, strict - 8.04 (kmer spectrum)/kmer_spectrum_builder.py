# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 11:22:00 2018

@author: Orlov
"""

'''
don't see any biological or computational need to write reverse complement of the found kmer
so just commenting corresponding part ot the code
'''
from Bio import Seq
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt


class kmer:
    sequence=''
    
    def __init__(self, seq):
        self.counter= 1
        self.sequence=seq
        
    def increase(self):
        self.counter+=1
    
    def show_info(self):
        print(self.sequence, 'occurs ', self.counter, 'times')
    
class kmer_spectrum_builder:
    
    def __init__(self, name, file):
        self.name=name
        self.filepath=file
        
    def analyse_file(self, k_length , quality_tres):
        self.kmer_length=k_length
        self.q_min=quality_tres
        self.kmer_dict={}

        for rec in SeqIO.parse(self.filepath, 'fastq'):
            curr_seq= str(rec.seq).upper()
            curr_qual=rec.letter_annotations['phred_quality']
            seq_len=len(curr_seq)
            for index in range(seq_len-self.kmer_length+1):
                start=index
                end=index+self.kmer_length
                if all(list(map(lambda x: x>self.q_min, curr_qual[start:end]))):
                    subseq = curr_seq[start:end]
                    ##rc_subseq=Seq.Seq(subseq).reverse_complement()
                    if subseq in self.kmer_dict.keys():
                        self.kmer_dict[subseq].increase()
                        ##self.kmer_dict[rc_subseq].increase()
                    else:
                        self.kmer_dict[subseq]=kmer(subseq)
                        ##self.kmer_dict[rc_subseq]=kmer(rc_subseq)
        
    def spectrum_data(self):
        self.trans_dict=defaultdict(int)
        self.x_freq, self.y_dist = [], []
        for key in self.kmer_dict.keys():
            self.trans_dict[self.kmer_dict[key].counter]+=1
        for key in sorted(self.trans_dict.keys()):
            self.x_freq.append(key)
            self.y_dist.append(self.trans_dict[key])
        self.trans_dict.clear()
        self.kmer_dict.clear()
        self.data_length=len(self.x_freq)
        
        ##assessing threshold
        step=int(self.x_freq[-1]/self.data_length)
        av_prev=self.y_dist[0]
        for i in range(self.data_length-step-1):
            av_cur=0
            for j in range(step):
                av_cur+=self.y_dist[i+j]
            av_cur=av_cur/step
            if av_cur*0.96>av_prev:
                self.noise_tres=self.x_freq[int(i+step/2)]
                break
            else:
                av_prev=av_cur
                    
            
    def draw_spectrum(self):
        ##setting axis
        ind=self.noise_tres
        while self.y_dist[ind]>0.01*self.y_dist[self.noise_tres]:
            ind+=1
        x_axis_max=self.x_freq[ind]
        self.pike_coord=self.y_dist.index(max(self.y_dist[self.noise_tres:self.data_length]))
        y_axis_max=1.1*self.y_dist[self.pike_coord]        
        plt.axis([-0.01, x_axis_max, -0.01, y_axis_max])
        plt.plot(self.x_freq, self.y_dist)
        plt.plot([self.noise_tres, self.noise_tres], [0, y_axis_max])
        
        
    def assess_genome_size(self):
        threshold=[0, self.noise_tres]
        c=self.x_freq[self.pike_coord]
        for start_ind in threshold:
            kmer_number=0
            for i in range(start_ind, self.data_length):
                kmer_number+=self.y_dist[i]*self.x_freq[i]
            if start_ind==0:
                self.g_size_w_noise=kmer_number/c
                print('Genome size without noise cut off', self.g_size_w_noise)
            else:
                self.g_size_corrected=kmer_number/c
                print('Corrected genome size with notoriuous noise cutted off', self.g_size_corrected)
        
        
    def clear_RAM(self):
        self.kmer_dict.clear()
        self.trans_dict.clear()


file_path='C:/Users/Orlov/Desktop/Bioinformatics institute/Homeworks/2nd sem/Python/soft - 1.04, strict - 8.04 (kmer spectrum)/test_kmer.fastq'

sp1=kmer_spectrum_builder('first_spectrum', file_path)
sp1.analyse_file(17, 0)

sp1.spectrum_data()

sp1.draw_spectrum()

sp1.assess_genome_size()



sp1.clear_RAM()

sp1.kmer_dict.clear()

