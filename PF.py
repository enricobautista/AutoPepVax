import subprocess
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from splinter import Browser
import math
import time
import pandas as pd
import os
import re
from mhcflurry import Class1AffinityPredictor
import warnings
import pickle


class peptide:
    def __init__(self, seq=None, wt_seq=None, cancer=None, mut=None, allele=None, immunogenicity=None, alt_immuno=None, antigenicity=None, halflife=None, toxicity=None, inf_gamma=None, allergenicity=None, isoelectric_point=None, instability_index=None, aliphatic_index=None, GRAVY_score=None, hydropathicity=None, MHCFlurry_binding=None, MHCFlurry_rank=None, MHCFlurry_wt_mut_rank=None, NetMHC_binding=None, NetMHC_rank=None, NetMHC_stab=None, score=None, id=None):
        self.seq = seq
        self.wt_seq = wt_seq
        self.cancer = cancer
        self.mut = mut
        self.allele = allele
        self.immunogenicity = immunogenicity
        self.alt_immuno = alt_immuno
        self.antigenicity = antigenicity
        self.halflife = halflife
        self.toxicity = toxicity
        self.inf_gamma = inf_gamma
        self.allergenicity = allergenicity
        self.isoelectric_point = isoelectric_point
        self.instability_index = instability_index
        self.aliphatic_index = aliphatic_index
        self.GRAVY_score = GRAVY_score
        self.hydropathicity = hydropathicity
        self.MHCFlurry_binding = MHCFlurry_binding
        self.MHCFlurry_rank = MHCFlurry_rank
        self.MHCFlurry_wt_mut_rank = MHCFlurry_wt_mut_rank
        self.NetMHC_binding = NetMHC_binding
        self.NetMHC_rank = NetMHC_rank
        self.NetMHC_stab = NetMHC_stab
        self.score = score
        self.id = id
    #this function allows us to print out the properties of each peptide
    def __str__(self):
        readout = 'Sequence: '+str(self.seq)+'\nWild Type Seqeunce: '+str(self.wt_seq)+'\nCancer: '+str(self.cancer)+'\nMutation: '+str(self.mut)+'\nAllele: '+str(self.allele)+'\nImmunogenicity: '+str(self.immunogenicity)+'\nAlternative Immunogenicity: '+str(self.alt_immuno)+'\nAntigenicity: '+str(self.antigenicity)+'\nHalf-Life: '+str(self.halflife)+' hours'+'\nToxicity: '+str(self.toxicity)+'\nINF_gamma: '+str(self.inf_gamma)+'\nAllergenicity: '+str(self.allergenicity)+'\nIsoelectric Point: '+str(self.isoelectric_point)+'\nInstability Index: '+str(self.instability_index)+'\nAliphatic Index: '+str(self.aliphatic_index)+'\nGRAVY Score: '+str(self.GRAVY_score)+'\nHydropathicity: '+str(self.hydropathicity)+'\nMHCFlurry Binding: '+str(self.MHCFlurry_binding)+'nM'+'\nMHCFlurry Rank: '+str(self.MHCFlurry_rank)+'\nMHCFlurry WT:MUT Rank: '+str(self.MHCFlurry_wt_mut_rank)+'\nNetMHC Binding: '+str(self.NetMHC_binding)+'nM'+'\nNetMHC Rank: '+str(self.NetMHC_rank)+'\nNetMHC Stability: '+str(self.NetMHC_stab)+'\nScore: '+str(self.score)+'\nID: '+str(self.id)+'\n\n'
        return readout

def to_dataframe(cancer_peptides):
    data = []
    columns = ['Sequence', 'Wild Type Seqeunce', 'Cancer', 'Mutation', 'Allele', 'Immunogenicity', 'Alternative Immunogenicity', 'Antigenicity', 'Half-life', 'Toxicity', 'INFgamma', 'Allergenicity', 'Isoelectric Point', 'Instability Index', 'Aliphatic Index', 'Gravy Score', 'Hydropathicity', 'MHCFlurry Binding', 'MHCFlurry Rank', 'MHCFlurry WT:MUT Rank', 'NetMHC Binding', 'NetMHC Rank', 'NetMHC Stability', 'Score', 'ID']
    for peptide in cancer_peptides:
        data. append([peptide.seq,
                      peptide.wt_seq,
                      peptide.cancer,
                      peptide.mut,
                      peptide.allele,
                      peptide.immunogenicity,
                      peptide.alt_immuno,
                      peptide.antigenicity,
                      peptide.halflife,
                      peptide.toxicity,
                      peptide.inf_gamma,
                      peptide.allergenicity,
                      peptide.isoelectric_point,
                      peptide.instability_index,
                      peptide.aliphatic_index,
                      peptide.GRAVY_score,
                      peptide.hydropathicity,
                      peptide.MHCFlurry_binding,
                      peptide.MHCFlurry_rank,
                      peptide.MHCFlurry_wt_mut_rank,
                      peptide.NetMHC_binding,
                      peptide.NetMHC_rank,
                      peptide.NetMHC_stab,
                      peptide.score,
                      peptide.id])
    return pd.DataFrame(data, columns = columns)

def to_peptide_list(csv_filename):
    peptide_list = []
    df = pd.read_csv(csv_filename)
    for i in range(len(df['Sequence'])):
        row = df.loc[i]
        peptide_list.append(peptide(*row[1:]))
    return peptide_list

def find_position(mis):
    if len(mis)==3:
        pos = int(mis[1])
    if len(mis)==4:
        pos = int(mis[1:3])
    if len(mis)==5:
        pos = int(mis[1:4])
    if len(mis)==6:
        pos = int(mis[1:5])
    return pos

def short_form_mutation(og_seq, missense, length):
    position = find_position(missense)
    amino_acid = missense[-1]
    og_amino_acid = missense[0]
    mut = og_seq[position-length:position-1] + amino_acid + og_seq[position:position+length-1]
    wt = og_seq[position-length:position+length-1]
    return mut, wt

def unique_sequences(cancer_name, cancer_peptides, FASTA=False):
    unique_seq = []
    for peptide in cancer_peptides:
        if peptide.seq not in unique_seq:
            unique_seq.append(peptide.seq)
    text_input = ''
    for seq in unique_seq:
        if FASTA:
            text_input += '>\n'
        text_input += seq
        text_input += '\n'
    sequences_filename = 'Data/'+cancer_name+'/'+cancer_name+' Sequences.txt'
    sequences_file = open(sequences_filename, 'w')
    sequences_file.write(text_input)
    sequences_file.close()
    return text_input

def remove_non_mutated(mutations, og_seq, CD8_csv_filename):
    epitopes = []
    for mis in mutations:
        seq,wt_seq = short_form_mutation(og_seq,mis,10)
        for i in range(len(seq)):
            if len(seq[i:i+10])==10:
                epitopes.append(seq[i:i+10])
            if len(seq[i:i+9])==9 and i!=0 and i+9!=len(seq):
                epitopes.append(seq[i:i+9])
    df = pd.read_csv(CD8_csv_filename)
    df = df[df['Sequence'].isin(epitopes)].reset_index(drop=True)
    df = df.drop(columns=['Unnamed: 0'], errors='ignore')
    df.to_csv(CD8_csv_filename)

def get_NetMHCI(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = []
    MHCI_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01']
    for epi_len in [9,10]:
        for missense in mutations:
            url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
            seq, wt_seq = short_form_mutation(og_seq, missense, epi_len)
            alleles = ''
            for i in range(len(MHCI_alleles)):
                if i!=0:
                    alleles+=','    
                alleles+=MHCI_alleles[i]
            length = str(epi_len)
            for i in range(len(MHCI_alleles)):
                if i!=0:
                    length+=','+str(epi_len)
            data = 'method=netmhcpan_ba&sequence_text='+seq+'&allele='+alleles+'&length='+length
            command = ['curl', '--data', data, url]
            result = subprocess.run(command, capture_output=True, text=True)
            result = str(result.stdout)
            result = result.split('\n')
            result = [row.split('\t') for row in result if 'HLA' in row]
            length = len(result[0])
            for row in result:
                if len(row)==length:
                    wt = wt_seq[int(row[2])-1:int(row[3])]
                    cancer_peptides.append(peptide(seq=row[5], wt_seq=wt, cancer=cancer_name, mut=missense, allele=row[0], NetMHC_rank=float(row[-1]), NetMHC_binding=float(row[-2])))
            time.sleep(1)
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD8_csv_filename)

def get_NetMHCII(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    lengths = [15,16,17,18]
    cancer_peptides = []
    MHCII_alleles = ['HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*04:05,HLA-DRB1*07:01,HLA-DRB1*08:02,HLA-DRB1*09:01,HLA-DRB1*11:01,HLA-DRB1*12:01,HLA-DRB1*13:02,HLA-DRB1*15:01,HLA-DRB3*01:01,HLA-DRB3*02:02,HLA-DRB4*01:01,HLA-DRB5*01:01','HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*05:01/DQB1*03:01,HLA-DQA1*03:01/DQB1*03:02,HLA-DQA1*04:01/DQB1*04:02,HLA-DQA1*01:01/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DPA1*02:01/DPB1*01:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DPA1*01:03/DPB1*04:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DPA1*02:01/DPB1*05:01,HLA-DPA1*02:01/DPB1*14:01']
    for missense in mutations:
        for l in lengths:
            for alleles in MHCII_alleles:
                url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
                seq, wt_seq = short_form_mutation(og_seq, missense, l)
                length = str(l)
                for i in range(len(alleles.split(','))):
                    if i!=0:
                        length+=','+length
                data = 'method=netmhciipan&sequence_text='+seq+'&allele='+alleles+'&length='+length
                with open('CD4_curl_input.txt', 'w') as file:
                    file.write(data)
                command = ['curl', '--data', '@CD4_curl_input.txt', url]
                result = subprocess.run(command, capture_output=True, text=True)
                result = str(result.stdout)
                result = result.split('\n')
                result = [row.split('\t') for row in result if 'HLA' in row]
                for row in result:
                    wt = wt_seq[int(row[2])-1:int(row[3])]
                    cancer_peptides.append(peptide(seq=row[6], wt_seq=wt, cancer=cancer_name, mut=missense, allele=row[0], NetMHC_rank=float(row[-1]), NetMHC_binding=float(row[-2])))
                time.sleep(1)
    os.remove('CD4_curl_input.txt')
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD4_csv_filename)
    
def get_MHCI_IEDB_immunogenicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = to_peptide_list(CD8_csv_filename)
    text_input = unique_sequences(cancer_name, cancer_peptides, False)
    url = 'http://tools.iedb.org/immunogenicity/'
    my_service = Service(executable_path=chromedriver_path)
    browser = Browser('chrome', service=my_service)
    time.sleep(0.5)
    browser.driver.maximize_window()
    time.sleep(0.5)
    browser.visit(url)
    time.sleep(0.5)
    browser.find_by_id('id_sequence_text').first.type(text_input)
    time.sleep(0.5)
    browser.find_by_name('submit').click()
    time.sleep(10)
    odd = browser.find_by_css('tr.odd')
    time.sleep(10)
    even = browser.find_by_css('tr.even')
    immunogenicity_data = []
    for i in range(len(even)):
        odd_temp = odd[i].value.split(' ')
        immunogenicity_data.append([odd_temp[0], float(odd_temp[2])])
        even_temp = even[i].value.split(' ')
        immunogenicity_data.append([even_temp[0], float(even_temp[2])])
    if len(odd)!=len(even):
        odd_temp = odd[-1].value.split(' ')
        immunogenicity_data.append([odd_temp[0], float(odd_temp[2])])
    browser.quit()
    for i in range(len(immunogenicity_data)):
        for peptide in cancer_peptides:
            if immunogenicity_data[i][0]==peptide.seq:
                peptide.immunogenicity = immunogenicity_data[i][1]
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD8_csv_filename)
    os.remove('Data/'+cancer_name+'/'+cancer_name+' Sequences.txt')

def get_MHCII_IEDB_immunogenicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, length=1000):
    length = length
    cancer_peptides = to_peptide_list(CD4_csv_filename)
    cancer_peptides = [cancer_peptides[i:i+length] for i in range(0, len(cancer_peptides), length)]
    for peptide_list in cancer_peptides:
        text_input = unique_sequences(cancer_name, peptide_list, False)
        url = 'http://tools.iedb.org/cd4episcore/'
        my_service = Service(executable_path=chromedriver_path)
        browser = Browser('chrome', service=my_service)
        time.sleep(0.5)
        browser.driver.maximize_window()
        time.sleep(0.5)
        browser.visit(url)
        time.sleep(0.5)
        browser.find_by_id('id_sequence_text').first.type(text_input)
        time.sleep(0.5)
        browser.find_by_xpath('//*[@id="id_threshold"]/option[10]').click()
        time.sleep(0.5)
        browser.find_by_name('submit').click()
        time.sleep(10)
        odd = browser.find_by_css('tr.odd')
        time.sleep(10)
        even = browser.find_by_css('tr.even')
        immunogenicity_data = []
        for i in range(len(even)):
            odd_temp = odd[i].value.split(' ')
            immunogenicity_data.append([odd_temp[2], float(odd_temp[6])])
            even_temp = even[i].value.split(' ')
            immunogenicity_data.append([even_temp[2], float(even_temp[6])])
        if len(odd)!=len(even):
            odd_temp = odd[-1].value.split(' ')
            immunogenicity_data.append([odd_temp[2], float(odd_temp[6])])
        browser.quit()
        for i in range(len(immunogenicity_data)):
            for peptide in peptide_list:
                if immunogenicity_data[i][0]==peptide.seq:
                    peptide.immunogenicity = immunogenicity_data[i][1]
    cancer_peptides = [element for sublist in cancer_peptides for element in sublist]
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD4_csv_filename)

def get_DeepImmuno(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = to_peptide_list(CD8_csv_filename)    
    df = pd.read_csv(CD8_csv_filename)
    df = df[['Sequence','Allele']]
    df = df.apply(lambda x: x.str.replace(':', ''))
    og_dir = os.getcwd()
    os.chdir('/Users/eric/Documents/Vaccine Programs/DeepImmuno/')
    df.to_csv('write/Input.csv', index=False, header=False)
    command = ['python3','deepimmuno-cnn.py','--mode','multiple','--intdir','write/Input.csv','--outdir',os.getcwd()+'/write']
    result = subprocess.run(command, capture_output=True, text=True)
    os.chdir(og_dir)
    df = pd.read_csv('DeepImmuno/write/deepimmuno-cnn-result.txt',sep='\t')
    for index, row in df.iterrows():
        for peptide in cancer_peptides:
            if peptide.seq==row['peptide'] and peptide.allele.replace(':', '')==row['HLA']:
                peptide.alt_immuno = row['immunogenicity']
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD8_csv_filename)

def get_antigenicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    text_input = unique_sequences(cancer_name, cancer_peptides, True)
    url = 'http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
    my_service = Service(executable_path=chromedriver_path)
    browser = Browser('chrome', service=my_service)
    time.sleep(0.5)
    browser.driver.maximize_window()
    browser.visit(url)
    time.sleep(0.5)
    file_input = browser.find_by_name('uploaded_file')
    time.sleep(0.5)
    file_path = 'Data/'+cancer_name+'/'+cancer_name+' Sequences.txt'
    absolute_path = os.path.abspath(file_path)
    file_input.type(absolute_path)
    time.sleep(0.5)
    browser.find_by_value('tumour').click()
    time.sleep(0.5)
    browser.find_by_value('Submit').click()
    time.sleep(10)
    seq = browser.find_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td')
    results = str(seq.text)
    browser.quit()
    results = results.replace('Model selected: tumour\nThreshold for this model: 0.5\n\n\n\n','')
    results = results.replace('Your Sequence:\n\n>\n\n','')
    results = results.replace('\n\n\nOverall Prediction for the Protective Antigen = ',' ')
    results = results.replace(' ( Probable ANTIGEN ).\n\n\n',' ')
    results = results.replace(' ( Probable NON-ANTIGEN ).\n\n\n',' ')
    results = results.replace(' ( Probable ANTIGEN ).','')
    results = results.replace(' ( Probable NON-ANTIGEN ).','')
    results = results.split(' ')
    antigenicity_list = []
    for i in range(0, len(results),2):
        antigenicity_list.append([results[i], float(results[i+1])])
    for i in range(len(antigenicity_list)):
        for peptide in cancer_peptides:
            if antigenicity_list[i][0]==peptide.seq:
                peptide.antigenicity = antigenicity_list[i][1]
    cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.antigenicity!=None]
    df = to_dataframe(cancer_peptides)
    if MHC==1:
        df.to_csv(CD8_csv_filename)
    if MHC==2:
        df.to_csv(CD4_csv_filename)
    os.remove('Data/'+cancer_name+'/'+cancer_name+' Sequences.txt')
    

def get_toxicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    length = 2500
    cancer_peptides = [cancer_peptides[i:i+length] for i in range(0, len(cancer_peptides), length)]
    for peptide_list in cancer_peptides:
        text_input = unique_sequences(cancer_name, peptide_list, True)
        url = 'https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php'
        my_service = Service(executable_path=chromedriver_path)
        browser = Browser('chrome', service=my_service)
        time.sleep(0.5)
        browser.driver.maximize_window()
        time.sleep(0.5)
        browser.visit(url)
        time.sleep(0.5)
        browser.find_by_name('seq').first.type(text_input)
        time.sleep(0.5)
        browser.find_by_value('Run Analysis!').click()
        element_xpath = '//*[@id="pagerTwo"]/td/select/option[8]'
        wait_time = 60000
        element = WebDriverWait(browser.driver, wait_time).until(EC.presence_of_element_located((By.XPATH, element_xpath)))
        browser.find_by_xpath('//*[@id="pagerTwo"]/td/select/option[8]').click()
        time.sleep(0.5)
        results = browser.find_by_xpath('//*[@id="tableTwo"]/tbody')
        time.sleep(0.5)
        result = results.text.split('\n')
        browser.quit()
        result[:] = [row.split(' ') for row in result]
        toxicity_data = []
        for row in result:
            toxicity_data.append([row[0],row[2]])
        for datum in toxicity_data:
            for peptide in peptide_list:
                if peptide.seq==datum[0]:
                    peptide.toxicity = datum[1]
    cancer_peptides = [element for sublist in cancer_peptides for element in sublist]
    df = to_dataframe(cancer_peptides)
    if MHC==1:
        df.to_csv(CD8_csv_filename)
    if MHC==2:
        df.to_csv(CD4_csv_filename)

def get_prot_parameters(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    unique_seq = []
    for peptide in cancer_peptides:
        if peptide.seq not in unique_seq:
            unique_seq.append(peptide.seq)
    url = 'https://web.expasy.org/protparam/'
    my_service = Service(executable_path=chromedriver_path)
    browser = Browser('chrome', service=my_service)
    time.sleep(0.5)
    browser.driver.maximize_window()
    time.sleep(0.5)
    ProtParam_output = []
    def get_ProtParam(string, end, results):
        start = results.find(string)+len(string)
        if end!=None:
            return results[start:results.find(end,start+1)]
        if end==None:
            return results[start:]
    for seq in unique_seq:
        browser.visit(url)
        time.sleep(0.5)
        browser.find_by_name('sequence').first.type(seq)
        time.sleep(0.5)
        browser.find_by_value('Compute parameters').click()
        time.sleep(0.5)
        parameters = []
        parameters.append(seq)
        results = browser.find_by_xpath('//*[@id="sib_body"]/pre[2]')
        results = str(results.text)
        hl = get_ProtParam('half-life is: ', ' ', results)
        if '>' in hl:
            hl = 20
        parameters.append(float(hl))
        pI = get_ProtParam('pI: ', '\n', results)
        parameters.append(float(pI))
        II = get_ProtParam('(II) is computed to be ', '\n', results)
        parameters.append(float(II))
        AI = get_ProtParam('Aliphatic index: ', '\n', results)
        parameters.append(float(AI))
        GS = get_ProtParam('(GRAVY): ', None, results)
        parameters.append(float(GS))
        ProtParam_output.append(parameters)
    browser.quit()
    for output in ProtParam_output:
        for peptide in cancer_peptides:
            if peptide.seq==output[0]:
                peptide.halflife = output[1]
                peptide.isoelectric_point = output[2]
                peptide.instability_index = output[3]
                peptide.aliphatic_index = output[4]
                peptide.GRAVY_score = output[5]
    df = to_dataframe(cancer_peptides)
    if MHC==1:
        df.to_csv(CD8_csv_filename)
    if MHC==2:
        df.to_csv(CD4_csv_filename)

def get_hydropathicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    kyte_doolittle_scale = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}
    for peptide in cancer_peptides:
        hydropathicity = 0
        for i,aa in enumerate(peptide.seq):
            if 2<i<len(peptide.seq)-1:
                hydropathicity+=kyte_doolittle_scale[aa]
        peptide.hydropathicity = hydropathicity
    df = to_dataframe(cancer_peptides)
    if MHC==1:
        df.to_csv(CD8_csv_filename)
    if MHC==2:
        df.to_csv(CD4_csv_filename)

def get_allergenicity(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    unique_seq = []
    for peptide in cancer_peptides:
        if peptide.seq not in unique_seq:
            unique_seq.append(peptide.seq)
    url = 'https://www.ddg-pharmfac.net/AllerTOP/'
    my_service = Service(executable_path=chromedriver_path)
    browser = Browser('chrome', service=my_service)
    time.sleep(0.5)
    browser.driver.maximize_window()
    time.sleep(0.5)
    for seq in unique_seq:
        browser.visit(url)
        time.sleep(0.5)
        browser.find_by_name('sequence').first.type(seq)
        time.sleep(0.5)
        browser.find_by_name('Submit').click()
        time.sleep(0.5)
        allergenicity = browser.find_by_xpath('//*[@id="box"]/h4[2]') 
        for peptide in cancer_peptides:
            if peptide.seq==seq:
                peptide.allergenicity = str(allergenicity.text)
    browser.quit()
    df = to_dataframe(cancer_peptides)
    if MHC==1:
        df.to_csv(CD8_csv_filename)
    if MHC==2:
        df.to_csv(CD4_csv_filename)

def get_NetMHC_stab(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = to_peptide_list(CD8_csv_filename)
    MHCI_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01']
    for epi_len in [9,10]:
        for missense in mutations:
            url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
            seq, wt_seq = short_form_mutation(og_seq, missense, epi_len)
            alleles = ''
            for i in range(len(MHCI_alleles)):
                if i!=0:
                    alleles+=','    
                alleles+=MHCI_alleles[i]
            length = str(epi_len)
            for i in range(len(MHCI_alleles)):
                if i!=0:
                    length+=','+str(epi_len)
            data = 'method=netmhcstabpan&sequence_text='+seq+'&allele='+alleles+'&length='+length
            command = ['curl', '--data', data, url]
            result = subprocess.run(command, capture_output=True, text=True)
            result = str(result.stdout)
            result = result.split('\n')
            result = [row.split('\t') for row in result if 'HLA' in row]
            length = len(result[0])
            for row in result:
                if len(row)==length:
                    seq=row[5]
                    allele=row[0]
                    for peptide in cancer_peptides:
                        if peptide.seq==seq and peptide.allele==allele:
                            peptide.NetMHC_stab = float(row[-2])
            time.sleep(1)
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD8_csv_filename)


def get_INFgamma(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = to_peptide_list(CD4_csv_filename)
    length = 10000
    cancer_peptides = [cancer_peptides[i:i+length] for i in range(0, len(cancer_peptides), length)]
    for peptide_list in cancer_peptides:
        text_input = unique_sequences(cancer_name, peptide_list, True)
        url = 'https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php'
        my_service = Service(executable_path=chromedriver_path)
        browser = Browser('chrome', service=my_service)
        time.sleep(0.5)
        browser.driver.maximize_window()
        browser.visit(url)
        time.sleep(0.5)
        browser.find_by_name('sequence').first.type(text_input)
        time.sleep(0.5)
        browser.find_by_value('Submit Peptides for Prediction').click()
        element_xpath = '//*[@id="example"]/tbody'
        wait_time = 60000
        element = WebDriverWait(browser.driver, wait_time).until(EC.presence_of_element_located((By.XPATH, element_xpath)))
        browser.find_by_xpath('//*[@id="example_length"]/label/select/option[4]').click()
        unique_seq = []
        for peptide in peptide_list:
            if peptide.seq not in unique_seq:
                unique_seq.append(peptide.seq)
        num_pages = math.ceil(len(unique_seq)/100)
        for i in range(num_pages):
            results = browser.find_by_xpath('//*[@id="example"]/tbody')
            results = results.text.split('\n')
            results[:] = [row.split(' ') for row in results]
            for row in results:
                for peptide in peptide_list:
                    if row[1]==peptide.seq:
                        peptide.inf_gamma = float(row[4])
            browser.find_by_xpath('//*[@id="example_next"]').click()
        browser.quit()
    cancer_peptides = [element for sublist in cancer_peptides for element in sublist]
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD4_csv_filename)
    

def get_MHCFlurry(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename):
    cancer_peptides = to_peptide_list(CD8_csv_filename)
    alleles_dict = {}
    for peptide in cancer_peptides:
        if peptide.seq in alleles_dict.keys():
            alleles_dict[peptide.seq].append(peptide.allele)
        if peptide.seq not in alleles_dict.keys():
            alleles_dict[peptide.seq] = [peptide.allele]
    alleles_wt_dict = {}
    for peptide in cancer_peptides:
        if peptide.wt_seq in alleles_wt_dict.keys():
            alleles_wt_dict[peptide.wt_seq].append(peptide.allele)
        if peptide.wt_seq not in alleles_wt_dict.keys():
            alleles_wt_dict[peptide.wt_seq] = [peptide.allele]
    predictor = Class1AffinityPredictor.load()
    for seq, alleles in alleles_dict.items():
        for allele in alleles:
            data = predictor.predict_to_dataframe(allele=allele, peptides=[seq])
            for peptide in cancer_peptides:
                if data['allele'][0]==peptide.allele and data['peptide'][0]==peptide.seq:
                    peptide.MHCFlurry_rank = data['prediction_percentile'][0]
                    peptide.MHCFlurry_binding = data['prediction'][0]
    for seq, alleles in alleles_wt_dict.items():
        for allele in alleles:
            data = predictor.predict_to_dataframe(allele=allele, peptides=[seq])
            for peptide in cancer_peptides:
                if data['allele'][0]==peptide.allele and data['peptide'][0]==peptide.wt_seq:
                    peptide.MHCFlurry_wt_mut_rank = data['prediction_percentile'][0]/peptide.MHCFlurry_rank
    df = to_dataframe(cancer_peptides)
    df.to_csv(CD8_csv_filename)

def get_population_coverage(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    df = df.sort_values('Sequence')
    df = df.reset_index(drop=True)
    df = df.drop(columns=['Unnamed: 0'], errors='ignore')
    df.to_csv(csv_filename)
    population_coverage = []
    unique_sequences = []
    for i in range(len(df['Allele'])-1):
        if df['Sequence'][i] not in unique_sequences:
            sequence_input = [df['Sequence'][i],df['Allele'][i]]
            unique_sequences.append(df['Sequence'][i])
            while(df['Sequence'][i]==df['Sequence'][i+1] and i<len(df['Allele'])-2):
                sequence_input[1]+=','+df['Allele'][i+1]
                i+=1
            population_coverage.append(sequence_input)
    url = 'http://tools.iedb.org/population/'
    my_service = Service(executable_path=chromedriver_path)
    browser = Browser('chrome', service=my_service)
    time.sleep(0.5)
    browser.driver.maximize_window()
    time.sleep(0.5)
    browser.visit(url)
    time.sleep(0.5)
    browser.find_by_name('number').first.clear()
    time.sleep(0.5)
    browser.find_by_name('number').first.type(len(population_coverage))
    time.sleep(0.5)
    browser.find_by_xpath('//*[@id="input-form"]/table[1]/tbody/tr[1]/td/button').click()
    time.sleep(0.5)
    browser.find_by_xpath('//*[@id="id_population_list"]/option[1]').click()
    time.sleep(0.5)
    browser.find_by_xpath('//*[@id="input-form"]/table[2]/tbody/tr[2]/td[2]/input').click()
    time.sleep(0.5)
    epitope_pre = 'id_epitope_'
    allele_pre = 'id_allele_'
    for i in range(len(population_coverage)):
        browser.find_by_id(epitope_pre+str(i)).type(population_coverage[i][0])
        time.sleep(0.5)
        browser.find_by_id(allele_pre+str(i)).type(population_coverage[i][1])
        time.sleep(0.5)
    browser.find_by_name('submit').click()
    time.sleep(5)
    cumulative_hist = browser.find_by_xpath('//*[@id="content"]/table[2]/tbody/tr[3]/td/img')
    time.sleep(2)
    browser.execute_script("window.scrollTo(0, document.body.scrollHeight);")
    time.sleep(5)
    img_path = cancer_name+'/'+cancer_name+' Population Coverage.png'
    absolute_path = os.path.abspath(img_path)
    cumulative_hist.screenshot(absolute_path)
    time.sleep(10)
    coverage = browser.find_by_xpath('//*[@id="content"]/table[2]/tbody/tr[3]/td/img')
    percent_coverage = browser.find_by_xpath('//*[@id="content"]/table[2]/tbody/tr[2]/td[2]')
    print('Percent Coverage: '+str(percent_coverage.text))
    average_hit = browser.find_by_xpath('//*[@id="content"]/table[2]/tbody/tr[2]/td[3]')
    print('Average Hit: '+str(average_hit.text))
    pc90 = browser.find_by_xpath('//*[@id="content"]/table[2]/tbody/tr[2]/td[4]')
    print('PC90: '+str(pc90.text))
    browser.find_by_name('submit').click()
    new_window_handle = browser.driver.window_handles[-1]
    browser.driver.switch_to.window(new_window_handle)
    hist_data = browser.find_by_xpath('//*[@id="content"]/table/tbody')
    hist_data = hist_data.text.split('\n')
    hist_data[:] = [row.split(' ') for row in hist_data]
    df = pd.DataFrame(hist_data, columns = ['Number of Epitopes','Percent','Cumutlative Percent'])
    if filtered==False:
        df.to_csv(cancer_name+'/'+cancer_name+' Histogram Data.csv')
    if filtered==True:
        df.to_csv(cancer_name+'/'+cancer_name+' Filtered Histogram Data.csv')
    browser.quit()

def allele_text_splitter(unique_alleles):
    allele_text = []
    if len(unique_alleles)>20:
        allele_string = unique_alleles[0]
        for allele in unique_alleles[1:20]:
            allele_string+=','+allele
        allele_text.append(allele_string)
        allele_string = unique_alleles[20]
        for allele in unique_alleles[21:]:
            allele_string+=','+allele
        allele_text.append(allele_string)
    if len(unique_alleles)<=20:
        if len(unique_alleles)==1:
            allele_string = unique_alleles[0]
            allele_text.append(allele_string)
        else:
            allele_string = unique_alleles[0]
            for allele in unique_alleles[1:]:
                allele_string+=','+allele
            allele_text.append(allele_string)
    return allele_text

def parameter_filter(cancer_name, chromedriver_path, og_seq, mutations, CD8_csv_filename, CD4_csv_filename, MHC=1, file_path=None, immunogenicity=False, antigenicity=False, halflife=False, instability_index=False, toxicity=False, allergenicity=False, inf_gamma=False, binding_affinity=False, id=False):
    if MHC==1:
        cancer_peptides = to_peptide_list(CD8_csv_filename)
    if MHC==2:
        cancer_peptides = to_peptide_list(CD4_csv_filename)
    if immunogenicity==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.immunogenicity>0.5]
    if antigenicity==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.antigenicity>0.4]
    if halflife==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.halflife>=1]
    if instability_index==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.instability_index<40]
    if toxicity==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.toxicity=='Non-Toxin']
    if allergenicity==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.allergenicity!='PROBABLE ALLERGEN']
    if inf_gamma==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.inf_gamma>0]
    if binding_affinity==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.NetMHC_binding<500]
    if id==True:
        cancer_peptides[:] = [peptide for peptide in cancer_peptides if peptide.id==1]
    df = to_dataframe(cancer_peptides)
    df.to_csv(file_path)

def rank_peptides(CD8_csv_filename, regression_model_file, classification_model_file):
    df = pd.read_csv(CD8_csv_filename)
    with open(regression_model_file, 'rb') as file:
        reg_model = pickle.load(file)
    columns = reg_model.feature_names_in_
    y_reg = reg_model.predict(df[columns])
    df['Score'] = y_reg
    with open(classification_model_file, 'rb') as file:
        class_model = pickle.load(file)
    columns = class_model.feature_names_in_
    y_class = class_model.predict(df[columns])
    df['ID'] = y_class
    df = df.sort_values(['ID','Score'],ascending=[False,False])
    df = df.reset_index(drop=True)
    df = df.drop(columns=['Unnamed: 0'], errors='ignore')
    df.to_csv(CD8_csv_filename)