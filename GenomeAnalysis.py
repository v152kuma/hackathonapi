from Bio import SeqIO
import pandas as pd
import OricFinder
from collections import Counter


def returnGenome(genomeId):
    covid_record = SeqIO.read(genomeId+".fasta", "fasta")
    covid_genome = covid_record.seq
    covid_genome_dict={}
    covid_genome_dict["genome"]=str(covid_genome)
    return covid_genome_dict

#Create one method for Oric

def returnOric(genomeId):
    covid_record = SeqIO.read(genomeId+".fasta", "fasta")
    covid_genome = covid_record.seq
    oric_dict = OricFinder.findOric(covid_genome, 5)
    oric_dict["genomeId"]=genomeId
    return oric_dict

#Create one method for transcription
def transcription(genomeId):
    covid_record = SeqIO.read(genomeId+".fasta", "fasta")
    covid_genome = covid_record.seq
    covid_mrna = covid_genome.transcribe()
    mrna_dict={"mrna":str(covid_mrna)}
    return mrna_dict

#Create one method for translation

def translation(genomeId):
    covid_record = SeqIO.read(genomeId + ".fasta", "fasta")
    covid_genome = covid_record.seq
    covid_mrna = covid_genome.transcribe()
    protienseq=covid_mrna.translate()
    protien_dict={"protien_seq":str(protienseq)}
    return protien_dict

#Create one method for amino acid
def amino_acids(genomeId):
    protien_seq=translation(genomeId)
    covid_aa=protien_seq.get("protien_seq").split("*")
    covid_aa_list = [str(i) for i in covid_aa]
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    amino_acid_df = pd.DataFrame({'polypeptides': covid_aa_list})
    return amino_acid_df

def top_ten_largest_amino_acids(genomeId):
    protien_seq=translation(genomeId)
    covid_aa=protien_seq.get("protien_seq").split("*")
    covid_aa_list = [str(i) for i in covid_aa]
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    amino_acid_df = pd.DataFrame({'polypeptides': covid_aa_list})
    amino_acid_df['count'] = amino_acid_df['polypeptides'].str.len()
    return amino_acid_df.nlargest(10,"count")

def top_amino_acids(genomeId):
    covid_record = SeqIO.read(genomeId + ".fasta", "fasta")
    covid_genome = covid_record.seq
    covid_mrna = covid_genome.transcribe()
    protienseq = covid_mrna.translate()
    aa_list = Counter(protienseq).most_common(10)
    aa_dict={}
    for i in range(len(aa_list)):
        if(aa_list[i][0]=='*'):
            continue
        aa_dict[aa_list[i][0]]=aa_list[i][1]
    return aa_dict


#print(amino_acid_df.nlargest(10,"count"))


#print(Counter(covid_protien).most_common(10))
