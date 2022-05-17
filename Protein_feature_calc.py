from Bio import SeqIO
import csv
from openpyxl import load_workbook
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
input_file = open ("Sequence.fasta","r")
for record in SeqIO.parse(input_file,"fasta"):
           my_sec= str(record.seq).rstrip('\\')
           analyse= ProteinAnalysis(my_sec)
           mol_weight = analyse.molecular_weight()
           count_amino= analyse.count_amino_acids()
           epsilon_prot = analyse.molar_extinction_coefficient()
           iso_point=analyse.isoelectric_point()
           ist_index=analyse.instability_index()
           aromati=analyse.aromaticity()
           gra_vy=analyse.gravy()
           flex=analyse.flexibility()
writer = pd.ExcelWriter('protein_feature_data.xlsx',engine='openpyxl')
wb= writer.book
df = pd.DataFrame({'Molecular_Weight':[mol_weight],
                    'Amino_Acid_Count':[count_amino],
                    'molar_extinction_coefficient':[epsilon_prot],
                    'isoelectric_point':[iso_point],
                    'instability_index':[ist_index],
                    'aromaticity':[aromati],
                    'Gravy':[gra_vy],
                    'Flexibility':[flex]})
df.to_excel(writer,index=False)
wb.save('protein_feature_data.xlsx')   
print("Success!")                 
