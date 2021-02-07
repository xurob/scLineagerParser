import pandas as pd
import os
import argparse
import subprocess
	
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", help="input directory: please add /")
parser.add_argument("-o", "--output", help="output path/to/outputdir/")
parser.add_argument("-ref", "--reference", help="Path/to/reffile")
parser.add_argument("-anno", "--annovar", help="Path/to/annovar/")


args = parser.parse_args()



file_dir = args.input
outpre = args.output
refpath = args.reference
annovar_path = args.annovar

temppath = str(outpre) + "temp/"
temppath_annovar_input = str(temppath) + "annovat_input/"
temppath_annovar_output = str(temppath) + "annovat_output/"
mutations_dir = outpre + "mutations/"

#create directories
if not os.path.exists(temppath):
		   os.makedirs(temppath)
if not os.path.exists(temppath_annovar_input):
		   os.makedirs(temppath_annovar_input)
if not os.path.exists(temppath_annovar_output):
		   os.makedirs(temppath_annovar_output)
if not os.path.exists(mutations_dir):
	   os.makedirs(mutations_dir)


# =============================================================================
# convert mgatk raw output into annovar input files
# =============================================================================
filelist = os.scandir(file_dir)
filenamelist = []
for file in filelist:
	if file.name.endswith(".A.txt") and file.name not in filenamelist:
		filenamelist.append(file.name.split(".A.txt")[0])


for filename in filenamelist:
	#make final dir for every cell; coverage already goes in; rest into temp
	
	cell_dir = mutations_dir + filename + "/"
	if not os.path.exists(cell_dir):
		   os.makedirs(cell_dir)
	
	
	
	#Parser thanks Paulinus
	Apath = file_dir + filename + ".A.txt"
	Cpath = file_dir + filename + ".C.txt"
	Gpath = file_dir + filename + ".G.txt"
	Tpath = file_dir + filename + ".T.txt"
	covPath = file_dir + filename + ".coverage.txt"
	REF = pd.read_table(refpath, sep= "\t", header = None)
	A = pd.read_table(Apath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
	T = pd.read_table(Tpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
	C = pd.read_table(Cpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
	G = pd.read_table(Gpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])

	ATCG = {'A': A, 'T': T, 'C': C, 'G': G}

	with open(temppath_annovar_input + filename + ".mutations.txt", "w") as V:
	    V.write(
	        "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "Ref" + "\t" + "Alt" + "\t" + "Normal_ref" + "\t" + "Normal_alt" + "\t" + "Tumor_ref" + "\t" + "Tumor_alt" + "\n")
	    for index, row in REF.iterrows():
	        start = int(row[0])
	        end = int(row[0])
	        ref_pos = int(row[0])
	        ref_base = row[1]
	        for Bkey in ATCG:
	            if ref_base == Bkey:
	                B = ATCG[Bkey]
	                if ref_pos in B['pos'].values:
	                    Blist = B[B['pos'] == ref_pos]
	                    normal_ref = int(Blist['fw']) + int(Blist['rv'])
	                    '''Warum lässt du den zwei mal die selben Werte auslesen? Der muss ja jedes mal das gesamte Objekt durchsuchen, dauert ewig'''
	                    # tumor_ref = int(B[B['pos'] == ref_pos]['fw']) + int(B[B['pos'] == ref_pos]['rv'])
	                    tumor_ref = normal_ref
	
	                    for Bkey2 in ATCG:
	                        if not Bkey2 == Bkey:
	                            B2 = ATCG[Bkey2]
	                            if ref_pos in B2['pos'].values:
	                                B2list = B2[B2['pos'] == ref_pos]
	                                normal_alt = int(B2list['fw']) + int(B2list['rv'])
	                                '''same here'''
	                                # tumor_alt = int(B2[B2['pos'] == ref_pos]['fw']) + int(B2[B2['pos'] == ref_pos]['rv'])
	                                tumor_alt = normal_alt
	                                V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(
	                                    ref_base) + "\t" + Bkey2 + "\t" + str(normal_ref) + "\t" + str(
	                                    normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n")
	    V.close()
	

	with open(cell_dir +  "coverage.txt", 'w') as V1:
		name = row[1]
		cov = pd.read_table(covPath,sep= ",", header = None)
		V1.write("Chr" + "\t" +"Pos" + "\t" +"Base" + "\t" + "Count" + "\n")
		for index, row in cov.iterrows():
			poscov = int(row[0])
			countcov = int(row[2])
			V1.write("2" + "\t" + str(poscov) +"\t" + "NA" +"\t" + str(countcov) + "\n" )
			
		V1.close()
	


# =============================================================================
# run annovar
# =============================================================================

with os.scandir(temppath_annovar_input)  as dir:
	for file in dir:
		if file.name.endswith(".mutations.txt"):
			
			subprocess.check_call("perl " + annovar_path+"table_annovar.pl " +file.path + " humandb/ -buildver hg19 -out "+ temppath_annovar_output + file.name.split(".txt")[0]+".annovar.txt -protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all -operation g,f,f,f,f,f -nastring .", shell=True)

'''
# =============================================================================
# convert annovar_output into sclineager_input
# =============================================================================


with os.scandir(temppath_annovar_out) as dir:
	for file in dir:
		annovar_out_file = file.path
		temp_file = temppath_annovar_input + file.name.split(".annovar.txt")[0]+".txt"
		
		
		annovar_table = pd.read_table(annovar_out_file, sep= "\t")
		temp_table = pd.read_table(temp_file, sep= "\t")
		
		
path = "/Users/robin/Desktop/SCLINEAGER/mutations/301/P301_2/tre123.txt"

file3 = pd.read_table(path, sep= "\t", header = None)
	

path = "/Users/robin/Desktop/annovar/myanno.hg19_multianno.txt"

file3 = pd.read_table(path, sep= "\t")

Chr	Start	End	Ref	Alt	Caller	Normal_ref	Normal_alt	Tumor_ref	Tumor_alt	Func.refGene	Gene.refGene	ExonicFunc.refGene	AAChange.refGene	SIFT_pred	Polyphen2_HVAR_pred	cosmic70	esp6500siv2_all	ExAC_ALL	X1000g2015aug_all

df1 = df[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'SIFT_pred', 'Polyphen2_HVAR_pred', 'cosmic70', 'esp6500siv2_all', 'ExAC_ALL', 'X1000g2015aug_all']

		 
with open("/Users/robin/Desktop/mycnsclin/Z.txt","w") as V:
			V.write("Chr"+"\t" + "Start" +"\t" +"End" +"\t" + "Ref" +"\t" + "Alt" +"\t" + "Normal_ref" +"\t" +"Normal_alt" +"\t" +"Tumor_ref" +"\t" + "Tumor_alt"+"\n")
			file3.drop(['SIFT_score', 'Polyphen2_HDIV_score'], axis=1)
				
	



columns = ['Chr','Start','End','Ref','Alt','Normal_ref', 'Normal_alt','Tumor_ref','Tumora_alt','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','SIFT_pred','Polyphen2_HVAR_pred','cosmic70','esp6500siv2_all','ExAC_ALL']

file2 = file3[columns]																																




	'''