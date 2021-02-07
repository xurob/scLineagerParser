import pandas as pd
import os
import argparse
import subprocess
	
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", help="input directory: please add /")
parser.add_argument("-o", "--output", help="output path/to/outputdir/")
parser.add_argument("-ref", "--reference", help="Path/to/reffile")
parser.add_argument("-anno", "--annovar", help="Path/to/annovar/")
parser.add_argument("-cor", "--correction", help="correction for different offset (int)")

args = parser.parse_args()



file_dir = args.input
outpre = args.output
refpath = args.reference
annovar_path = args.annovar
correction = args.correction

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
	                                V.write("2" + "\t" + str(start+correction) + "\t" + str(end+correction) + "\t" + str(
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
			
			subprocess.check_call("perl " +"table_annovar.pl " +file.path + " humandb/ -buildver hg19 -out "+ temppath_annovar_output + file.name.split(".txt")[0]+" -protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all -operation g,f,f,f,f,f -nastring . --nopolish", shell=True)


# =============================================================================
# convert annovar_output into sclineager_input
# =============================================================================

temppath_annovar_output = "/Users/robin/Documents/projects/scLineagerParser/out/temp/annovat_output/"
temppath_annovar_input = "/Users/robin/Documents/projects/scLineagerParser/out/temp/annovat_input/"
mutations_dir = "/Users/robin/Documents/projects/scLineagerParser/out/mutations/"
filelist_anno_out = []
with os.scandir(temppath_annovar_output) as dir:
	for file in dir:
		if file.name.endswith(".hg19_multianno.txt"):
			filelist_anno_out.append(file.name.split(".hg19_multianno.txt")[0])
			
			
	for filename in filelist_anno_out:
		annovar_out_file = temppath_annovar_output + filename + ".hg19_multianno.txt"
		temp_file = temppath_annovar_input + filename + ".txt"
		
		
		annovar_table = pd.read_table(annovar_out_file, header = [0,1], sep= "\t")
		temp_table = pd.read_table(temp_file,header = [0], sep= "\t")
		with open(mutations_dir + filename.split(".mutations")[0] + "/" + filename.split(".mutations")[0] + "_mutations_hg19.txt", 'w') as V3:
			V3.write("Chr" + "\t" + "Start" + "\t" + "End" + "\t"+ "Ref" +
					 "\t" + "Normal_ref" + "\t" + "Normal_alt" + "\t"+ "Tumor_ref" +
					 "\t"+ "Tumor_alt" + "\t" + "Func.refGene" + "\t"+ "Gene.refGene" 
					 + "\t" + "ExonicFunc.refGene" + "\t"+ "AAChange.refGene"
					 + "\t" + "SIFT_pred" + "\t"+ "Polyphen2_HVAR_pred"
					 + "\t" + "cosmic70" + "\t"+ "esp6500siv2_all"
					 + "\t" + "ExAC_ALL" + "\t"+ "1000g2015aug_all"
					 + "\n")
			for index, row in annovar_table.iterrows():
				if index >=0:
					V3.write("2" + "\t" + str(row['Start'][0])+ "\t" + str(row['End'][0])+ "\t" + str(row['Ref'][0])+ "\t" + str(row['Alt'][0])
					  + "\t" + str(temp_table.loc[index,'Normal_ref']) + "\t" + str(temp_table.loc[index,'Normal_alt']) + "\t" 
					  + str(temp_table.loc[index,'Normal_ref']) + "\t" + str(temp_table.loc[index,'Normal_alt']) + "\t" +
					  str(row['Func.refGene'][0]) + "\t" + str(row['Gene.refGene'][0]) + "\t" + str(row['ExonicFunc.refGene'][0]) + "\t" +
					  str(row['AAChange.refGene'][0]) + "\t" + str(row['SIFT_pred'][0]) + "\t" + str(row['Polyphen2_HVAR_pred'][0]) + "\t" + str(row['cosmic70'][0])
					  + "\t" + str(row['esp6500siv2_all'][0]) + "\t" + str(row['ExAC_ALL'][0]) + "\t" + str(row['1000g2015aug_all'][0]) + "\n" 
					   )
		