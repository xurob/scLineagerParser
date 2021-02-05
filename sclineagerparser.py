import pandas as pd
import os

			
			
Apath = "/Users/robin/Desktop/mycnsclin/A1_rmdup_MYCN.rf.A.txt"
Cpath = "/Users/robin/Desktop/mycnsclin/A1_rmdup_MYCN.rf.C.txt"
Gpath = "/Users/robin/Desktop/mycnsclin/A1_rmdup_MYCN.rf.G.txt"
Tpath = "/Users/robin/Desktop/mycnsclin/A1_rmdup_MYCN.rf.T.txt"
refpath = "/Users/robin/Desktop/mycnsclin/ref.txt"
covPath = "/Users/robin/Desktop/mycnsclin/A1_rmdup_MYCN.rf.coverage.txt"
REF = pd.read_table(refpath, sep= "\t", header = None)
A = pd.read_table(Apath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
T = pd.read_table(Tpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
C = pd.read_table(Cpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])
G = pd.read_table(Gpath, sep= ",", header = None, names=["pos", "name", "fw", "rv"])


path = "/Users/robin/Desktop/mycnsclin/"
filelist = []
f = os.listdir(path)
for element in f:
	t = {"A.txt","G.txt","T.txt","C.txt","coverage.txt"}
	g = element.split(t)
	if not g[0] in filelist:
		filelist.append(g[0])
		
for el in filelist:
	Apath = 
	Cpath =
	Gpath =
	Tpath =
	covPath = 



with open("/Users/robin/Desktop/mycnsclin/Z.txt","w") as V:
			V.write("Chr"+"\t" + "Start" +"\t" +"End" +"\t" + "Ref" +"\t" + "Alt" +"\t" + "Normal_ref" +"\t" +"Normal_alt" +"\t" +"Tumor_ref" +"\t" + "Tumor_alt"+"\n")
			for index, row in REF.iterrows():
				start = int(row[0])
				end = int(row[0])
				ref_pos = int(row[0])
				ref_base = row[1]
				if(ref_base == "A" and ref_pos in A['pos'].values):
					normal_ref = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
					tumor_ref = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
					if ref_pos in T['pos'].values:
						normal_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						tumor_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "T" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in C['pos'].values:
						normal_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
						tumor_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "C" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in G['pos'].values:
						normal_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
						tumor_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "G" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
			
				if(ref_base == "T" and ref_pos in T['pos'].values):
						normal_ref = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						tumor_ref = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						if ref_pos in A['pos'].values:
							normal_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
							tumor_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
							V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "A" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
						if ref_pos in C['pos'].values:
							normal_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
							tumor_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
							V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "C" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
						if ref_pos in G['pos'].values:
							normal_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
							tumor_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
							V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "G" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
				if(ref_base == "C" and ref_pos in C['pos'].values):
					normal_ref = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
					tumor_ref = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
					if ref_pos in T['pos'].values:
						normal_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						tumor_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "T" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in A['pos'].values:
						normal_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
						tumor_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "A" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in G['pos'].values:
						normal_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
						tumor_alt = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "G" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
			
				if(ref_base == "G" and ref_pos in G['pos'].values):
					normal_ref = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
					tumor_ref = int(G[G['pos']==ref_pos]['fw']) + int(G[G['pos']==ref_pos]['rv'])
					if ref_pos in T['pos'].values:
						normal_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						tumor_alt = int(T[T['pos']==ref_pos]['fw']) + int(T[T['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "T" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in C['pos'].values:
						normal_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
						tumor_alt = int(C[C['pos']==ref_pos]['fw']) + int(C[C['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "C" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )
					if ref_pos in A['pos'].values:
						normal_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
						tumor_alt = int(A[A['pos']==ref_pos]['fw']) + int(A[A['pos']==ref_pos]['rv'])
						V.write("2" + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref_base) + "\t" + "A" + "\t" + str(normal_ref) + "\t" + str(normal_alt) + "\t" + str(tumor_ref) + "\t" + str(tumor_alt) + "\n" )

			V.close()


	
with open("/Users/robin/Desktop/mycnsclin/ZZ.txt", 'w') as V1:
	name = row[1]
	cov = pd.read_table(covPath,sep= ",", header = None)
	V1.write("Pos" + "\t" + "Count" + "\n")
	for index, row in cov.iterrows():
		poscov = int(row[0])
		countcov = int(row[2])
		V1.write(str(poscov) + "\t" + str(countcov) + "\n" )
		
	V1.close()
	