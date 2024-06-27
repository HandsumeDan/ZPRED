from pymol import cmd

def saveAsPdb(arg1,arg2):
	#print arg1
	cmd.load(arg1);
	fNm=arg1.split(".")
	#print fNm
	pdbFile=fNm[0] + ".pdb"
	cmd.save(arg2)

cmd.extend ('saveAsPdb',saveAsPdb)
