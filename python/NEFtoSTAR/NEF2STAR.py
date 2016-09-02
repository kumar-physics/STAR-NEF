'''
Created on Aug 30, 2016

@author: kumaran
'''
import sys
from scipy.constants.constants import atm, atmosphere
sys.path.append('/home/kumaran/git/PyNMRSTAR')
import bmrb,csv
import json
from Bio.PDB import *
import ntpath
import re


class NEFtoSTAR(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.IUPAC={'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'H', 'HA', 'HB2', 'HB3', 'HD2'],
                    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE21', 'HE22'],
                    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'],
                    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13'],
                    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3'],
                    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'],
                    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
                    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'H', 'HA', 'HB2', 'HB3', 'HD21', 'HD22'],
                    'GLY': ['N', 'CA', 'C', 'O', 'H', 'HA2', 'HA3'],
                    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2'],
                    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', 'HB2', 'HB3', 'HG', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
                    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE', 'HH11', 'HH12', 'HH21', 'HH22'],
                    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'],
                    'ALA': ['N', 'CA', 'C', 'O', 'CB', 'H1', 'HA', 'HB1', 'HB2', 'HB3', 'H2', 'H3'],
                    'VAL': ['N', 'CA', 'C', "O'", 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', "O''"],
                    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE2'],
                    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
                    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1', 'HE2', 'HE3']}

    
    def get_atm_list(self,res,nefAtom):
        try:
            atms=self.IUPAC[res]
            alist=[]
            try:
                refatm=re.findall(r'(\S+)([XY])([%*])|(\S+)([%*])|(\S+)([XY])',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\d\d+'%(refatm[0]))
                    alist=[i for i in atms if re.search(pattern, i)]
                elif set==[3,4]:
                    if refatm[4]=="%":
                        pattern=re.compile(r'%s\d+'%(refatm[3]))
                    elif refatm[4]=="*":
                        pattern=re.compile(r'%s\S+'%(refatm[3]))
                    else:
                        print "something wrong"
                    alist=[i for i in atms if re.search(pattern, i)]
                elif set==[5,6]:
                    pattern=re.compile(r'%s\d+'%(refatm[5]))
                    alist=[i for i in atms if re.search(pattern, i)]
                else:
                    print "Wrong regular expression"
            except IndexError:
                print nefAtom
            if len(alist)==0:
                if nefAtom in atms:
                    alist.append(nefAtom)
        except KeyError:
            print "Residue not found"
            alist=[]
        return alist
        
        
    def generate_lookup_table(self):
        for res in self.bmrb_standard.keys():
            for atm in self.bmrb_standard[res]:
                for i in range(1,len(atm)+1):
                    s=[j for j in self.bmrb_standard[res] if j[:i]==atm[:i]]
                    if len(s)==1:
                        print res,atm,atm
                    elif len(s)==2:
                        print res,"%sX"%(atm[:i]),",".join(s)
                        print res,"%sY"%(atm[:i]),",".join(s)
                        print res,"%s%%"%(atm[:i]),",".join(s)
                    elif len(s)==3:
                        print res,"%s%%"%(atm[:i]),",".join(s)
                    else:
                        print res,"%s*"%(atm[:i]),",".join(s)
                        
        
    def get_standards(self):
        pdb_parser=PDBParser(QUIET=True)
        aa_names=pdb_parser.get_structure('aa_name','/home/kumaran/nef/aa_normal_20.pdb')
        self.bmrb_standard={}
        for m in aa_names:
            for c in m:
                for r in c:
                    atm=[]
                    for a in r:
                        atm.append(a.name)
                    self.bmrb_standard[r.resname]=atm
        
        
    def read_map_file(self,mapfile):
        with open(mapfile,'rb') as csvfile:
            spamreader = csv.reader(csvfile,delimiter=',')
            map_dat=[]
            for r in spamreader:
                if r.count('') != 10:
                    map_dat.append(r)
        self.map=map(list,zip(*map_dat))


if __name__=="__main__":
    p=NEFtoSTAR()
    for i in range(10):
        res=raw_input("Enter res\n")
        nfatm=raw_input("Enter nef atom\n")
        print p.get_atm_list(res, nfatm)
    #p.generate_lookup_table()