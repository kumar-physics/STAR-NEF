'''
Created on Jun 8, 2016

@author: kumaran
'''
import sys
sys.path.append('/home/kumaran/git/PyNMRSTAR')
import bmrb,csv
import json
from Bio.PDB import *
import ntpath

class StarToNef(object):
    
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        #self.starData=bmrb.entry.fromDatabase(entryid)
        #bmrb.enableNEFDefaults()
        #bmrb.enable_nef_defaults()
        self.nefData=bmrb.Entry.from_scratch('test01')
        self.read_map_file('/home/kumaran/nef2start.csv')
        self.get_standards()
        
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
     
    
    def nef_sf_category(self,star_sf_category):
        try:
            star_sf_cat=self.map[0][self.map[4].index(star_sf_category)]
        except ValueError:
            star_sf_cat="MISSING"
        if star_sf_cat=="":star_sf_cat="MISSING"
        return star_sf_cat
    
    def star_sf_category(self,nef_sf_category):
        try:
            nef_sf_cat=self.map[4][self.map[0].index(nef_sf_category)]
        except ValueError:
            nef_sf_cat="MISSING"
        if nef_sf_cat=="":nef_sf_cat="MISSING"
        return nef_sf_cat
    
    def nef_tag_category(self,star_tag_category):
        try:
            star_tag_cat=self.map[1][self.map[5].index(star_tag_category)]
        except ValueError:
            star_tag_cat="MISSING"
        if star_tag_cat=="":star_tag_cat="MISSING"
        return star_tag_cat
    
    def star_tag_category(self,nef_tag_category):
        try:
            nef_tag_cat=self.map[5][self.map[1].index(nef_tag_category)]
        except ValueError:
            nef_tag_cat="MISSING"
        if nef_tag_cat=="":nef_tag_cat="MISSING"
        return nef_tag_cat
        
    def nef_tag(self,star_tag):
        try:
            star_tag=self.map[2][self.map[6].index(star_tag)]
        except ValueError:
            star_tag="MISSING"
        if star_tag=="":star_tag="MISSING"
        return star_tag
    
    def star_tag(self,nef_tag):
        try:
            nef_tag=self.map[6][self.map[2].index(nef_tag)]
        except ValueError:
            nef_tag="MISSING"
        if nef_tag=="":nef_tag="MISSING"
        return nef_tag
    
    def nmrstar_to_nef(self,star_file):
        self.starData=bmrb.Entry.from_file(star_file)
        self.nefData=bmrb.Entry.from_scratch(self.starData.entry_id)
        for saveframe in self.starData:
            sf=bmrb.Saveframe.from_scratch(saveframe.name)
            for tag in saveframe.tags:
                s_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                n_tag=self.nef_tag(s_tag)
                if n_tag!="MISSING" and n_tag!="":
                    if n_tag.split(".")[1]=="sf_category":
                        nsf_name=self.nef_sf_category(tag[1])
                        #print n_tag,nsf_name
                        sf.add_tag(n_tag,nsf_name)
                    else:
                        sf.add_tag(n_tag,tag[1])
            for loop in saveframe:
                ll=bmrb.Loop.from_scratch()
                missing_col=[]
                for coln in loop.columns:
                    sl_tag="%s.%s"%(loop.category,coln)
                    nl_tag=self.nef_tag(sl_tag)
                    if nl_tag!="MISSING" and nl_tag!="":
                        ll.add_column(nl_tag)
                    else:
                        missing_col.append(loop.columns.index(coln))
                    
                if loop.category=="_Atom_chem_shift":
                    res_pos=loop.columns.index('Auth_comp_ID')
                    atm_pos=loop.columns.index('Auth_atom_ID')
                if loop.category=="_Gen_dist_constraint":
                    res_pos1=loop.columns.index('Auth_comp_ID_1')
                    atm_pos1=loop.columns.index('Auth_atom_ID_1')
                    res_pos2=loop.columns.index('Auth_comp_ID_2')
                    atm_pos2=loop.columns.index('Auth_atom_ID_2')
                
                    
                
                for dat in loop.data:
                    dat2=dat
                    
                    if loop.category=="_Atom_chem_shift":
                        res=dat[res_pos]
                        atm=dat[atm_pos]
                        amb_id=loop.columns.index("Ambiguity_code")
                        
                        if dat2[amb_id]!='1':
                            amb_atm=[i for i in self.bmrb_standard[res] if atm[:-1] in i]
                            if (len(amb_atm)==2 and res!="TYR" and res!="PHE") or  (len(amb_atm)==2 and (res=="TYR" or res=="PHE") and atm[:-1]=="HB"):
                                if amb_atm.index(atm)==0:
                                    dat2[atm_pos]="%sX"%(atm[:-1])
                                elif amb_atm.index(atm)==1:
                                    dat2[atm_pos]="%sY"%(atm[:-1])
                                else:
                                    print "Something wrong1"
                                ll.add_data(dat2[:-1])
                            elif len(amb_atm)==3 or (len(amb_atm)==2 and (res=="TYR" or res=="PHE")):
                                dat2[atm_pos]="%s%%"%(atm[:-1])
                                if dat2[:-1] not in ll.data:
                                    ll.add_data(dat2[:-1])
                            else:
                                print "Something wrong2"
                        else:
                            ll.add_data(dat2[:-1])
                    elif loop.category=="_Gen_dist_constraint":
                        res1=dat[res_pos1]
                        atm1=dat[atm_pos1]
                        res2=dat[res_pos2]
                        atm2=dat[atm_pos2]
                        logic_col=loop.columns.index("Member_logic_code")
                        if len(atm1[:-1])>1:
                            amb_atm1=[i for i in self.bmrb_standard[res1] if atm1[:-1] in i]
                        else:
                            amb_atm1=[]
                        if len(atm2[:-1])>1:
                            amb_atm2=[i for i in self.bmrb_standard[res2] if atm2[:-1] in i]
                        else:
                            amb_atm2=[]
                        if dat2[logic_col]=="OR":
                            if (len(amb_atm1)==2 and res1!="TYR" and res1!="PHE") or  (len(amb_atm1)==2 and (res1=="TYR" or res1=="PHE") and atm1[:-1]=="HB"):
                                if amb_atm1.index(atm1)==0:
                                    dat2[atm_pos1]="%sX"%(atm1[:-1])
                                elif amb_atm1.index(atm1)==1:
                                    dat2[atm_pos1]="%sY"%(atm1[:-1])
                                else:
                                    print "Something wrong1"
                                
                            if (len(amb_atm2)==2 and res2!="TYR" and res2!="PHE") or  (len(amb_atm2)==2 and (res2=="TYR" or res2=="PHE") and atm2[:-1]=="HB"):
                                if amb_atm2.index(atm2)==0:
                                    dat2[atm_pos2]="%sX"%(atm2[:-1])
                                elif amb_atm2.index(atm2)==1:
                                    dat2[atm_pos2]="%sY"%(atm2[:-1])
                                else:
                                    print "Something wrong1"
                                    
                            
                            if len(amb_atm1)==3 or (len(amb_atm1)==2 and (res1=="TYR" or res1=="PHE")):
                                dat2[atm_pos1]="%s%%"%(atm1[:-1])
                                
                            if len(amb_atm2)==3 or (len(amb_atm2)==2 and (res2=="TYR" or res2=="PHE")):
                                dat2[atm_pos2]="%s%%"%(atm2[:-1])
                                
                            if len(amb_atm1)==3 or len(amb_atm2)==3 or len(amb_atm1)==2 or len(amb_atm2)==2:
                                for i in missing_col: del(dat2[i])
                                if dat2[:] not in ll.data:
                                    ll.add_data(dat2[:])
                            else:
                                print "no way"
                        else:
                            for i in missing_col: del(dat2[i])
                            ll.add_data(dat2[:])
                    else:
                        for i in missing_col: del(dat2[i])
                        ll.add_data(dat2[:])
                        #ll.add_data(dat2[:])
                        
                        
                sf.add_loop(ll)
            self.nefData.add_saveframe(sf)
        print self.nefData
                
    
    
    def nef_to_nmrstar(self,nef_file):
        self.nefData=bmrb.Entry.from_file(nef_file)
        self.starData=bmrb.Entry.from_scratch(self.nefData.entry_id)
        for saveframe in self.nefData:
            #ssf_name=self.star_sf_category(saveframe.name)
            sf=bmrb.Saveframe.from_scratch(saveframe.name)
            for tag in saveframe.tags:
                n_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                s_tag=self.star_tag(n_tag)
                if s_tag!="MISSING" and s_tag!="":
                    if s_tag.split(".")[1]=="Sf_category":
                        ssf_name=self.star_sf_category(tag[1])
                        sf.add_tag(s_tag,ssf_name)
                    else:
                        sf.add_tag(s_tag,tag[1])
            for loop in saveframe:
                ll=bmrb.Loop.from_scratch()
                missing_col=[]
                for coln in loop.columns:
                    nl_tag="%s.%s"%(loop.category,coln)
                    sl_tag=self.star_tag(nl_tag)
                    if sl_tag!="MISSING" and sl_tag!="":
                        ll.add_column(sl_tag)
                    else:
                        missing_col.append(loop.columns.index(coln))
                if loop.category=="_nef_chemical_shift":
                    ll.add_column('_Atom_chem_shift.Ambiguity_code')
                        
                if loop.category=="_nef_chemical_shift":
                    res_pos=loop.columns.index('residue_type')
                    atm_pos=loop.columns.index('atom_name')
                if loop.category=="_nef_distance_restraint":
                    ll.add_column('_Gen_dist_constraint.Member_logic_code')
                if loop.category=="_nef_distance_restraint":
                    res_pos1=loop.columns.index('residue_type_1')
                    atm_pos1=loop.columns.index('atom_name_1')
                    res_pos2=loop.columns.index('residue_type_2')
                    atm_pos2=loop.columns.index('atom_name_2')
                for dat in loop.data:
                    dat2=dat
                    if loop.category=="_nef_distance_restraint":
                        dat2.append(".")
                        
                        if "X" in dat[atm_pos2] or "Y" in dat[atm_pos2] or "%" in dat[atm_pos2] or "X" in dat[atm_pos1] or "Y" in dat[atm_pos1] or "%"  in dat[atm_pos1]:
                            if ("X" in dat[atm_pos1] or "Y" in dat[atm_pos1] or "%" in dat[atm_pos1]):
                                atm1=dat2[atm_pos1]
                                atm1=atm1.replace("X","")
                                atm1=atm1.replace("Y","")
                                atm1=atm1.replace("%","")
                                atm_list1=[i for i in self.bmrb_standard[dat[res_pos1]] if atm1 in i]
                            if ("X" in dat[atm_pos2] or "Y" in dat[atm_pos2] or "%" in dat[atm_pos2]):
                                atm2=dat2[atm_pos2]
                                atm2=atm2.replace("X","")
                                atm2=atm2.replace("Y","")
                                atm2=atm2.replace("%","")
                                atm_list2=[i for i in self.bmrb_standard[dat[res_pos2]] if atm2 in i]
                            stereo_atms=["X","Y"]
                            # % -
                            if "%" in dat2[atm_pos1] and "X" not in dat2[atm_pos1] and "Y" not in dat2[atm_pos1]  and "%" not in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                for aa1 in atm_list1:
                                    dat2[atm_pos1]=aa1
                                    dat2[-1]="OR"
                                    ll.add_data(dat2[:])
                            # - %
                            elif "%" not in dat2[atm_pos1] and "X" not in dat2[atm_pos1] and "Y" not in dat2[atm_pos1]  and "%"  in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                for aa2 in atm_list2:
                                    dat2[atm_pos2]=aa2
                                    dat2[-1]="OR"
                                    ll.add_data(dat2[:])
                            # XY - 
                            elif "%" not in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1])  and "%" not in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                for s_atm in stereo_atms:
                                    s_pos1=dat2[atm_pos1].find(s_atm)
                                    if s_pos1>=0: dat2[atm_pos1]=atm_list1[stereo_atms.index(s_atm)]
                                dat2[-1]="OR"
                                ll.add_data(dat2[:])
                            # - XY
                            elif "%" not in dat2[atm_pos1] and "X"  not in dat2[atm_pos1] and "Y"  not in dat2[atm_pos1]  and "%" not in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                for s_atm in stereo_atms:
                                    s_pos2=dat2[atm_pos2].find(s_atm)
                                    if s_pos2>=0: dat2[atm_pos2]=atm_list2[stereo_atms.index(s_atm)]
                                dat2[-1]="OR"
                                ll.add_data(dat2[:])   
                            # % %               
                            elif "%" in dat2[atm_pos1] and "X" not in dat2[atm_pos1] and "Y" not in dat2[atm_pos1]  and "%" in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                for aa1 in atm_list1:
                                    for aa2 in atm_list2:
                                        dat2[atm_pos1]=aa1
                                        dat2[atm_pos2]=aa2
                                        dat2[-1]="OR"
                                        ll.add_data(dat2[:])
                            # XY XY       
                            elif "%" not in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1])  and "%" not in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                for s_atm in stereo_atms:
                                    s_pos1=dat2[atm_pos1].find(s_atm)
                                    if s_pos1>=0: dat2[atm_pos1]=atm_list1[stereo_atms.index(s_atm)]
                                    s_pos2=dat2[atm_pos2].find(s_atm)
                                    if s_pos2>=0: dat2[atm_pos2]=atm_list2[stereo_atms.index(s_atm)]
                                dat2[-1]="OR"
                                ll.add_data(dat2[:])
                            # XY% -
                            elif "%" in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1])  and "%" not in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                if "X" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="1":
                                            dat2[atm_pos1]=aa1
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="2":
                                            dat2[atm_pos1]=aa1
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                            # - XY%
                            elif "%" not in dat2[atm_pos1] and "X"  not in dat2[atm_pos1] and  "Y"  not in dat2[atm_pos1]   and "%"  in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                if "X" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="1":
                                            dat2[atm_pos2]=aa2
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="2":
                                            dat2[atm_pos2]=aa2
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                            # XY XY%
                            elif "%" not in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1]) and "%"  in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                for s_atm in stereo_atms:
                                    s_pos1=dat2[atm_pos1].find(s_atm)
                                    if s_pos1>=0: dat2[atm_pos1]=atm_list1[stereo_atms.index(s_atm)]
                                if "X" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="1":
                                            dat2[atm_pos2]=aa2
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="2":
                                            dat2[atm_pos2]=aa2
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                                    
                            #XY% XY
                            elif "%"  in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1]) and "%" not in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                for s_atm in stereo_atms:
                                    s_pos2=dat2[atm_pos2].find(s_atm)
                                    if s_pos2>=0: dat2[atm_pos2]=atm_list2[stereo_atms.index(s_atm)]
                                if "X" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="1":
                                            dat2[atm_pos1]=aa1
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="2":
                                            dat2[atm_pos1]=aa1
                                            dat2[-1]="OR"
                                            ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                            
                            # % XY
                            elif "%" in dat2[atm_pos1] and "X" not in dat2[atm_pos1] and "Y" not in dat2[atm_pos1] and "%" not in dat2[atm_pos2] and ("X" in dat2[atm_pos2] or "Y" in dat2[atm_pos2]):
                                for s_atm in stereo_atms:
                                    s_pos2=dat2[atm_pos2].find(s_atm)
                                    if s_pos2>=0: dat2[atm_pos2]=atm_list2[stereo_atms.index(s_atm)]
                                for aa1 in atm_list1:
                                    dat2[atm_pos1]=aa1
                                    dat2[-1]="OR"
                                    ll.add_data(dat2[:])
                            # XY %
                            elif "%" not in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1]) and "%"  in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                for s_atm in stereo_atms:
                                    s_pos1=dat2[atm_pos1].find(s_atm)
                                    if s_pos1>=0: dat2[atm_pos1]=atm_list1[stereo_atms.index(s_atm)]
                                for aa2 in atm_list2:
                                    dat2[atm_pos2]=aa2
                                    dat2[-1]="OR"
                                    ll.add_data(dat2[:])
                            #XY% %
                            elif "%" in dat2[atm_pos1] and ("X" in dat2[atm_pos1] or "Y" in dat2[atm_pos1]) and "%" in dat2[atm_pos2] and "X" not in dat2[atm_pos2] and "Y" not in dat2[atm_pos2]:
                                if "X" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="1":
                                            dat2[atm_pos1]=aa1
                                            for aa2 in atm_list2:
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos1]:
                                    for aa1 in atm_list1:
                                        s_pos1=dat2[atm_pos1].find("X")
                                        if aa1[s_pos1]=="2":
                                            dat2[atm_pos1]=aa1
                                            for aa2 in atm_list2:
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                                
                            # % XY%
                            elif "%" in dat2[atm_pos1] and "X" not in dat2[atm_pos1] and "Y" not in dat2[atm_pos1] and "%" in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                if "X" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="1":
                                            dat2[atm_pos2]=aa2
                                            for aa1 in atm_list1:
                                                dat2[atm_pos1]=aa1
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos2]:
                                    for aa2 in atm_list2:
                                        s_pos2=dat2[atm_pos2].find("X")
                                        if aa2[s_pos2]=="2":
                                            dat2[atm_pos2]=aa2
                                            for aa1 in atm_list1:
                                                dat2[atm_pos1]=aa1
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                                
                                    
                                
                            
                            # XY% XY%
                            elif "%" in dat2[atm_pos1] and ("X"  in dat2[atm_pos1] or "Y"  in dat2[atm_pos1])  and "%"  in dat2[atm_pos2] and ("X"  in dat2[atm_pos2] or "Y"  in dat2[atm_pos2]):
                                #print dat2
                                if "X" in dat2[atm_pos1] and "X" in dat2[atm_pos2]:
                                    s_pos1=dat2[atm_pos1].find("X")
                                    s_pos2=dat2[atm_pos2].find("X")
                                    for aa1 in atm_list1:
                                        for aa2 in atm_list2:
                                            if aa1[s_pos1]=="1" and aa2[s_pos2]=="1":
                                                dat2[atm_pos1]=aa1
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos1] and "X" in dat2[atm_pos2]:
                                    s_pos1=dat2[atm_pos1].find("Y")
                                    s_pos2=dat2[atm_pos2].find("X")
                                    for aa1 in atm_list1:
                                        for aa2 in atm_list2:
                                            if aa1[s_pos1]=="2" and aa2[s_pos2]=="1":
                                                dat2[atm_pos1]=aa1
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                elif "X" in dat2[atm_pos1] and "Y" in dat2[atm_pos2]:
                                    s_pos1=dat2[atm_pos1].find("X")
                                    s_pos2=dat2[atm_pos2].find("Y")
                                    for aa1 in atm_list1:
                                        for aa2 in atm_list2:
                                            if aa1[s_pos1]=="1" and aa2[s_pos2]=="2":
                                                dat2[atm_pos1]=aa1
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos1] and "Y" in dat2[atm_pos2]:
                                    s_pos1=dat2[atm_pos1].find("Y")
                                    s_pos2=dat2[atm_pos2].find("Y")
                                    for aa1 in atm_list1:
                                        for aa2 in atm_list2:
                                            if aa1[s_pos1]=="2" and aa2[s_pos2]=="2":
                                                dat2[atm_pos1]=aa1
                                                dat2[atm_pos2]=aa2
                                                dat2[-1]="OR"
                                                ll.add_data(dat2[:])
                                else:
                                    print "Failed",dat2
                                #print dat2
                            else:
                                print "Failed",dat2
                        else:
                            ll.add_data(dat2[:])
                    
                    elif loop.category=="_nef_chemical_shift":
                        dat2.append("1")
                        if "X" in dat[atm_pos] or "Y" in dat[atm_pos] or "%" in dat[atm_pos]:
                            atm=dat2[atm_pos]
                            atm=atm.replace("X","")
                            atm=atm.replace("Y","")
                            atm=atm.replace("%","")
                            #print atm
                            atm_list=[i for i in self.bmrb_standard[dat[res_pos]] if atm in i]
                            if (dat[res_pos]=="PHE" or dat[res_pos]=="TYR") and ("HD" in atm_list[0] or "HE" in atm_list[0]):
                                amb_code='3'
                            else:
                                amb_code='2'
                            if "%" in  dat2[atm_pos] and "X" not in dat2[atm_pos] and "Y" not in dat2[atm_pos]:
                                for aa in atm_list:
                                    dat2[atm_pos]=aa
                                    dat2[-1]=amb_code
                                    ll.add_data(dat2[:])
                                    #print ll
                            elif "%" not in  dat2[atm_pos] and ("X" in dat2[atm_pos] or "Y"  in dat2[atm_pos]):
                                if "X" in dat2[atm_pos]:
                                    dat2[atm_pos]=atm_list[0]
                                elif "Y" in dat2[atm_pos]:
                                    dat2[atm_pos]=atm_list[1]
                                else:
                                    dat2[atm_pos]="XXX"
                                dat2[-1]=amb_code
                                ll.add_data(dat2[:])
                            elif "%"  in  dat2[atm_pos] and ("X" in dat2[atm_pos] or "Y"  in dat2[atm_pos]):
                                if "X" in dat2[atm_pos]:
                                    X_pos=dat2[atm_pos].find("X")
                                    for aa in atm_list:
                                        if aa[X_pos]=="1":
                                            dat2[atm_pos]=aa
                                            dat2[-1]=amb_code
                                            ll.add_data(dat2[:])
                                elif "Y" in dat2[atm_pos]:
                                    Y_pos=dat2[atm_pos].find("Y")
                                    for aa in atm_list:
                                        if aa[Y_pos]=="2":
                                            dat2[atm_pos]=aa
                                            dat2[-1]=amb_code
                                            ll.add_data(dat2[:])
                                else:
                                    dat2[atm_pos]="XXX"
                                    dat2[-1]=amb_code
                                    ll.add_data(dat2[:])
                                            
                                        
                                
                            else:
                                ll.add_data(dat2[:])
                        else:
                            ll.add_data(dat2[:])
                    else:
                        for i in missing_col: del(dat2[i])
                        ll.add_data(dat2[:])
                #print ll
                
                    
                    
                    
                sf.add_loop(ll) 
            self.starData.add_saveframe(sf)
        (file_path,file_name)=ntpath.split(nef_file)
        if file_path=="": file_path="."
        out_file_name="%s.str"%(file_name.split(".nef")[0])
        outfile=file_path+"/"+out_file_name
        with open(outfile,'w') as strfile:
            strfile.write(str(self.starData))

            
    
    def read_nef_file(self,nef_file):
        self.starData=bmrb.Entry.from_file(nef_file)
        print self.starData.bmrb_id
        for saveframe in self.starData:
            print "saveframe",saveframe.name,self.star_sf_category(saveframe.name)
            for t in saveframe.tags:
                sf_tag="%s.%s"%(saveframe.tag_prefix,t[0])
                print "saveframeTag"+","+sf_tag+","+self.star_tag(sf_tag)
            for loop in saveframe:
                for coln in loop.columns:
                    loop_tag="%s.%s"%(loop.category,coln)
                    print "loopTag(Columns)"+","+loop_tag+","+self.star_tag(loop_tag)
            print "================End of Saveframe==============="
    
    def read_star_file(self,star_file):
        self.starData=bmrb.Entry.from_database(star_file)
        for saveframe in self.starData:
            print saveframe.name,self.nef_sf_category(saveframe.name)
            for t in saveframe.tags:
                sf_tag="%s.%s"%(saveframe.tag_prefix,t[0])
                print t[0],self.nef_tag(sf_tag)
                for loop in saveframe:
                    for coln in loop.columns:
                        loop_tag="%s.%s"%(loop.category,coln)
                        print loop_tag,self.nef_tag(loop_tag)
    
    def test(self):
        for saveframe in self.starData:
            print "===================Saveframe==========================="
            try:
                star_sf_cat=self.map[4][self.map[0].index(saveframe.name)]
            except ValueError:
                star_sf_cat="MISSING in STAR"
            try:
                star_sf_tag=self.map[5][self.map[1].index(saveframe.tag_prefix[1:])]
            except ValueError:
                star_sf_tag="MISSING in STAR"
            print saveframe.name,"-----",star_sf_cat
            print saveframe.tag_prefix[1:],"-----",star_sf_tag
            sf=bmrb.Saveframe.from_scratch(star_sf_cat,tag_prefix=star_sf_tag)
            
            print "------------------SaveframeTags------------------------"
            sf_tags=[i[0] for i in saveframe.tags]
            sf_values=[i[1] for i in saveframe.tags]
            print sf_tags
            print sf_values
            sf_dat=[]
            for i in sf_tags:
                nef_tag="%s.%s"%(saveframe.tag_prefix,i)
                try:
                    star_tag=self.map[6][self.map[2].index(nef_tag)]
                except ValueError:
                    star_tag="MISSING.MISSING_in_STAR_%d"%(sf_tags.index(i))
                print nef_tag,"-----",star_tag
                if star_tag=="":star_tag="MISSING.MISSING_in_STAR_%d"%(sf_tags.index(i))
                sf_dat.append((star_tag.split(".")[1],sf_values[sf_tags.index(i)]))
            sf.add_tags(sf_dat)
            print "--------------------Loops-------------------------------"
            for loop in saveframe:
                
                #print loop.category,loop.columns
                print loop.get_tag()
                print loop.data
                ll=bmrb.Loop.from_scratch()
                lt=[]
                cat=""
                for i in loop.columns:
                    loop_tag="%s.%s"%(loop.category,i)
                    try:
                        star_loop_tag=self.map[6][self.map[2].index(loop_tag)]
                        print star_loop_tag
                        cat=star_loop_tag.split(".")[0]
                        print "Cat",cat
                    except ValueError:
                        star_loop_tag="_Missing%s.MISSING_in_STAR_%d"%(loop.category,loop.columns.index(i))
                    if star_loop_tag=="":star_loop_tag="_Missing%s.MISSING_in_STAR_%d"%(loop.category,loop.columns.index(i))
                    print loop_tag,"-----",star_loop_tag
                    lt.append(star_loop_tag.split(".")[0])
                if cat=="":
                    lt.category="Missiong%s"%(loop.category)
                else:
                    lt.category=cat
                ll.add_column(lt)
                for i in loop.data:
                    ll.add_data(i)
                sf.add_loop(ll)
            self.nefData.add_saveframe(sf)
           
    def read_map_file(self,mapfile):
        with open(mapfile,'rb') as csvfile:
            spamreader = csv.reader(csvfile,delimiter=',')
            map_dat=[]
            for r in spamreader:
                if r.count('') != 10:
                    map_dat.append(r)
        self.map=map(list,zip(*map_dat))
        
    
                
        
if __name__=="__main__":
    #p=StarToNef('15060')
    fname=sys.argv[1]
    p=StarToNef()
    p.nmrstar_to_nef(fname)
    #p.nef_to_nmrstar(fname)
    #p.nef_to_nmrstar('/home/kumaran/nef/CCPN_H1GI.nef')
    #p.nef_to_nmrstar('/home/kumaran/nef/CCPN_2l9r_Paris_155.nef')
    #p.nef_to_nmrstar('/home/kumaran/nef/CCPN_2lci_Piscataway_179.nef')
    #p.read_nef_file('/home/kumaran/nef/CCPN_2l9r_Paris_155.nef')
    #p.read_nef_file('/home/kumaran/nef/CCPN_2lci_Piscataway_179.nef')