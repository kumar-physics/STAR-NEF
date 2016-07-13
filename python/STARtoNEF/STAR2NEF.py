'''
Created on Jun 8, 2016

@author: kumaran
'''
import sys
from Tix import ROW
sys.path.append('/home/kumaran/git/PyNMRSTAR')
import bmrb,csv
import json
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
    
    def read_nef_file(self,nef_file):
        self.starData=bmrb.Entry.from_file(nef_file)
        for saveframe in self.starData:
            #print "saveframe",saveframe.name,self.star_sf_category(saveframe.name)
            for t in saveframe.tags:
                sf_tag="%s.%s"%(saveframe.tag_prefix,t[0])
                print "saveframeTag",t[0],self.star_tag(sf_tag)
            for loop in saveframe:
                for coln in loop.columns:
                    loop_tag="%s.%s"%(loop.category,coln)
                    print "loop_tag",loop_tag,self.star_tag(loop_tag)
    
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
    p=StarToNef()
    #p.read_nef_file('/home/kumaran/nef/CCPN_H1GI.nef')
    p.read_nef_file('/home/kumaran/nef/CCPN_2l9r_Paris_155.nef')
    #p.read_nef_file('/home/kumaran/nef/CCPN_2lci_Piscataway_179.nef')