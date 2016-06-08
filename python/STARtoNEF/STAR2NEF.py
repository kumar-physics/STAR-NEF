'''
Created on Jun 8, 2016

@author: kumaran
'''
import sys
sys.path.append('/home/kumaran/git/PyNMRSTAR')
import bmrb

class StarToNef(object):
    '''
    classdocs
    '''


    def __init__(self, entryid):
        '''
        Constructor
        '''
        self.starData=bmrb.entry.fromDatabase(entryid)
        print self.starData.getTags()
        
if __name__=="__main__":
    p=StarToNef('15060')