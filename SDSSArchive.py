#Define an archive specific to the Pan Stars-1 medium field Sn surveys.
# This can just help me to expedite analysis that want to focus just on those fields.
# There are 10 fields, which are (arbtrarily) numbers 0-9.  They can be specifically requested.

from AstronomicalParameterArchive import AstronomicalParameterArchive
from SNDataArchive import SNDataArchive
from RawSNDataStorer import RawSNDataStorer

class SDSSArchive:

    def __init__(self):
        self.ps_survey = 'SDSS'

        #These are the on-sky regions that contain a given set of observations.
        #  They are not actual regions of observation; just regions in which all observations
        #  from a given field are contained.  
        self.fields = {0:[0.0,7.0,-2.0,2.0],   1:[7.0,14.0,-2.0,2.0],
                       2:[14.0,21.0,-2.0,2.0], 3:[21.0,28.0,-2.0,2.0],
                       4:[28.0,35.0,-2.0,2.0], 5:[35.0,42.0,-2.0,2.0],
                       6:[42.0,49.0,-2.0,2.0], 7:[49.0,56.0,-2.0,2.0],
                       8:[300.0,307.0,-2.0,2.0], 9:[307.0,314.0,-2.0,2.0],
                       10:[314.0,321.0,-2.0,2.0], 11:[321.0,328.0,-2.0,2.0],
                       12:[335.0,342.0,-2.0,2.0], 13:[342.0,349.0,-2.0,2.0],
                       14:[349.0,360.0,-2.0,2.0]}
        
        self.fields = {0:[310,417, -2.0, 2.0]}
        
