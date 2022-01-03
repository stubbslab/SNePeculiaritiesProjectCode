def linef(file1,line1):
    co=0
    with open(file1, 'r') as inF:
        for line in inF:
            if line1 in line:
                return co
            co=co+1
            
