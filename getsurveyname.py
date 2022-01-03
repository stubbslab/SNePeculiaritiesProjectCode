
def fromsurveyid(sid):
    lines = open('survey.def').readlines()
    for l in lines:
        try:
            if int(l.split()[2]) == int(sid):
                return l.split()[1]
        except:
            pass
    print('could not find surveyid in survey.def')
    raise
