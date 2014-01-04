
rhoinp = '''\
crystal {basename}-{np_margin}-{margin}-rho.cube     
load {basename}-{np_margin}-{margin}-rho.cube     

yt
'''
mepinp = '''\
crystal {basename}-{np_margin}-{margin}-mep.cube     
load {basename}-{np_margin}-{margin}-mep.cube     
load {basename}-{np_margin}-{margin}-rho.cube     

yt
'''
