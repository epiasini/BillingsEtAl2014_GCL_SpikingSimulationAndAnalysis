import numpy as np
from matplotlib import pyplot as plt
import networkx
from networkx.algorithms import bipartite

structures_dir = '/home/ucbtepi/code/network/data'
tissue_c = []
random_c = []
gd_range = range(1,21)
for gd in gd_range:
    print ' '
    print gd
    tissue_model = networkx.read_graphml(structures_dir+'/gd{0}/cr0/gd{0}_cr0.graphml'.format(gd))
    random_bipartite_model = networkx.read_graphml(structures_dir+'/gd{0}/cr1/gd{0}_cr1.graphml'.format(gd))
    for k,g in enumerate([tissue_model, random_bipartite_model]):
        mf_nodes = list(set(n for n,d in g.nodes(data=True) if d['bipartite']==0))
        grc_nodes = list(set(n for n,d in g.nodes(data=True) if d['bipartite']==1))
        print(['tissue', 'random'][k])
        c = bipartite.robins_alexander_clustering(g)
        print("robins-alexander clustering: {}".format(c))
        if k==0:
            tissue_c.append(c)
        else:
            random_c.append(c)

fig, ax = plt.subplots()
ax.plot(gd_range, tissue_c, linewidth=1.5, color='k', label='tissue model')
ax.plot(gd_range, random_c, linewidth=1.5, color='r', label='random bipartite graph')
ax.set_xlabel('dendrites')
ax.set_ylabel('Robins-Alexander clustering')
ax.legend(loc='best')
plt.show()
