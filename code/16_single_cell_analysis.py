import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import matplotlib as mpl
#####################################
cell_age=pd.read_csv("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/cell_age_61841.csv",sep="\t",header=0)
############重命名行############
cell_age.rename(index=cell_age['Unnamed: 0'],inplace=True)
########################################
scvis_map = pd.read_csv("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/cell_scvis_61841.csv",sep="\t",header=0)
scvis_map.columns = ['cell','z_coordinate_0','z_coordinate_1','week']
xmin=-17
xmax=17
X, Y = np.mgrid[xmin:xmax:341j, xmin:xmax:341j]
positions = np.vstack([X.ravel(), Y.ravel()])
#path_week='/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/cell_density/ple_30_38_loca.txt'
values = np.vstack([scvis_map['z_coordinate_0'], scvis_map['z_coordinate_1']])
kernel = stats.gaussian_kde(values)
h=kernel.factor


#week_loca=pd.read_csv(path_week,sep="\t",header=0)
week=["W6","W12","W18","W24","W30","W38"]
for w in week:
    print(w)
    scvis_map_m=scvis_map[scvis_map.index.str.contains(w)]
    week_m=scvis_map[scvis_map.index.str.contains(w)]
    #values = np.vstack([week_loca['X0'],week_loca['X1']])
    values = np.vstack([scvis_map_m['z_coordinate_0'], scvis_map_m['z_coordinate_1']])
    #Z = np.reshape(kernel(positions).T, X.shape)
    marker=pd.read_csv("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/gene_density/ple_nonple_mean_47_24_61841_7.1.txt",sep="\t",header=0)
    marker_m=marker[marker.index.str.contains(w)]
    gene1=stats.gaussian_kde(values,weights=marker_m['nonple_mean'])
    gene1.set_bandwidth(bw_method=h)
    Z = np.reshape(gene1(positions).T, X.shape)
    fig, ax = plt.subplots()
    ax.plot(scvis_map['z_coordinate_0'], scvis_map['z_coordinate_1'], 'k.', markersize=0.1,zorder=10)
    ax.imshow(np.rot90(Z), cmap=mpl.cm.OrRd, extent=[xmin, xmax, xmin, xmax],zorder=1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([xmin, xmax])
    nonple_week='/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/gene_density/6_38figure/week'+str(w)+'nonple.png'
    plt.savefig(nonple_week)
    plt.draw()
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.5)
    cmap =  mpl.cm.OrRd
    norm = mpl.colors.Normalize(vmin=Z.min(), vmax=Z.max())
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
              orientation='vertical', label='')
    nonple_colorbar='/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/gene_density/6_38figure/week'+str(w)+'nonple_bar.pdf'
    plt.savefig(nonple_colorbar)
    plt.draw()
    gene1=stats.gaussian_kde(values,weights=marker_m['ple_mean'])
    gene1.set_bandwidth(bw_method=h)
    Z = np.reshape(gene1(positions).T, X.shape)
    fig, ax = plt.subplots()
    ax.plot(scvis_map['z_coordinate_0'], scvis_map['z_coordinate_1'], 'k.', markersize=0.1,zorder=10)
    #ax.imshow(np.rot90(Z), cmap=plt.cm.Spectral_r, extent=[xmin, xmax, xmin, xmax],zorder=1)
    #ax.scatter(X, Y)
    ax.imshow(np.rot90(Z), cmap=mpl.cm.OrRd,extent=[xmin, xmax, xmin, xmax],zorder=1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([xmin, xmax])
    #plt.colorbar()
    ple_week='/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/gene_density/6_38figure/week'+str(w)+'ple.png'
    plt.savefig(ple_week)
    plt.draw()
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.5)
    cmap =  mpl.cm.OrRd
    norm = mpl.colors.Normalize(vmin=Z.min(), vmax=Z.max())
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
              orientation='vertical', label='')
    ple_colorbar='/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/count_table/human_retinal_organoid_f49b7/paper_scvis/30-38week/gene_density/6_38figure/week'+str(w)+'ple_bar.pdf'
    plt.savefig(ple_colorbar)
    plt.draw()

