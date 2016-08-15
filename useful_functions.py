import numpy as np
from matplotlib import pyplot as plt
import os
import scipy.cluster.hierarchy as H
from copy import deepcopy

def get_ratio_mother_sisters_volume(mother, vals, lin_tree, inv_lin_tree, around=5):
    """ From a mother cell id (right before a division),
        computes the average value for a given attribute `vals` accros `around` timepoints
        for the mother cell and the daughter cells
        Return the ratio of the sum of the average of attribute 
        of the two sisters over the average of the attribute for the mother.
        Args:
            mother(int): id of the mother
            vals(dict): dictionary that maps a cell id to a numerical value
            around(int): time window to consider in time-points
        Returns:
            (int): The ration of the sum of the two sister cell values
                   over the mother cell value.
    """
    before = [mother]
    for i in range(around):
        before.append(inv_lin_tree.get(before[-1], before[-1]))
    after1 = [lin_tree[mother][0]]
    after2 = [lin_tree[mother][1]]
    for i in range(around):
        after1.append(lin_tree.get(after1[-1], [after1[-1]])[0])
        after2.append(lin_tree.get(after2[-1], [after2[-1]])[0])
    return ((np.mean([vals.get(c, 0) for c in after1]) + np.mean([vals.get(c, 0) for c in after2])) / 
            np.mean([vals.get(c, 0) for c in before]))

def get_symetric_volumes(couple, vals, lin_tree):
    """ Computes the average values of a given couple of cells
        throughout their life-span of the attribute `vals`
        Args:
            couple(list): list of two cells (usually symetrical cells but not necessary)
            vals(dict): dictionary that maps a cell id to a numerical value
        Returns:
            ((int), (int)): The average for the two cells of 
                            the given attribute throughout their life-spans
    """
    c1 = [couple[0]]
    c2 = [couple[1]]
    while len(lin_tree.get(c1[-1], []))==1 or len(lin_tree.get(c2[-1], []))==1:
        if len(lin_tree.get(c1[-1], []))==1:
            c1.append(lin_tree[c1[-1]][0])
        if len(lin_tree.get(c2[-1], []))==1:
            c2.append(lin_tree[c2[-1]][0])
    return np.mean([vals.get(c) for c in c1]), np.mean([vals.get(c) for c in c2])

def get_values_around_division(c, vals, time_around, lin_tree, inv_lin_tree):
    """Compute the metric for a cell c around its division. 
    Before the division the metric used is the one of the mother,
    after the division it is the one of the clone of the mother.
    Args:
        c(int): cell id from the lineage tree
        vals(dict): dictionary that maps a cell id to a numerical value
        time_around(int): the time to consider before and after the division
    
    Returns:
        4 lists:
            The metric for the mother from *time_around* before division to the division`
            The metric for the daughter 1 from the division to *time_around* after it
            The metric for the daughter 2 from the division to *time_around* after it
    """
    mother_life = [c]
    for i in range(time_around):
        mother_life = [inv_lin_tree.get(mother_life[0], mother_life[0])] + mother_life
    sister_life = [lin_tree[c][0]]
    sister_life2 = [lin_tree[c][1]]
    for i in range(time_around):
        sister_life = sister_life + [lin_tree.get(sister_life[-1], [sister_life[-1]])[0]]
        sister_life2 = sister_life2 + [lin_tree.get(sister_life2[-1], [sister_life2[-1]])[0]]
    return ([vals.get(c, vals[lin_tree[c][0]]) for c in mother_life],
            [vals.get(c, None) for c in sister_life], 
            [vals.get(c, None) for c in sister_life2])

def plot_around_div(vals, lin_tree, inv_lin_tree, interest_cells, names, col=6, row=3, ylim=(), x_label='', y_label='', 
                    around=20, title='', saving_folder=''):
    """ For the cell cycles from 7 to 10, plots the evolution of a given attribute
        `vals` around the division events.
        Args:
            vals(dict): dictionary that maps a cell id to a numerical value
            col(int): #columns in the plot
            row(int): #rows in the plot
            ylim(tuple(int, int)): fixed ylim of the plot
            x_label(string): label of the x axis
            y_label(string): label of the y axis
            around(int): number of time-points around the divisions to consider
            title(string): title of the plot
            saving_folder(string): path to the directory where to save the plot
    """
    fig = plt.figure(figsize=(10, 8))
    i = 1
    for z_c in range(7, 11):
        ax = fig.add_subplot(row, col, i)
        whole = []
        whole_sis = []
        interest_cells_f = [c for c in interest_cells if int(names.get(c).split('.')[0][1:]) == z_c]
        for c in interest_cells_f:
            Yb, Yfs1, Yfs2 = get_values_around_division(c, vals, around, lin_tree, inv_lin_tree)
            Xb = range(-len(Yb)*2+2, 2, 2)
            Xe = range(0, max(len(Yfs1), len(Yfs2))*2, 2)
            Y = Yb
            Y1 = Yfs1
            ax.plot(Xe, Y1, '-', alpha=.2, color='k')
            ax.plot(Xb, Y, '-', alpha=.2, color='k')
            whole.append(list(Y))
            Y2 = Yfs2
            ax.plot(Xe, Y2, '-', alpha=.2, color='k')
            whole_sis.append(Y1)
            whole_sis.append(Y2)

        ax.set_xlim(2-len(Xb)*2, (len(Xe)-1)*2)            
        if i%row==1 or col==1:
            ax.set_ylabel(y_label, fontsize=30)
        else:
            ax.set_yticks([])
        if (i>=col*row-1 and col!=1) or (col==1 and i==row):
            ax.set_xlabel(x_label, fontsize=30)
        else:
            ax.set_xticks([])
        ax.text(.95,.9,'%d cell cycle'%z_c,
            horizontalalignment='right',
            transform=ax.transAxes, fontsize=20,
            bbox={'facecolor':'white', 'pad':10, 'linewidth':3})
        if i==1:
            ax.plot([0, 0], [ylim[0], ylim[1]], 'k--', label='Moment of division')            
            whole = np.array(whole)
            ax.plot(Xb, np.median(whole, axis=0), 'r-', lw=3, label='Distribution median')
            ax.plot(Xb, np.percentile(whole, 25, axis=0), 'r--', lw=2)
            ax.plot(Xb, np.percentile(whole, 75, axis=0), 'r--', lw=2)
            
            X=range(0, len(whole_sis[0])*2, 2)
            whole_sis = np.array(whole_sis)
            ax.plot(Xe, np.median(whole_sis, axis=0), 'r-', lw=3)
            ax.plot(Xe, np.percentile(whole_sis, 25, axis=0), 'r--', lw=2)
            ax.plot(Xe, np.percentile(whole_sis, 75, axis=0), 'r--', lw=2)
            ax.legend(fontsize=15, loc='lower left')
        else:
            ax.plot([0, 0], [ylim[0], ylim[1]], 'k--')      
            whole = np.array(whole)
            ax.plot(Xb, np.median(whole, axis=0), 'r-', lw=3)
            ax.plot(Xb, np.percentile(whole, 25, axis=0), 'r--', lw=2)
            ax.plot(Xb, np.percentile(whole, 75, axis=0), 'r--', lw=2)
            
            X=range(0, len(whole_sis[0])*2, 2)
            whole_sis = np.array(whole_sis)
            ax.plot(Xe, np.median(whole_sis, axis=0), 'r-', lw=3)
            ax.plot(Xe, np.percentile(whole_sis, 25, axis=0), 'r--', lw=2)
            ax.plot(Xe, np.percentile(whole_sis, 75, axis=0), 'r--', lw=2)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_ylim(ylim)

                
        i += 1
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.02)
    if saving_folder!='':
        newpath=saving_folder
        if not os.path.exists(newpath): os.makedirs(newpath)
        fig.savefig(newpath+y_label+'_around_division.pdf')

def get_x_y(distribution, bins, range=(0, 0.8), normed=False):
    """ From a distribution computes an histogram and 
        returns the x-y position of the bins
        Args:
            distribution(list): list of numerical values
            bins(int): number of bins to compute
            range(float, float): range of the histogram
            normed(bool): whether the histogram should be normed or not
        Returns:
            X(list): x positions of the bins
            Y(list): y positions of the bins
    """
    vals=np.histogram(distribution, bins=bins, range=range, normed=normed)
    X=np.array(vals[1])
    X=(X[1:]+X[:-1])/2.
    return X, np.array(vals[0])

def build_distance_matrix(tree_distances, size=64):
    """ Build a size x size matrix of distances between lineage trees
        Args:
            tree_distance(dict): dictionary that maps a couple of 
                                 cell ids onto their pairwise distance
            size(int): size of the squared matrix
        Returns:
            X(np.array): size x size np.array where X[i, j] is the distance between i and j
            corres(dict): a dictionary that maps a cell id onto an index in X
            corres_inv(dict): a dictionary that maps an index in X onto a cell id
    """
    corres = {}
    index = 0
    X = np.zeros((size, size))
    for c1, c2 in tree_distances.keys():
        if c1<c2:
            if not corres.has_key(c1):
                corres[c1] = index
                index += 1
            if not corres.has_key(c2):
                corres[c2] = index
                index += 1

    corres_inv={v:k for k, v in corres.iteritems()}
    for i in range(max(corres.values())+1):
        for j in range(i+1, max(corres.values())+1):
            c1, c2=corres_inv[i], corres_inv[j]
            if tree_distances.has_key((c1, c2)):
                X[i, j]=tree_distances[(c1, c2)]
                X[j, i]=tree_distances[(c1, c2)]
            else:
                X[i, j]=tree_distances[(c2, c1)]
                X[j, i]=tree_distances[(c2, c1)]
    return X, corres, corres_inv

def format_name(n):
    """ Function that format a name n removing the unnecessary '0'
        while keeping its original lenght by adding spaces at the end
        Args:
            n(string): cell name formated as follow: <a/b>#.####<*/_>
        Returns:
            (string): the formated name
    """
    size = len(n)
    first_part = n.split('.')[0]+'.'
    second_part = str(int(n.split('.')[1][:-1]))
    last_part = n.split('.')[1][-1]
    all_name = first_part + second_part + last_part
    #print n, size - len(all_name)
    for i in range(size - len(all_name)):
        all_name+=' '
    return all_name
    
def perform_and_plot_dendrogram(D1, corres, corres_inv, f, ColorMap, names,
                                method='ward', saving_folder='', distance_name='',
                                color_threshold=None):
    """ Given a distance matrix `D1` performs a hierarchical linkage
        and plot the corresponding dendrogram
        Args:
            D1(np.array): a n x n symmetrical matrix
            corres(dict): a dictionary that maps a cell id onto an index in D1
            corres_inv(dict): a dictionary that maps an index in D1 onto a cell id
            f(dict): fate map dictionary
            saving_folder(string): path to the directory where to save the plot
            distance_name(string): distance name to include in the title
            color_threshold(float): threshold to cut the tree and assign different
                                    colors to different resulting clusters
    """
    CMapFates=ColorMap(list(set(f.values()))+['Undetermined'], 'rainbow')
    Y = H.linkage(D1, method=method)
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(111)
    tmp = H.dendrogram(Y, leaf_rotation=.5, color_threshold=color_threshold)
    labels = []
    for c in tmp['leaves']:
        n = ''
        for s in f.get(corres_inv[c], 'Undetermined').split(' '):
            n += s[0]
        labels.append(format_name(names[corres_inv[c]])+' '+f.get(corres_inv[c], 'Undetermined'))
        
    ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=30)
    for i in plt.gca().get_xticklabels():
        i.set_color(CMapFates(i.get_text()[9:]))
        i.set_weight('bold')
        i.set_fontsize(30)
    fig.subplots_adjust(left=0.04, bottom=0.6, right=0.96, top=0.96)
    ax.tick_params(axis='both', which='major', labelsize=22)
    if saving_folder!='':
        fig.savefig(saving_folder+'dendogram.pdf')
def get_min_vol(k, vals, lin_tree, size=10, decal=5):
    """ Given a cell id return the average attribute of its two daughter that is minimum
        (starting after `decal` time points for `size` time points)
        Args:
            k(int): a cell id that divide next time point
            vals(dict): dictionary that maps a cell id to a numerical value
            size(int): number of time point over which to average
            decal(int): number of time points not to consider after the division
        Returns:
            (float): the average attribute of its two daughter that is minimum
    """
    c1=lin_tree[k][0]
    c2=lin_tree[k][1]
    out=[]
    for i in range(decal):
        c1=lin_tree.get(c1, [c1])[0]
        c2=lin_tree.get(c2, [c2])[0]
    for i in range(size):
        out.append([vals[c1], vals[c2]])
        c1=lin_tree.get(c1, [c1])[0]
        c2=lin_tree.get(c2, [c2])[0]
    return min(np.mean(out, axis=0))

def get_max_vol(k, vals, lin_tree, size=10, decal=5):
    """ Given a cell id return the average attribute of its two daughter that is maximum
        (starting after `decal` time points for `size` time points)
        Args:
            k(int): a cell id that divide next time point
            vals(dict): dictionary that maps a cell id to a numerical value
            size(int): number of time point over which to average
            decal(int): number of time points not to consider after the division
        Returns:
            (float): the average attribute of its two daughter that is maximum
    """
    c1=lin_tree[k][0]
    c2=lin_tree[k][1]
    out=[]
    for i in range(decal):
        c1=lin_tree.get(c1, [c1])[0]
        c2=lin_tree.get(c2, [c2])[0]
    for i in range(size):
        out.append([vals[c1], vals[c2]])
        c1=lin_tree.get(c1, [c1])[0]
        c2=lin_tree.get(c2, [c2])[0]
    return max(np.mean(out, axis=0))

def short_name(n):
    """ Shorten a cell name
        Args:
            n(string): cell name
        Returns:
            (string): shorted cell name
    """
    return n.split('.')[0]+'.'+str(int(n.split('.')[1][:-1]))

def get_mother_name(n):
    """ From a cell name `n` gives the name that gave rise to this cell name
        Args:
            n(string): cell name
        Returns:
            (string): Name of the mother
    """
    letter=n[0]
    z_c=str(int(n.split('.')[0][1:])-1)
    num='%04d'%np.ceil(float(n.split('.')[1][:-1])/2)
    end=n[-1]
    return letter+z_c+'.'+num+end

def get_sister_name(n):
    """ From a cell name `n` gives the sister name
        Args:
            n(string): cell name
        Returns:
            (string): Name of the sister
    """
    return n.split('.')[0]+'.'+'%04d'%(int(n.split('.')[1][:-1])+1)+n[-1]

def generate_32_cell_stage(lin_tree, names, inv_lin_tree, sim_tree, sim_nv, vol):
    """ From a given lineage tree, names, lineage tree distances and volumes
        build a new lineage tree similar of the initial one with added the
        mother cells of the first cells of the initial lineage tree
        Args:
            lin_tree(dict): lineage tree
            names(dict): names dictionary
            inv_lin_tree(dict): inv_lin_tree dictionary
            sim_tree(dict): sim_tree dictionary
            sim_nv(dict): sim_nv dictionary
            vol(dict): vol dictionary
        Returns:
            new_lin_tree(dict): new lin_tree
            new_names(dict): new names
            new_vol(dict): new vol
            new_sim_nv(dict): new sim_nv
    """
    sim_tree = {(2*10**4+k[0], 2*10**4+k[1]):v for k, v in sim_tree}
    new_lin_tree = deepcopy(lin_tree)
    new_names = deepcopy(names)
    new_vol = deepcopy(vol)
    new_sim_nv = deepcopy(sim_nv)
    cells = [k for k in lin_tree.keys() if k/10**4==1]
    i = 1
    for c in cells:
        if int(names[c].split('.')[1][:-1])%2==1:
            s_n = get_sister_name(names[c])
            for ci in cells:
                if names[ci] == s_n:
                    if not inv_lin_tree.has_key(c):
                        new_lin_tree[i] = [c, ci]
                    if sim_tree.has_key((new_lin_tree[c][0], new_lin_tree[ci][0])):
                        new_sim_nv[i] = sim_tree[(new_lin_tree[c][0], new_lin_tree[ci][0])]
                    else:
                        new_sim_nv[i] = sim_tree[(new_lin_tree[ci][0], new_lin_tree[c][0])]
                    new_names[i] = get_mother_name(names[c])
                    new_vol[i] = vol[c] + vol[ci]
                    i += 1
    return new_lin_tree, new_names, new_vol, new_sim_nv
