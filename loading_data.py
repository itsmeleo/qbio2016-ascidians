def loading_data(path):
    cell_fate2={
        'A line Endoderm': (["a7.0001", "a7.0002", "a7.0005"], 1),
        'B line Endoderm':(["b7.0001", "b7.0002", "b9.0034"], 2),

        'Germline': (["b7.0006"], 3),

        'Mesoderm Notochord 1':(["a7.0003", 'a7.0007'], 4),
        'Mesoderm Notochord 2':(["b8.0006"], 5),
        'Mesoderm Trunk Lateral Cell':(['a7.0006'], 6),
        'Mesoderm Trunk Ventral Cell':(['b7.0005'], 7),
        'Mesoderm Muscle 1':(['b7.0004', 'b7.0008'], 8),
        'Mesoderm Muscle 2':(['a9.0031', 'b9.0033'], 9),
        'Mesoderm Mesenchyme':(['b7.0007', 'b8.0005'], 10),

        'Neural Plate Tail Lateral':(['a8.0015', 'a9.0032'], 11),
        'Neural Plate Tail Ventral':(['a9.0015', 'a9.0013'], 11),
        'Neural Plate Tail Dorsal':(['b8.0019'], 11),
        'Neural Plate Head Ventral': (['a9.0016', 'a9.0014', 'a7.0009', 'a7.0010'],11),
        'Neural Plate Head Dorsal': (['a7.0013'],11),
    
        'Epidermis Head': (['a7.0011', 'a7.0012', 'a7.0014', 'a7.0015', 'a7.0016'], 14),
        'Epidermis Tail': (['b8.0020', 'b8.0018', 'b7.0011', 'b7.0012', 'b7.0013', 'b7.0014', 'b7.0015', 'b7.0016'], 15)
    }


    cell_fate3={
        'Endoderm Head':(["a7.0001", "a7.0002", "a7.0005", "b7.0001", "b8.0003", "b9.0007"], 1),
        'Endodermal Strand 1':(["b9.0008"], 2),
        'Endodermal Strand 2':(["b9.0034"], 2),
        
        'Germline': (["b7.0006"], 3),
        
        'Mesoderm Notochord 1':(["a7.0003", 'a7.0007'], 4),
        'Mesoderm Notochord 2':(["b8.0006"], 5),
        'Mesoderm Trunk Lateral Cell':(['a7.0006'], 6),
        'Mesoderm Trunk Ventral Cell':(['b7.0005'], 7),
        'Mesoderm Muscle 1':(['b7.0004', 'b7.0008'], 8),
        'Mesoderm Muscle 2':(['a9.0031', 'b9.0033'], 9),
        'Mesoderm Mesenchyme':(['b7.0007', 'b8.0005'], 10),
        
        'Neural Plate Tail Lateral':(['a8.0015', 'a9.0032'], 11),
        'Neural Plate Tail Ventral':(['a9.0015', 'a9.0013'], 11),
        'Neural Plate Tail Dorsal':(['b8.0019'], 11),
        'Neural Plate Head Ventral': (['a9.0016', 'a9.0014', 'a7.0009', 'a7.0010'],11),
        'Neural Plate Head Dorsal': (['a7.0013'],11),
    
        'Epidermis Head': (['a7.0011', 'a7.0012', 'a7.0014', 'a7.0015', 'a7.0016'], 14),
        'Epidermis Tail': (['b8.0020', 'b8.0018', 'b7.0011', 'b7.0012', 'b7.0013', 'b7.0014', 'b7.0015', 'b7.0016'], 15)
    }

    cell_fate={
        'Endoderm': (["a7.0001", "a7.0002", "a7.0005", "b7.0001", "b7.0002", "b9.0034"], 7),
        'germ line': (["b7.0006"], 3),
        'Epidermis': (['a7.0011', 'a7.0012', 'a7.0014', 'a7.0015', 'a7.0016', 'b8.0020', 'b8.0018', 'b7.0011', 'b7.0012', 'b7.0013', 'b7.0014', 'b7.0015', 'b7.0016'], 2),
        'Mesoderm': (["a7.0003", 'a7.0006','a7.0007', "b7.0003", 'b7.0005', 'b7.0004', 'b7.0008', 'a9.0031', 'b9.0033', 'b7.0007'], 6),
        'Nervous system': (["a7.0004", "a7.0009", "a7.0010", 'a7.0013', 'a8.0015', 'a9.0032', 'b8.0019'], 5)
    }

    import numpy as np

    # TOP_1=(151*2, 237*2, 46*2)
    # BOTTOM_1=(247*2, 208*2, 260*2)

    # TOP_2=(156*2, 230*2, 50*2)
    # BOTTOM_2=(261*2, 234*2, 272*2)

    # TOP_3=(159*2, 277*2, 15*2)
    # BOTTOM_3=(239*2, 174*2, 270*2)

    # TOP_4=(175*2, 321*2, 36*2)
    # BOTTOM_4=(226*2, 53*2, 235*2)

    # TOP_5=(117*2, 320*2, 33*2)
    # BOTTOM_5=(234*2, 43*2, 241*2)

    # pos_AP={1  :(TOP_1, BOTTOM_1),
    #         60 :(TOP_2, BOTTOM_2),
    #         120:(TOP_3, BOTTOM_3),
    #         180:(TOP_4, BOTTOM_4),
    #         192:(TOP_5, BOTTOM_5)}

    # tmp_top=np.array([np.linspace(TOP_1[0], TOP_2[0], 60-1), np.linspace(TOP_1[1], TOP_2[1], 60-1), np.linspace(TOP_1[2], TOP_2[2], 60-1)])
    # tmp_bottom=np.array([np.linspace(BOTTOM_1[0], BOTTOM_2[0], 60-1), np.linspace(BOTTOM_1[1], BOTTOM_2[1], 60-1), np.linspace(BOTTOM_1[2], BOTTOM_2[2], 60-1)])
    # pos_AP.update({t:(tmp_top[:,t-1], tmp_bottom[:,t-1]) for t in range(1, 60)})

    # tmp_top=np.array([np.linspace(TOP_2[0], TOP_3[0], 120-60), np.linspace(TOP_2[1], TOP_3[1], 120-60), np.linspace(TOP_2[2], TOP_3[2], 120-60)])
    # tmp_bottom=np.array([np.linspace(BOTTOM_2[0], BOTTOM_3[0], 120-60), np.linspace(BOTTOM_2[1], BOTTOM_3[1], 120-60), np.linspace(BOTTOM_2[2], BOTTOM_3[2], 120-60)])
    # pos_AP.update({t:(tmp_top[:,t-60], tmp_bottom[:,t-60]) for t in range(60, 120)})

    # tmp_top=np.array([np.linspace(TOP_3[0], TOP_4[0], 180-120), np.linspace(TOP_3[1], TOP_4[1], 180-120), np.linspace(TOP_3[2], TOP_4[2], 180-120)])
    # tmp_bottom=np.array([np.linspace(BOTTOM_3[0], BOTTOM_4[0], 180-120), np.linspace(BOTTOM_3[1], BOTTOM_4[1], 180-120), np.linspace(BOTTOM_3[2], BOTTOM_4[2], 180-120)])
    # pos_AP.update({t:(tmp_top[:,t-120], tmp_bottom[:,t-120]) for t in range(120, 180)})

    # tmp_top=np.array([np.linspace(TOP_4[0], TOP_5[0], 192-180), np.linspace(TOP_4[1], TOP_5[1], 192-180), np.linspace(TOP_4[2], TOP_5[2], 192-180)])
    # tmp_bottom=np.array([np.linspace(BOTTOM_4[0], BOTTOM_5[0], 192-180), np.linspace(BOTTOM_4[1], BOTTOM_5[1], 192-180), np.linspace(BOTTOM_4[2], BOTTOM_5[2], 192-180)])
    # pos_AP.update({t:(tmp_top[:,t-180], tmp_bottom[:,t-180]) for t in range(180, 192)})

    class ColorMap(object):

        def __init__(self, values, cmap='jet', default='w'):
            values = sorted(values)
            if None in values: values.remove(None)        
            if all([isinstance(i, int) for i in values]):
                self.values = range(values[0], values[-1]+1)
                self.cmap = plt.cm.get_cmap(cmap, len(self.values))
            elif all([isinstance(i, float) for i in values]):
                self.values = None   
                self.cmap = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=values[0], vmax=values[-1]), cmap=cmap)
            else:
                self.values = values
                self.cmap = plt.cm.get_cmap(cmap, len(self.values))
            self.default = default
            
        def __call__(self, value):
            if value is None:
                return self.default
            elif self.values is None:
                return self.cmap.cmap(self.cmap.norm(value))
            else:
                if not value in self.values:
                    raise ValueError('`value` parameter')
                return self.cmap(self.values.index(value))

    import cPickle as pkl
    import numpy as np
    from matplotlib import pyplot as plt

    f=open(path+'final_properties.pkl')
    lin_tree, properties=pkl.load(f)
    f.close()
    names=properties['Names'][0]
    inv_names={v:k for k, v in names.iteritems()}
    fates=properties['fate'][0]
    fates2=properties['fate2'][0]
    vol=properties['volumes_information'][0]
    inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }
    surf_ex=properties['cell_cell_surface_information'][0]
    surfaces=properties['surface'][0]
    bary=properties['barycenter'][0]
    fates_set=list(set(fates2.values()))
    fates_set.append('ext')
    fates_set.append('undeter')
    CMapFates2=ColorMap(fates_set, 'jet')
    CMapTime=ColorMap(range(1, 193), 'jet')
    fates2[1]='ext'
    fates2[0]='undeter'
    col={'Epidermis':'c',
         'Endoderm':'r',
         'Mesoderm':'g',
         'Nervous system':'b'
         }
    cell_col={1:'b',
             2:'g',
             3:'r',
             4:'c',
             5:'m',
             6:'y',
             0:'k'}

    def find_color(name, cell_fate):
        for k, v in cell_fate.iteritems():
            if name[:-1].replace(' ', '') in v[0]:
                return k
        return ''
    nodes=list(set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values])))
    fates2={}
    for t in range(1, 193, 1):
        for n in nodes:
            if t*10**4<n<(t+1)*10**4: 
                if fates2.get(inv_lin_tree.get(n, ''), '')!='':
                    fates2[n]=fates2[inv_lin_tree[n]]
                else:
                    col=find_color(names.get(n, ''), cell_fate2)
                    if col!='':
                        fates2[n]=col

    fates2[1]='ext'
    fates2[0]='undeter'
    fates3={}
    for t in range(1, 193, 1):
        for n in nodes:
            if t*10**4<n<(t+1)*10**4: 
                if fates3.get(inv_lin_tree.get(n, ''), '')!='':
                    fates3[n]=fates3[inv_lin_tree[n]]
                else:
                    col=find_color(names.get(n, ''), cell_fate3)
                    if col!='':
                        fates3[n]=col
    fates3[1]='ext'
    fates3[0]='undeter'
    properties.pop('axis_n')
    properties.pop('Errors')
    properties.pop('vals_clone_15')
    # properties.pop('surface')
    properties.pop('vals_fate')
    properties.pop('barycenter')
    # properties.pop('Names')
    properties.pop('neighbourhood')
    # properties.pop('Sim no vol')
    properties.pop('axis_mother')
    properties.pop('axis_clone')
    # properties.pop('fate')
    properties.pop('Cum_diff')
    properties.pop('vals_n')
    # properties.pop('h_mins_information')
    properties.pop('vals_mother')
    properties.pop('axis_fate')
    # properties.pop('cell_cell_surface_information')
    # properties.pop('volumes_information')
    properties.pop('axis_clone_15')
    properties.pop('surf_fate')
    # properties.pop('fate3')
    # properties.pop('fate2')
    # properties.pop('sigmas_information')
    properties.pop('vals_clone')
    properties['compactness'] = properties.pop('compacity')
    cells=[k for k in lin_tree.keys() if k/10**4==1]

    cells=[k for k in lin_tree.keys() if k/10**4==1]
    for c in cells:
        vol[c]=vol[lin_tree[c][0]]
    
    return lin_tree, fates, fates2, fates3, vol, inv_lin_tree, surf_ex, surfaces, properties['Names'][0], properties, ColorMap
