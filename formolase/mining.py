import matplotlib.pyplot as plt

def get_uniprot_metadata(accessions, batch_size = 300):

    import requests
    from io import StringIO
    
    cols = 'id,entry name,reviewed,protein names,genes,organism,organism-id,fragment,annotation score,sequence'

    data = ''
    
    for istart in range(0,len(accessions),batch_size):
    
        seqids = accessions[istart:(istart+batch_size)]

        query = ('https://www.uniprot.org/uniprot/?query=%s&columns=%s&format=tab' %
                 ('+OR+'.join(['accession:%s' % s for s in seqids]), cols))

        r = requests.get(query)
        
        # remove header line
        lines = r.content.decode().split('\n')
        header = lines[0]
        data += '\n'.join(lines[1:])

    df = pd.read_table(StringIO(header+'\n'+data))
    df.columns = ['accession','alt-accession','status','description','genes','species',
              'taxon-id','fragment','annotation-score','protein-sequence']
    return df
    
    
def clustered_link_color_func(linkmat=None, clus=None, default_color='#AAAAAA', cmap = plt.cm.prism, randomize_colors = False):
    if linkmat is None:
        return lambda x: default_color

    leaf_cmap = cmap(np.linspace(0,1,len(set(clus))))
    if randomize_colors:
        leaf_cmap = np.random.permutation(leaf_cmap)
    leaf_cmap_hex = ['#%02x%02x%02x' % (int(c[0]*255),int(c[1]*255),int(c[2]*255)) for c in leaf_cmap]
    link_cols = {}
    for i, i12 in enumerate(linkmat[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(linkmat) else leaf_cmap_hex[clus[x]-1] for x in i12)
        link_cols[i+1+len(linkmat)] = c1 if np.all(c1 == c2) else default_color
    return lambda x: link_cols[x]

def clustergram(mat, linkmat1, clus=None, linkmat2=None, dendro_linewidth = None,
                cm_heatmap=plt.cm.inferno_r, cm_dendrogram=plt.cm.prism,
                randomize_dendro_colors = True,
                figsize=(8,8), heatmap_xlabel='', heatmap_ylabel='', colorbar_label=''):

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    
    if dendro_linewidth is None:
        dendro_linewidth = plt.rcParamsDefault['lines.linewidth']

    if linkmat2 is None:
        fig, axs = plt.subplots(nrows=3,ncols=4,figsize=figsize, 
                        gridspec_kw = {'width_ratios':[.5, 3, 0.01, 0.1], 'height_ratios':[3,0.1,0.1]})
        ax_heatmap = axs[0,1]
        ax_left_dgram = axs[0,0]
        ax_colorbar = axs[2,1]
        ax_annot = axs[0,3]
    else:
        fig, axs = plt.subplots(nrows=4,ncols=4,figsize=figsize, 
                        gridspec_kw = {'width_ratios':[.5, 3, 0.01, 0.1], 'height_ratios':[.5,3,0.1,0.1]})
        ax_heatmap = axs[1,1]
        ax_left_dgram = axs[1,0]
        ax_top_dgram = axs[0,1]
        ax_colorbar = axs[3,1]
        ax_annot = axs[1,3]
        
    plt.subplots_adjust(wspace=0,hspace=0)

    lcf = (clustered_link_color_func(linkmat1, clus, default_color='#808080', cmap=cm_dendrogram, 
                                     randomize_colors = randomize_dendro_colors)
           if clus is not None else lambda x: '#808080')
           
    with plt.rc_context({'lines.linewidth': dendro_linewidth}):
        dgram1 = dendrogram(linkmat1,
                         orientation="left",
                         get_leaves=True,
                         no_labels=True,
                         #labels=ids,
                         distance_sort="ascending",
                         ax=ax_left_dgram,
                         link_color_func = lcf,
                        )
    idx1 = [int(x) for x in dgram1['ivl']]
    tmp = mat[idx1,:]

    if linkmat2 is not None:
        with plt.rc_context({'lines.linewidth': dendro_linewidth}):
            dgram2 = dendrogram(linkmat2,
                             orientation="top",
                             get_leaves=True,
                             no_labels=True,
                             #labels=ids,
                             distance_sort="ascending",
                             ax=ax_top_dgram,
                             link_color_func = lambda x: '#808080',
                            )
        idx2 = [int(x) for x in dgram2['ivl']]
        tmp = tmp[:,idx2]
        ax_top_dgram.set_xlabel(heatmap_xlabel);
        ax_top_dgram.xaxis.set_label_position('top') 
    else:
        ax_heatmap.set_xlabel(heatmap_xlabel);
        ax_heatmap.xaxis.set_label_position('top') 

    img = ax_heatmap.imshow(tmp,
                 aspect="auto",origin="lower",
                 cmap=cm_heatmap,
                 vmin=np.percentile(mat,5),
                 vmax=np.percentile(mat,95),
                 interpolation="none"
                )

    ax_left_dgram.set_ylabel(heatmap_ylabel);

    for ax in axs.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    cbar = plt.colorbar(mappable=img, cax=ax_colorbar, orientation='horizontal')
    cbar.ax.set_xlabel(colorbar_label);
    
    return fig, axs, idx1, ax_annot

def clustergram_sideticks(ax, idx, idxvals=None, colors = None,
                          colormap=plt.cm.Dark2, labels=None, linewidth = None):
    from formolase.plotting import deduplicate_legend_labels

    if idxvals is None:
        idxvals = set(idx)-set([0])

    if colors is None:
        colors = colormap(np.linspace(0,1,len(idxvals)))

    if labels is None:
        labels = np.arange(len(idxvals))
    
    for i,c,lbl in zip(idxvals,colors,labels):
        pos = np.where(idx == i)[0]
        if len(pos) == 0:
            continue
        tmp = pos.reshape(len(pos),1)
        y = np.hstack((tmp,tmp))
        ax.plot([0,1],y.T,'-',c=c, label = lbl, linewidth = linewidth);

    ax.set_ylim([-.5,len(idx)-.5]);
    ax.legend(*deduplicate_legend_labels(ax), loc='center left',bbox_to_anchor=(1.01,0.5))

def load_clustalo_dist_matrix(filename):
    dist = pd.read_table(filename,delim_whitespace=True,skiprows=1,header=None)
    dist.index = dist[0]
    dist.index.name = 'accession'
    dist = dist.drop(0,axis=1)
    dist.columns = dist.index
    return dist

def get_subcluster(dist, nclus, iclus):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    lm = linkage(squareform(dist.values), method='ward')
    clus = fcluster(lm,nclus,criterion='maxclust')
    if type(iclus) is list:
        idx = np.any([clus==i for i in iclus],axis=0)
    else:
        idx = clus==iclus
    dist = dist.loc[idx,idx]
    return dist

def diverse_subsample(distmat, n):
    import numpy as np
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, fcluster
    
    lm = linkage(squareform(distmat.values), method = "ward", metric = "euclidean")
    clus = fcluster(lm,criterion="maxclust",t=n)

    idx_rep = []
    for iclus in range(min(clus),max(clus)+1):
        idx = np.array(range(len(distmat)))[clus==iclus]
        dmsubset = distmat.loc[clus==iclus,clus==iclus]
        idx_rep.append(idx[np.argmin(np.median(dmsubset,axis=1))])
    return list(distmat.index[idx_rep])

import json
import pandas as pd
import numpy as np

def get_uniprot_seqs(pidlist, ntries=3, delay=0.1):
    import requests
    from StringIO import StringIO
    from Bio import SeqIO
    import time
    host = 'http://www.uniprot.org/uniprot/'
    host2 = 'http://www.uniprot.org/uniparc/'
    seqlist = []
    for pid in pidlist:
        itry = 0
        found = False
        while itry < ntries and not found:
            req = requests.get(host+pid+'.fasta')
            if req.status_code==200:
                seqlist.append(list(SeqIO.parse(StringIO(req.content),'fasta'))[0])
                found = True
                continue
            req = requests.get(host2+pid+'.fasta')
            if req.status_code==200:
                seqlist.append(list(SeqIO.parse(StringIO(req.content),'fasta'))[0])
                found = True
                continue
            time.sleep(delay)
            itry += 1
        if not found:
            print('failed: '+pid)
    return seqlist

def hmmer_json_to_dataframe(fn):
    from collections import OrderedDict
    data = json.load(open(fn))
    reclist = []
    for dat in data['results']['hits']:
        rec = OrderedDict()
        rec['accession'] = dat['acc']
        rec['alt-accession'] = dat['acc2']
        rec['description'] = dat['desc']
        rec['kingdom'] = dat['kg']
        rec['species'] = dat['species']
        rec['taxon-id'] = dat['taxid']
        rec['e-value'] = dat['evalue']
        reclist.append(rec)
    return pd.DataFrame.from_records(reclist)

def load_ebihmmer_df(full_length_fasta_file, json_file):
    """Loads search results from HMMER EBI server (http://www.ebi.ac.uk/Tools/hmmer/search/phmmer) into a 
    DataFrame. Inputs are the "Full-length FASTA" and "JSON" files downloaded from an EBI HMMER search.
    
    Parameters
    ----------
    fn-prefix : str
        File name prefix for the full-length fasta file and json metadata file.
    
    Returns
    -------
    :class:`pandas.DataFrame`
        Table containing each sequence and associated metadata (UniProt ID, species, taxon ID, etc.)
    """
    from Bio import SeqIO
    import gzip
    import os.path

    if full_length_fasta_file[-3:] == '.gz': 
        fin = gzip.open(full_length_fasta_file)
    else: 
        fin = open(full_length_fasta_file)
        
    seqlist = list(SeqIO.parse(fin,'fasta'))
    df = pd.DataFrame.from_items([
            ('accession',[sr.id for sr in seqlist]),
            ('protein-sequence',[str(sr.seq) for sr in seqlist])
        ])
    df = df.merge(hmmer_json_to_dataframe(json_file),on='accession')

    tmp = df['accession'].copy()
    df['accession'] = df['alt-accession'].copy()
    df['alt-accession'] = tmp

    return df

def filter_ebihmmer_df(df, **kwargs):
    filter_seq_df(df, **kwargs)

def filter_seq_df(df, seq_field='protein-sequence', name_field='accession', desc_field='description'):
    """Filters out sequences typically considered bad or useless from a table of sequences.
    Assumes that sequences are formatted as if they were loaded by `load_ebihmmer_df`."""
    
    print('%d total sequences.' % len(df))
    
    df_filt = df.drop_duplicates(seq_field)
    print('Dropped duplicate sequences. %d remaining.' % len(df_filt))
    
    df_filt = df_filt.drop_duplicates(name_field)
    print('Dropped duplicate names (isoforms). %d remaining.' % len(df_filt))
    
    df_filt = df_filt[~df_filt[seq_field].str.contains('X')]
    print('Dropped ambiguous sequences (containing \'X\' residues). %d remaining.' % len(df_filt))
    
    df_filt = df_filt[df_filt[seq_field].str.match('^M')]
    print('Dropped no-START sequences. %d remaining.' % len(df_filt))

    df_filt = df_filt[~df_filt[desc_field].str.match('.*Fragment.*')]
    print('Dropped annotated fragments. %d remaining.' % len(df_filt))

    print('Reset DataFrame index (renumber rows).')
    return df_filt.reset_index(drop=True)


def load_uniref_xml(fn):
    import xml.etree.ElementTree as et
    import pandas as pd
    tree = et.parse(fn)
    reclist = []
    for entry in tree.findall('.//{http://uniprot.org/uniref}dbReference'):
        rec = {}
        rec['Name'] = entry.attrib['id']
        for prop in entry:
            if prop.attrib['type']=='UniProtKB accession':
                rec['UniProtKB'] = prop.attrib['value']
            if prop.attrib['type']=='UniParc ID':
                rec['UniParc'] = prop.attrib['value']
            if prop.attrib['type']=='protein name':
                rec['Description'] = prop.attrib['value']
            if prop.attrib['type']=='source organism':
                rec['Species'] = prop.attrib['value']
            if prop.attrib['type']=='NCBI taxonomy':
                rec['TaxID'] = prop.attrib['value']
        reclist.append(rec)
    return pd.DataFrame.from_records(
        reclist,
        columns=['Name','UniProtKB','UniParc','Description','Species','TaxID']
    )



