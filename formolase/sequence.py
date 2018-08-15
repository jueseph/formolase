def df_to_fasta(df, filename, name_field='accession', seq_field = 'protein-sequence', description_field = None):

    with open(filename,'w') as f:
        for i,row in df.iterrows():
            print('>%s' % row[name_field]+((' %s' % row[description_field]) if description_field is not None else ''), file=f)
            print(row[seq_field],file=f)

def contains_restriction_sites(enzymes,sequence):
    import re
    found = [re.search(e.site,sequence,flags=re.I)!=None for e in enzymes]
    foundenzymes = [e for e,f in zip(enzymes,found) if f]
    return (any(found), foundenzymes)

def remove_start_stop(df, seq_field='protein-sequence', mut_field='protein-mutations'):
    idx = df[seq_field].str.contains('^M')
    df.loc[idx,seq_field] = [s[1:] for s in df.loc[idx,seq_field]]
    if mut_field is not None:
        df.loc[idx,mut_field] = df.loc[idx,mut_field].apply(lambda x: 'd1' if x=='wt' else 'd1 '+x)
    n1 = sum(idx)
        
    idx = df[seq_field].str.contains('\*$')
    df.loc[idx,seq_field] = [s[:-1] for s in df.loc[idx,seq_field]]
    n2 = sum(idx)
    
    print('Removed %d STARTs and %d STOPs among %d sequences total.' % (n1,n2,len(df)))
    
    return df

def truncate_phobius_sp(df, seq_field='protein-sequence', mut_field='protein-mutations'):
    import re
    for i,row in df.iterrows():
        if '/' in row['PREDICTION']:
            altstart = int(re.search('(\d*)$',row['PREDICTION'].split('/')[0]).groups()[0]) # 1-indexed
            newseq = 'G'+df.loc[i,seq_field][altstart:]
            df.loc[i,seq_field] = newseq
            if mut_field is not None:
                mutstr = 'd1-%d %s%dG%s' % (altstart, 
                                            newseq[1],
                                            altstart+1,
                                            newseq[1]) 

                currmutstr = df.loc[i,mut_field]
                df.loc[i,mut_field] = mutstr if currmutstr=='wt' else mutstr+' '+currmutstr
    return df
