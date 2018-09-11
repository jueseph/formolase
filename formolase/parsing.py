import pandas as pd
import numpy as np

def well2coord(well):
    well = well.upper()
    row = 'ABCDEFGH'.index(well[0])
    col = int(well[1:])-1
    return (row,col)

def coord2well(coord):
    r,c = coord
    return 'ABCDEFG'[r]+str(c+1)

def flatten_96w_plate(df, value_name, well_order='row major'):
    if well_order == 'row major':
        wells = [(a,b) for a in df.index for b in df.columns]
    elif well_order == 'column major':
        wells = [(a,b) for b in df.columns for a in df.index]
    from collections import OrderedDict
    return pd.DataFrame.from_dict(OrderedDict([
        ('well',['%s%02d' % (a,int(b)) for (a,b) in wells]),
        (value_name,[df.loc[a,b] for (a,b) in wells])
    ]))

def load_96w_metadata(filename, well_order='row major'):

    xl = pd.ExcelFile(filename)

    metadata = None

    for sheet in xl.sheet_names:
        tmp = xl.parse(sheet, header=None)

        num_pad_cols = 12 - len(tmp.columns)

        tmp = pd.concat([tmp,pd.DataFrame(np.full((len(tmp),num_pad_cols),np.nan))],axis=1)

        tmp.columns = range(1,13)
        tmp.index = list('ABCDEFGH')

        tmp = flatten_96w_plate(tmp,sheet, well_order=well_order)

        if metadata is None:
            metadata = tmp
        else:
            metadata = metadata.merge(tmp, on = 'well', how = 'outer')

    return metadata

def load_spectramax_plate(filename, label=False):
    dat = pd.read_table(filename,skiprows=3, header = None)
    dat = dat.iloc[:8,2:14]
    dat.columns = range(1,13)
    dat.index = list('ABCDEFGH')
    if label:
        dat = flatten_96w_plate(dat, label)
    return dat

def load_spectramax_timecourse(filename, wells=None, metadatafile=None):

    df = pd.read_table(filename,skiprows=3,skipfooter=1,header=None, engine='python')

    rowlist = []

    for iplate in range(int(df.shape[0]/9)):
        time = df.iloc[iplate*9,0]
        temp = df.iloc[iplate*9,1]
        dfplate = df.iloc[iplate*9:(iplate*9+8),2:14]
        dfplate.columns = range(1,13)
        dfplate.index = list('ABCDEFGH')
        dfplate = flatten_96w_plate(dfplate,time).T
        dfplate.columns = dfplate.loc['well']
        dfplate = dfplate.drop('well')
        dfplate.insert(0,'Temperature',[temp])
        dfplate.insert(0,'Time',[time])
        dfplate.columns.name = None
        rowlist.append(dfplate)

    df2 = pd.concat(rowlist).reset_index(drop=True)
            
    def convert_time_to_minutes(string):
        tokens = string.split(':')
        return int(tokens[0])*60 + int(tokens[1])
    
    df2.insert(0,'Seconds',df2['Time'].apply(convert_time_to_minutes))
    df2 = df2.drop('Time',axis=1)
    
    if wells is not None:
        df2 = df2[['Seconds','Temperature']+wells]
    
    if metadatafile is not None:
        meta = load_96w_metadata(metadatafile)
        meta = meta.T
        meta.columns = meta.loc['well']
        meta.columns.name = None
        meta = meta.drop('well')
        meta.insert(0,'Temperature',np.nan)
        meta.insert(0,'Seconds',np.nan)
        df2 = meta.append(df2, sort=False)
        df2 = df2.loc[:,meta.iloc[0] != 'blank']
    
    return df2.T

