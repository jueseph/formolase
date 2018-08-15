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
    return pd.DataFrame.from_items([
        ('well',['%s%02d' % (a,int(b)) for (a,b) in wells]),
        (value_name,[df.loc[a,b] for (a,b) in wells])
    ])

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
