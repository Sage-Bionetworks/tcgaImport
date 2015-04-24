import synapseclient
from synapseclient import Table
import pandas as pd
import synapseHelpers
from multiprocessing.dummy import  Pool

FILEQUERY = ("select * from file where benefactorId=='syn2812961' "
                               "and fileType!='clinicalMatrix'"
                               "and fileType!='maf'")
#TABLEID='syn3505659'  #New
TABLEID='syn3281840' #Original
syn = synapseclient.login()


def findUpdates(files, tableId=TABLEID):
    """Compares the lastmodified date in the table and file anntoations.
    
    files is a table of file medata including the column lastModified.
    """
    table = syn.tableQuery('SELECT distinct id, versionNumber FROM %s' %tableId)
    df = table.asDataFrame().set_index('id')
    toUpdate=[]
    for i, f in files.iterrows():
        if f.id not in df.index:
            toUpdate.append(i)
        else:
            if df.ix[f.id,'versionNumber'] < f.versionNumber:
                toUpdate.append(i)
    return files.ix[toUpdate,:]


def deleteAffectedRows(updatedFiles, tableId=TABLEID):
    """Given a list of file Entities.  Extracts all rows in the table
    that match the id column in the updatedFiles data frame and
    deletes them from Synapse."""

    #Build query to find rows to delete for update
    idSet = "('" + "','".join(updatedFiles.id) +"')"
    queryStr = "select id from %s where id in %s" %(tableId, idSet)
    #Query and delete
    rowSet = syn.tableQuery(queryStr).asRowSet()
    if len(rowSet['rows'])>0:
        print 'DELETING', len(rowSet['rows']), 'ROWS'
        syn.delete(rowSet)


def countAndUpdateTable(input, path='out/', tableId=TABLEID):
    i, fileMeta = input
    if fileMeta.fileType=='ttl':
        return []
    print 'updating table:%s' %tableId, i, fileMeta.id, fileMeta['name'], fileMeta['basename']
    ent = syn.get(fileMeta.id, downloadLocation=path)
    if fileMeta.fileType =='bed5':
        data = pd.read_csv(ent.path, sep='\t', header=None)
        nFeatures = 0
        samples = list(set(data[3].dropna()))
    else: #All other fileTypes
        data = pd.read_csv(ent.path, sep='\t', index_col=0)
        nFeatures, nSamples = data.shape
        samples = data.columns
    metadata = pd.DataFrame([fileMeta]*len(samples))
    metadata['nFeatures'] = nFeatures
    metadata['samples'] = samples
    metadata['patient_barcode'] = [x[:12] for x in metadata.samples]
    metadata.drop(['projectId', 'tissue', u'md5', u'assembly'], axis=1, inplace=True)
    metadata.nFeatures = metadata.nFeatures.astype('int')

    #Update rows in table
    print 'adding', metadata.shape[0]
    t = syn.store(Table(tableId, metadata))
    return metadata



if __name__ == '__main__':
    files = synapseHelpers.query2df(syn.chunkedQuery(FILEQUERY), savedSynapseFields=('id', 'name', 'versionNumber'))
    updatedFiles= findUpdates(files)
    print 'NEED TO UPDATE:', updatedFiles.shape[0], 'FILES'
    deleteAffectedRows(updatedFiles)
    dfs = map(countAndUpdateTable, updatedFiles.iterrows())
