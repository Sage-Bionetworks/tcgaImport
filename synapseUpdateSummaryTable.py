import argparse
import synapseclient
from synapseclient import Table
import pandas as pd
import synapseHelpers
from multiprocessing.dummy import  Pool

FILEQUERY = ("select * from file where benefactorId=='%s' "
                               "and fileType!='clinicalMatrix' "
                               "and fileType!='maf' "
                               "and fileType!='ttl'")
syn = synapseclient.login()


def findUpdates(files, tableId):
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


def deleteAffectedRows(updatedFiles, tableId):
    """Given a list of file Entities.  Extracts all rows in the table
    that match the id column in the updatedFiles data frame and
    deletes them from Synapse."""
    if len(updatedFiles)==0:
        return
    
    #Build query to find rows to delete for update
    idSet = "('" + "','".join(updatedFiles.id) +"')"
    queryStr = "select id from %s where id in %s" %(tableId, idSet)
    #Query and delete
    rowSet = syn.tableQuery(queryStr).asRowSet()
    if len(rowSet['rows'])>0:
        print 'DELETING', len(rowSet['rows']), 'ROWS'
        syn.delete(rowSet)


def countAndUpdateTable(input, tableId):
    i, fileMeta = input
    print 'updating table:%s' %tableId, 'with file %s(%s)' %(fileMeta['name'], fileMeta.id), fileMeta['basename']
    ent = syn.get(fileMeta.id)
    if fileMeta.fileType =='bed5':
        data = pd.read_csv(ent.path, sep='\t')
        nFeatures = 0
        samples = list(set(data.Sample.dropna()))
    else: #All other fileTypes
        data = pd.read_csv(ent.path, sep='\t', index_col=0)
        nFeatures, nSamples = data.shape
        samples = data.columns
    metadata = pd.DataFrame([fileMeta]*len(samples))
    metadata['nFeatures'] = nFeatures
    metadata['samples'] = samples
    metadata['patient_barcode'] = [x[:12] for x in metadata.samples]
    metadata.drop(['tissue', u'md5', u'assembly'], axis=1, inplace=True)
    metadata.nFeatures = metadata.nFeatures.astype('int')
    cols = syn.tableQuery('select * from %s limit 1' %args.tableId).asDataFrame().columns

    #Update rows in table
    print 'adding', metadata.shape[0]
    t = syn.store(Table(tableId, metadata[cols]))
    return metadata



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Updates a Synapse Table with '
                                                  'sample sampling information'))
    parser.add_argument('-t', '--table',  dest='tableId', default='syn3281840',
            help='Table where results are stored (e.g. syn3281840) ')
    parser.add_argument('-p', '--project',  dest='projectId', default='syn2812961',
            help='Project (benefactorId) where output files are stored. (e.g. syn2812961)')
    args = parser.parse_args()



    files = synapseHelpers.query2df(syn.chunkedQuery(FILEQUERY % args.projectId), savedSynapseFields=('id', 'name', 'versionNumber'))
    updatedFiles= findUpdates(files, args.tableId)
    print 'NEED TO UPDATE:', updatedFiles.shape[0], 'FILES'
    deleteAffectedRows(updatedFiles, args.tableId)
    dfs = [countAndUpdateTable(row, tableId=args.tableId) for row in updatedFiles.iterrows()]
