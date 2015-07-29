import synapseclient
from synapseclient import File
import synapseHelpers
import pandas as pd
import hashlib
from synapseHelpers import query2df, thisCodeInSynapse
#from multiprocessing.dummy  import Pool

QUERY_STR = "select * from file where benefactorId=='syn2812961' and acronym=='PANCAN'"
WHITELISTID = 'syn4551248'

def isUptodate(name, files):
    q = syn.chunkedQuery("select id from file where name=='%s' and parentId=='%s'" %(name, 'syn4557014'))
    try:
        id =  q.next()['file.id']
    except StopIteration:
        print 'not found',
        return False
    activity = syn.getProvenance(id)
    used = set(['%s.%s' % (x['reference']['targetId'], x['reference']['targetVersionNumber']) for x in activity['used'] if x['wasExecuted']==False])
    currentVersions = set(['%s.%s' % (x.id, x.versionNumber) for x in files])
    return currentVersions==used

#mp = Pool(8)
syn = synapseclient.login(silent=True)

whitelistEntity = syn.get(WHITELISTID)
whitelist = pd.read_csv(whitelistEntity.path, sep='\t')
inputFiles = synapseHelpers.query2df(syn.chunkedQuery(QUERY_STR))

code=synapseHelpers.thisCodeInSynapse(parentId='syn1774100')
for i, row in inputFiles.iterrows():
    print row.id, row['name'] 
    inputFileEntity = syn.get(row.id)
    outFileName = row['name'][:-4]+'_whitelisted'+row['name'][-4:]
    
    # Get the platform type for the file from the filename
    platform = row['platform']

    toRemove = whitelist.ix[whitelist.Do_not_use & 
        (whitelist.platform == platform), 'aliquot_barcode'].tolist()

    if isUptodate(outFileName, [whitelistEntity, inputFileEntity]):
        print ' is up to date'
        continue
    if row.fileType =='bed5':  #Do the filtering for bed files
        df = pd.read_csv(inputFileEntity.path, sep='\t', header=None)
        print df.shape,
        idx = ~df[3].isin(toRemove)
        df = df[idx]
        print '->', df.shape
        df.to_csv('/gluster/home/lomberg/tcgaImport/out/'+outFileName, sep='\t', header=False, index=False)
        nFeatures = 0
        nSamples = len(set(df[3]))
    else: #All other fileTypes
        df = pd.read_csv(inputFileEntity.path, sep='\t', index_col=0)
        print df.shape,
        colsToKeep = [col for col in df.columns if col.startswith('TCGA') and col.split('.')[0] not in toRemove]
        df = df.ix[:, colsToKeep]
        print '->', df.shape
        df.to_csv('/gluster/home/lomberg/tcgaImport/out/'+outFileName, sep='\t')
        nFeatures , nSamples = df.shape

    annots = syn.getAnnotations(row.id)
    del annots['etag']
    del annots['uri']
    del annots['id']
    annots['nSamples'] = nSamples
    annots['nFeatures'] = nFeatures 

    f = File('/gluster/home/lomberg/tcgaImport/out/'+outFileName,  parent='syn4557014', annotations=annots)
    f = syn.store(f, used=[inputFileEntity, whitelistEntity], executed=code)
