import synapseclient
from synapseclient import File
import synapseHelpers
import pandas as pd
import hashlib
from synapseHelpers import query2df, thisCodeInSynapse
#from multiprocessing.dummy  import Pool

QUERY_STR = "select * from file where benefactorId=='syn2812961' and acronym=='PANCAN'"
WHITELISTID = 'syn4551248'

def getFileIdFromName(name):
    q = syn.chunkedQuery("select id from file where name=='%s' and parentId=='%s'" %(name, 'syn4557014'))
    id =  q.next()['file.id']
    return id

def isUptodate(name, files, toRemove, platform):
    try:
        activity = syn.getProvenance(getFileIdFromName(name))
    except StopIteration:
        print 'File not found'
        return False
    used = set([(x['reference']['targetId'], x['reference']['targetVersionNumber']) for x in activity['used'] if x['wasExecuted']==False])
    oldWhitelist = id, version = [i for i in used if i[0]==WHITELISTID][0]
    used = set(['%s.%s' %i for i in used if i[0]!=WHITELISTID])
    currentVersions = set(['%s.%s' % (x.id, x.versionNumber) for x in files])
    # * If upstream data files changed return False
    if currentVersions!=used:
        return False
    #else check if the whitelisting is different for this specific platform
    oldToRemove = getChangeSet(oldWhitelist[1], platform)
    return oldToRemove==toRemove


def getChangeSet(version, platform):
    """Extracts the old whitelist id and version of used and filters the changes down
    to a specific platform."""
    old_whitelist  = syn.get(WHITELISTID, version=version)
    whitelist = pd.read_csv(whitelistEntity.path, sep='\t')
    oldToRemove = set(whitelist.ix[whitelist.Do_not_use & (whitelist.platform==platform), 
                                'aliquot_barcode'])
    return oldToRemove
    

#mp = Pool(8)
syn = synapseclient.login(silent=True)

whitelistEntity = syn.get(WHITELISTID)
whitelist = pd.read_csv(whitelistEntity.path, sep='\t')
inputFiles = synapseHelpers.query2df(syn.chunkedQuery(QUERY_STR))

code=synapseHelpers.thisCodeInSynapse(parentId='syn1774100')
for i, row in inputFiles.iterrows():
    print row.id, row['name'],
    inputFileEntity = syn.get(row.id)
    outFileName = row['name'][:-4]+'_whitelisted'+row['name'][-4:]
    
    toRemove = set(whitelist.ix[whitelist.Do_not_use & (whitelist.platform == row['platform']), 
                                'aliquot_barcode'])

    if isUptodate(outFileName, [inputFileEntity], toRemove, row['platform']):
        print ' is up to date - but update provenance'
        e = syn.get(getFileIdFromName(outFileName), downloadFile=False)
        syn.store(e, used=[inputFileEntity, whitelistEntity], executed=code)
        continue
    if row.fileType =='bed5':  #Do the filtering for bed files
        df = pd.read_csv(inputFileEntity.path, sep='\t')
        print df.shape,
        idx = ~df.Sample.isin(toRemove)
        df = df[idx]
        print '->', df.shape
        df.to_csv('/gluster/home/lomberg/tcgaImport/out/'+outFileName, sep='\t', index=False)
        nFeatures = 0
        nSamples = len(set(df.Sample))
    else: #All other fileTypes
        df = pd.read_csv(inputFileEntity.path, sep='\t', index_col=0)
        print df.shape,
        colsToKeep = [col for col in df.columns if (col.startswith('TCGA') and 
                                                    (col.split('.')[0] not in toRemove) and
                                                    ('.' not in col))]
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
