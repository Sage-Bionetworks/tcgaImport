import synapseclient
import pandas as pd
import hashlib
from synapseHelpers import query2df, thisCodeInSynapse
from multiprocessing.dummy  import Pool

PARENT_ID = 'syn4301332'
PLATFORMS = [('MDA_RPPA_Core', 'RPPA', 'mdanderson.org_PANCAN_MDA_RPPA_Core.RPPA.tsv'),
             ('IlluminaGA_RNASeqV2', 'isoformExp', 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.isoformExp.tsv'),
             ('IlluminaGA_RNASeqV2', 'geneExp', 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.geneExp.tsv'),
             ('IlluminaHiSeq_RNASeqV2', 'isoformExp', 'unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.isoformExp.tsv'),
             ('IlluminaHiSeq_RNASeqV2', 'geneExp', 'unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'),
             ('IlluminaGA_miRNASeq','miRNAExp', 'bcgsc.ca_PANCAN_IlluminaGA_miRNASeq.miRNAExp.tsv'), 
             ('IlluminaHiSeq_miRNASeq', 'miRNAExp', 'bcgsc.ca_PANCAN_IlluminaHiSeq_miRNASeq.miRNAExp.tsv'),
             ('HumanMethylation27','betaValue', 'jhu-usc.edu_PANCAN_HumanMethylation27.betaValue.tsv'),
             ('HumanMethylation450', 'betaValue', 'jhu-usc.edu_PANCAN_HumanMethylation450.betaValue.tsv'),
             #These are bed files
             ('IlluminaHiSeq_DNASeqC', 'cna', 'hms.harvard.edu_PANCAN_IlluminaHiSeq_DNASeqC.cna.bed'),
             ('Genome_Wide_SNP_6', 'cna', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna.bed'),
             ('Genome_Wide_SNP_6', 'cna_nocnv', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_nocnv.bed'), 
             ('Genome_Wide_SNP_6', 'cna_nocnv_probecount', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_nocnv_probecount.bed'),
             ('Genome_Wide_SNP_6', 'cna_probecount', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_probecount.bed')]
             #MSI,
             #Maf


QUERY_STR = "select * from file where benefactorId=='syn2812961'"

def isUptodate(name, files):
    q = syn.chunkedQuery("select id from file where name=='%s' and parentId=='%s'" %(name, PARENT_ID))
    try:
        id =  q.next()['file.id']
    except StopIteration:
        return False
    activity = syn.getProvenance(id)
    used = set(['%s.%s' % (x['reference']['targetId'], x['reference']['targetVersionNumber']) for x in activity['used'] if x['wasExecuted']==False])
    currentVersions = set(['%s.%s' % (x.id, x.versionNumber) for x in files])
    return currentVersions==used




mp = Pool(8)
syn = synapseclient.login(silent=True)
allFiles =  query2df(syn.chunkedQuery(QUERY_STR))

for platform, dataSubType, name in PLATFORMS:
    print platform, dataSubType,
    filteredMeta = allFiles[(allFiles.platform==platform) & (allFiles.dataSubType==dataSubType) & (allFiles.acronym!='PANCAN')]
    files = mp.map(syn.get, filteredMeta.id)
    if isUptodate(name, files):
        print ' is up to date'
        continue
    
    if list(set(filteredMeta.fileType))[0] =='bed5':
        dfs = mp.map(lambda f: pd.read_csv(f.path, sep='\t'), files)
        df = pd.concat(dfs, axis=0)
        df.to_csv('/gluster/home/lomberg/tcgaImport/out/'+name, sep='\t', index=False)
        nSamples = len(set(df.Sample))
        nFeatures = 0
    else: #All other fileTypes
        dfs = mp.map(lambda f: pd.read_csv(f.path, sep='\t', index_col=0), files)
        df = pd.concat(dfs, axis=1)
        df.to_csv('/gluster/home/lomberg/tcgaImport/out/'+name, sep='\t')
        nFeatures, nSamples = df.shape
    print 'Created', name, df.shape
    #Add file to Synapse
    entity = synapseclient.File('/gluster/home/lomberg/tcgaImport/out/'+name, parentId=PARENT_ID)
    #Set annotations
    entity.platform = platform
    entity.dataSubType = dataSubType
    entity.acronym='PANCAN'
    entity.dataProducer='TCGA'
    entity.disease='cancer'
    entity.center =  list(set(filteredMeta.center))
    entity.centerTitle = list(set(filteredMeta.centerTitle))
    entity.fileType = list(set(filteredMeta.fileType))
    entity.platformTitle = list(set(filteredMeta.platformTitle))
    entity.nSamples = nSamples 
    entity.nFeatures = nFeatures 
    entity = syn.store(entity, used = files, executed=thisCodeInSynapse(parentId='syn1774100'))
