import synapseclient
import pandas as pd
import hashlib
from synapseHelpers import query2df, thisCodeInSynapse
from multiprocessing.dummy  import Pool
import argparse

PLATFORMS = [('MDA_RPPA_Core', 'RPPA', 'mdanderson.org_PANCAN_MDA_RPPA_Core.RPPA.tsv'),
             ('IlluminaGA_RNASeqV2', 'isoformExp', 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.isoformExp.tsv'),
             ('IlluminaGA_RNASeqV2', 'geneExp', 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.geneExp.tsv'),
             ('IlluminaHiSeq_RNASeqV2', 'isoformExp', 'unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.isoformExp.tsv'),
             ('IlluminaHiSeq_RNASeqV2', 'geneExp', 'unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'),
             ('IlluminaGA_miRNASeq','miRNAExp', 'bcgsc.ca_PANCAN_IlluminaGA_miRNASeq.miRNAExp.tsv'), 
             ('IlluminaHiSeq_miRNASeq', 'miRNAExp', 'bcgsc.ca_PANCAN_IlluminaHiSeq_miRNASeq.miRNAExp.tsv'),
             ('HumanMethylation27','betaValue', 'jhu-usc.edu_PANCAN_HumanMethylation27.betaValue.tsv'),
             ('HumanMethylation450', 'betaValue', 'jhu-usc.edu_PANCAN_HumanMethylation450.betaValue.tsv'),
             #These are bed and seg files
             ('IlluminaHiSeq_DNASeqC', 'cna', 'hms.harvard.edu_PANCAN_IlluminaHiSeq_DNASeqC.cna.bed'),
             ('Genome_Wide_SNP_6', 'cna', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna.seg'),
             ('Genome_Wide_SNP_6', 'cna_nocnv', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_nocnv.seg'), 
             ('Genome_Wide_SNP_6', 'cna_nocnv_probecount', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_nocnv_probecount.seg'),
             ('Genome_Wide_SNP_6', 'cna_probecount', 'broad.mit.edu_PANCAN_Genome_Wide_SNP_6.hg19.cna_probecount.seg')]
             #MSI,
             #Maf

# Generate string of unique platforms from PLATFORMS array.
availPlatforms = '\n'.join(set([x[0] for x in PLATFORMS]))

# Argument parser to allow user to indicate with synapse project to merge files from, which project to load the merged file into, 
# and an optional platform argument if the user only has files of a single platform type. 
parser = argparse.ArgumentParser()
parser.add_argument('benefactorId',help='ID of synapse project to merge files from.')
parser.add_argument('parentId',help='ID of synapse project to add merged file to.')
parser.add_argument('-p','--platform',help='If merging files of single platform type, add platform name after option. Available platforms' + '\n' +  availPlatforms,type=str)
args = parser.parse_args()

if args.platform:
	PLATFORMS = [x for x in PLATFORMS if x[0] == args.platform]	

QUERY_STR = "select * from file where benefactorId==" + ("'{0}'").format(args.benefactorId)

def isUptodate(name, files):
    q = syn.chunkedQuery("select id from file where name=='%s' and parentId=='%s'" %(name, args.parentId))
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

    if list(set(filteredMeta.fileType))[0] =='seg' or list(set(filteredMeta.fileType))[0] == 'bed5':
        dfs = mp.map(lambda f: pd.read_csv(f.path, sep='\t'), files)
        df = pd.concat(dfs, axis=0)
        df.to_csv('/Users/duncan/sageBio/tcgaImport/out/'+name, sep='\t', index=False)
        nSamples = len(set(df.Sample))
        nFeatures = 0
    else: #All other fileTypes
        dfs = mp.map(lambda f: pd.read_csv(f.path, sep='\t', index_col=0), files)
        df = pd.concat(dfs, axis=1)
        df.to_csv('/Users/duncan/sageBio/tcgaImport/out/'+name, sep='\t')
        nFeatures, nSamples = df.shape

    print 'Created', name, df.shape

    #Add file to Synapse
    entity = synapseclient.File('/Users/duncan/sageBio/tcgaImport/out/'+name, parentId=args.parentId)

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
    entity = syn.store(entity, used = files, executed=thisCodeInSynapse(parentId='syn5016984'))
