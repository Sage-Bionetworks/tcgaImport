#!/usr/bin/env python

import synapseclient
import tcgaImport
from argparse import ArgumentParser
import json
import sys
import logging
import time

#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def getAllSynapseFiles(projectId):
    """Extracts a list of files in Synapse creating a name to id dictionary."""
    results = syn.chunkedQuery("select id, name from entity where benefactorId=='%s'" % projectId)
    return {r['entity.name']:r['entity.id'] for r in results}


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("basename", nargs="?", default=None)
    parser.add_argument("-o", "--output", help="Output File", default=None)
    parser.add_argument("--project", help="Project", default=None)

    args = parser.parse_args()

    if args.basename is not None:
        basename_list = [args.basename]
    else:
        basename_list = []
        for plat in tcgaImport.platform_list():
            logger.info("Queueing: %s" % (plat))
            basename_list += tcgaImport.archive_list(plat)
        basename_list.extend(tcgaImport.clinicnal_archive_list())
        basename_list.extend(tcgaImport.mutation_archive_list())

    handle = open(args.output, "w") if args.output is not None else  sys.stdout
    syn = synapseclient.login()
    allSynapseFiles = getAllSynapseFiles(args.project)

    for basename in basename_list:
        logger.info("Checking: %s" % (basename))
        basename_platform_alias = tcgaImport.get_basename_platform(basename)
        logger.info(basename_platform_alias)
        conf = tcgaImport.getBaseBuildConf(basename, basename_platform_alias, "./")
        request = conf.buildRequest()
        
        for subtypeName, subtypeData in tcgaImport.tcgaConfig[basename_platform_alias].dataSubTypes.items():
            filename = subtypeData['nameGen'](basename)
            entity_id = allSynapseFiles.get(filename)
            if entity_id is None:
                handle.write("MISSING: %s\n" % (basename)) 
                #We have already notified that basename is missing skip rest of subtypes
                break
            else:
                prov = syn.getProvenance(entity_id)
                #Filter out only used not executed items
                prov = [x['url'] for x in prov['used'] if not x.get('wasExecuted', False)]
                req_elems = [x['url'] for x in request['provenance']['used']]
                if set(prov)==set(req_elems):
                    handle.write("READY: %s\n" % (basename)) 
                else:
                    logger.info("Not all used archives Found: %s" %', '.join(set(req_elems) - set(prov)))
                    handle.write("UPDATE: %s\n" % (basename))
                    #We have already notified that basename needs to be updated skip rest of subtypes
                    break
    handle.close()
