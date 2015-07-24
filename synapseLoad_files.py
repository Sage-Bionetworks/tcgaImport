#!/usr/bin/env python
import os
import sys
import json
import re
from glob import glob
from  multiprocessing.dummy import Pool
import synapseclient
from synapseHelpers import query2df
from synapseclient import Activity, File, Folder
import hashlib
import logging
from argparse import ArgumentParser

PROJECTHIERARCHYQUERY =  "select id, name, parentId from file where benefactorId=='%s'"


def find_child(syn, project, name):
    query = "select id from entity where parentId=='%s' and name=='%s'" % (project, name)
    for res in syn.query(query)['results']:
        return res['entity.id']
    return None
    
def getParentFolder(syn, project, meta):
    fid = find_child(syn, project, meta['annotations']['acronym'])
    if fid is None:
        folder = syn.store(Folder(name=meta['annotations']['acronym'], parentId=project))
        fid = folder.id
    pid = find_child(syn, fid, meta['annotations']['platform'])
    if pid is None:
        folder = syn.store(Folder(name=meta['annotations']['platform'], parentId=fid))
        pid = folder.id
    return pid

def loadOneSample(a):
    """Goes through a single json annotation file a and:
        1) Finds the parent Folder where to store the file (or makes directories)
        2) Fetches the md5 of any existing file and compares
        3) If new or different md5 upload file.
    """
    logging.debug( "Loading:" + a )
    with open(a) as handle:
        meta = json.load(handle)
    dpath = re.sub(r'.json$', '', a)
    #Skip the rest of the loop if data file is empty or we are not doing the current acronyms
    if os.stat(dpath).st_size==0 or (args.acronym != meta['annotations']['acronym'] and args.acronym is not None):
        return 

    parentId= getParentFolder(syn, args.project, meta)
    #Determine if we are updating an existing file and if we should update based on md5
    query = "select id from entity where parentId=='%s' and name=='%s'" % (parentId, meta['name'])
    res = list(syn.chunkedQuery(query))
    if len(res) != 0:
        tmp_ent = syn.get(res[0]['entity.id'], downloadFile=False)
        upload = (tmp_ent.md5 != meta['annotations']['md5'])
        logging.debug( "\tFound: %s and upload (MD5 %s match)" %(tmp_ent.id, 'DOESN\'T' if upload else 'does'))
    else:
        logging.debug("\tNot found:" + meta['name'])
        upload = True
    #Prepare the entity for upload
    if upload and not args.push:
        logging.info( "\tWILL UPLOAD: %s" %meta['name'])
    if upload and args.push: 
        entity = File(dpath, name=meta['name'], parentId=parentId, annotations=meta['annotations'])
        if 'provenance' in meta:
            #Fix labels for urls
            for u in meta['provenance']['used']:
                if 'name' not in u and 'url' in u:
                    u['name'] = u['url']
            prov = Activity(data=meta['provenance'])
            prov.executed('https://github.com/Sage-Bionetworks/tcgaImport')

        else:
            prov=None
        logging.debug('\tUploading:%s' %entity.name)
        entity = syn.store(entity, activity=prov)
        logging.debug('\tCreated/Updated: **** %s ****' %entity.id)



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("src", help="Scan directory", default=None)
    parser.add_argument("--project", help="Project", default=None)
    parser.add_argument("--skip-md5", help="Skip MD5", action="store_true", default=False)
    parser.add_argument("--push", help="Push", action="store_true", default=False)
    parser.add_argument("--acronym", help="Limit to one Acronym", default=None)

    args = parser.parse_args()
    syn = synapseclient.login()

    logging.basicConfig(level=logging.INFO)

    #Get a list of the files and folders in Synapse (will be updated when folder added)
    #hierarchy = query2df(syn.query(PROJECTHIERARCHYQUERY %args.project))

    p=Pool(5)
    p.map(loadOneSample, glob(os.path.join( args.src, "*.json")))
        

