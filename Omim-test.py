import json
import urllib
from urllib import request
from urllib import parse
import neo4j
from neo4j import GraphDatabase, ServiceUnavailable, basic_auth
import datetime
from datetime import date
import time
import logging, logging.config
import itertools
from itertools import islice
import pandas
import jmespath
import xml.etree.ElementTree as ET
import re
  

def find_OMIM_articles(OMIMNumber):
        """fetch articles and return a map"""        
        params = urllib.parse.urlencode({'mimNumber': OMIMNumber, 'include':"all", 'format': 'json', 'apiKey':""})
        data = request.urlopen("https://api.omim.org/api/entry?%s" % params).read().decode()
        return json.loads(data)

omim_reference = find_OMIM_articles("300100")
import jmespath
pmids = jmespath.search("omim.entryList[0].entry.referenceList[*].reference.pubmedID",omim_reference)
refNumbers = jmespath.search("omim.entryList[0].entry.referenceList[*].reference.[referenceNumber,pubmedID]",omim_reference)
textSections = jmespath.search("omim.entryList[0].entry.textSectionList[*].textSection",omim_reference)

references = {}
for t in textSections:
        refs = re.findall("({.*?}|{.*?:.*?})",t['textSectionContent'])
        refs
        if (refs is not None):
                #sectionReferenced = [(str.split(ref,":")[0])[1:] for ref in refs if (ref is not None and len(str.split(ref,":"))>1) and not re.search(",",str.split(ref,":")[0])]
                sectionReferenced = []
                for ref in refs:
                        splitRef= ref[1:].split(":")
                        if (len(splitRef)>1):
                                if (not re.search(",",splitRef[0])):
                                        sectionReferenced.append(splitRef[0])
                                else:
                                        multiples =  splitRef[0].split(",")
                                        for r in multiples:
                                                sectionReferenced.append(r)

        references[t['textSectionName']] = list(set(sectionReferenced))

articleString = {}
for refNumber,pmid in refNumbers:
  if (pmid is not None):
    #print (refNumber, pmid)
    tsections = ['OMIM:reference']
    for idx, sectionName in enumerate(references):
      if len(set([str(refNumber)]) & set(references[sectionName]))>0:
        tsections.append("OMIM:" + sectionName)
    articleString[pmid] = tsections
print (articleString)

