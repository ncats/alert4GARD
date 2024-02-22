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
import jmespath
import re


def get_GARD_diseases_list():
  """Returns list of the GARD Diseases"""
  with GraphDatabase.driver("bolt://disease.ncats.io:80") as driver:
    try:
      with driver.session() as session:
        cypher_query = '''
        match p = (d:DATA)-[]-(m:S_GARD) return distinct d.gard_id, d.name
        '''
        results = session.run(cypher_query, parameters={})
        myData = {}
        for record in results:
          myData[record['d.gard_id']] = record['d.name']
              
      return myData
    except ServiceUnavailable as e:
            print(e)

def find_mondo_articles(MondoID):
  """fetch mondo articles and return a map"""
  # https://api.monarchinitiative.org/api/bioentity/disease/MONDO:0007846/publications

  data = request.urlopen("https://api.monarchinitiative.org/api/bioentity/disease/" + MondoID + "/publications").read().decode()

  return json.loads(data)


def get_Mondo_and_GARD_diseases_list():
  """Returns list of the Mondo and GARD Diseases"""

  with GraphDatabase.driver("bolt://disease.ncats.io:80") as driver:
    try:
      with driver.session() as session:
        cypher_query = '''
        match p = (d:DATA)-[:PAYLOAD]->(m:S_MONDO)<-[r:R_equivalentClass]-(n:S_GARD)<-[:PAYLOAD]-(d1:DATA) 
        RETURN d.notation as MONDO_ID, d1.gard_id as GARD_ID
        '''
        results = session.run(cypher_query, parameters={})
        myData = {}
        for record in results:
          myData[record['MONDO_ID']] = record['GARD_ID']

      return myData
    except ServiceUnavailable as e:
        print(e)


def get_OMIM_and_GARD_diseases_list():
  """Returns list of the OMIM and GARD Diseases"""

  with GraphDatabase.driver("bolt://disease.ncats.io:80") as driver:
    try:
      with driver.session() as session:
          cypher_query = '''
          match (d:DATA)-[]-(m:S_GARD)
          where exists(d.gard_id)
          unwind d.xrefs as x
          return distinct d.gard_id, d.name, x
          '''
          results = session.run(cypher_query, parameters={})
          myData = {}
          for record in results:
            if (record['x'].startswith("OMIM")):
              #myData[record['x']] = record['d.gard_id']
              myData[record['d.gard_id']] = record['x']
      return myData
    except ServiceUnavailable as e:
      print(e)

def find_articles(keyword, mindate, maxdate):
  """fetch articles and return a map"""
  params = urllib.parse.urlencode({'db': 'pubmed', 'term': keyword, 'mindate':mindate, 'maxdate':maxdate, 'retmode': 'json', 'retmax':"1000"})
  data = request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?%s" % params).read().decode()
  return json.loads(data)

def fetch_abstract(pubmedID): 
  """fetch abstract for an article"""
  ebiUrl = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:" + str(pubmedID) + "&resulttype=core&format=json"
  ebiData = request.urlopen(ebiUrl).read().decode()
  return json.loads(ebiData)

def fetch_pubtator_annotations(pubmedId):
  """fetch annotations from pubtator"""
  pubtatorUrl = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmids=" + str(pubmedId)
  pubData = request.urlopen(pubtatorUrl).read().decode()
  return json.loads(pubData)

def fetch_pmc_fulltext_xml(pmcId):
  """fetch full text xml from pmc"""
  pmcUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + str(pmcId)
  pmcData = request.urlopen(pmcUrl).read().decode()
  return json.loads(pmcData)

def fetch_pmc_fulltext_json(pubmedId):
  """fetch full text json from pmc"""
  #pmcUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + pmcId
  pmcUrl = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/" + str(pubmedId) + "/unicode" 
  pmcData = request.urlopen(pmcUrl).read().decode()
  return json.loads(pmcData)

def find_OMIM_articles(OMIMNumber):
        """fetch articles and return a map"""        
        params = urllib.parse.urlencode({'mimNumber': OMIMNumber, 'include':"all", 'format': 'json', 'apiKey':"wr20AKemSkWnoFR8b-vSkA"})
        data = request.urlopen("https://api.omim.org/api/entry?%s" % params).read().decode()
        return json.loads(data)

def get_omim_refs(omim_reference):
  """return a dictionary of pubmed article ids and a list of Omim text sections where that article was referenced"""
  """{8415729: ['OMIM:reference', 'OMIM:mapping'], 9537923: ['OMIM:reference'],...}"""

  pmids = jmespath.search("omim.entryList[0].entry.referenceList[*].reference.pubmedID",omim_reference)
  refNumbers = jmespath.search("omim.entryList[0].entry.referenceList[*].reference.[referenceNumber,pubmedID]",omim_reference)
  textSections = jmespath.search("omim.entryList[0].entry.textSectionList[*].textSection",omim_reference)

  references = {}
  for t in textSections:
    refs = re.findall("({.*?}|{.*?:.*?})",t['textSectionContent'])
    refs
    if (refs is not None):
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
      tsections = ['OMIM:reference']
      for idx, sectionName in enumerate(references):
        if len(set([str(refNumber)]) & set(references[sectionName]))>0:
          tsections.append("OMIM:" + sectionName)
      articleString[pmid] = tsections
  return articleString

def create_gard_disease(session, gard_id, name):
  """Create or merge a GARD Disease node with ID and Name"""
  
  return session.run("MERGE (d:Disease {gard_id:$gard_id, name:$name}) RETURN id(d)", 
    gard_id=gard_id, name=name).single().value()

def get_gard_disease_node(session, gard_id):
  """Fetch a GARD Disease node with ID and Name"""
  
  return session.run("MATCH (d:Disease {gard_id:$gard_id}) RETURN id(d)", 
    gard_id=gard_id).single().value()


def create_article(tx, abstractDataRel, disease_node, search_source):
  """Create the node for an article"""

  create_article_query = '''
  MATCH (d:Disease) WHERE id(d)=$id
  MERGE (n:Article {pubmed_id:$pubmed_id, doi:$doi, title:$title, abstractText:$abstractText, 
          authorString:$authorString, affiliation:$affiliation, firstPublicationDate:$firstPublicationDate, citedByCount:$citedByCount,
          isOpenAccess:$isOpenAccess, inEPMC:$inEPMC, inPMC:$inPMC, hasPDF:$hasPDF})
  MERGE (d)-[r:MENTIONED_IN]->(n)
  MERGE (n)-[a:APPEARS_IN]->(s:Source {source:$source})
  ON CREATE SET n.created=datetime()
  ON MATCH SET n.updated=datetime()
  RETURN id(n)
  '''
  params={
    "id":disease_node,
    "pubmed_id":abstractDataRel['pmid'] if 'pmid' in abstractDataRel else '',
    "source":abstractDataRel['source'] if 'source' in abstractDataRel else '',
    "doi":abstractDataRel['doi'] if 'doi' in abstractDataRel else '',
    "title":abstractDataRel['title'] if 'title' in abstractDataRel else '',
    "abstractText":abstractDataRel['abstractText'] if 'abstractText' in abstractDataRel else '',
    "authorString":abstractDataRel['authorString'] if 'authorString' in abstractDataRel else '',
    "affiliation":abstractDataRel['affiliation'] if 'affiliation' in abstractDataRel else '',
    "firstPublicationDate":abstractDataRel['firstPublicationDate'] if 'firstPublicationDate' in abstractDataRel else '',
    "isOpenAccess":abstractDataRel['isOpenAccess'] if 'isOpenAccess' in abstractDataRel else '',
    "inEPMC":abstractDataRel['inEPMC'] if 'inEPMC' in abstractDataRel else '',
    "inPMC":abstractDataRel['inPMC'] if 'inPMC' in abstractDataRel else '',
    "hasPDF":abstractDataRel['hasPDF'] if 'hasPDF' in abstractDataRel else '',
    "citedByCount":abstractDataRel['citedByCount'] if 'citedByCount' in abstractDataRel else ''
    }
  return tx.run(create_article_query,**params ).single().value()

def update_article_references(session, pubmedID, search_source):
  """update the article references"""

  match_article_query = '''MATCH (a:Article) WHERE pubmed_id = $pubmed_id RETURN  a.references'''
  params={
    "pubmed_id":pubmedID
    }
  txout = session.run(match_article_query,**params )
  if (txout):
    references = txout.value()
    updatedRef = set(references) + set(search_source)
    update_article_query = ''''MATCH(a:Article {pubmed_id:$pubmed_id}) 
      SET a.references=$references
      return a'''
    params={
      "pubmed_id":pubmedID,
      "references":updatedRef
      }
    session.run(update_article_query,**params )



def create_authors(tx, abstractDataRel, article_node, disease_node):
  create_author_query = '''
  MATCH (a:Article), (d:Disease) WHERE id(a) = $article_id and id(d) = $disease_node
  MERGE (p:Person {fullName:$fullName, firstName:$firstName, lastName:$lastName, affiliation:$affiliation})
  MERGE (p) - [r:WROTE] -> (a)
  ON CREATE SET p.created=datetime()
  ON MATCH SET p.updated=datetime()
  MERGE (p)-[s:STUDIES]->(d)
  '''
  for author in abstractDataRel['authorList']['author']:
    tx.run(create_author_query, parameters={
      "article_id":article_node,
      "disease_node":disease_node,
      "fullName": author['fullName'] if 'fullName' in author else '',
      "firstName": author['firstName'] if 'firstName' in author else '',
      "lastName": author['lastName'] if 'lastName' in author else '',
      "affiliation": author['affiliation'] if 'affiliation' in author else '',
    })

def create_journal(tx, abstractDataRel, article_node):
  create_journal_query = '''
  MATCH (a:Article) WHERE id(a) = $article_id
  MERGE (j:Journal{title:$title,medlineAbbreviation:$medlineAbbreviation,essn:$essn,issn:$issn,nlmid:$nlmid})
  MERGE (ji:JournalInfo{issue:$issue, volume:$volume, journalIssueId:$journalIssueId,
          dateOfPublication:$dateOfPublication, monthOfPublication:$monthOfPublication,yearOfPublication:$yearOfPublication,
          printPublicationDate:$printPublicationDate})
  MERGE (a)-[:APPEARS_IN]->(ji)
  MERGE (ji)-[:CONTENT_OF]->(j)
  ON CREATE SET j.created=datetime()
  '''
  tx.run(create_journal_query, parameters={
    "article_id":article_node,
    "issue": abstractDataRel['issue'] if 'issue' in abstractDataRel else '',
    "volume": abstractDataRel['volume'] if 'volume' in abstractDataRel else '',
    "journalIssueId": abstractDataRel['journalIssueId'] if 'journalIssueId' in abstractDataRel else '',
    "dateOfPublication": abstractDataRel['dateOfPublication'] if 'dateOfPublication' in abstractDataRel else '',
    "monthOfPublication": abstractDataRel['monthOfPublication'] if 'monthOfPublication' in abstractDataRel else '',
    "yearOfPublication": abstractDataRel['yearOfPublication'] if 'yearOfPublication' in abstractDataRel else '',
    "printPublicationDate": abstractDataRel['printPublicationDate'] if 'printPublicationDate' in abstractDataRel else '',
    "title": abstractDataRel['journal']['title'] if 'title' in abstractDataRel['journal'] else '',
    "medlineAbbreviation": abstractDataRel['journal']['medlineAbbreviation'] if 'medlineAbbreviation' in abstractDataRel['journal'] else '',
    "essn": abstractDataRel['journal']['essn'] if 'essn' in abstractDataRel['journal'] else '',
    "issn": abstractDataRel['journal']['issn'] if 'issn' in abstractDataRel['journal'] else '',
    "nlmid": abstractDataRel['journal']['nlmid'] if 'nlmid' in abstractDataRel['journal'] else ''
  })

def create_keywords(tx, abstractDataRel, article_node):
  create_keyword_query = '''
  MATCH (a:Article), (d:Disease) WHERE id(a) = $article_id
  MERGE (k:Keyword {keyword:$keyword}) 
  MERGE (k)- [r:for] -> (a)
  ON CREATE SET k.created=datetime()
  '''
  for keyword in abstractDataRel:
    if keyword:
      tx.run(create_keyword_query, parameters={
        "article_id":article_node,                        
        "keyword": keyword
      })

def create_fullTextUrls(tx, abstractDataRel, article_node):
  create_fullTextUrls_query = '''
  MATCH (a:Article) WHERE id(a) = $article_id
  MERGE (u:fullTextUrl {availability:$availability, availabilityCode:$availabilityCode, documentStyle:$documentStyle,site:$site,url:$url})
  MERGE (u) - [r:CONTENT_FOR] -> (a)
  ON CREATE SET u.created=datetime()
  ON MATCH SET u.updated=datetime()
  '''
  for fullTextUrl in abstractDataRel:
    tx.run(create_fullTextUrls_query, parameters={
      "article_id":article_node,
      "availability": fullTextUrl['availability'] if 'availability' in fullTextUrl else '',
      "availabilityCode": fullTextUrl['availabilityCode'] if 'availabilityCode' in fullTextUrl else '',
      "documentStyle": fullTextUrl['documentStyle'] if 'documentStyle' in fullTextUrl else '',
      "site": fullTextUrl['site'] if 'site' in fullTextUrl else '',
      "url": fullTextUrl['url'] if 'url' in fullTextUrl else ''
    })

def create_meshHeadings(tx, abstractDataRel, article_node):
  create_meshHeadings_query = '''
  MATCH (a:Article) WHERE id(a) = $article_id
  MERGE (m:MeshTerm {majorTopic_YN:$majorTopic_YN, descriptorName:$descriptorName}) 
  MERGE (m) - [r:MESH_TERM_FOR] -> (a)
  ON CREATE SET m.created=datetime()
  ON MATCH SET m.updated=datetime()
  RETURN id(m)
  '''
  create_meshQualifiers_query = '''
  MATCH (m:MeshTerm) WHERE id(m) = $meshHeading_id
  MERGE (mq:MeshQualifier {abbreviation:$abbreviation, qualifierName:$qualifierName,majorTopic_YN:$majorTopic_YN}) 
  MERGE (mq) - [r:MESH_QUALIFIER_FOR] -> (m)
  ON CREATE SET mq.created=datetime()
  '''
  for meshHeading in abstractDataRel:
    parameters={
      "article_id":article_node,
      "majorTopic_YN": meshHeading['majorTopic_YN'] if 'majorTopic_YN' in meshHeading else '',
      "descriptorName": meshHeading['descriptorName'] if 'descriptorName' in meshHeading else ''
    }
    txout = tx.run(create_meshHeadings_query,**parameters).single()
    if (txout):
      meshHeadingId = txout.value()
      if ('meshQualifierList' in meshHeading and 'meshQualifier' in meshHeading['meshQualifierList']):
        for meshQualifier in meshHeading['meshQualifierList']['meshQualifier']:
          tx.run(create_meshQualifiers_query, parameters={
          "meshHeading_id":meshHeadingId,
          "abbreviation": meshQualifier['abbreviation'] if 'abbreviation' in meshQualifier else '',
          "qualifierName": meshQualifier['qualifierName'] if 'qualifierName' in meshQualifier else '',
          "majorTopic_YN": meshQualifier['majorTopic_YN'] if 'majorTopic_YN' in meshQualifier else '',
          })


def create_chemicals(tx, abstractDataRel, article_node):
  create_chemicals_query = '''
  MATCH (a:Article) WHERE id(a) = $article_id
  MERGE (u:Chemical {name:$name, registryNumber:$registryNumber}) - [r:CHEMICAL_FOR] -> (a)
  ON CREATE SET u.created=datetime()
  ON MATCH SET u.updated=datetime()
  '''
  for chemical in abstractDataRel:
    tx.run(create_chemicals_query, parameters={
      "article_id":article_node,
      "name": chemical['name'] if 'name' in chemical else '',
      "registryNumber": chemical['registryNumber'] if 'registryNumber' in chemical else '',
    })

def create_annotations(tx, pubtatorData, article_node):
  create_passages_query = '''
  MATCH(a:Article) WHERE id(a) = $article_id
  MERGE (p:PubtatorPassage {type:$type, text:$text}) - [r:passage_for] -> (a)
  ON CREATE SET p.created=datetime()
  ON MATCH SET p.updated=datetime()
  RETURN id(p)
  '''
  create_annotations_query = '''
  MATCH(pp:PubtatorPassage) WHERE id(pp) = $pubtator_passage_id
  MERGE (pa:PubtatorAnnotations {infons_identifier:$infons_identifier, infons_type:$infons_type, text:$text, locations_offset:$locations_offset, locations_length:$locations_length}) - [r:annotations_for] -> (pp)
  ON CREATE SET pa.created=datetime()
  ON MATCH SET pa.updated=datetime()
  '''
  for passage in pubtatorData['passages']:
    txout = tx.run(create_passages_query, parameters={
      "article_id":article_node,
      "type": passage['infons']['type'] if 'type' in passage['infons'] else '',
      "text": passage['text'] if 'text' in passage else '',
    }).single()
    passageId = txout.value()
    for annotation in passage['annotations']:
      parameters={
        "pubtator_passage_id":passageId,
        "infons_identifier": annotation['infons']['identifier'] if ('identifier' in annotation['infons'] and annotation['infons']['identifier'])  else '',
        "infons_type": annotation['infons']['type'] if ('type' in annotation['infons'] and annotation['infons']['type']) else '',
        "text": annotation['text'] if 'text' in annotation else '',
        "locations_offset": annotation['locations'][0]['offset'] if 'offset' in annotation['locations'][0] else '',
        "locations_length": annotation['locations'][0]['length'] if 'length' in annotation['locations'][0] else '',
      }
      txout = tx.run(create_annotations_query, **parameters)

def save_all_data(abstractData, disease_node, pubmedID, search_source, session):
  if ('resultList' in abstractData and 'result' in abstractData['resultList'] and len(abstractData['resultList']['result']) > 0):
    
    tx = session.begin_transaction()

    logging.info(f'Invoking create_article')
    article_node = create_article(tx, abstractData['resultList']['result'][0], disease_node, search_source)
    
    if ('meshHeadingList' in abstractData['resultList']['result'][0] and 'meshHeading' in abstractData['resultList']['result'][0]['meshHeadingList']):
      logging.info(f'Invoking create_meshHeading')
      create_meshHeadings(tx, abstractData['resultList']['result'][0]['meshHeadingList']['meshHeading'], article_node)

    if ('authorList' in abstractData['resultList']['result'][0] and 'author' in abstractData['resultList']['result'][0]['authorList']):
      logging.info(f'Invoking create_authors')
      create_authors(tx, abstractData['resultList']['result'][0], article_node, disease_node)
    
    if ('journalInfo' in abstractData['resultList']['result'][0]):
      logging.info(f'Invoking create_journal')
      create_journal(tx, abstractData['resultList']['result'][0]['journalInfo'], article_node)

    if ('keywordList' in abstractData['resultList']['result'][0] and 'keyword' in abstractData['resultList']['result'][0]['keywordList']):
      logging.info(f'Invoking create_keywords')
      create_keywords(tx, abstractData['resultList']['result'][0]['keywordList']['keyword'], article_node)
    
    if ('fullTextUrlList' in abstractData['resultList']['result'][0] and 'fullTextUrl' in abstractData['resultList']['result'][0]['fullTextUrlList']):
      logging.info(f'Invoking create_fullTextUrls')
      create_fullTextUrls(tx, abstractData['resultList']['result'][0]['fullTextUrlList']['fullTextUrl'], article_node)

    if ('chemicalList' in abstractData['resultList']['result'][0] and 'chemical' in abstractData['resultList']['result'][0]['chemicalList']):
      logging.info(f'Invoking create_chemical')
      create_chemicals(tx, abstractData['resultList']['result'][0]['chemicalList']['chemical'], article_node)

    tx.commit()

    #begin another transaction save article annotations
    tx = session.begin_transaction()
    logging.info(f'Invoking create annotations')
    try:
      create_annotations(tx, fetch_pubtator_annotations(pubmedID),article_node)
    except Exception as e:
      logging.error(f'Exception creating annotations for article {pubmedID}:  {e}')
    finally:
      tx.commit()

    #begin another transaction save article full text
    # tx = session.begin_transaction()
    # logging.info(f'Invoking create full text')
    # try:
    #         create_fulltext(tx, fetch_pmc_fulltext(pubmedID),article_node)
    # except Exception as e:
    #         logging.error(f'Exception creating full text for article {pubmedID}:  {e}')
    # finally:
    #         tx.commit()

def create_fulltext(tx, pmcData, article_node):
  create_passages_query = '''
  MATCH(a:Article) WHERE id(a) = $article_id
  MERGE (p:PMCText {type:$type, title:$title, text:$text}) - [r:passage_for] -> (a)
  ON CREATE SET p.created=datetime()
  ON MATCH SET p.updated=datetime()
  RETURN id(p)
  '''
  create_annotations_query = '''
  MATCH(pp:PubtatorPassage) WHERE id(pp) = $pubtator_passage_id
  MERGE (pa:PubtatorAnnotations {infons_identifier:$infons_identifier, infons_type:$infons_type, text:$text, locations_offset:$locations_offset, locations_length:$locations_length}) - [r:annotations_for] -> (pp)
  ON CREATE SET pa.created=datetime()
  ON MATCH SET pa.updated=datetime()
  '''
  for passage in pmcData['passages']:
    txout = tx.run(create_passages_query, parameters={
      "article_id":article_node,
      "type": passage['infons']['type'] if 'type' in passage['infons'] else '',
      "text": passage['text'] if 'text' in passage else '',
    }).single()
    passageId = txout.value()
    for annotation in passage['annotations']:
      parameters={
        "pubtator_passage_id":passageId,
        "infons_identifier": annotation['infons']['identifier'] if ('identifier' in annotation['infons'] and annotation['infons']['identifier'])  else '',
        "infons_type": annotation['infons']['type'] if ('type' in annotation['infons'] and annotation['infons']['type']) else '',
        "text": annotation['text'] if 'text' in annotation else '',
        "locations_offset": annotation['locations'][0]['offset'] if 'offset' in annotation['locations'][0] else '',
        "locations_length": annotation['locations'][0]['length'] if 'length' in annotation['locations'][0] else '',
      }
      txout = tx.run(create_annotations_query, **parameters)

def save_articles_pubmed(mindate, maxdate):
  """Find and save article data for GARD diseases between from mindate to maxdate"""

  with GraphDatabase.driver("bolt://3.128.119.166:7687",auth=("neo4j", "back up Generic Plains 27")) as driver:
    try:
      with driver.session() as session:
        results = get_GARD_diseases_list()
        results = dict(itertools.islice(results.items(),1000,3000))
        for idx, gard_id in enumerate(results):
          time.sleep(0.34)
          logging.warning(f'{idx}: Invoking find_articles({results[gard_id]},{mindate},{maxdate})')
          try:
            pubmedIDs = find_articles(results[gard_id],mindate,maxdate)
          except urllib.request.URLError as e:
            logging.error(f'Exception when finding articles: {e}')
            continue

          logging.debug('Invoking create_gard_disease')
          disease_node = create_gard_disease(session, gard_id, results[gard_id])
          
          for pubmedID in pubmedIDs['esearchresult']['idlist']:
            try:
              res = session.run("match(a:Article{pubmed_id:$pmid}) set a.pubmed_evidence=TRUE return a", parameters={"pmid":pubmedID})
              matching_articles = len(list(res))
              if (matching_articles):
                logging.warning(f'pubmed evidence article already exists: {pubmedID}')
              else:
                logging.warning(f'Fetch pubmed abstract({pubmedID})')
                abstractData = fetch_abstract(pubmedID)
                search_source = 'pubmed_evidence'
                save_all_data(abstractData, disease_node, pubmedID, search_source, session)
            except urllib.request.URLError as e:
              logging.error(f'Exception when finding articles: {e}')
            continue
    except ServiceUnavailable as e:
        logging.error(f'Exception: {e}')
        pass
                  
def save_articles_mondo():
  """Find and save article data from mondo for GARD diseases"""

  with GraphDatabase.driver("bolt://3.128.119.166:7687",auth=("neo4j", "back up Generic Plains 27")) as driver:
    try:
      with driver.session() as session:          
        #check for mondo data
        mondo_gard_list = get_Mondo_and_GARD_diseases_list()
        for idx, mondo_id in enumerate(mondo_gard_list):
          time.sleep(0.34)
          #logging.warning(f'{idx}: Invoking find_articles({results[gard_id]},{mindate},{maxdate})')
          try:
            pubmedIDs = find_mondo_articles(mondo_id)
            if ('associations' in pubmedIDs and len(pubmedIDs['associations']) > 0):
              for association in pubmedIDs['associations']:
                if('publications' in association and len(association['publications']) > 0 
                  and association['publications'][0] and association['publications'][0]['id']):
                  pubmedid = association['publications'][0]['id']
                  pubmedid = str.split(pubmedid,":")[1]
                  res = session.run("match(a:Article{pubmed_id:$pmid}) set a.mondo_evidence=TRUE return a", parameters={"pmid":pubmedid})
                  matching_articles = len(list(res))
                  if (matching_articles):
                    logging.warning(f'Mondo evidence article already exists: {pubmedid}')
                  else:
                    logging.warning(f'Adding Mondo evidence article : {pubmedid}')
                    logging.warning(f'Fetch mondo abstract({pubmedid})')
                    abstractData = fetch_abstract(pubmedid)
                    search_source = 'mondo_evidence'
                    save_all_data(abstractData, disease_node, pubmedid, search_source, session)
          except urllib.request.URLError as e:
            logging.error(f'Exception when finding articles: {e}')
          continue
    except ServiceUnavailable as e:
        logging.error(f'Exception: {e}')
        pass

def save_articles_omim():
  """Find and save article data from OMIM for GARD diseases"""
  with GraphDatabase.driver("bolt://3.128.119.166:7687",auth=("neo4j", "back up Generic Plains 27")) as driver:
    try:
      with driver.session() as session:          
        #check for omim updates
        OMIM_and_gard_diseases = get_OMIM_and_GARD_diseases_list()
        for idx, gard_id in enumerate(OMIM_and_gard_diseases):
          time.sleep(0.34)
          try:
            oidNumber = str.split(OMIM_and_gard_diseases[gard_id],sep=":")[1]
            omim_result = find_OMIM_articles(oidNumber)
            omim_pmids = jmespath.search("omim.entryList[0].entry.referenceList[*].reference.pubmedID",omim_result)
            omim_refs = get_omim_refs(omim_result)

            logging.debug('Invoking get_gard_disease_node')
            disease_node = get_gard_disease_node(session, gard_id)
    
            for pubmedID in omim_pmids:
              try:
                res = session.run("MATCH (a:Article) WHERE a.pubmed_id='$pmid' RETURN a", parameters={"pmid":pubmedID})
                matching_articles = len(list(res.value()))
                search_source = omim_refs[pubmedID]
                if (matching_articles):
                  update_article_references(session, pubmedID, search_source)
                  logging.warning(f'omim evidence article already exists: {pubmedID}')
                else:
                  logging.warning(f'Fetch pubmed abstract({pubmedID})')
                  abstractData = fetch_abstract(pubmedID)
                  save_all_data(abstractData, disease_node, pubmedID, search_source, session)
              except urllib.request.URLError as e:
                logging.error(f'Exception when finding articles: {e}')
              continue

          except urllib.request.URLError as e:
            logging.error(f'Exception when finding articles: {e}')
          continue
    except ServiceUnavailable as e:
        logging.error(f'Exception: {e}')
        pass

def save_initial_articles():
  """Write initial article data for the GARD Diseases"""
  mindate = '2000/01/01'
  today = date.today().strftime("%Y/%m/%d")
  logging.info(f'Started save_new_artilcles. Mindate: {mindate}  Maxdate: {today}')
  save_articles_pubmed(mindate, today)
                
def save_new_articles():
  yesterday = (date.today() - datetime.timedelta(days=1)).strftime("%Y/%m/%d")
  today = date.today().strftime("%Y/%m/%d")
  logging.info(f'Started save_new_artilcles. Mindate: {yesterday}  Maxdate: {today}')

  save_articles_pubmed(yesterday, today)

def main():
# logging.config.fileConfig('logging.conf')

# create logger
# logger = logging.getLogger('simpleExample')
  logging.info('Started')
  #save_initial_articles()
  #save_new_articles()
  save_articles_omim()
  logging.info('Finished')


if __name__ == '__main__':
  main()

#%% 
print("hello world")

