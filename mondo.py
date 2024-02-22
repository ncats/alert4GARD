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


def save_mondo_data():
with GraphDatabase.driver("bolt://3.128.119.166:7687",auth=("neo4j", "back up Generic Plains 27")) as driver:
    with driver.session() as session:
      results = get_Mondo_and_GARD_diseases_list()
      #results = dict(itertools.islice(results.items(),1,2))
      for idx, mondo_id in enumerate(results):                        
        time.sleep(0.34)
        #logging.warning(f'{idx}: Invoking find_articles({results[gard_id]},{mindate},{maxdate})')
        try:
          pubmedIDs = find_mondo_articles(mondo_id)
          if ('associations' in pubmedIDs and len(pubmedIDs['associations']) > 0):
            for association in pubmedIDs['associations']:
              if('publications' in association and association['publications'][0] and association['publications'][0]['id']):
                pubmedid = association['publications'][0]['id']
                pubmedid = str.split(pubmedid,":")[1]
                res = session.run("match(a:Article) where a.pubmed_id = $pmid return a", parameters={"pmid":pubmedid})
                matching_articles = len(res)
                if (matching_articles<1):
                  #we need to add the article to the list


            except urllib.request.URLError as e:
              logging.error(f'Exception when finding articles: {e}')
              continue

            logging.debug('Invoking create_gard_disease')
            disease_node = create_mondo_node(session, gard_id, results[gard_id])


def main():
  save_mondo_data()


if __name__ == '__main__':
  main()