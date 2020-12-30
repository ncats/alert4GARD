# alerts4GARD

### Key Notebook
- *mondo_and_pubmed_evidence_comparison.ipynb* : Various tools to compare the mondo articles and evidence against Pubmed terms api and code to build a dataframe of mondo, gard and omim comparisons

### Python files
- *invoke-esearch.py*: Code to fetch, parse and load pubmed, omim and mondo articles into alerts neo4j
- *gene_review.py*: Code to fetch, parse and load the gene review books into alerts neo4j
- *mondo.py*: Code to fetch, parse and load the mondo data into alerts neo4j
- *PMC.py*: Code to fetch, parse and load the full pmc text into alerts neo4j

### Settings files
- *settings.ini*: settings for connecting to the neo4j instances

### Documentation
- *NIH Alerts System.drawio*: architecture diagram
