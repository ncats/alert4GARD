# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import numpy as np

#graph network packages
import obonet
import networkx as nx

#API requests packages
import requests
import json
from pandas.io.json import json_normalize

# For Analysis
import sklearn.metrics as metrics

# %% [markdown]
# # HPO Annotations File
# 

# %%
def get_hpo_annotations_and_clean():
    """ This function loads data from the HPO annotations file:  
    http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
    
    It then uses the inner function **get_references_one_per_row function** to ensure that references,
    which are the sources indicated for the HPO code, are each represented in their own row

    Input: None
    Output: Dataset from HPO with thier annotations for each disease
    """
    
    def get_references_one_per_row (dataset):

        semicolon = dataset['reference'].str.split(";", n = 0, expand = True)
        semicolon_data = pd.concat([dataset['reference'], semicolon], axis=1)

        semi_done = pd.DataFrame()
        for i in range(1,len(semicolon_data.columns)):
            one = semicolon_data.iloc[:,[0,i]]
            one.columns = ['reference', 'one_reference'] 

            semi_done = semi_done.append(one)

        semi_done = semi_done[(semi_done['one_reference'] != '')]
        semi_done = semi_done[~semi_done['one_reference'].isnull()]

    #Repeating first step with comma.  Should be one function
        semicolon_comma = semi_done['one_reference'].str.split(",", n = 0, expand = True)

        semicolon_comma_data = pd.concat([semi_done['reference'], semicolon_comma], axis=1)

        semi_comma_done = pd.DataFrame()
        for i in range(1,len(semicolon_comma_data.columns)):
            one = semicolon_comma_data.iloc[:,[0,i]]
            one.columns = ['reference', 'one_reference'] 

            semi_comma_done = semi_comma_done.append(one)

        semi_comma_done = semi_comma_done[(semi_comma_done['one_reference'] != '')]
        semi_comma_done = semi_comma_done[~semi_comma_done['one_reference'].isnull()]

        final_data = semi_comma_done.drop_duplicates()
        return final_data
    
    hpo_file = pd.read_table('http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab')

    corrected_references = get_references_one_per_row(hpo_file)

    final_data = hpo_file.merge(corrected_references, how = 'left', on = 'reference')
    final_data['reference'] = final_data['one_reference']

    final_data = final_data.drop(['one_reference'], axis=1)

#ADDING THESE COLUMNS into the cleaning process
    final_data = final_data.rename(columns = {'HPO-ID': 'hpo'})
    final_data['uniqueid'] = final_data['hpo'] + final_data['reference']


    return final_data


# %%
def hpo_annotations_pmid_and_disease(hpo_annotations_dataset):
    """  This function selects all of the pubmed ids from the hpo annotations dataset.
    """

    hpo_annotations_dataset['text_id'] = hpo_annotations_dataset['reference']
    
    #We only want the pubmed articles and thier related diseases
    direct_annotations = hpo_annotations_dataset[hpo_annotations_dataset['text_id'].str.contains('PMID:', na = False)]
    direct_annotations['text_id'] = direct_annotations['text_id'].str.replace('PMID:', '')


    each_unique_pair = direct_annotations[['disease-name', 'text_id']].drop_duplicates()

    #Some pubmeds reference more than one disease.  We want to make sure we remove
    #those pubmeds
    checking_singilarity = each_unique_pair.groupby(['text_id']).count().reset_index()
    two_count = checking_singilarity[checking_singilarity['disease-name'] > 1]
    
    final_list_for_diseases = each_unique_pair[~each_unique_pair['text_id'].isin(two_count['text_id'])]

    
    return final_list_for_diseases

# %% [markdown]
# # Functions for HPO Graph (OBO) File

# %%
def hpo_hierarchy_graph_load():
    """Loads all data from the graph/obo file located here:

    https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
    
    Input: None

    Output: HPO IDs and their hierarchy
    """
    url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
    graph = obonet.read_obo(url)

    return graph


# %%
def graph_phenotypic_abnormality(graph_network):
    """This creates a list of HP codes that are children of the phenotypic abnormality.
        It checks if each HPO code in the graph is a child of HP:0000118, which is phenotypic abnormality.

        Input: HPO graph_network
        Output: List of hpo codes that are under phenotypic abnormality
        
    """
    phenotypic_list = []

    for i in range(len(list(graph_network.nodes))):
    
        if  nx.has_path(graph_network, list(graph_network.nodes)[i],'HP:0000118') == True:
        
            phenotypic_list.append(list(graph_network.nodes)[i])
    
    return phenotypic_list


# %%
def get_all_external_references_for_hpo_codes(graph_network):
    """ This function pulls all of the other codes associated with each hpo code

    Input: HPO graph network
    Output: DataFrame with hpo codes and associated external references
    """
    
    a = nx.get_node_attributes(graph_network, 'xref')
    a_df = pd.DataFrame.from_dict(a,orient='index')
    a_df['HPO'] = a_df.index

    all_references = pd.DataFrame()
    
    for i in range(0,(len(a_df.columns)-1)):
        one = a_df.iloc[:,[12,i]]
        one.columns = ['hpo', 'x_ref'] 
        
        all_references = all_references.append(one)

    all_references = all_references[(all_references['x_ref'] != '')]
    all_references = all_references[~all_references['x_ref'].isnull()]

    all_references = all_references.rename_axis(None, axis = 1)

    return all_references


# %%
def get_hpos_for_umls_code(graph_network):
    """Metamap is coded under UMLS concept codes.  We need to convert these into HPO code.
    There are external references for HPOs codes in the graph network of HPO phenotypes.
    
    This function gets all of the umls codes for hpos that are phenotypic abnormalities.
    
    Child Function: get_all_external_references_for_hpo_codes
    """

    all_references = get_all_external_references_for_hpo_codes(graph_network)
    umls_references = all_references[all_references['x_ref'].str.contains('UMLS')]
    umls_references['x_ref'] = umls_references['x_ref'].str.replace('UMLS:', '') 
    
    # this matches the naming convension for metamap
    umls_references = umls_references.rename(columns={'x_ref': 'conceptId'})

    return umls_references

# %% [markdown]
# ## Graph Functions that Support Comparison

# %%
def get_alt_ids(hpo_code, graph_network):

    a = 'alt_id' in graph_network.node[hpo_code]

    if a == True:
        one_time = graph_network.node[hpo_code]['alt_id']
        one_time_dataset = pd.DataFrame(one_time)
        one_time_dataset['hpo'] = hpo_code
        one_time_dataset.columns = ['alt_hpo', 'hpo']

    else:
        one_time_dataset = pd.DataFrame()

    return one_time_dataset


# %%
def graph_alternate_direct_ids(List_of_HPOs, graph_network):
    """This gets all alternate HPO codes for each given HPO.  
    There can be many HPO codes being used for the same phenotype issues.
    
    Child function:
        - get_alt_ids - this supplies the alternate HPO codes for one given HPO code"""
    
    final_data = pd.DataFrame()

    for i in range(0, len(List_of_HPOs)):
        ex_ref = get_alt_ids(List_of_HPOs[i], graph_network)
        final_data = final_data.append(ex_ref)

    return final_data
    


# %%
def get_child_parent(hpo_code, graph_network):
    """ This function takes one hpo code and returns the parent and child hpos.
    """

        # Internal Function to get the children of an hpo code
    def get_child(hpo_code, graph_network):
        children = list(graph_network.predecessors(hpo_code))

        one = pd.DataFrame(children, columns=['related_hpo'])
        one['relationship'] = 'Child'
        one['hpo'] = hpo_code

        return one

        # Internal Function to get the parent of an hpo code
    def get_parent(hpo_code, graph_network):
        parent = list(graph_network.successors(hpo_code))

        one = pd.DataFrame(parent, columns=['related_hpo'])
        one['relationship'] = 'Parent'
        one['hpo'] = hpo_code

        return one 

    child = get_child(hpo_code, graph_network)
    parent = get_parent(hpo_code, graph_network)

    full_child_parent = parent.append(child)

    full_child_parent = full_child_parent[['related_hpo', 'hpo', 'relationship']]

    return full_child_parent


# %%
def get_child_parent_and_alternatives(hpo_code, graph_network):
    """ Combines the:
        1. get_child_parent function, with
        2. graph alternate_direct_ids function

        To get a full list of all parent and child and thier alternative ids for one hpo code.
    """
    #parent and children
    full_child_parent = get_child_parent(hpo_code, graph_network)
    
    #we now check the parent and children for alternate ids
    alt_id_check_list = list(full_child_parent['related_hpo'])
    alternate_ids = graph_alternate_direct_ids(alt_id_check_list, graph_network)
    
    #this is to check if there are alternate ids for any of the relatives
    if len(alternate_ids) > 0:
        alternate_ids.columns = ['all_hpo_options', 'related_hpo']

        merging = full_child_parent.merge(alternate_ids, left_on = 'related_hpo', right_on = 'related_hpo' )

        full_child_parent['all_hpo_options'] = full_child_parent['related_hpo']


        final_data = full_child_parent.append(merging)
    
    else:
        final_data = full_child_parent
        final_data['all_hpo_options'] = final_data['related_hpo']
    
    final_data = final_data[['all_hpo_options', 'hpo', 'relationship']]
    final_data.columns = ['related_hpos_with_alternates', 'hpo', 'relationship']

    return final_data  


# %%
def graph_parent_child(dataset, graph_network):
    """Want to get the HPO codes for each child and parent of a given HPO code. 
    We also want to get the alternate ids so that we have all of the possible HPO codes that are children and parents.

    Child Formula:
     - get_child_parent_and_alternatives
    """
     
    final_data = pd.DataFrame()
    
    for i in range(0,len(dataset)):
        one_row = dataset.iloc[i,:]


        ex_ref = get_child_parent_and_alternatives(one_row['hpo'], graph_network)
        ex_ref['text_id'] = one_row['text_id']
        final_data = final_data.append(ex_ref)
    
    if len(final_data) == 0:
        final_data = pd.DataFrame(columns=['alt_id', 'hpo', 'relationship', 'text_id'])
    else:
        final_data.columns = ['alt_id', 'hpo', 'relationship', 'text_id'] 
        
    return final_data

# %% [markdown]
# # API Calls and Cleaning for Mondo and MetaMap
# 
# ## Mondo Annotation

# %%
def get_one_mondo_annotation(text, text_id,  min_word_length = 4, longest_only = 'true', include_abbreviation = 'false', include_acronym = 'false',
                            include_numbers = 'false'):

    """Creates a call to the mondo anotator.  it allows for the same modifications to the mondo annotator as the api itself.
    
        Requires input of:
            - Text to be examined
            - ID for the text, most commonly PubMedId
    """

    #needs to be a string to add to the call to the api
    min_word_length = str(min_word_length)

    call = 'https://api.monarchinitiative.org/api/nlp/annotate/entities/?content=' + text + '&min_length=' + min_word_length + '&longest_only=' + longest_only + '&include_abbreviation=' + include_abbreviation + '&include_acronym=' + include_acronym + '&include_numbers='+ include_numbers
    
    call = str(call)
    
    r = requests.get(call)
    
    # Checking to make sure the request went through to the api sucessfully
    if r.status_code == 200:
        
        x = r.json()
        data_call_df = json_normalize(x['spans'], record_path=['token'])
        one = data_call_df.reset_index()
        one['text_id'] = text_id    
        one = one[['id', 'category', 'terms', 'text_id']].astype(str)
        one = one.drop_duplicates()
    
    else:
        one = pd.DataFrame(columns = [ 'id', 'category', 'terms', 'text_id'])
        
    return one


# %%
def get_all_mondo_annotations(dataset, min_word_length_all = 4, longest_only_all = 'true', 
    include_abbreviation_all ='false', include_acronym_all = 'false',include_numbers_all = 'false'):

    """ Loops through a dataset with text and IDs. 
        Requires a dataset with:
        - A "text" column
        - A "text_id" column to specify the id for the article

    """
   
    mondo_long = pd.DataFrame()

    for i in range(0,len(dataset)):

        one_abstract = dataset.iloc[i, :]

        #getting mondo_long
        one_mondo_long = get_one_mondo_annotation(one_abstract['text'], one_abstract['text_id'],  
                                                            min_word_length = min_word_length_all, longest_only = longest_only_all, 
        include_abbreviation = include_abbreviation_all, include_acronym = include_acronym_all,
                                                            include_numbers = include_numbers_all)

        mondo_long = mondo_long.append(one_mondo_long)


    mondo_long.rename(columns = {'id':'annotation_id'}, inplace = True)   

    return mondo_long


# %%
def cleaning_mondo_to_hpo(mondo_data_frame):
    """Makes sure HPO code is being referenced.  
        Can be used with one or many annotations.
    """

    mondo_hpo = mondo_data_frame[mondo_data_frame['annotation_id'].str.contains('HP:', na = False)]
    #mondo_hpo = mondo_hpo[mondo_hpo['id'].isin(phenotypic_abnormality_hpo)]

    mondo_hpo.rename(columns = {'annotation_id':'hpo'}, inplace = True) 
    return mondo_hpo

# %% [markdown]
# ## Metamap

# %%
def get_one_metamap_annotation(text, text_id):
    """The metamap api cannont take more than ~1000(?) characters, so we break the text into one sentance each,
    then run each sentance through the metamap annotator.
    """
    metamap_test_data = pd.DataFrame()

   # MetaMap has a short limit on characters that can be put into the annotator through the api
    abstract_sentaces = text.split('. ')
   
    for abstract_sentace_part in range(len(abstract_sentaces)): # send each sentance through the annotator one by one
        
        call = 'https://knowledge.ncats.io/ks/umls/metamap/' + abstract_sentaces[abstract_sentace_part]
        call = str(call)
        r = requests.get(call)
        
        if r.status_code == 200: # checks if api call is successful
            x = r.json()
            metamap_test_data_one_example = pd.DataFrame()
            for utterance_len in range(len(x['utteranceList'])):  # Each  annotation is nested, so we need to loop through to get all of them
                
                for pcmlist_len in range(len(x['utteranceList'][utterance_len]['pcmlist'])):
                    
                    if "mappingList" in x['utteranceList'][utterance_len]['pcmlist'][pcmlist_len]:
                        
                        for mapping in range(len(x['utteranceList'][utterance_len]['pcmlist'][pcmlist_len]['mappingList'])):
                            #print('mapping')
                            for evList1 in range(len(x['utteranceList'][utterance_len]['pcmlist'][pcmlist_len]['mappingList'][mapping]['evList'])):
                                
                                one1 = json_normalize(x['utteranceList'][utterance_len]['pcmlist'][pcmlist_len]['mappingList'][mapping]['evList'][evList1])
                                metamap_test_data_one_example = metamap_test_data_one_example.append(one1)
            metamap_test_data = metamap_test_data.append(metamap_test_data_one_example)
    
    metamap_test_data['text_id'] = text_id
   

    return metamap_test_data


# %%
def get_all_metamap_annotations(dataset):
    """ Runs a loop through a full dataset of text and text_ids through the get_one_metamap_annotation function."""
    metamap = pd.DataFrame()

    for i in range(0,len(dataset)):

        one_abstract = dataset.iloc[i, :]

        #getting mondo_long
        one_metamap = get_one_metamap_annotation(one_abstract['text'], one_abstract['text_id'])

        metamap = metamap.append(one_metamap)
    
    return metamap
#removing duplicate values


# %%
def cleaning_metamap_adding_hpo(metamap_data_frame, graph_network):
    """ Function requires a metamap data frame AND the graph network for HPO. This graph network is required
        To get all of the umls codes from metamap returned as hpos

        This function:
        1. Gets the hpo codes for the umls returned in the metamap.
        2. Drops duplicate values in the data frame
        3. Attached the hpo code to the metamap data frame and removes umls codes that do not have an hpo equivilant.

        Uses functions: get_hpos_for_umls_code
    """

    umls_to_hpo = get_hpos_for_umls_code(graph_network)
    
    metamap_data_frame['semanticTypes'] = metamap_data_frame['semanticTypes'].astype(str)
    #drop_dups = metamap_data_frame.drop_duplicates(['conceptId', 'pmid'])
    drop_dups = metamap_data_frame.drop_duplicates(['conceptId', 'text_id', 'semanticTypes' ])
    
    metamap_hpo = drop_dups.merge(umls_to_hpo, how = 'left', on = 'conceptId')
    metamap_hpo = metamap_hpo[~metamap_hpo['hpo'].isnull()]

    #changing text_id to 
    metamap_hpo['text_id'] = metamap_hpo['text_id'].astype(str)

    return metamap_hpo

# %% [markdown]
# # Comparison of Known Annotations and Data to Test

# %%
def one_annotation_direct_matching(one_hpo_annotation, one_annotation_to_test, graph_network):

    """ 
    In this function, we take one set of known hpo annotations and one set of annotations to test against this known set. For this test set, we check for hpo matches directly and with alternate ids.

    Input: known_annotations dataset (with at least unqiueid and hpo columns); dataset of annotations to test (with at least hpo and text_id columns); graph network for hpo; graph network
    Output: dataset with exact matches, annotations to test that do not have matches, known hpo annotations that do not have a match.

    Child functions: 
    1. graph_alternate_direct_ids
    """

    #Testing up up the list
    data_to_test = one_annotation_to_test[['hpo', 'text_id']]
    data_to_test['hpo'] = data_to_test['hpo'].astype(str)
    data_to_test['exact_match'] = 1
    data_to_test['alt_id'] = 0
    data_to_test['test_set_annotations_with_no_match'] = 0


    # We want to determine which of our test set hpos are an exact match for the test data
    exact_match = one_hpo_annotation.merge(data_to_test, how = 'left', on =['hpo'])#, 'pmid'])#, 'uniqueid'])
    exact_match = exact_match[exact_match['exact_match'] == 1]

    ##############
    #TESTING ALTERNATE ID
    ##############
    #create the new annotations and test data sets without the hpos that are part of the exact match
    alternate_id = one_hpo_annotation[~one_hpo_annotation['uniqueid'].isin(exact_match['uniqueid'])]
    test_group_for_alt_id = data_to_test[~data_to_test['hpo'].isin(exact_match['hpo'])]

    alt_id = graph_alternate_direct_ids(test_group_for_alt_id['hpo'].tolist(), graph_network)

    # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        exact_alt_match = pd.DataFrame()
    else:
        alt_id['exact_match'] = 0
        alt_id['alt_id'] = 1
        alt_id['test_set_annotations_with_no_match'] = 0
        alt_id = alt_id.merge(test_group_for_alt_id[['hpo', 'text_id']], how = 'left', on = 'hpo') # adding back in the text_id
        alt_id_final = alt_id[['alt_hpo','exact_match', 'alt_id', 'test_set_annotations_with_no_match', 'text_id']]

        exact_alt_match = alternate_id.merge(alt_id_final , how = 'left', left_on = 'hpo', right_on = 'alt_hpo')
        exact_alt_match = exact_alt_match[exact_alt_match['alt_id'] == 1]
        exact_alt_match = exact_alt_match.drop('alt_hpo' , axis='columns')

    ##################
    #NO MATCH
    ###################
    #we still want to capture hpos from the test group that don't match any of the known hpos
    no_match = data_to_test[~data_to_test['hpo'].isin(exact_match['hpo'])]

        # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        none = 1
        #skip this step if there is no alternate ids
    else:
        no_match = no_match[~no_match['hpo'].isin(exact_alt_match['hpo'])]
   
    no_match['exact_match'] = 0
    no_match['alt_id'] = 0
    no_match['test_set_annotations_with_no_match'] = 1
    #test_parent_child_final = test_parent_child[['alt_hpo_code','exact_match', 'alt_id', 'test_set_annotations_with_no_match', 'text_id']]

    ##############################################
    #COMBINING ALL OF THE VALUES FROM TESTS ABOVE
    ###############################################

    combined_findings = exact_match
    combined_findings = combined_findings.append(exact_alt_match)
    combined_findings = combined_findings.append(no_match)

    #combining direct and alternate id since they both represent exact matches
    combined_findings['direct_alt'] = combined_findings['exact_match'] + combined_findings['alt_id']
    combined_findings['exact_match'] = np.where(combined_findings['direct_alt'] >= 1, 1, 0)

    combined_findings = combined_findings[['uniqueid', 'hpo', 'text_id', 'exact_match', 'test_set_annotations_with_no_match']]
    combined_findings['known_annotations_with_no_match'] = 0
  
    #####################
    ### ORGANIZING Known annotations that were not selected

    remaining_from_one_hpo_annotation = one_hpo_annotation[~one_hpo_annotation['uniqueid'].isin(combined_findings['uniqueid'])]
    remaining_from_one_hpo_annotation['exact_match'] = 0
    remaining_from_one_hpo_annotation['test_set_annotations_with_no_match'] = 0
    remaining_from_one_hpo_annotation['known_annotations_with_no_match'] = 1
    remaining_from_one_hpo_annotation['text_id'] =  np.NaN

    remaining_from_one_hpo_annotation = remaining_from_one_hpo_annotation[['uniqueid', 'hpo', 'text_id', 'exact_match', 'test_set_annotations_with_no_match','known_annotations_with_no_match']]

    final_data = combined_findings.append(remaining_from_one_hpo_annotation)
    
    return final_data


# %%
#one_hpo_annotation = hpo_annotations
#one_annotation_to_test = one_mondo_annotation
#graph_network = g
def one_annotation_matching_relatives(one_hpo_annotation, one_annotation_to_test, graph_network):

    """ 
    In this function, we take one set of known hpo annotations and one set of annotations to test against this known set. For this test set, we check for hpo matches directly, with alternate ids, and finally, with any parent, child relationships (along with thier alternate ids).

    Input: known_annotations dataset (with at least unqiueid and hpo columns); dataset of annotations to test (with at least hpo and text_id columns); graph network for hpo; graph network
    Output: dataset with exact matches, matches to relatives, annotations to test that do not have matches, known hpo annotations that do not have a match.

    Child functions: 
    1. graph_alternate_direct_ids
    2. graph_parent_child
    """

    #Testing up up the list
    data_to_test = one_annotation_to_test[['hpo', 'text_id']]
    data_to_test['hpo'] = data_to_test['hpo'].astype(str)
    data_to_test['exact_match'] = 1
    data_to_test['alt_id'] = 0
    data_to_test['relative_match'] = 0
    data_to_test['test_set_annotations_with_no_match'] = 0


    # We want to determine which of our test set hpos are an exact match for the test data
    exact_match = one_hpo_annotation.merge(data_to_test, how = 'left', on =['hpo'])#, 'pmid'])#, 'uniqueid'])
    exact_match = exact_match[exact_match['exact_match'] == 1]

    ##############
    #TESTING ALTERNATE ID
    ##############
    #create the new annotations and test data sets without the hpos that are part of the exact match
    alternate_id = one_hpo_annotation[~one_hpo_annotation['uniqueid'].isin(exact_match['uniqueid'])]
    test_group_for_alt_id = data_to_test[~data_to_test['hpo'].isin(exact_match['hpo'])]

    alt_id = graph_alternate_direct_ids(test_group_for_alt_id['hpo'].tolist(), graph_network)
    
    # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        exact_alt_match = pd.DataFrame()
    else:
        alt_id['exact_match'] = 0
        alt_id['alt_id'] = 1
        alt_id['test_set_annotations_with_no_match'] = 0
        alt_id = alt_id.merge(test_group_for_alt_id[['hpo', 'text_id']], how = 'left', on = 'hpo') # adding back in the text_id
        alt_id_final = alt_id[['alt_hpo','exact_match', 'alt_id', 'test_set_annotations_with_no_match', 'text_id']]

        exact_alt_match = alternate_id.merge(alt_id_final , how = 'left', left_on = 'hpo', right_on = 'alt_hpo')
        exact_alt_match = exact_alt_match[exact_alt_match['alt_id'] == 1]
        exact_alt_match = exact_alt_match.drop('alt_hpo' , axis='columns')


    ####################
    ## ADD IN TESTING OF RELATIVES
    #################
    #now we need to remove any hpos that have been found in the exact match or alternate id process
    relatives_id = one_hpo_annotation[~one_hpo_annotation['uniqueid'].isin(exact_match['uniqueid'])]

    # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        none = 1
        #skip this step if there is no alternate ids
    else:
        relatives_id = relatives_id[~relatives_id['uniqueid'].isin(exact_alt_match['uniqueid'])]

    #remove hpos from the test group that have been found already
    test_group_for_relatives = data_to_test[~data_to_test['hpo'].isin(exact_match['hpo'])]
    
        # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        none = 1
        #skip this step if there is no alternate ids
    else:
        test_group_for_relatives = test_group_for_relatives[~test_group_for_relatives['hpo'].isin(exact_alt_match['hpo'])]
    

    #This code checks the test group hpos against children and parent hpos. It also checks against the alternative ids
    test_parent_child = graph_parent_child(test_group_for_relatives, graph_network)
    
    test_parent_child.rename(columns={'alt_id':'alt_hpo_code', 'hpo': 'original_hpo'}, inplace=True)
    test_parent_child['exact_match'] = 0
    test_parent_child['alt_id'] = 0
    test_parent_child['relative_match'] = 1
    test_parent_child['test_set_annotations_with_no_match'] = 0
    test_parent_child_final = test_parent_child[['alt_hpo_code','exact_match', 'alt_id', 'relative_match', 'test_set_annotations_with_no_match','text_id','original_hpo']]

    # Combine the two datasets to test 
    exact_relative_match = relatives_id.merge(test_parent_child_final, how = 'left', left_on = 'hpo', right_on = 'alt_hpo_code' )
    exact_relative_match = exact_relative_match[exact_relative_match['relative_match'] == 1]
    exact_relative_match = exact_relative_match.drop('alt_hpo_code' , axis='columns')

    ##################
    #NO MATCH
    ###################
    #we still want to capture hpos from the test group that don't match any of the known hpos
    no_match = data_to_test[~data_to_test['hpo'].isin(exact_match['hpo'])]
    
    no_match = no_match[~no_match['hpo'].isin(exact_relative_match['original_hpo'])]

    # Sometimes there are no alternate ids
    if len(alt_id) == 0:
        none = 1
        #skip this step if there is no alternate ids
    else:
        no_match = no_match[~no_match['hpo'].isin(exact_alt_match['hpo'])]


    no_match['exact_match'] = 0
    no_match['alt_id'] = 0
    no_match['relative_match'] = 0
    no_match['test_set_annotations_with_no_match'] = 1
    #test_parent_child_final = test_parent_child[['alt_hpo_code','exact_match', 'alt_id', 'relative_match', 'test_set_annotations_with_no_match', 'text_id']]

    ##############################################
    #COMBINING ALL OF THE VALUES FROM TESTS ABOVE
    ###############################################

    combined_findings = exact_match
    combined_findings = combined_findings.append(exact_alt_match)
    combined_findings = combined_findings.append(exact_relative_match)
    combined_findings = combined_findings.append(no_match)

    #combining direct and alternate id since they both represent exact matches
    combined_findings['direct_alt'] = combined_findings['exact_match'] + combined_findings['alt_id']
    combined_findings['exact_match'] = np.where(combined_findings['direct_alt'] >= 1, 1, 0)

    combined_findings = combined_findings[['uniqueid', 'hpo', 'original_hpo', 'text_id', 'exact_match', 'relative_match', 'test_set_annotations_with_no_match']]
    combined_findings['known_annotations_with_no_match'] = 0

    #####################
    ### ORGANIZING Known annotations that were not selected

    remaining_from_one_hpo_annotation = one_hpo_annotation[~one_hpo_annotation['uniqueid'].isin(combined_findings['uniqueid'])]
    remaining_from_one_hpo_annotation['exact_match'] = 0
    remaining_from_one_hpo_annotation['relative_match'] = 0
    remaining_from_one_hpo_annotation['test_set_annotations_with_no_match'] = 0
    remaining_from_one_hpo_annotation['known_annotations_with_no_match'] = 1
    remaining_from_one_hpo_annotation['text_id'] =  np.NaN
    remaining_from_one_hpo_annotation['original_hpo'] =  np.NaN

    remaining_from_one_hpo_annotation = remaining_from_one_hpo_annotation[['uniqueid', 'hpo','original_hpo', 'text_id', 'exact_match', 'relative_match', 'test_set_annotations_with_no_match','known_annotations_with_no_match']]

    final_data = combined_findings.append(remaining_from_one_hpo_annotation)
    #combining the exact and relative match

    #final_data['exact_and_relative_match'] = np.where((final_data['exact_match'] == 1) | (final_data['relative_match'] == 1), 1, 0)
    #final_data_df = final_data[['']]

    return final_data

# %% [markdown]
# # Scoring Process for Dataset

# %%
def scoring(dataset):

     """
     This function scores the datasets produced in comparing the test data to known annotators.
     """

# Just indicated whether relatives were included
     if 'relative_match' in dataset:
          ScoringType = 'Relative'
     else :
          ScoringType = 'Direct'

     #We need to change the acutally false depending on whether we're counting exact matches or relatives are included
     dataset['predicted_true'] = np.where(~ pd.isna(dataset['text_id']), 1, 0)
     dataset['actually_true'] = np.where(~ pd.isna(dataset['uniqueid']), 1, 0)
     #dataset['actually_false'] = np.where(pd.isna(dataset['uniqueid']), 1, 0)

     #confusion matrix
     dataset['true_positive'] = np.where((dataset['predicted_true'] == 1) & (dataset['actually_true'] == 1), 1, 0) 
     dataset['false_positive'] = np.where((dataset['predicted_true'] == 1) & (dataset['actually_true'] == 0), 1, 0)
     dataset['false_negative'] = np.where((dataset['predicted_true'] == 0) & (dataset['actually_true'] == 1), 1, 0)
     dataset['true_negative'] = 0 #we don't have any true negatives since our annotators are not predicting what is NOT there.
     
     true_positive = dataset['true_positive'].sum()
     false_positive = dataset['false_positive'].sum()
     false_negative = dataset['false_negative'].sum()
     true_negative = dataset['true_negative'].sum()


     #basic statistics
     number_of_annotations_to_test = dataset['predicted_true'].sum()
     known_annotations_to_test_against = dataset['actually_true'].sum()
     accurately_predicted = dataset['true_positive'].sum()

     #Measures
     #accuracy = metrics.accuracy_score(dataset['actually_true'], dataset['predicted_true'])
     precision = metrics.precision_score(dataset['actually_true'], dataset['predicted_true'])
     recall = metrics.recall_score(dataset['actually_true'], dataset['predicted_true'])
     f1_score = metrics.f1_score(dataset['actually_true'], dataset['predicted_true'])



     final_data = {'ScoringType': ScoringType, 'Annotations_to_Test': [number_of_annotations_to_test],
                    'Known_Annotations': [known_annotations_to_test_against],
                    'Accurately_Predicted': [accurately_predicted],
                    #'Accuracy': [accuracy],
                    'Precision': [precision],
                    'Recall': [recall],
                    'F1_Score': [f1_score],
                    'True_Positive': true_positive,
                    'False_Positive': false_positive,
                     'False_Negative': false_negative,
                    'True_Negative': true_negative}

     final_data_df = pd.DataFrame(final_data)
 

     return final_data_df


# %%



