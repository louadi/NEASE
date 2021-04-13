import numpy as np 
import networkx as nx
import scipy.stats as stats
import pandas as pd
import csv
import pickle
import statsmodels.api as sm
import plotly.graph_objects as go




pd.options.mode.chained_assignment = None  # default='warn'

'''

import os

def load_obj(data_folder):
    data_folder=os.path.join(os.path.dirname(__file__), data_folder)
    with open(data_folder + '.pkl', 'rb') as f:
        return pickle.load(f)
        
        
# Databases
database_mapping={}
Pathways={}
# Join graph
network={}
# The PPI
PPI={}

# Human
database_mapping['Human']= pd.read_csv(os.path.join(os.path.dirname(__file__),"data/database/Human.csv"))
Pathways['Human']= pd.read_pickle(os.path.join(os.path.dirname(__file__),"data/pathways/pathways_human"))
network['Human']=load_obj('data/network/graph_human')
PPI['Human']=load_obj('data/network/PPI_Human')





'''



#load the data
def load_obj(data_folder):
    with open(data_folder + '.pkl', 'rb') as f:
        return pickle.load(f)

    
    
# Databases
database_mapping={}
Pathways={}
# Join graph
network={}
# The PPI
PPI={}

database_mapping['Human']= pd.read_pickle("data/database/Human.pkl")
Pathways['Human']= pd.read_pickle("data/pathways/pathways_human")
network['Human']=load_obj('data/network/graph_human')
PPI['Human']=load_obj('data/network/PPI_Human')









##########################################################################################################################
# Functions to process AS tools outputs
import re
import collections

pd.set_option('display.max_colwidth', 1000)




# DIGGER Exon-level link: 
"""Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interactions and domain-domain interactions
into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set 
to visually explore their interactions. The following modes of DIGGER can be used interchangeably:
"""
DIGGER='https://exbio.wzw.tum.de/digger/ID/exon/'






def splitDataFrameList(df,target_column):

        ''' 
        Efficiently split Pandas Dataframe cells containing lists into multiple rows,
        duplicating the other column's values.
        Original code from https://gist.github.com/jlln/338b4b0b55bd6984f883

        df = dataframe to split,
        target_column = the column containing the values to split
        separator = the symbol used to perform the split
        returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
        The values in the other columns are duplicated across the newly divided rows.
        '''

        def splitListToRows(row,row_accumulator,target_column):
            split_row = row[target_column]
            for s in split_row:
                new_row = row.to_dict()
                new_row[target_column] = s
                row_accumulator.append(new_row)
        new_rows = []
        df.apply(splitListToRows,axis=1,args = (new_rows,target_column))
        new_df = pd.DataFrame(new_rows)
        return new_df

    
    
    
# Standard data processing 
# Gene ID     Start       End        Delta  (optionally) 


def process_standard (data,
                      mapping,
                      min_delta):
    
        # verify the input data format
        columns=data.columns

        # make sure you have at least two columns
        if len(columns)<3:
                raise ValueError('Make sure your table have at least 3 columns:    Gene ensembl ID    EXON START    EXON END    dPSI (optional)')

        else:
                genes=list(data[columns[0]])

                # Verify the gene IDs
                if not  all(x.startswith('ENS') for x in genes):
                     raise ValueError(' Could not recognize Ensembl gene ID. Please make sure that the first column corresponds to gene IDs.')

                # Verify the start and end
                try:
                    data[columns[1]]=data[columns[1]].astype(int)
                    data[columns[2]]=data[columns[2]].astype(int)
                except:
                    raise ValueError('Could not find exons coordinates. Please make sure that the second column corresponds to the exon start and the third to the exon end (hg38).')


        # map to domains by calculating the overlap of exon coordinate and domain
        mapping_tb=pd.merge(data, mapping,   left_on=columns[0], right_on='Gene stable ID').drop_duplicates()    
        
        
        # get all coding genes affected by splicing
        # Only genes with Pfam domain will be considred here
        # NCBI id used here for the network visualization
        spliced_genes=list(mapping_tb['NCBI gene ID'].unique())

        if len(mapping_tb)==0:
            return []
        
        try:
                #try to get the delta PSI from the user input
                mapping_tb['max_change']=mapping_tb[columns[3]].astype(float)
                
                # significance filter (delta PSI)
                mapping_tb=mapping_tb[mapping_tb['max_change'].abs()>=min_delta]
        except:
                # No delta PSI given
                print('Delta PSI column was not found. Proceeding with all events (no filtering)')
                mapping_tb['max_change']='-'
                
                
    
        
        # any overlap is considered
        mapping_tb['overl']=mapping_tb[["Genomic coding start", "new_start"]].max(axis=1) <= mapping_tb[["Genomic coding end", "new_end"]].min(axis=1)
        mapping_tb=mapping_tb[mapping_tb['overl']].drop_duplicates()
        
        
        
        mapping_tb=mapping_tb[['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID','max_change']]
        #mapping_tb=mapping_tb.sort_values(['max_change'].abs(), ascending=False)
        
        mapping_tb=mapping_tb.reindex(mapping_tb['max_change'].abs().sort_values(ascending=False).index)
        
        mapping_tb=mapping_tb[mapping_tb['NCBI gene ID'].notnull()]
        mapping_tb['NCBI gene ID']=mapping_tb['NCBI gene ID'].astype('int').astype('str')
        #mapping_tb=mapping_tb.drop_duplicates(['Gene name','NCBI gene ID','Gene stable ID','Pfam ID'],keep= 'first')
        
        
        
        #mapping_tb=mapping_tb.groupby(['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID']).max()['max_change']

        
        return mapping_tb,spliced_genes






# Majiq output
def process_MAJIQ(data,
                  mapping,
                  Majiq_confidence=0.95, 
                  min_delta=0.05 ):
    
        # extract exon skipping events:
        data=data[ data['ES']==True]
        
        
        # helper functions:
        print('Processing MAJIQ format...')
        f = lambda x: [abs(float(y)) for y in x.split(';')]
        l=lambda x:  any( y >= min_delta for y in x )
        data=data
        # get Delta PSI values for each junction
        data["delta"] = data['E(dPSI) per LSV junction'].apply(f)
        data["P(|dPSI|>=0.20)"] = data['P(|dPSI|>=0.20) per LSV junction'].apply(f)

        data["delta_sif"] = data['delta'].apply(l)
        data=data[ data['delta_sif']]

        #filter for significant of diff. AS events
        # only keep diff used junction with confidence higher than 'Majiq_confidence' (for instance: 0.95)

        l=lambda x:  any( y >=Majiq_confidence for y in x )
        data["delta_sif"] = data['P(|dPSI|>=0.20)'].apply(l)
        data=data[ data['delta_sif']]

        z= lambda x: x.split(':')[1]
        data['Gene ID'] = data['Gene ID'].apply(z)
    
    
        junctions=list(data['Junctions coords'])
        confidence=list(data['P(|dPSI|>=0.20)'])   
        junc_confid=dict(zip(junctions, confidence))
        
        
        complexx=0
        jun_to_target={}
        jun_to_source={}
        for j in junc_confid.keys():
            #find the source of the junction
            x=re.findall(r"[\w']+", j)
            source=list(set([l for l in x if x.count(l) == len(x)/2]))

            #complex events : skip for now
            # TO DO later
            if len(source)==0:      complexx+=1

            #source found
            #search for all possible targets
            else: 
                targets=[ y for y in x if y not in source ]


                #check if we have correct confidence for every target
                if len(targets)==len(junc_confid[j]):

                    #filter low confident diff used event
                    targets=[ x  for x,y in zip(targets, junc_confid[j]) if y>= 0.95]

                    #save
                    jun_to_target[j]=targets
                    jun_to_source[j]=int(source[0])
        
        
        t=lambda x:  jun_to_target[x] if x in jun_to_target.keys() else False
        data["targets"] = data['Junctions coords'].apply(t)
        
        
        g=lambda x:  jun_to_source[x] if x in jun_to_source.keys() else False
        data["source"] = data['Junctions coords'].apply(g)
        
        
        data=data[data['targets']!=False]
        
        data=splitDataFrameList(data,'targets')
        
        
        mapping_tb=[]
        
        # check mapped exons
        if len(data)==0:
                    print('None of MAJIQ junctions maps to annotated exons (Ensembl Exons). Try to use the standard input instead.')

        else:
                data=data  [ [ 'Gene ID','E(dPSI) per LSV junction' , 'Junctions coords','P(|dPSI|>=0.20) per LSV junction','delta','targets','source' ]]

                data['targets']=data['targets'].astype(int)


                #Map exons to domain
                mapping_tb=pd.merge(mapping, data,  left_on='Genomic coding start', right_on='targets')
                
                # get all coding genes affected by splicing
                # Only genes with Pfam domain will be considred here
                # NCBI id used here for the network visualization
                spliced_genes=list(mapping_tb['NCBI gene ID'].unique())


                #check if the mapping based on coordinate match the gene ID provided from Majiq
                mapping_tb=mapping_tb[mapping_tb['Gene ID']==mapping_tb['Gene stable ID']]

                m= lambda x: max(list(x))
                mapping_tb['max_change']=mapping_tb['delta'].apply(m)

                mapping_tb=mapping_tb[['Exon stable ID','Gene name',
                               'NCBI gene ID','Gene stable ID','Genomic coding start','Genomic coding end','max_change','Pfam ID','source','targets']]


                mapping_tb=mapping_tb.sort_values(['max_change'], ascending=False)
                mapping_tb=mapping_tb[mapping_tb['NCBI gene ID'].notnull()]
                mapping_tb['NCBI gene ID']=mapping_tb['NCBI gene ID'].astype('int').astype('str')
                mapping_tb=mapping_tb.drop_duplicates()

                # an extra filtering step
                # make sure that the source belong to an annotated exon 
                #mapping_tb=pd.merge(mapping[['Genomic coding end']], mapping_tb,  left_on='Genomic coding end', right_on='source')




                mapping_tb=mapping_tb[['Gene name','NCBI gene ID','Gene stable ID','Exon stable ID','Pfam ID','max_change']].drop_duplicates()


                print('MAJIQ output converted successfully to NEASE format.')


        return mapping_tb,spliced_genes
    




##########################################################################################################################
# Proccess functions


# Functions to convert IDs 
def Entrez_to_name(gene,mapping):
    try:
        return mapping[mapping['NCBI gene ID']==int(gene)]['Gene name'].unique()[0]
    
    except :
        return id

    
       
def exons_to_edges(mapped,G):
        # check if domains have known interactions/binding:
        mapped['domain']=mapped['NCBI gene ID']+'/'+mapped['Pfam ID']
        mapped['Interacting domain']=mapped['domain'].apply(lambda x: G.has_node(x))
        mapped=mapped.rename(columns={"max_change": "dPSI",
                                      "domain": "Domain ID"}).reset_index(drop=True) 
    
        mapped['Visualization link']=''
        mapped.loc[mapped['Interacting domain'],['Visualization link']]=DIGGER+mapped['Exon stable ID']
        return mapped
    
    

def affected_edges(data,Join,mapping):

            # get domains with DDIs
            interacting_domains=data[data['Interacting domain']]



            # Identify binding of affected domains = Edges in the PPI

            t=lambda node: [  x for x in list(set([x.split('/')[0] for x in 
                                            [n for n in Join.neighbors(node)] ])) ]
            interacting_domains['Affected binding (NCBI)']=interacting_domains['Domain ID'].apply(t)

            #Convert IDs to names
            c=lambda x: [ Entrez_to_name(gene,mapping) for gene in list(set(x))]
            interacting_domains['Affected binding']=interacting_domains['Affected binding (NCBI)'].apply(c)

            # count number of affected PPI for every domain
            count=lambda x: len(x)
            interacting_domains['Number of affected interactions']=interacting_domains['Affected binding'].apply(count)

    
            return interacting_domains
    
    
def gene_to_edges(data):
    
    #For every gene get all edges
        gene_edges={}
        for gene in data['NCBI gene ID'].unique():
            edges=data[data['NCBI gene ID']==gene]['Affected binding (NCBI)']
            edges=[item for sublist in edges for item in sublist]
            
            gene_edges[gene]=list(set(edges))
        return gene_edges
    
    
def pathway_enrichment(g2edges,paths, mapping,organism):
    # General enrichment analysis
    
    pathway_genes=[]
    # Totat degree of structural network for human (pre-computer)
    # For statistical test: edge enrichment
    # TO DO for mouse
    if organism=='Human':
        n=52467
    
    # number of effected edges 
    affected_edges=len([item for sublist in g2edges.values() for item in sublist])



    
    # for every path :
    path_name=[]
    path_id=[]
    source=[]
    genes=[]
    p_values=[]
    score=[]
    
    for path in list(paths['pathway']):
        # count of affected edges connected to the pathway
        # specific to that pathway list p
        # initialise the variable for every path 
        connected=0
        genes_tmp=[]
        gene_count=0
        
        try:
            # get path total degree "p" and gene list
            p=int(paths [paths['pathway']==path]['Degree in the structural PPI'])
            path_genes=list(paths [paths['pathway']==path]['entrez_gene_ids'])[0]
            
        except:
            pass
        
        for gene in g2edges:
                # for every affected gene
                # count affected gene edges connected to the 
                # specific to the gene and to the pathway list
                tmp=len([x for x in g2edges[gene] if x in path_genes ])

                if tmp>0:
                    # gene with edges connected to the pathway
                    
                    gene_count=gene_count+1
                    # increment for path edges
                    connected=connected+tmp
                    
                    # add gene to the gene list of the pathway
                    genes_tmp.append(Entrez_to_name(gene,mapping) +" ("+str(tmp)+")")
        
        
        #  affected edges not connected to tha pathway
        not_connected=affected_edges-connected
        
        # Join function is slow can be optimized later 
        if genes_tmp==[]: genes_tmp= ''
        else: genes_tmp=(", ").join(genes_tmp)
            
        # Run hypergeometric test on affected edges
        _,p_value_tmp=edge_enrich(connected ,not_connected , p, n)
        
    

        p_values.append(p_value_tmp)
        
        #compute combined score
        #score.append( (1+(gene_count/2)) * -np.log(p_value_tmp) )
        path_name.append(path)
        source.append(list(paths[paths['pathway']==path]['source'])[0])   
        path_id.append(list(paths[paths['pathway']==path]['external_id'])[0]) 
        genes.append(genes_tmp)

        
    # save results
    Enrichment = pd.DataFrame(list(zip(path_id, path_name,source, genes,p_values)), 
                              columns= ['Pathway ID', 'Pathway name', 'Source', 'Spliced genes (number of interactions affecting the pathway)', "p_value"] )
    

    
    
    return Enrichment.sort_values(['p_value'], ascending=True)



def single_path_enrich(path_id,Pathways,g2edges,mapping,organism):
        

        
        # Totat degree of structural network for human (pre-computer)
        # For statistical test: edge enrichment
        # TO DO for mouse
        if organism=='Human':
            n=52467
            
        p=int(Pathways [Pathways['external_id']==path_id]['Degree in the structural PPI'])
        path_genes=list(Pathways [Pathways['external_id']==path_id]['entrez_gene_ids'])[0]
        
        #collect:
        spliced_genes=[]
        spliced_genes_entrez=[]
        gene_association=[]
        num=[]
        affected_edges=[]
        affected_edges_entrez=[]
        p_val=[]
        
        # graph to save affected edges
        G = nx.Graph()
        
        for g in g2edges:
            
            # affected edges of the gene g
            affected=g2edges[g]
            
            # edges connected to the pathway
            edges=[x for x in affected if x in path_genes ]
            a=len(edges)
            
            # Not connected
            b=len(affected)-a

            if a!=0:
                # calculate gene specific p_value:
                _,p_value=edge_enrich(a,b,p,n)
                
                # Save results
                spliced_genes.append(Entrez_to_name(g,mapping))
                spliced_genes_entrez.append(g)
                gene_association.append(g in path_genes)
                num.append(str(a)+'/'+str(a+b))
                affected_edges.append((',').join([ Entrez_to_name(x,mapping) for x in edges]))
                affected_edges_entrez.append((',').join(edges))
                p_val.append(p_value)
                
                # save affected edges
                G.add_edges_from( [ (g,x) for x in edges ]  )
                
        Enrichment = pd.DataFrame(list(zip(spliced_genes,spliced_genes_entrez, gene_association, num,p_val, affected_edges,affected_edges_entrez)), 
                          columns= ['Spliced genes','NCBI gene ID','Gene knwon to be in the pathway','Percentage of edges associated to the pathway', 'p_value', 'Affected binding (edges)','Affected binding (NCBI)'] )
        
        
        
        
        
        return Enrichment,G
    
    
    
    
    
    
    
####################################################################################################################################################################

# stats functions
# Statistical test
# Edge enrichment test

def edge_enrich(a,b,p,n):
    #function to calculate P value
    #test if affected edges by AS are significally enriched in a pathway
    # fisher exact test
    # a+b affected adges, with a the one linked to pathway p
    #p total degree of pathway p
    #n total edges in the ppi used.
    
    #linked to pathway but not affected
    c=p-a
    
    # background of test: not linked to p and not affected edges
    d=(2*n)-p-b
    
    # retun oddsratio and pvalue from fisher exact test
    return   stats.fisher_exact([[a, b], [c, d]], alternative='greater')









####################################################################################################################################################################

# Main Code

class run(object):

    
    
    
    def __init__(self, 
                 data ,
                 organism,
                 input_type='Standard',
                 min_delta=0.05,
                 Majiq_confidence=0.95):
        
        
        
        """
            data: dataframe with list of (diff.) splicing events or junctions. 
                
                Standard input: Ensemble gene ID    Start of exon  End of exon
                for external tools, Please change the  input_type to "MAJIQ",...
                #TO DO
         """
        self.data=[]
        self.organism=organism
        
        if organism!='Human' and organism!='Mouse':
            print('Error: Please choose one of the  supported  organism: "Human" and "Mouse".')
        
     

            
        else:
        
            # TO DO
            
            #Open the Join graph and databases of the selected organism:
            Join=network[organism]
            self.mapping=database_mapping[organism]
            self.path=Pathways[organism]
            self.ppi=PPI[organism]
            
            self.data=[]
            if input_type=='MAJIQ':

                # Processing Majiq output
                try:
                        self.data,self.spliced_genes=process_MAJIQ(data,self.mapping, Majiq_confidence, min_delta )
                        if len(self.data)==0:
                            print('Found no overlap with protein domains.')
                except:
                        print('Could not recognize MAJIQ format. Please make sure your table matches MAJIQ output or use the standard format.')
            
            
            
            elif input_type=='Standard':
                    
                try:
                    self.data,self.spliced_genes=process_standard(data,self.mapping,min_delta )
                    if len(self.data)==0:
                        print('Found no overlap with protein domains.')
                        print('Make sure that the genomic coordinates of the exons correspond to the human genome build hg38 (GRCh38).')

                    
                except:
                        print('Could not recognize the standard format. Please make sure your table matches the standard format.')
                        print('Gene ensembl ID          EXON START        EXON END          dPSI (optional)')
                        print('Make sure that the genomic coordinates of the exons correspond to the human genome build hg38 (GRCh38).')
                    
            
            
        


            self.data=self.data.drop_duplicates(['Gene name','NCBI gene ID','Gene stable ID','Pfam ID'],keep= 'first')


            if len(self.data)==0:
                print('process canceled...')


            else : 
                # check interaction of the domains
                self.data=exons_to_edges(self.data,Join)
                print('\n\t\tData Summary')
                print('**************************************************')

                print(str(len(self.data['Pfam ID'].unique()))+' protein domains are affected by AS.\n'
                      + str(len(self.data[self.data['Interacting domain']]['Pfam ID'].unique()))+' of the affected domains have known interactions.' ) 


                # Identify binding of affected domains = Edges in the PPI
                self.interacting_domains=affected_edges(self.data,Join,self.mapping)



                #get all edges of a gene
                self.g2edges=gene_to_edges(self.interacting_domains)
                
                print(str(len([item for sublist in self.g2edges.values() for item in sublist]))+' protein interactions/binding affected.')

                
                # Runing Enrichment analysis
                
                print('\n**************************************************')
                print('Running enrichment analysis...')
                
                self.supported_database=  list(self.path['source'].unique())
                self.enrichment=pathway_enrichment(self.g2edges,self.path, self.mapping,organism).reset_index(drop=True)
                print('NEASE enrichment done.')
                
    def get_domains(self):
        
        """
            Display the list of AS events in NEASE format.
         """
                
        if len(self.data)==0 :
            print('Processing failed')

            
        elif self.organism=="Mouse":
                #no visualization available for mouse in DIGGER
                return self.data.drop(columns=['Domain ID','Visualization link'])
        else:
            
            #DIGGER visualization available for Human
            return self.data.drop(columns=[ 'Domain ID'])
    
    
    def get_edges(self):
        
        """   
            Display affected interactions from AS. 
        """
        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.')
        else:
            edges=self.interacting_domains[['Gene name','NCBI gene ID','dPSI','Pfam ID','Number of affected interactions','Affected binding','Affected binding (NCBI)']]
            a=lambda x: ",".join(x)
            edges['Affected binding']=edges['Affected binding'].apply(a)
            edges['Affected binding (NCBI)']=edges['Affected binding (NCBI)'].apply(a)
            edges=edges.drop_duplicates()
            edges=edges.sort_values('Number of affected interactions', ascending=False)

            return edges.reset_index(drop=True)
    
          
    
    def enrich(self, database=  ['PharmGKB',
                                 'HumanCyc',
                                 'Wikipathways',
                                 'Reactome',
                                 'KEGG',
                                 'SMPDB',
                                 'Signalink',
                                 'NetPath',
                                 'EHMN',
                                 'INOH',
                                 'BioCarta',
                                 'PID'], cutoff=0.05 ):
        
        """ 
        Run enrichement analysis
        database: List of gene set sources for enrichment.
        
        """
        
        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.') 
        
        else:
            # Check if user input matches the available databases
            database=[ x for x in database if x in self.supported_database]

            if len(database)==0: 
                print('Please select a supported pathway database as an argument. ')
                print('supported databases for '+self.organism+' are :',[ x for x in  self.supported_database],".")
                print('\n')
            else:


                enrich_results=self.enrichment[self.enrichment['Source'].isin(database)]
                
                # Correct for multiple testing
                # fdr_bh : Benjamini/Hochberg (non-negative)
                #enrich_results['adj p_value']=sm.stats.multipletests(list(enrich_results['p_value']),method='fdr_bh',alpha=0.01)[1]
                enrich_results['adj p_value']=sm.stats.fdrcorrection(list(enrich_results['p_value']),alpha=0.01)[1]
                
                

                print('NEASE enrichment for the pathway databases:\n',[ x for x in database])
                num=len(enrich_results[enrich_results['adj p_value']<=cutoff])
                if num==0:
                    print('No enrichment found with the cutoff '+str(cutoff)+'.')
                else:
                    print("Found "+str(num)+" enriched pathways after multiple testing correction.\n")

                return enrich_results.sort_values(['p_value']).reset_index(drop=True)

    
    def path_analysis(self,path_id):
        
        '''
        Run enrichment analysis on a specific pathway with details of impact of AS.
        '''
            
        if len(self.data)==0:
            print('Processing failed')
        elif len(self.interacting_domains)==0:
            print('No affected edges identified.') 
            
            
        path_info=self.enrichment[self.enrichment['Pathway ID']==path_id]
        
        if len(path_info)==0:
            print('No pathway with the given ID found.')
        
        else:
            path_name=list(path_info['Pathway name'])[0]
            print('Enrichment of the pathway: '+path_name+'.\n')
            print('Overall p_value: ',list(path_info['p_value'])[0])
            print('\n')
            # run enrichment
            enrich,_=single_path_enrich(path_id,self.path,self.g2edges,self.mapping,self.organism)
            
            if len(enrich)==0:
                print('No enrichment or genes found for the selected pathway.')
                return
            else:
                return enrich.sort_values(['p_value']).reset_index(drop=True)
            
            


    def Vis_path(self,path_id, k=0.8):

            '''
               Visualize the network module of a specific pathway.

                    K: float (default=None))is a parameter to be tuned by the user:
                    Position nodes using Fruchterman-Reingold force-directed algorithm.
                    Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                    Increase this value to move nodes farther apart.
                    Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html
            '''


            enrich,affected_graph=single_path_enrich(path_id,self.path,self.g2edges,self.mapping,self.organism)

            if len(enrich)==0:
                return

            else:

                # Get genes of the pathway (Entrez IDs)
                path_genes=list(self.path[ self.path['external_id']==path_id ]['entrez_gene_ids'])[0]
                
                significant= list(enrich [enrich ['p_value']<=0.05] ['NCBI gene ID'].unique())
                
                graph_data= extract_subnetwork(path_genes,
                                               self.ppi,
                                               list(enrich['NCBI gene ID'].unique()),
                                               self.spliced_genes,
                                               k,
                                               self.mapping,
                                               affected_graph,
                                               significant)
                
                
                
                
                path_info=self.enrichment[self.enrichment['Pathway ID']==path_id]
                path_name=list(path_info['Pathway name'])[0]
                
                fig = go.Figure(data=graph_data,
                                 layout=go.Layout(
                                    title='<br>'+path_name,
                                    titlefont_size=16,
                                    showlegend=False,
                                    hovermode='closest',
                                    margin=dict(b=20,l=5,r=5,t=40),
                                    annotations=[ dict(
                                        text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> p/</a>",
                                        showarrow=False,
                                        xref="paper", yref="paper",
                                        x=0.005, y=-0.002 ) ],
                                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                                    )



                fig.write_html(path_name+'.html', auto_open=True)
                return
                




                
                
###################################################################################################################################

# Visualization Function

def extract_subnetwork(path_genes,
                       ppi,
                       affected_genes,
                       all_spliced_genes,
                       k,
                       mapping,
                       affected_graph,
                       significant):

    
        # Affected_genes: genes with lost/gained interaction in the pathway
        # all_spliced_genes: all genes affected with splicing
        
        

        # Extract the pathway module for the complete PPI
        # We would like to visualize the pathway with affected edges:

        G=ppi.subgraph(path_genes+affected_genes)
        G = nx.Graph(G)
        G.remove_nodes_from(list(nx.isolates(G)))

        
        
        #Position nodes using Fruchterman-Reingold force-directed algorithm.
        pos = nx.spring_layout(G, k=k, iterations=100)
        
        
        # Prepare the visualization
        for n, p in pos.items():
            G.nodes[n]['pos'] = p
            
        
        # Define node and edges in plott
        node_trace = go.Scatter(x=[],
                                y=[],
                                text=[],
                                mode='markers+text',
                                hoverinfo='text',
                                textposition='top center',
                                marker=dict(
                                reversescale=True,
                                    color=[],
                                    size=[],
                                    line=dict(width=0)))
        
        
        edge_trace = go.Scatter(x=[],
                                y=[],
                                mode='lines',
                                line=dict(width=0.7,
                                          color='#888'),
                                hoverinfo='none')
        
        colored_edge_trace = go.Scatter(x=[],
                                y=[],
                                mode='lines',
                                line=dict(width=4,
                                          color='red'),
                                hoverinfo='none')

        for node in G.nodes():
            x, y = G.nodes[node]['pos']
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_info=Entrez_to_name(node,mapping)
            node_trace['text']+=tuple([node_info])
            
            if not node in path_genes:
                #Node not in pathway
                color='orange'
                if node in significant:
                    size=50
                else:size=30

            elif int(node) in all_spliced_genes:
                # spliced and part of the pathway
                color='red'
                if node in significant:
                    size=50
                else:size=30
                
            else:
                # other pathway nodes
                color='#888'
                if node in significant:
                    size=50
                else:size=20
                
            node_trace['marker']['color']+=tuple([color])
            node_trace['marker']['size']+=tuple([size])
            
            
        
        for edge in G.edges():
            
                x0, y0 = G.nodes[edge[0]]['pos']
                x1, y1 = G.nodes[edge[1]]['pos']
                
                

                # Check if edge is affected
                if affected_graph.has_edge(*edge):
                        colored_edge_trace['x'] += tuple([x0, x1, None])
                        colored_edge_trace['y'] += tuple([y0, y1, None])
                    
                    
                else:
                        edge_trace['x'] += tuple([x0, x1, None])
                        edge_trace['y'] += tuple([y0, y1, None])


            
        
            
        

        return [ colored_edge_trace,node_trace,edge_trace]
            
        
        
