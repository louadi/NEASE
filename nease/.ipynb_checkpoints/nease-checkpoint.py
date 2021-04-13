
from .functions import *
import statsmodels.api as sm



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
            
            


    def Vis_path(self,path_id,file='', k=0.8):

            '''
               Visualize the network module of a specific pathway.
            path_id:    -  str: An unique pathway id. Run  enrich() to get all pathways and their ids.
            K:          -  Float (default=None))is a parameter to be tuned by the user:
                            Position nodes using Fruchterman-Reingold force-directed algorithm.
                            Optimal distance between nodes. If None the distance is set to 1/sqrt(n) where n is the number of nodes. 
                            Increase this value to move nodes farther apart.
                        Link: networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html.
            file         - A string representing a local file path.
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
                                        text="<br> The large nodes have significant p_value (affecting the pathway).<br> 🔴 Spliced gene and known to be part of the patwhay.<br> 🟠 Spliced gene but not known to be in the pathway.",
                                        showarrow=False,
                                        font=dict(size=16),
                                        xref="paper", yref="paper",
                                        x=0.005, y=-0.002 ) ],
                                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))


                fig.write_html(os.path.join(os.path.dirname(file),path_name+'.html'), auto_open=True)
                
                return
                