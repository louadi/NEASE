from .process import *
import scipy.stats as stats
import plotly.graph_objects as go

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np




# main functions for nease
      
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
    
    
def pathway_enrichment(g2edges,paths, mapping,organism,p_value_cutoff):
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

                    # increment for path edges
                    connected=connected+tmp
                    
                    # add gene to the gene list of the pathway
                    genes_tmp.append(Entrez_to_name(gene,mapping) +" ("+str(tmp)+")")
                    
                    
                    # gene specific test
                    _,p_gene=edge_enrich(tmp ,len(g2edges[gene])-tmp , p, n)
                        
                    if p_gene<=0.05:
                        # gene with edges siginifically connected to the pathway
                    
                        gene_count=gene_count+1
                    
        
        
        #  affected edges not connected to tha pathway
        not_connected=affected_edges-connected
        
        # Join function is slow can be optimized later 
        if genes_tmp==[]: genes_tmp= ''
        else: genes_tmp=(", ").join(genes_tmp)
            
        # Run hypergeometric test on affected edges
        _,p_value_tmp=edge_enrich(connected ,not_connected , p, n)
        
    

        p_values.append(p_value_tmp)
        
        #compute combined score
        if p_value_tmp<p_value_cutoff:
            s= -(np.sqrt(gene_count) * np.log10(p_value_tmp))
        else: s=0
            
        score.append( s  )
        path_name.append(path)
        source.append(list(paths[paths['pathway']==path]['source'])[0])   
        path_id.append(list(paths[paths['pathway']==path]['external_id'])[0]) 
        genes.append(genes_tmp)

        
    # save results
    Enrichment = pd.DataFrame(list(zip(path_id, path_name,source, genes,p_values,score)), 
                              columns= ['Pathway ID', 'Pathway name', 'Source', 'Spliced genes (number of interactions affecting the pathway)', "p_value",'Nease score'] )
    

    
    
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
                          columns= ['Spliced genes','NCBI gene ID','Gene is known to be in the pathway','Percentage of edges associated to the pathway', 'p_value', 'Affected binding (edges)','Affected binding (NCBI)'] )
        
        
        
        
        
        return Enrichment,G

    
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
        #pos = nx.kamada_kawai_layout(G)
        

        
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
                                #textfont_size = 20,
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
            
        
        

        
        
# plot domains stats

def stats_domains(affecting_percentage,
                  binding_percentage, 
                  file_path):
    
        # from https://matplotlib.org/stable/gallery/pie_and_polar_charts/bar_of_pie.html#sphx-glr-gallery-pie-and-polar-charts-bar-of-pie-py
        # make figure and assign axis objects
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
        fig.subplots_adjust(wspace=0)

        # pie chart parameters
        ratios = [affecting_percentage, 1-affecting_percentage, ]
        labels = ['Affecting a domain', 'Not affecting a domain']
        explode = [0.1, 0]
        # rotate so that first wedge is split by the x-axis
        angle = -180 * ratios[0]
        ax1.pie(ratios, autopct='%1.1f%%', startangle=angle,
                labels=labels, explode=explode, shadow=True,)
        ax1.set_title("Genes with AS affecting protein domains")
        # bar chart parameters

        xpos = 0
        bottom = 0
        ratios = [binding_percentage, 1-binding_percentage]
        width = .2
        colors = ['#66b3ff','#ff9999',]

        for j in range(len(ratios)):
            height = ratios[j]
            ax2.bar(xpos, height, width, bottom=bottom, color=colors[j])
            ypos = bottom + ax2.patches[j].get_height() / 2
            bottom += height
            ax2.text(xpos, ypos, "%d%%" % (ax2.patches[j].get_height() * 100),
                     ha='center')

        ax2.set_title('Domains with known interactions')
        ax2.legend(('Yes', 'No'))
        ax2.axis('off')
        ax2.set_xlim(- 2.5 * width, 2.5 * width)

        # use ConnectionPatch to draw lines between the two plots
        # get the wedge data
        theta1, theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
        center, r = ax1.patches[0].center, ax1.patches[0].r
        bar_height = sum([item.get_height() for item in ax2.patches])

        # draw top connecting line
        x = r * np.cos(np.pi / 180 * theta2) + center[0]
        y = r * np.sin(np.pi / 180 * theta2) + center[1]
        con = ConnectionPatch(xyA=(-width / 2, bar_height), coordsA=ax2.transData,
                              xyB=(x, y), coordsB=ax1.transData)
        con.set_color([0, 0, 0])
        con.set_linewidth(4)
        ax2.add_artist(con)

        # draw bottom connecting line
        x = r * np.cos(np.pi / 180 * theta1) + center[0]
        y = r * np.sin(np.pi / 180 * theta1) + center[1]
        con = ConnectionPatch(xyA=(-width / 2, 0), coordsA=ax2.transData,
                              xyB=(x, y), coordsB=ax1.transData)
        con.set_color([0, 0, 0])
        ax2.add_artist(con)
        con.set_linewidth(4)
        
        file_path=os.path.join(os.path.dirname(file_path),'domains_stats.pdf')

        plt.savefig(file_path,format='pdf',bbox_inches='tight')
        plt.show()
        
        # Save the figure


            
        return