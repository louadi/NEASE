import pandas as pd
import pickle
import os


pd.options.mode.chained_assignment = None  # default='warn'


# Helper function

def load_obj(data_folder):
    data_folder=os.path.join(os.path.dirname(__file__), data_folder)
    with open(data_folder + '.pkl', 'rb') as f:
        return pickle.load(f)
   

# Save data as dictionary        
# Databases
database_mapping={}
Pathways={}
# Join graph
network={}
# The PPI
PPI={}




# Human
database_mapping['Human']= pd.read_pickle(os.path.join(os.path.dirname(__file__),"data/database/Human.pkl"))
Pathways['Human']= pd.read_pickle(os.path.join(os.path.dirname(__file__),"data/pathways/pathways_human"))

network['Human']=load_obj(os.path.join(os.path.dirname(__file__),'data/network/graph_human'))
PPI['Human']=load_obj(os.path.join(os.path.dirname(__file__),'data/network/PPI_Human'))


