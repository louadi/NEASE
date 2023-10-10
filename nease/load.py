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

elm={}
elm_interactions={}

pdb={}

here=os.path.dirname(__file__)


# Human
database_mapping['Human']= pd.read_pickle(os.path.join(here,"data/database/Human"))
Pathways['Human']= pd.read_pickle(os.path.join(here,"data/pathways/pathways_human"))
network['Human']=load_obj(os.path.join(here,'data/network/graph_human'))
PPI['Human']=load_obj(os.path.join(here,'data/network/PPI_Human'))
elm['Human']= pd.read_pickle(os.path.join(here,"data/database/elm"))
elm_interactions['Human']= pd.read_pickle(os.path.join(here,"data/database/ELM_interactions"))
pdb['Human']= pd.read_pickle(os.path.join(here,"data/database/pdb"))
non_coding=load_obj(os.path.join(here,'data/database/non_coding'))

database_mapping['Mouse']= pd.read_pickle(os.path.join(here,"data/database/Mouse"))
Pathways['Mouse']= pd.read_pickle(os.path.join(here,"data/pathways/pathways_mouse"))
network['Mouse']=load_obj(os.path.join(here,'data/network/graph_mouse'))
PPI['Mouse']=load_obj(os.path.join(here,'data/network/PPI_mouse'))
elm['Mouse']= pd.read_pickle(os.path.join(here,"data/database/elm_mouse"))
elm_interactions['Mouse']= pd.read_pickle(os.path.join(here,"data/database/ELM_interactions_mouse"))
pdb['Mouse']= pd.read_pickle(os.path.join(here,"data/database/pdb_mouse"))

