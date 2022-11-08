import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import AllChem
import pickle
from stmol import showmol
import py3Dmol
import requests
from streamlit_lottie import st_lottie

#==========================================================================================================================================
#==========================================================================================================================================

# load the model with pickle
predicting_model = pickle.load(open('website/stack_model_pickle.pickle', 'rb'))

# fuction to generate the descriptors
@st.cache
def generate_descriptors(*smiles):
  descriptors = pd.DataFrame(columns=["MolWt","MolLogP","MolMR","HeavyAtomCount","NumHAcceptors","NumHDonors","NumHeteroatoms",
                 "NumRotatableBonds","NumValenceElectrons","NumAromaticRings","NumSaturatedRings",
                 "NumAliphaticRings","RingCount","TPSA","LabuteASA","BalabanJ","BertzCT"])
  
  for smile in smiles: 
    if Chem.MolFromSmiles(smile) is None:
      pass

    else:
        desc_MolWt = Descriptors.MolWt(Chem.MolFromSmiles(smile))
        desc_MolLogP = Descriptors.MolLogP(Chem.MolFromSmiles(smile))
        desc_MolMR = Descriptors.MolMR(Chem.MolFromSmiles(smile))
        desc_HeavyAtomCount = Descriptors.HeavyAtomCount(Chem.MolFromSmiles(smile))
        desc_NumHAcceptors = Descriptors.NumHAcceptors(Chem.MolFromSmiles(smile))
        desc_NumHDonors = Descriptors.NumHDonors(Chem.MolFromSmiles(smile))
        desc_NumHeteroatoms = Descriptors.NumHeteroatoms(Chem.MolFromSmiles(smile))
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(Chem.MolFromSmiles(smile))
        desc_NumValenceElectrons = Descriptors.NumValenceElectrons(Chem.MolFromSmiles(smile))           
        desc_NumAromaticRings = Descriptors.NumAromaticRings(Chem.MolFromSmiles(smile))      
        desc_NumSaturatedRings = Descriptors.NumSaturatedRings(Chem.MolFromSmiles(smile))      
        desc_NumAliphaticRings = Descriptors.NumAliphaticRings(Chem.MolFromSmiles(smile))
        desc_RingCount = Descriptors.RingCount(Chem.MolFromSmiles(smile))
        desc_TPSA = Descriptors.TPSA(Chem.MolFromSmiles(smile))
        desc_LabuteASA = Descriptors.LabuteASA(Chem.MolFromSmiles(smile))       
        desc_BalabanJ = Descriptors.BalabanJ(Chem.MolFromSmiles(smile))
        desc_BertzCT = Descriptors.BertzCT(Chem.MolFromSmiles(smile))
                
        row = [desc_MolWt,desc_MolLogP,desc_MolMR,desc_HeavyAtomCount,desc_NumHAcceptors,desc_NumHDonors,desc_NumHeteroatoms,
                        desc_NumRotatableBonds,desc_NumValenceElectrons,desc_NumAromaticRings,desc_NumSaturatedRings,
                        desc_NumAliphaticRings,desc_RingCount,desc_TPSA,desc_LabuteASA,desc_BalabanJ,desc_BertzCT]
        descriptors.loc[len(descriptors)] = row

  return descriptors


# function to predict the solubility
@st.cache
def SolPred(*smiles):
  solpred_data = pd.DataFrame(columns=['SMILES', 'Chemical_Formula', 
  'Predicited_LogS', 'Predicited_Solubility (g/l)','MolWt', 'NumRotatableBonds', 'NumHAcceptors', 'NumHDonors', 
  'NumHeteroatoms', 'NumValenceElectrons', 'NumAromaticRings',
  'NumAliphaticRings', 'RingCount'
  ])
  for smile in smiles: 
    if Chem.MolFromSmiles(smile) is None:
      solpred_data.loc[len(solpred_data)] = [smile] + [None for i in range(len(solpred_data.columns)-1)]
    else:
      solpred_data.loc[len(solpred_data)] = [
                            smile, 
                            CalcMolFormula(Chem.MolFromSmiles(smile)), 
                            round(predicting_model.predict(generate_descriptors(smile))[0], 2),
                            (10 ** (predicting_model.predict(generate_descriptors(smile))[0])) * generate_descriptors(smile)['MolWt'][0],
                            generate_descriptors(smile)['MolWt'][0],
                            generate_descriptors(smile)['NumRotatableBonds'][0],
                            generate_descriptors(smile)['NumHAcceptors'][0],
                            generate_descriptors(smile)['NumHDonors'][0],
                            generate_descriptors(smile)['NumHeteroatoms'][0],
                            generate_descriptors(smile)['NumValenceElectrons'][0],
                            generate_descriptors(smile)['NumAromaticRings'][0],
                            generate_descriptors(smile)['NumAliphaticRings'][0],
                            generate_descriptors(smile)['RingCount'][0],
                            ]       

  return solpred_data

#==========================================================================================================================================
#==========================================================================================================================================
# Streamlit page configuration
st.set_page_config(
    page_title="SolPred",
    page_icon="website/chemical.png",
    layout="wide",
    initial_sidebar_state='collapsed',
)
# remove menu and the footer
st.markdown("""
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>""", unsafe_allow_html=True)


#==========================================================================================================================================
#==========================================================================================================================================
# The sidebar configuration
sidebar = st.sidebar
sidebar.title("S O L P R E D")
sidebar.subheader('Input SMILES string')
sidebar.caption('In case of multiple SMILES string, please separate them with a comma')


smiles_input = [smile.strip() for smile in sidebar.text_input(label='SMILES', label_visibility='collapsed', value='CN1C=NC2=C1C(=O)N(C)C(=O)N2C').split(',')]
smiles_input = list(filter(None, smiles_input))


predict_button_sidebar = sidebar.button('Predict')
sidebar.caption('Check the box below to view the 3D structure of the molecule')
show_mol_button = st.sidebar.checkbox('Show 3D structure')


sidebar.subheader('You can also upload a CSV file')
sidebar.caption('The file should have a column named "SMILES"')
upload_button = sidebar.file_uploader("Upload CSV", label_visibility='collapsed', type=['csv', 'txt'])
st.sidebar.markdown('''
---
Created with ❤️ by [Shalash](https://twitter.com/__shalash__).

You can find the source code on my [GitHub](https://github.com/Shalash96/SolPred)
''')


#==========================================================================================================================================
#==========================================================================================================================================
# The input section configuration

if 'predict_button_sidebar' not in st.session_state:
  st.session_state.predict_button_sidebar = False

if predict_button_sidebar or st.session_state.predict_button_sidebar:
  st.session_state.predict_button_sidebar = True
  if len(smiles_input) > 0:
    st.header('Predicted data of the input SMILES')
    predicted_values_from_input = SolPred(*smiles_input)
    st.dataframe(predicted_values_from_input)
    download_button = st.download_button("Download as CSV", data=predicted_values_from_input.to_csv(index=False).encode('utf-8'), file_name='SolPred(@__shalash__).csv', mime='text/csv')
  
  else:
    st.error('Please enter a valid SMILES string in the text input in the sidebar to predict the solubility or to view the 3D structure of the molecule')
#==========================================================================================================================================
#==========================================================================================================================================
# The molecule viewer configuration

# the right column configuration

if 'show_mol_button' not in st.session_state:
  st.session_state.show_mol_button = False

if show_mol_button or st.session_state.show_mol_button:
  st.session_state.show_mol_button = True

left , right = st.columns(2)
with right:
  if show_mol_button:
    st.session_state.show_mol_button = True
    if len(smiles_input) > 0:
      mol_number = st.number_input('Index of the molecule to view', min_value=0, max_value=len(smiles_input)-1, value=0, step=1, key='mol_index')
      style = st.selectbox('The prefered style to show the molecule',['stick','line','sphere', 'cross'])
      color_bg = st.color_picker('Background color', '#FFFFFF')
      spin = st.checkbox('Spin', value=False)
      
    else:
      style = None
      color_bg = None
      mol_number = None
      spin = None

#==========================================================================================================================================
# the left column configuration

with left:
  if show_mol_button:
    st.session_state.show_mol_button = True
    def mol_viewer(smile, style=style, spin=spin):
      mol = Chem.MolFromSmiles(smile)
      mol = Chem.AddHs(mol)
      AllChem.EmbedMolecule(mol)
      mblock = Chem.MolToMolBlock(mol)

      xyzview = py3Dmol.view()
      xyzview.addModel(mblock,'mol')
      xyzview.setStyle({style:{}})
      xyzview.setBackgroundColor(color_bg)
      if spin:
        xyzview.spin()

      xyzview.zoomTo()
      # xyzview.png()
      showmol(xyzview,height=500,width=350)
    try:
      mol_viewer(smiles_input[mol_number])
    except:
      pass


#==========================================================================================================================================
#==========================================================================================================================================
# The upload section configuration
    
if 'upload_button' not in st.session_state:
  st.session_state.upload_button = False
if upload_button or st.session_state.upload_button:
  st.header('Predicted data of the uploaded SMILES')
  st.session_state.upload_button = True
  try:
    uploaded_file = pd.read_csv(upload_button)
    if 'SMILES' in uploaded_file.columns:
      smiles_from_upload = uploaded_file['SMILES'].tolist()
      predicted_values_from_upload = SolPred(*smiles_from_upload)
      download_button = st.download_button("Download as CSV", data=predicted_values_from_upload.to_csv(index=False).encode('utf-8'), file_name='SolPred(@__shalash__).csv', mime='text/csv')
      st.dataframe(predicted_values_from_upload.style.format(precision=2))

    else:
      st.error('SMILES column not found. Please check your CSV file and make sure that smiles are in a column named "SMILES".')
  except:
    st.error('The file uploaded has been removed')

#==========================================================================================================================================
#==========================================================================================================================================
# The start up configuration

if not st.session_state.predict_button_sidebar and not st.session_state.show_mol_button and not st.session_state.upload_button:
  url = 'https://assets9.lottiefiles.com/packages/lf20_spxbu9ky.json'

  def load_lottieurl(url):
    data = requests.get(url)
    if data.status_code == 200:
      return data.json()

  st.markdown("""
  
  # Welcome to SolPred
  The name derived from the word "Solubility" and "Prediction"
  ## It is a web app to predict the solubility of molecules using machine learning models. It is easy to use just use the SMILES string.

  ### How to use it?
  First open the sidebar from the left side by clicking on the small arrow at the top of the page then you will be able to:
  1. Enter the SMILES string of the molecule in the text input in the sidebar then click on the "Predict" button. There is a molecule already entered as an example.
  2. You can also upload a CSV file containing the SMILES strings of the molecules in a column named "SMILES". 
  3. You can also view the 3D structure of the molecule by clicking on the "Show 3D structure" button in the sidebar.
  4. Enjoy using the app and don't forget to share it with your friends and colleagues. All feedbacks are highly appreciated.

  """)

  res_lottie = load_lottieurl(url)
  st_lottie(res_lottie, quality='high', height=400, speed=1, key='lottie')
