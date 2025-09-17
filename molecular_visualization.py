# Import libraries
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import py3Dmol
import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
import pandas as pd
import os

# Step 1: Define the dataset (common electrolyte solvents)
smiles_data = [
    {'name': 'EC', 'smiles': 'O=C1OCCO1'},
    {'name': 'PC', 'smiles': 'CC1COC(=O)O1'},
    {'name': 'DMC', 'smiles': 'COC(=O)OC'},
    {'name': 'DEC', 'smiles': 'CCOC(=O)OCC'},
    {'name': 'EMC', 'smiles': 'CCOC(=O)OC'},
    {'name': 'VC', 'smiles': 'C1=COC(=O)O1'},
    {'name': 'FEC', 'smiles': 'FC1COC(=O)O1'},
    {'name': 'DME', 'smiles': 'COCCOC'},
    {'name': 'Diglyme', 'smiles': 'COCCOCCOC'},
    {'name': 'THF', 'smiles': 'C1CCOC1'},
    {'name': 'DOL', 'smiles': 'C1COCO1'},
    {'name': 'GBL', 'smiles': 'C1CC(=O)OC1'},
    {'name': 'AN', 'smiles': 'CC#N'},
    {'name': 'SL', 'smiles': 'C1CCS(=O)(=O)C1'},
    {'name': 'EA', 'smiles': 'CCOC(=O)C'}
]

# Step 2: Compute properties
molecules = []
mw = []
for entry in smiles_data:
    mol = Chem.MolFromSmiles(entry['smiles'])
    if mol is None:
        print(f'Invalid SMILES for {entry["name"]}, skipping.')
        continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coords
    AllChem.MMFFOptimizeMolecule(mol)  # Optimize geometry
    mw.append(rdMolDescriptors.CalcExactMolWt(mol))  # Molecular weight
    molecules.append(mol)
    entry['mol'] = mol
    entry['mw'] = mw[-1]
    print(f'{entry["name"]}: MW = {mw[-1]:.2f}')

# Save data to CSV
df = pd.DataFrame(smiles_data, columns=['name', 'smiles', 'mw'])
df.to_csv('electrolyte_properties.csv', index=False)
print('Data saved to electrolyte_properties.csv')

# Step 3: Visualize molecules in 3D (save as HTML)
view = py3Dmol.view(width=400, height=400)
for mol in molecules[:5]:  # Show first 5 for demo
    mb = Chem.MolToMolBlock(mol)
    view.addModel(mb, 'mol')
view.setStyle({'stick': {}})
view.zoomTo()
with open('3d_molecules.html', 'w') as f:
    f.write(view._make_html())
print('3D visualization saved to 3d_molecules.html (open in browser)')

# Save 2D images for all molecules
img = Draw.MolsToGridImage(molecules, molsPerRow=5, subImgSize=(200, 200))
img.save('molecule_2d_grid.png')
print('2D grid image saved to molecule_2d_grid.png')

# Step 4: Correlate properties (using MW only for now)
plt.figure(figsize=(8, 6))
plt.scatter(range(len(mw)), mw, color='blue')
for i, name in enumerate([d['name'] for d in smiles_data if d['mol'] is not None]):
    plt.annotate(name, (i, mw[i]), fontsize=8)
plt.xlabel('Molecule Index')
plt.ylabel('Molecular Weight (g/mol)')
plt.title('Molecular Weight of Electrolytes')
plt.grid(True)
plt.savefig('mw_plot.png')
plt.close()

# Step 5: Streamlit Interface
st.title('Molecular Visualization and Property Correlation Tool')
st.write('This tool visualizes and correlates properties of electrolyte molecules for lithium-ion battery research.')

# Display dataset
st.subheader('Electrolyte Properties')
st.write(df)

# Display 2D grid image
st.subheader('2D Molecular Structures')
st.image('molecule_2d_grid.png', caption='2D Grid of Electrolyte Molecules', use_container_width=True)

# Display MW plot
st.subheader('Molecular Weight Plot')
st.image('mw_plot.png', caption='Molecular Weight vs. Index', use_container_width=True)

# Link to 3D visualization
st.subheader('3D Molecular Visualization')
st.write('Click [here](3d_molecules.html) to view the 3D structures in your browser (requires opening the local file).')
if os.path.exists('3d_molecules.html'):
    st.success('3D file generated successfully!')
else:
    st.error('3D file generation failed.')

# Add interactivity (optional simple feature)
st.subheader('Add a New Molecule')
new_smiles = st.text_input('Enter SMILES string (e.g., O=C1OCCO1 for EC)')
if st.button('Add and Visualize'):
    new_mol = Chem.MolFromSmiles(new_smiles)
    if new_mol:
        new_mol = Chem.AddHs(new_mol)
        AllChem.EmbedMolecule(new_mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(new_mol)
        new_mw = rdMolDescriptors.CalcExactMolWt(new_mol)
        st.write(f'Molecular Weight: {new_mw:.2f}')
        img = Draw.MolToImage(new_mol)
        st.image(img, caption=f'2D Structure of {new_smiles}', use_container_width=True)
    else:
        st.error('Invalid SMILES string!')

# Save updated data if new molecule added
if 'new_mol' in locals() and new_mol:
    new_entry = {'name': f'Custom_{len(smiles_data)}', 'smiles': new_smiles, 'mw': new_mw}
    smiles_data.append(new_entry)
    df = pd.DataFrame(smiles_data, columns=['name', 'smiles', 'mw'])
    df.to_csv('electrolyte_properties.csv', index=False)
    st.success('New molecule added to dataset!')