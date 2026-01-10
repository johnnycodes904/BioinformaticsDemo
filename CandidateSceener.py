import pandas as pd
import io
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski


csv_data = """Name,Class,SMILES
Aspirin,Painkiller,CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol,Painkiller,CC(=O)NC1=CC=C(O)C=C1
Ibuprofen,Anti-inflammatory,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
Amoxicillin,Antibiotic,CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C
Atorvastatin,Cholesterol,CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3NC(=O)C4=CC=CC=C4)C(=O)NC5=CC=CC=C5
Vancomycin,Antibiotic,CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)CN)C)O)Cl)CO)O)O)(C)N
Warfarin,Blood Thinner,CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O
Metformin,Diabetes,CN(C)C(=N)NC(=N)N
Cyclosporine,Immunosuppressant,CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C
Digoxin,Heart Failure,CC1C(C(CC(O1)OC2C(OC(CC2O)OC3C(OC(CC3O)OC4CCC5(C(C4)CCC6C5CCC7(C6(CCC7C8=CC(=O)OC8)O)C)O)C)C)O)O"""

#load this into Pandas
df = pd.read_csv(io.StringIO(csv_data))

#Iterate through SMILES strings and have RDKit do the physics
def get_molecular_specs(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return pd.Series({
                'Weight': Descriptors.ExactMolWt(mol),
                'LogP': Descriptors.MolLogP(mol),
                'H_Donors': Lipinski.NumHDonors(mol),
                'H_Acceptors': Lipinski.NumHAcceptors(mol)
            })
    except Exception as e:
        print(f"Error processing {smiles}: {e}")
        pass
    return pd.Series([None, None, None, None])
    
#Appy this function and display
df_specs = df['SMILES'].apply(get_molecular_specs)
final_df = pd.concat([df, df_specs], axis=1)

print(final_df[['Name', 'Weight', 'LogP', 'H_Donors']])


#Defining Lipinski Rule filter
rule_weights = final_df['Weight'] <= 500
rule_logp = final_df['LogP'] <= 5
rule_donors = final_df['H_Donors'] <= 5
rule_acceptors = final_df['H_Acceptors'] <= 10

#Creating boolean column which is TRUE if all condiditions/rules above are met
final_df['Passes_Lipinski'] = rule_weights & rule_logp & rule_donors & rule_acceptors

#Break list into winners and losers categories
winners = final_df[final_df['Passes_Lipinski'] == True]
losers = final_df[final_df['Passes_Lipinski'] == False]

print(f"Passed: {len(winners)}")
print(f"Failed: {len(losers)}")

#show the losers and which specific rule(s) failed
print("\n--- The Failures and Reason ---")
print(losers[['Name', 'Weight', 'LogP', 'Passes_Lipinski']])



#DO THE VIZ NOW

#Setup the chart/figure
plt.figure(figsize=(10, 6))

#Plot the "passing" candidates in GREEN
plt.scatter(
    winners['Weight'],
    winners['LogP'],
    color='green',
    label='Passes Rule of 5',
    s=100,
    edgecolor='black'
)

#Plot the "failing" candidates in RED
plt.scatter(
    losers['Weight'],
    losers['LogP'],
    color='red',
    label='Fails Rule of 5',
    s=100,
    edgecolor='black'
)

#labels for data points
for i, row in final_df.iterrows():
    plt.text(
        row['Weight'] + 15,
        row['LogP'],
        row['Name'],
        fontsize=9
    )

#Draw the boundaries
plt.axvline(x=500, color='gray', linestyle='--', label='Max Weight (500)')
plt.axhline(y=5, color='gray', linestyle='--', label='Max LogP (5)')

#Labels and Title
plt.xlabel('Molecular Weight (Daltons)')
plt.ylabel('LogP (Lipophilicity)')
plt.title('Drug Screening: Lipinski Rule of 5 Analysis')
plt.legend()
plt.grid(True, alpha=0.3)

plt.show()