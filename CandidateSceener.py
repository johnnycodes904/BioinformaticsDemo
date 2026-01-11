import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

#Config
INPUT_FILE = 'drugs.csv'
WEIGHT_LIMIT = 500
LOGP_LIMIT = 5

#Iterate through SMILES and calculate the physical specs of each molecule
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
        print(f"Error processing {smiles} -> {e}")
    return pd.Series([None, None, None, None])

#read from drugs.csv
def run_analysis():
    print(f"Loading data from {INPUT_FILE}...")
    try:
        df = pd.read_csv(INPUT_FILE)
    except FileNotFoundError:
        print(f"Error: Could not find {INPUT_FILE}. Is it in the same directory as this script?")
        return

#Appy this function and display
    print("Calculating, please wait....")
    df_specs = df['SMILES'].apply(get_molecular_specs)
    df = pd.concat([df, df_specs], axis=1)

#Apply Lipinski logic gates
#Create a boolean mask for each rule
    passes_filter = (
        (df['Weight'] <= WEIGHT_LIMIT) &
        (df['LogP'] <= LOGP_LIMIT) &
        (df['H_Donors'] <= 5) &
        (df['H_Acceptors'] <= 10)
    )
    df['Passes_Lipinski'] = passes_filter


#Break list into winners and losers categories
    winners = df[df['Passes_Lipinski'] == True]
    losers = df[df['Passes_Lipinski'] == False]

#Output results
    print(f"\nAnalysis Complete:")
    print(f"-> Total Compounds: {len(df)}")
    print(f"-> Qualified Candidates: {len(winners)}")
    print(f"-> Rejected: {len(losers)}")

    print("\n--- The Winners (Oral Candidates) ---")
    print(winners[['Name', 'Weight', 'LogP']])

    print("\n--- The Outliers (Injectable Candidates/Failed) ---")
    print(losers[['Name', 'Weight', 'LogP']])

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

#Draw the boundaries
    plt.axvline(x=WEIGHT_LIMIT, color='gray', linestyle='--', alpha=0.5, label='Max Weight')
    plt.axhline(y=LOGP_LIMIT, color='gray', linestyle='--', alpha=0.5, label='Max LogP')

#labels for data points offset slightly (+15)
    for i, row in df.iterrows():
        plt.text(
            row['Weight'] + 15,
            row['LogP'],
            row['Name'],
            fontsize=9
        )

#Labels and Title
    plt.xlabel('Molecular Weight (Daltons)')
    plt.ylabel('LogP (Lipophilicity)')
    plt.title('Drug Screening: Lipinski Rule of 5 Analysis')
    plt.legend()
    plt.grid(True, alpha=0.3)

    print("\nDisplaying plot...")
    plt.show()

if __name__ == "__main__":
    run_analysis()