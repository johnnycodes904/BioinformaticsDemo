Drug Candidate Screener (Lipinski Filter)
Automated Screening of Chemical Compounds for Oral Bioavailability

Project Status: Done
Complete Tech Stack: Python, Pandas, RDKit, Matplotlib

About The Project

In the pharmaceutical industry, finding a potential drug molecule is only half the battle. The molecule must also be "deliverable";  can it survive the human digestive system and enter the bloodstream?

This project implements a Computational Drug Discovery pipeline that screens raw chemical data against Lipinski's Rule of 5, the industry standard for oral drug properties. It acts as a digital filter, accepting thousands of chemical structures (SMILES strings) and rejecting those that are too heavy, too oily, or too chemically reactive to be effective as pills.
The "Hardware" Logic (Lipinski's Rule)

Just as PC components must fit specific form factors (ATX, ITX, etc.), drug molecules must fit specific biological constraints to be orally active. This script calculates "specs" for every molecule:

    Molecular Weight: Must be ≤ 500 Daltons (Molecules that are too heavy have poor absorption).

    LogP (Lipophilicity): Must be ≤ 5 (Molecules that are too oily get stuck in fat. Conversely, too watery are flushed out).

    Hydrogen Bond Donors: Must be ≤ 5.

    Hydrogen Bond Acceptors: Must be ≤ 10.

How It Works

    Ingestion: Loads a dataset of chemical compounds via their SMILES strings (Simplified Molecular Input Line Entry System).

    Calculation: Uses RDKit (Cheminformatics library) to compute the physicochemical properties of each molecule in 3D space.

    Filtration: Applies a Boolean mask to the dataframe, separating "Pass" and "Fail" candidates.

    Visualization: Generates a scatter plot mapping Molecular Weight vs. LogP, to identify clusters of viable drug candidates.

Key Insights & Analysis

The script correctly identifies the split between "small molecule" drugs and specialized treatments.

    The "Successes": Common household drugs like Aspirin, Tylenol, and Ibuprofen cluster in the bottom-left quadrant (representing Low Weight and Low LogP), validating why they work oral tablets.

    The "Failures": The script correctly flagged Vancomycin and Cyclosporine as failing the Lipinski rules.

        Bio-Insight: This is a True Negative. Vancomycin is a massive molecule (~1449 Da) used for severe infections. Because it violates these rules, it cannot be absorbed by the gut and must be administered via IV injection. This demonstrates that "failing" the filter only predicts the administration method, not necessarily the drug's potency.
