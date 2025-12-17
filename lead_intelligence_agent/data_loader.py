import pandas as pd
import os
from pubmed_agent import search_pubmed_for_lead
from scorer import calculate_score

DATA_DIR = os.path.join(os.path.dirname(__file__), 'sample_data')
LEADS_FILE = os.path.join(DATA_DIR, 'leads_input.csv')
FUNDING_FILE = os.path.join(DATA_DIR, 'funding_data.csv')

def load_and_enrich_data():
    """
    Orchestrates the data pipeline:
    1. Load Leads and Funding data
    2. Merge on Company
    3. Enrich with PubMed API (Scientific Intent)
    4. Calculate Probability Scores
    5. Rank and Sort
    """
    print("Loading data...")
    leads_df = pd.read_csv(LEADS_FILE)
    funding_df = pd.read_csv(FUNDING_FILE)
    
    # Merge Funding Data
    # Assuming exact company name match for simplicity in this prototype
    merged_df = pd.merge(leads_df, funding_df, on='Company', how='left')
    
    # Fill missing funding with "Unknown" or handle gracefully
    merged_df['Funding_Stage'] = merged_df['Funding_Stage'].fillna('Unknown')
    
    print("Enriching with PubMed data (this may take a moment)...")
    # Apply PubMed Enrichment
    # We use apply but to be mindful of API rate limits/performance in a demo, 
    # we iterate.
    
    pubmed_counts = []
    pubmed_scores = []
    
    for index, row in merged_df.iterrows():
        name = row['Name']
        print(f"Querying PubMed for: {name}")
        count, score = search_pubmed_for_lead(name)
        pubmed_counts.append(count)
        pubmed_scores.append(score)
        
    merged_df['PubMed_Papers_Count'] = pubmed_counts
    merged_df['Scientific_Intent_Score'] = pubmed_scores
    
    # Calculate Score
    print("Calculating scores...")
    merged_df['Probability_Score'] = merged_df.apply(calculate_score, axis=1)
    
    # Generate Simulated Email
    # Format: first.last@company.com (simplified)
    merged_df['Email'] = merged_df.apply(lambda x: f"{x['Name'].split()[0].lower()}.{x['Name'].split()[-1].lower()}@{x['Company'].replace(' ', '').lower()}.com", axis=1)
    
    # Rank by Score Descending
    merged_df = merged_df.sort_values(by='Probability_Score', ascending=False)
    
    # Add Rank Column
    merged_df.reset_index(drop=True, inplace=True)
    merged_df.index += 1
    merged_df.insert(0, 'Rank', merged_df.index)
    
    return merged_df

if __name__ == "__main__":
    df = load_and_enrich_data()
    print(df.head())
    df.to_csv("processed_leads.csv", index=False)
