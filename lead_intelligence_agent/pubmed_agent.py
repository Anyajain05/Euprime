import time
from Bio import Entrez
from datetime import datetime, timedelta

# Set email for Entrez (Required by NCBI)
Entrez.email = "demo_user@leadintelligence.com"

KEYWORDS = [
    "Drug-Induced Liver Injury",
    "Hepatic Toxicity",
    "3D Cell Culture",
    "Organ-on-Chip",
    "In-vitro Toxicology"
]

def search_pubmed_for_lead(name):
    """
    Searches PubMed for publications by the author containing specific keywords
    within the last 2 years.
    Returns: (count, score)
    """
    # Clean name: remove "Dr. " and split
    clean_name = name.replace("Dr. ", "").strip()
    
    # Construct query: (Name[Author]) AND (Keyword1 OR Keyword2 ...) AND date range
    # Note: Precise author matching is hard without affiliation, relying on keywords for relevance context.
    
    keyword_query = " OR ".join([f'"{k}"' for k in KEYWORDS])
    
    # Date range: Last 24 months
    # PubMed date format YYYY/MM/DD
    two_years_ago = datetime.now() - timedelta(days=730)
    date_str = two_years_ago.strftime("%Y/%m/%d")
    
    query = f'("{clean_name}"[Author]) AND ({keyword_query}) AND ("{date_str}"[Date - Publication] : "3000"[Date - Publication])'
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        
        # Scoring Rule: "PubMed paper in last 2 years" -> +40
        # We can implement a binary score or scaled. The prompt assumes a flat +40 if present?
        # "Scientific Intent | PubMed paper in last 2 years | +40"
        score = 40 if count > 0 else 0
        
        # Sleep to be nice to API
        time.sleep(0.34) # Max 3 req/sec without API key
        
        return count, score
        
    except Exception as e:
        print(f"Error searching for {name}: {e}")
        return 0, 0

if __name__ == "__main__":
    # Test
    print(search_pubmed_for_lead("Dr. Sarah Chen"))
