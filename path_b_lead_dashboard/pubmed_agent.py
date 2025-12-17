import time
from Bio import Entrez
from datetime import datetime, timedelta

# Email required by NCBI for API usage
Entrez.email = "demo_user@leadintelligence.com"

KEYWORDS = [
    "Drug-Induced Liver Injury",
    "Hepatic Toxicity",
    "3D Cell Culture",
    "Organ-on-Chip",
    "Hepatic spheroids",
    "Investigative Toxicology"
]

def get_scientific_enrichment(name):
    """
    Queries PubMed for publications by the author containing specific keywords
    within the last 24 months.
    
    Args:
        name (str): Author Name (e.g., "Dr. Sarah Miller")
    
    Returns:
        dict: {
            "pubmed_paper_count": int,
            "scientific_intent": bool
        }
    """
    clean_name = name.replace("Dr. ", "").strip()
    
    # Construct Query
    # (Name[Author]) AND (Keyword1 OR Keyword2 ...) AND date range
    keyword_query = " OR ".join([f'"{k}"' for k in KEYWORDS])
    
    # Date logic: Last 24 months (730 days)
    two_years_ago = datetime.now() - timedelta(days=730)
    date_str = two_years_ago.strftime("%Y/%m/%d")
    
    query = f'("{clean_name}"[Author]) AND ({keyword_query}) AND ("{date_str}"[Date - Publication] : "3000"[Date - Publication])'
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        
        # Be nice to the API
        time.sleep(0.34)
        
        return {
            "pubmed_paper_count": count,
            "scientific_intent": count > 0
        }
        
    except Exception as e:
        print(f"Error searching PubMed for {name}: {e}")
        return {
            "pubmed_paper_count": 0,
            "scientific_intent": False
        }

if __name__ == "__main__":
    print(get_scientific_enrichment("Dr. Sarah Miller"))
