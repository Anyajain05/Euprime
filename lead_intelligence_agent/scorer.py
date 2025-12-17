def calculate_score(row):
    """
    Calculates the probability score based on the strict rules provided.
    
    Signals:
    1. Role Fit: Title contains Toxicology / Safety / Hepatic (+30)
    2. Company Intent: Series A/B or IPO (+20)
    3. Technographic: Mentions in-vitro / NAMs (+15)
    4. Location: Bio hub HQ (+10)
    5. Scientific Intent: PubMed paper in last 2 years (+40)
    
    Returns: Raw Score and Normalized Score (0-100)
    """
    score = 0
    signals = []
    
    title = str(row.get('Title', '')).lower()
    company_hq = str(row.get('Company_HQ', '')).lower()
    funding = str(row.get('Funding_Stage', '')).lower()
    pubmed_score = row.get('Scientific_Intent_Score', 0)
    
    # 1. Role Fit (+30)
    if any(keyword in title for keyword in ['toxicology', 'safety', 'hepatic']):
        score += 30
        signals.append("Role Fit")

    # 2. Company Intent (+20)
    # "Series A", "Series B", "IPO"
    # Matches "series a", "series b", "ipo"
    if funding in ['series a', 'series b', 'ipo']:
        score += 20
        signals.append("Company Funding")

    # 3. Technographic (+15)
    # Mentions in-vitro / NAMs (Checking Title for this prototype as a proxy)
    if any(keyword in title for keyword in ['in-vitro', 'nams', '3d', 'organ-on-chip']):
        score += 15
        signals.append("Technographic")

    # 4. Location (+10)
    # Bio hub HQ: Boston / Cambridge, Bay Area, Basel, UK Golden Triangle
    bio_hubs = ['boston / cambridge', 'bay area', 'basel', 'uk golden triangle']
    if company_hq in bio_hubs:
        score += 10
        signals.append("Location Match")
        
    # 5. Scientific Intent (+40)
    # Passed from PubMed Agent
    if pubmed_score > 0:
        score += 40
        signals.append("Scientific Intent")
        
    # Normalize to 0-100
    # Max possible score = 30 + 20 + 15 + 10 + 40 = 115
    # Let's cap at 100 or just normalize standardly: (score / 115) * 100
    # Or just return raw score if it's treated as probability?
    # Prompt says "Normalize final score to 0–100".
    
    max_score = 115
    normalized_score = min(int((score / max_score) * 100), 100)
    
    return normalized_score
