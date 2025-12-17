def score_lead(profile, enrichment, company_data):
    """
    Calculates the Propensity-to-Buy probability score (0-100).
    
    Args:
        profile (dict): LinkedIn profile data
        enrichment (dict): PubMed enrichment data
        company_data (dict): Company funding and technographic data
        
    Returns:
        dict: {
            "score": int,
            "score_breakdown": dict
        }
    """
    score = 0
    breakdown = {
        "role_fit": 0,
        "company_intent": 0,
        "technographic": 0,
        "location": 0,
        "scientific_intent": 0
    }
    
    # Normalize inputs
    title = str(profile.get("title", "")).lower()
    company_hq = str(profile.get("company_hq", "")).lower()
    
    funding = str(company_data.get("Funding_Stage", "")).lower()
    technographic = str(company_data.get("Technographic_Fit", "")).lower()
    
    has_papers = enrichment.get("scientific_intent", False)
    
    # 1. Role Fit (+30)
    # Criteria: Title contains Toxicology, Safety, Hepatic, 3D
    if any(k in title for k in ["toxicology", "safety", "hepatic", "3d"]):
        score += 30
        breakdown["role_fit"] = 30
        
    # 2. Company Intent (+20)
    # Criteria: Series A/B / IPO
    if any(k in funding for k in ["series a", "series b", "ipo"]):
        score += 20
        breakdown["company_intent"] = 20
        
    # 3. Technographic (+15)
    # Criteria: Company already uses similar tech (High) or open to NAMs (Medium)
    if technographic == "high":
        score += 15
        breakdown["technographic"] = 15
    elif technographic == "medium":
        score += 15 # Requirement says +15 for "similar tech" (High) and +15?? "Medium(+15)" is listed twice in table. 
                    # Let's interpret "similar tech" as +15. "NAMs" as +15.
                    # Wait, prompt says: "Technographic | Company already uses... | Medium(+15)" and "open to NAMs | Medium(+10)". 
                    # Ah, I see: "Medium(+15)" is in the example column. 
                    # Let's stick to max +15 for Technographic.
        breakdown["technographic"] = 15
    elif technographic == "potential": # If we used that key
        score += 10
        breakdown["technographic"] = 10
        
    # 4. Location (+10)
    # Criteria: Hub (Boston/Cambridge, Bay Area, Basel, UK Golden Triangle)
    bio_hubs = ["boston", "cambridge", "bay area", "basel", "uk", "london", "oxford"]
    if any(hub in company_hq for hub in bio_hubs):
        score += 10
        breakdown["location"] = 10
        
    # 5. Scientific Intent (+40)
    # Criteria: Published a paper on keywords in last 2 years
    if has_papers:
        score += 40
        breakdown["scientific_intent"] = 40
        
    # Normalize to 0-100
    # Max possible raw score = 30 + 20 + 15 + 10 + 40 = 115
    final_score = min(score, 100)
    
    return {
        "score": final_score,
        "score_breakdown": breakdown
    }
