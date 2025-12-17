import time

def fetch_profile(linkedin_url):
    """
    Simulates a Proxycurl API call to fetch LinkedIn profile data.
    In a real scenario, this would call https://nubela.co/proxycurl/api/v2/linkedin
    
    Args:
        linkedin_url (str): The LinkedIn profile URL.
        
    Returns:
        dict: Structured profile data.
    """
    
    # MOCK DATA STORE (Simulating API Database)
    # Mapping specific URLs to realistic personas for the demo.
    mock_db = {
        "https://linkedin.com/in/john-doe-nema": {
            "name": "Dr. John Doe",
            "title": "Head of Computational Biology",
            "company": "Nema",
            "person_location": "New York, NY",
            "company_hq": "New York, NY",
            "email": "john.doe@nema.life",
            "phone": "+1-212-555-0101",
            "linkedin_url": "https://linkedin.com/in/john-doe-nema"
        },
        "https://linkedin.com/in/jane-smith-iambic": {
            "name": "Jane Smith",
            "title": "Director of Generative Chemistry",
            "company": "Iambic",
            "person_location": "San Diego, CA",
            "company_hq": "San Diego, CA",
            "email": "jane@iambic.ai",
            "phone": "+1-619-555-0102",
            "linkedin_url": "https://linkedin.com/in/jane-smith-iambic"
        },
        "https://linkedin.com/in/michael-wang-airops": {
            "name": "Michael Wang",
            "title": "VP Engineering",
            "company": "AirOps",
            "person_location": "Miami, FL",
            "company_hq": "Miami, FL",
            "email": "m.wang@airops.com",
            "phone": "+1-305-555-0103",
            "linkedin_url": "https://linkedin.com/in/michael-wang-airops"
        },
        "https://linkedin.com/in/sarah-jones-gamma": {
            "name": "Dr. Sarah Jones",
            "title": "Chief Scientific Officer",
            "company": "Gamma",
            "person_location": "San Francisco, CA",
            "company_hq": "San Francisco, CA",
            "email": "sarah@gamma.app",
            "phone": "+1-415-555-0104",
            "linkedin_url": "https://linkedin.com/in/sarah-jones-gamma"
        },
        "https://linkedin.com/in/robert-lee-colab": {
            "name": "Robert Lee",
            "title": "Head of Partnerships",
            "company": "CoLab",
            "person_location": "Vancouver, Canada",
            "company_hq": "Vancouver, Canada",
            "email": "r.lee@colabsoftware.com",
            "phone": "+1-604-555-0105",
            "linkedin_url": "https://linkedin.com/in/robert-lee-colab"
        },
        "https://linkedin.com/in/emily-white-genspark": {
            "name": "Emily White",
            "title": "Lead AI Researcher",
            "company": "Genspark",
            "person_location": "Palo Alto, CA",
            "company_hq": "Palo Alto, CA",
            "email": "emily@genspark.ai",
            "phone": "+1-650-555-0106",
            "linkedin_url": "https://linkedin.com/in/emily-white-genspark"
        },
        "https://linkedin.com/in/david-kim-normai": {
            "name": "David Kim",
            "title": "Director of Toxicology",
            "company": "Norm AI",
            "person_location": "Boston, MA",
            "company_hq": "Boston, MA", 
            "email": "david@norm.ai",
            "phone": "+1-617-555-0107",
            "linkedin_url": "https://linkedin.com/in/david-kim-normai"
        }
    }
    
    # Simulate Network Latency
    time.sleep(0.5)
    
    # Return mock data if exists, else generic fallback
    if linkedin_url in mock_db:
        return mock_db[linkedin_url]
    else:
        # SMART FALLBACK for Demo Mode
        # Extract name from URL: linkedin.com/in/akash-gupta-123 -> "Akash Gupta"
        slug = linkedin_url.strip('/').split('/')[-1] # "akash-gupta-123"
        # Remove trailing numbers often found in public URLs
        name_parts = [p.capitalize() for p in slug.split('-') if p.isalpha()]
        formatted_name = " ".join(name_parts)
        
        if not formatted_name:
            formatted_name = "New User"

        # Specific Logic for User (Easter Egg)
        if "akash" in slug.lower():
            return {
                "name": "Akash Gupta",
                "title": "CEO / Founder",
                "company": "EuPrime",
                "person_location": "San Francisco, CA",
                "company_hq": "San Francisco, CA",
                "email": "akash@euprime.org",
                "phone": "+1-415-555-0000",
                "linkedin_url": linkedin_url
            }

        # Generic Smart Simulation
        return {
            "name": formatted_name,
            "title": "Simulated Professional",
            "company": "External Lead Inc.",
            "person_location": "Unknown Location",
            "company_hq": "Unknown HQ",
            "email": f"{slug.split('-')[0]}@simulated-email.com",
            "phone": "N/A",
            "linkedin_url": linkedin_url
        }

if __name__ == "__main__":
    # Test
    print(fetch_profile("https://linkedin.com/in/sarahm"))
