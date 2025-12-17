import streamlit as st
import pandas as pd
import os
from linkedin_agent import fetch_profile
from pubmed_agent import get_scientific_enrichment
from scorer import score_lead

# Configuration
st.set_page_config(
    page_title="3D In-Vitro Lead Qualification Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Dark Theme Polish
st.markdown("""
<style>
    [data-testid="stHeader"] {
        background-color: rgba(0,0,0,0);
    }
    .main .block-container {
        padding-top: 2rem;
    }
    footer {
        visibility: hidden;
    }
    .footer-text {
        text-align: center;
        color: #888;
        font-size: 0.8rem;
        margin-top: 50px;
    }
</style>
""", unsafe_allow_html=True)

# Constants
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
URLS_FILE = os.path.join(DATA_DIR, 'linkedin_urls.txt')
FUNDING_FILE = os.path.join(DATA_DIR, 'funding_data.csv')

# Load Funding Data
def load_company_data():
    if os.path.exists(FUNDING_FILE):
        return pd.read_csv(FUNDING_FILE)
    return pd.DataFrame(columns=['Company', 'Domain', 'Amount', 'Round', 'Investors', 'Technographic_Fit'])

COMPANY_DF = load_company_data()

def get_company_info(company_name):
    """Retrieves funding and technographic data."""
    match = COMPANY_DF[COMPANY_DF['Company'].str.lower() == company_name.lower()]
    if not match.empty:
        return match.iloc[0].to_dict()
    # Default fallback
    return {
        "Domain": "N/A",
        "Amount": "Unknown",
        "Round": "Unknown",
        "Investors": "Unknown",
        "Technographic_Fit": "Low"
    }

def process_leads(urls):
    """
    Orchestrates the fetching, enrichment, and scoring pipeline.
    """
    results = []
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total = len(urls)
    
    for i, url in enumerate(urls):
        if not url.strip():
            continue
            
        status_text.text(f"Processing ({i+1}/{total}): {url}")
        
        # 1. Identification
        profile = fetch_profile(url.strip())
        
        # 2. Scientific Enrichment
        enrichment = get_scientific_enrichment(profile['name'])
        
        # 3. Business Enrichment
        company_info = get_company_info(profile['company'])
        
        # 4. Scoring
        score_data = score_lead(profile, enrichment, company_info)
        
        # Logic for Work Mode
        p_loc = profile['person_location'].lower()
        c_hq = profile['company_hq'].lower()
        is_remote = False
        
        if "remote" in p_loc:
            is_remote = True
        elif not any(part.strip() in c_hq for part in p_loc.split(',')):
            if "bay area" in c_hq and "san francisco" in p_loc:
                is_remote = False
            elif "uk" in c_hq and ("london" in p_loc or "cambridge" in p_loc):
                is_remote = False
            elif "cambridge" in c_hq and "boston" in p_loc:
                is_remote = False
            elif c_hq in p_loc:
                 is_remote = False
            else:
                 is_remote = True

        work_mode = "Remote" if is_remote else "Onsite"

        # Assembler Row - RICH FORMAT matching Google Sheet + Person Data
        row = {
            "rank": 0, # Placeholder
            "probability_score": score_data['score'],
            "name": profile['name'],
            "title": profile['title'],
            "company": profile['company'],
            # New Rich Data Columns matching User's Sheet
            "domain": company_info.get("Domain", "N/A"),
            "funding_amount": company_info.get("Amount", "N/A"),
            "funding_round": company_info.get("Round", "Unknown"),
            "investors": company_info.get("Investors", "N/A"),
            # End New Columns
            "person_location": profile['person_location'],
            "company_hq": profile['company_hq'],
            "work_mode": work_mode,
            "pubmed_papers": enrichment['pubmed_paper_count'],
            "email": profile['email'],
            "phone": profile.get("phone", "N/A"),
            "linkedin_url": profile['linkedin_url']
        }
        results.append(row)
        progress_bar.progress((i + 1) / total)
        
    status_text.empty()
    progress_bar.empty()
    
    return pd.DataFrame(results)

def main():
    st.title("3D In-Vitro Lead Qualification Dashboard")
    st.markdown("This dashboard identifies, enriches, and ranks life-science professionals based on their probability of working with 3D in-vitro models.")
    
    # Sidebar: Match the screenshot
    st.sidebar.markdown("### 🔍 Filters")
    
    # Search
    search_query = st.sidebar.text_input("Search (Name, Title, Company, Location)")
    
    # Score Slider
    min_score = st.sidebar.slider("Minimum Probability Score", 0, 100, 0)
    
    # Data Loader (Expander to keep UI clean like screenshot)
    with st.sidebar.expander("⚙️ Data Source", expanded=True):
        default_urls = ""
        if os.path.exists(URLS_FILE):
            with open(URLS_FILE, 'r') as f:
                default_urls = f.read()
        
        input_urls = st.text_area("LinkedIn URLs", value=default_urls, height=100)
        fetch_btn = st.button("Fetch & Rank Leads", type="primary")

    # Processing Logic
    if fetch_btn:
        with st.spinner("Processing leads..."):
            urls = input_urls.strip().split('\n')
            df = process_leads(urls)
            
            # Ranking
            df = df.sort_values(by='probability_score', ascending=False)
            df.reset_index(drop=True, inplace=True)
            df.index += 1
            df['rank'] = df.index
            
            st.session_state['leads_data'] = df

    # Main Table Display
    if 'leads_data' in st.session_state:
        df = st.session_state['leads_data']
        
        # Application of Sidebar Filters
        filtered_df = df.copy()
        
        # 1. Min Score
        filtered_df = filtered_df[filtered_df['probability_score'] >= min_score]
        
        # 2. Search Query
        if search_query:
            filtered_df = filtered_df[filtered_df.apply(lambda row: search_query.lower() in str(row).lower(), axis=1)]
            
        # Display Table - CLEAN VIEW (Matching User's "Demo" screenshot)
        # We hide the extensive financial details here to keep it clean like the requested screenshot
        display_cols = [
            'rank', 'probability_score', 'name', 'title', 'company', 
            'person_location', 'company_hq', 'work_mode', 
            'funding_stage', 'pubmed_papers', 'email', 'phone', 'linkedin_url'
        ]
        
        # Ensure columns exist before selecting
        clean_display_cols = [c for c in display_cols if c in filtered_df.columns]
        
        st.dataframe(
            filtered_df[clean_display_cols],
            use_container_width=True,
            column_config={
                "linkedin_url": st.column_config.LinkColumn("linkedin"),
                "probability_score": st.column_config.ProgressColumn(
                    "probability_score",
                    format="%d",
                    min_value=0,
                    max_value=100
                ),
                "email": st.column_config.LinkColumn("email", display_text="Email")
            },
            hide_index=True
        )
        
        # Prepare Data for Export (Matching Google Sheet Headers - RICH DATA)
        export_df = filtered_df.copy()
        
        # Rename columns to match User's Spreadsheet + Lead Info
        export_cols = {
            "company": "Company",
            "domain": "Domain",
            "linkedin_url": "Linkedin",
            "funding_amount": "Amount (USD)",
            "funding_round": "Round",
            "investors": "Investors",
            "name": "Lead Name",
            "title": "Lead Title",
            "probability_score": "Propensity Score",
            "email": "Email",
            "phone": "Phone",
            "person_location": "Location",
            "work_mode": "Work Mode"
        }
        
        # Select and Rename
        available_cols = [c for c in export_cols.keys() if c in export_df.columns]
        export_df = export_df[available_cols].rename(columns=export_cols)
        
        st.success(f"Ready to download {len(export_df)} rows.")

        # CSV Download Button
        # Using utf-8-sig for better Excel compatibility
        csv = export_df.to_csv(index=False).encode('utf-8-sig')
        
        st.download_button(
            label="📥 Download Qualified Leads (Excel/CSV)",
            data=csv,
            file_name='lead_intelligence_output.csv',
            mime='application/csv',
            key='download_csv_v3'
        )

    # Footer
    st.markdown("<div class='footer-text'>Demo Web Agent • Identification → Enrichment → Ranking • Built with Streamlit</div>", unsafe_allow_html=True)

if __name__ == "__main__":
    main()
