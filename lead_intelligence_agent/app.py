import streamlit as st
import pandas as pd
from data_loader import load_and_enrich_data

st.set_page_config(
    page_title="3D In-Vitro Lead Intelligence Engine",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS for "Premium" feel
st.markdown("""
<style>
    .reportview-container {
        background: #0e1117;
    }
    .main-header {
        font-family: 'Inter', sans-serif;
        color: #f0f2f6;
    }
    .metric-card {
        background-color: #262730;
        padding: 20px;
        border-radius: 10px;
        border: 1px solid #41444b;
    }
    div[data-testid="stDataFrame"] {
        border: 1px solid #41444b;
        border-radius: 10px;
    }
</style>
""", unsafe_allow_html=True)

@st.cache_data
def get_data():
    return load_and_enrich_data()

def color_probability(val):
    color = 'red'
    if val > 70:
        color = '#00cc66' # Green
    elif val >= 40:
        color = '#ffcc00' # Yellow
    else:
        color = '#ff3333' # Red
    return f'color: {color}; font-weight: bold;'

def main():
    st.title("🧬 3D In-Vitro Lead Intelligence Engine")
    st.markdown("""
    **Prototype Description**: This web agent simulates identifying and ranking high-probability business development leads 
    for 3D in-vitro therapy models. It fuses signals from LinkedIn (simulated), Funding Data, and PubMed scientific publications.
    """)
    
    with st.spinner('Ingesting data and running AI enrichment pipeline...'):
        df = get_data()

    # Sidebar Filters
    st.sidebar.header("🔍 Filters")
    
    # Location Filter
    locations = ['All'] + sorted(df['Person_Location'].unique().tolist())
    selected_loc = st.sidebar.selectbox("Person Location", locations)
    
    # Company Filter
    companies = ['All'] + sorted(df['Company'].unique().tolist())
    selected_company = st.sidebar.selectbox("Company", companies)
    
    # Keyword Filter
    keyword_filter = st.sidebar.text_input("Filter by Keyword (Oncology, Toxicology, etc.)")
    
    # Filtering Logic
    filtered_df = df.copy()
    
    if selected_loc != 'All':
        filtered_df = filtered_df[filtered_df['Person_Location'] == selected_loc]
        
    if selected_company != 'All':
        filtered_df = filtered_df[filtered_df['Company'] == selected_company]
        
    if keyword_filter:
        filtered_df = filtered_df[
            filtered_df.apply(lambda row: keyword_filter.lower() in str(row).lower(), axis=1)
        ]

    # Metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Leads Identified", len(df))
    with col2:
        high_prob = len(df[df['Probability_Score'] > 70])
        st.metric("High Probability (>70)", high_prob)
    with col3:
        avg_score = int(df['Probability_Score'].mean())
        st.metric("Avg Portfolio Score", avg_score)
    with col4:
        papers = df['PubMed_Papers_Count'].sum()
        st.metric("Total Papers Found", papers)

    st.markdown("### 🎯 Ranked Lead Opportunities")

    # Apply Styling
    # We want to display specific columns
    display_cols = [
        'Rank', 'Probability_Score', 'Name', 'Title', 'Company', 
        'Person_Location', 'Company_HQ', 'Funding_Stage', 
        'PubMed_Papers_Count', 'Email', 'LinkedIn_URL'
    ]
    
    # Pandas Styler
    styler = filtered_df[display_cols].style.\
        applymap(color_probability, subset=['Probability_Score'])

    st.dataframe(
        styler,
        use_container_width=True,
        column_config={
            "LinkedIn_URL": st.column_config.LinkColumn("LinkedIn Profile"),
            "Probability_Score": st.column_config.ProgressColumn(
                "Score",
                help="Probability Score 0-100",
                format="%d",
                min_value=0,
                max_value=100,
            ),
        },
        hide_index=True
    )

    # CSV Export
    csv = filtered_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Export Filtered Leads to CSV",
        data=csv,
        file_name='lead_intelligence_export.csv',
        mime='text/csv',
    )
    
    st.markdown("---")
    st.caption("3D In-Vitro Lead Intelligence Engine | Prototype v1.0")

if __name__ == "__main__":
    main()
