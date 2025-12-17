# 3D In-Vitro Lead Intelligence Engine

## Business Problem
Business Development (BD) teams in the 3D in-vitro / organ-on-chip space face a "needle in the haystack" problem. They spend hours manually searching LinkedIn for toxicologists and safety pharmacologists, but fail to prioritize them efficiently.

A generic "Director of Toxicology" at a massive pharma company might be a low-quality lead if they focus purely on in-vivo methods. Conversely, a "Senior Scientist" at a specific biotech who just published a paper on "Liver-on-Chip" is a **high-value target**.

## System Design
This prototype automates the qualification process by fusing three distinct signals:
1.  **Role & Profile (LinkedIn)**: Does their title match decision-maker criteria?
2.  **Scientific Intent (PubMed)**: Have they published on relevant topics (e.g., "Hepatic Toxicity", "3D Cell Culture") in the last 24 months?
3.  **Business Intent (Funding)**: Is their company well-capitalized (Series A/B, IPO) and located in a key Bio Hub?

The **Probability Scoring Engine** synthesizes these signals into a single score (0-100), allowing BD teams to focus on the top 10% of leads that matter.

## Data Sources
-   **LinkedIn Data**: Simulated for this prototype via `leads_input.csv` (15-20 realistic profiles).
-   **Funding Data**: Simulated via `funding_data.csv`.
-   **PubMed Data**: **REAL TIME**. The agent queries the PubMed API (via Biopython) to find actual publication counts for the simulated names to demonstrate the logic.

## Usage
### Prerequisites
- Python 3.8+
- Internet connection (for PubMed API)

### Installation
1.  Navigate to the project directory:
    ```bash
    cd lead_intelligence_agent
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Running the App
Run the Streamlit dashboard:
```bash
streamlit run app.py
```

## Limitations
-   **LinkedIn Scraping**: Direct scraping is against TOS and strictly avoided here. Input is CSV-based.
-   **Identity Resolution**: Name matching on PubMed is simple (Name + Keywords) for demo purposes. Production systems would require stricter author disambiguation.

---
**Lead Intelligence Engine** - *Stop searching. Start closing.*
