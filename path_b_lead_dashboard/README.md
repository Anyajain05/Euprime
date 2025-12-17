# 3D In-Vitro Lead Intelligence Dashboard

## Problem Statement
Biotech Business Development (BD) teams struggle to identify the "right" scientists to sell to. Cold calling generic titles on LinkedIn has a low conversion rate. To identify a high-value lead, you need to know:
1.  **Role Fit**: Are they a toxicologist or safety scientist?
2.  **Scientific Intent**: Are they *actually* using 3D models/in-vitro methods (proven by publications)?
3.  **Business Intent**: Does their company have the funding (Series A/B/IPO) to pay?

## Solution
This dashboard acts as an intelligent layer on top of LinkedIn. It takes a list of profiles, enriches them with real-time scientific data (PubMed) and business data (Funding), and calculates a **Propensity-to-Buy Score (0-100)**.

## Data Sources
1.  **LinkedIn Profile Data**: Fetched via Agent (Simulated Proxycurl API for this demo).
2.  **Scientific Data**: Real-time query to **PubMed (NCBI)** via Biopython.
3.  **Funding Data**: Internal lookup CSV mapping companies to funding stages.

## Lead Scoring Logic
The **Propensity-to-Buying Score** is calculated as follows:

| Signal | Condition | Score |
| :--- | :--- | :--- |
| **Scientific Intent** | Published relevant paper (last 2 yr) | **+40** |
| **Role Fit** | Title: Toxicology / Safety / Hepatic | **+30** |
| **Company Intent** | Funding: Series A/B or IPO | **+20** |
| **Technographic** | Interest in In-vitro/NAMs (inferred) | **+15** |
| **Location** | Bio Hub (Boston, Cambridge, Basel, etc.) | **+10** |

*Score is normalized to max 100.*

## How to Run
1.  **Install Requirements**:
    ```bash
    pip install -r requirements.txt
    ```
2.  **Run the App**:
    ```bash
    streamlit run app.py
    ```
3.  **Use the Dashboard**:
    -   Click **"Fetch & Rank Leads"** to process the sample list.
    -   Watch the system enrich each profile LIVE.
    -   Filter by "Boston" or "High Probability".
    -   Download the CSV for your CRM.

*Note: This is a demo-ready prototype focusing on business logic.*
