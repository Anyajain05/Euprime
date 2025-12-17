🧬 3D In-Vitro Lead Intelligence Engine

A Business Development–Focused Web Agent Prototype

📌 Overview

This project is a working prototype of a web-based lead intelligence system designed to help biotech business development teams identify, enrich, and prioritize scientists and decision-makers who are most likely to work with 3D in-vitro models for therapy design.

Instead of generating generic lead lists, the system applies scientific relevance, business readiness, and role-based intent to rank leads by their propensity to buy.

The prototype is built as a Python + Streamlit dashboard and is intentionally demo-scale, focusing on clarity of logic and decision-making, not production scraping.

🎯 Business Problem

Biotech BD teams face three major challenges:

Too many potential contacts (LinkedIn, conferences, papers)

No clear way to know who actually needs the technology

Significant time wasted on low-probability outreach

This tool answers a simple BD question:

“Who should we speak to first, and why?”

🧠 Solution Approach

The system follows a 3-stage pipeline, directly aligned with the assignment description.

🔹 Stage 1: Identification (Input)

The tool ingests a list of target profiles that would typically come from:

LinkedIn Search / Sales Navigator

Conference attendee lists

For this prototype, profiles are provided as a structured input dataset containing:

Name

Job Title

Company

Person Location

Company HQ

LinkedIn URL

This represents the initial universe of potential leads.

🔹 Stage 2: Enrichment (Data Gathering)

Each identified profile is enriched with additional signals to assess intent and feasibility.

1. Scientific Intent (PubMed)

Queries PubMed for recent publications (last 24 months)

Keywords include:

Drug-Induced Liver Injury

Hepatic Toxicity

3D Cell Culture

Organ-on-Chip

Indicates whether the person is actively working on relevant problems

2. Business & Funding Intent

Company funding stage is mapped (Seed, Series A, Series B, etc.)

Indicates budget availability for adopting premium technologies

3. Location Intelligence

Distinguishes between:

Person’s location (remote / home)

Company HQ

Flags major biotech hubs:

Boston / Cambridge

Bay Area

Basel

UK Golden Triangle

🔹 Stage 3: Ranking (Probability Engine)

A transparent rule-based scoring engine assigns each lead a Propensity-to-Buy score (0–100).

Scoring Logic
Signal Category	Criteria	Weight
Role Fit	Title contains Toxicology, Safety, Hepatic, 3D	+30
Scientific Intent	Published relevant paper in last 2 years	+40
Company Intent	Series A / Series B / IPO funded	+20
Technographic	Uses / open to in-vitro or NAMs	+15
Location	HQ in major biotech hub	+10

This ensures:

Senior, funded, scientifically active profiles rank highest

Junior or unfunded profiles naturally rank lower

📊 Output: Lead Generation Dashboard

The final output is a dynamic Streamlit dashboard showing a ranked table of leads.

Dashboard Features

Probability-based ranking

Search and filter (e.g., “Boston”, “Toxicology”)

Clear separation of:

Person location

Company HQ

Exportable to CSV / Excel

Table Columns
Rank
Probability Score
Name
Title
Company
Person Location
Company HQ
Funding Stage
PubMed Paper Count
Email (simulated)
LinkedIn URL


This allows BD teams to:

Decide who to call first

Plan in-person meetings vs remote outreach

Focus efforts on highest-ROI leads

🛠️ Tech Stack

Python 3

Streamlit

Pandas

Requests

Biopython (PubMed / Entrez)

⚠️ Important Notes & Limitations

Direct LinkedIn scraping is not performed

LinkedIn APIs are restricted

The prototype simulates LinkedIn data or uses provider-style logic (e.g., Proxycurl)

This is a demo-scale implementation

Designed to show logic, not scale

The architecture is intentionally modular

Can be extended to real APIs and larger datasets

▶️ How to Run
pip install -r requirements.txt
streamlit run app.py


Then open the local URL shown in the terminal.

🚀 Why This Matters

This prototype demonstrates:

Business development thinking

Scientific signal interpretation

Clear prioritization logic

Practical dashboard delivery

It shows how a BD team can move from “who are all the scientists?” to “who should we speak to today?”

👤 Author

Anya Jain
