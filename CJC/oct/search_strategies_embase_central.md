# Updated search: Embase and CENTRAL translations of the PubMed strategy

These are translations of PubMed search #7 (search.docx). Like #7, they combine the population/procedure block with the drug-pair block, use no RCT filter, and are limited to records added from 1 September 2022 to the search date.

## Embase (Ovid, Embase 1974 to present, via McGill Library)

```
1  exp acute coronary syndrome/ or exp heart infarction/ or exp heart muscle ischemia/ or unstable angina pectoris/
2  ("acute coronary syndrome*" or ((myocardial or heart) adj infarction*) or "acute mi" or ((myocardial or heart muscle) adj2 isch?emi*) or unstable angina* or STEMI or NSTEMI).ti,ab,kw.
3  exp percutaneous coronary intervention/ or (((percutaneous coronary or heart muscle) adj (intervention* or revasculari#ation*)) or PCI).ti,ab,kw.
4  exp stent/ or stent*.ti,ab,kw.
5  exp coronary artery bypass graft/ or "coronary artery bypass*".ti,ab,kw.
6  or/1-5
7  ((clopidogrel and ticagrelor) or (clopidogrel and prasugrel) or (prasugrel and ticagrelor)).ti,ab,kw.
8  6 and 7
9  limit 8 to dc=20220901-20261020
```

Notes:
- If Ovid rejects `dc=` on your platform, use `limit 8 to em=202236-202643` (Embase entry week) or `limit 8 to yr="2022-Current"` and record which limit you used.
- Optional, if line 9 returns several hundred records, you can add an RCT filter to reduce screening:
  `10  9 and (randomized controlled trial/ or controlled clinical trial/ or random*.ti,ab,kw. or trial.ti,ab,kw.)`
  The PubMed search had no such filter, so for strict parallelism keep line 9 and report it in Supplementary S1.
- Export: select all → Export → format **RIS**, fields **Complete Reference** (including abstract and DOI). Save the file as `embase_update.ris`.

## Cochrane CENTRAL (Cochrane Library, Advanced search → Search manager)

```
#1  MeSH descriptor: [Acute Coronary Syndrome] explode all trees
#2  MeSH descriptor: [Myocardial Infarction] explode all trees
#3  MeSH descriptor: [Myocardial Ischemia] explode all trees
#4  MeSH descriptor: [Angina, Unstable] explode all trees
#5  ("acute coronary syndrome*" OR (myocardial NEXT infarction*) OR (heart NEXT infarction*) OR "acute MI" OR ((myocardial OR "heart muscle") NEAR/2 isch*emi*) OR "unstable angina*" OR STEMI OR NSTEMI):ti,ab,kw
#6  MeSH descriptor: [Percutaneous Coronary Intervention] explode all trees
#7  (("percutaneous coronary" NEXT (intervention* OR revasculari*ation*)) OR PCI):ti,ab,kw
#8  MeSH descriptor: [Stents] explode all trees
#9  stent*:ti,ab,kw
#10 MeSH descriptor: [Coronary Artery Bypass] explode all trees
#11 "coronary artery bypass*":ti,ab,kw
#12 {OR #1-#11}
#13 ((clopidogrel AND ticagrelor) OR (clopidogrel AND prasugrel) OR (prasugrel AND ticagrelor)):ti,ab,kw
#14 #12 AND #13
```

Notes:
- After running #14, click the **Trials** tab. Then set **Date added to CENTRAL: 01/09/2022 – 20/10/2026** using "Limits / Custom range" in the Search manager.
- CENTRAL includes trial-registry records (ClinicalTrials.gov and WHO ICTRP). Keep them in the export; the merge script reports them separately.
- Export: select all → **Export selected** → format **RIS (with abstract)**. Save the file as `central_update.ris`.

## Then

Put both `.ris` files in `CJC/oct`, together with the PubMed screening log that is already there, and run `search_merge_dedup.R`. It does three things:

1. It de-duplicates across PubMed, Embase and CENTRAL by DOI, then PMID, then normalized title.
2. Records already screened in PubMed inherit that decision.
3. Any genuinely new records are written to `search_update_new_to_screen.csv` for screening. Once their decisions are entered, `prisma_flow_update.R` picks up all the counts automatically.
