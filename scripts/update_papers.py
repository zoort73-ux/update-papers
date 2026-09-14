#!/usr/bin/env python3
"""
Update data/yonsei.json, data/ilsan.json and data/meta.json by searching
PubMed (NCBI E-utilities) directly.

Queries reproduced from the original manual research:

  Yonsei:
    "Institute of Kidney Disease Research"[Affiliation] AND Yonsei[Affiliation]

  Ilsan (broad query, then filtered — see below):
    "National Health Insurance Service Ilsan Hospital"[Affiliation]
    AND "Internal Medicine"[Affiliation]

The Ilsan broad query can match papers where the two phrases appear in
DIFFERENT authors' affiliations (false positive: an Ilsan Hospital author
from an unrelated department, e.g. Radiology, co-authoring with an
unrelated author whose affiliation happens to contain "Internal
Medicine"). So after fetching full records we keep only articles where
BOTH phrases co-occur in a SINGLE author's affiliation string.

Only articles with year >= 2020 are kept (a small number of PubMed date
hits have an epub/print year mismatch vs. the pdat filter).

Run with no arguments. Designed to run unattended from GitHub Actions,
but safe to run locally too:

    python3 scripts/update_papers.py
"""

import json
import os
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()  # optional, raises rate limit if set

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "data")

MONTHS = {
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06",
    "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
}

QUERIES = {
    "yonsei": '"Institute of Kidney Disease Research"[Affiliation] AND Yonsei[Affiliation]',
    "ilsan": '"National Health Insurance Service Ilsan Hospital"[Affiliation] AND "Internal Medicine"[Affiliation]',
}


def http_get(url, params, retries=5):
    if API_KEY:
        params = {**params, "api_key": API_KEY}
    qs = urllib.parse.urlencode(params)
    full_url = f"{url}?{qs}"
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(full_url, headers={"User-Agent": "paper-archive-updater/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read()
        except Exception as e:
            last_err = e
            wait = 2 * (attempt + 1)
            print(f"  request failed ({e}); retrying in {wait}s...", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"Failed after {retries} retries: {last_err}")


def esearch_all_pmids(query):
    """Return the full list of PMIDs matching `query`, paginating via retstart."""
    pmids = []
    retstart = 0
    retmax = 200
    while True:
        raw = http_get(f"{EUTILS}/esearch.fcgi", {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retstart": retstart,
            "retmax": retmax,
            "sort": "pub+date",
        })
        data = json.loads(raw)
        result = data.get("esearchresult", {})
        ids = result.get("idlist", [])
        pmids.extend(ids)
        count = int(result.get("count", "0"))
        retstart += retmax
        # be polite to NCBI
        time.sleep(0.34 if not API_KEY else 0.11)
        if retstart >= count or not ids:
            break
    return pmids


def _text(el):
    return "".join(el.itertext()).strip() if el is not None else ""


def parse_article(article_el):
    """Parse one <PubmedArticle> element into our record schema (+ raw author affiliations)."""
    medline = article_el.find("MedlineCitation")
    art = medline.find("Article")
    pmid = _text(medline.find("PMID"))

    title = _text(art.find("ArticleTitle")).rstrip()

    journal_el = art.find("Journal")
    journal = ""
    if journal_el is not None:
        iso = journal_el.find("ISOAbbreviation")
        title_el = journal_el.find("Title")
        journal = _text(iso) or _text(title_el)

    # Date: prefer ArticleDate (epub), fall back to Journal/JournalIssue/PubDate
    year = month = None
    art_date = art.find("ArticleDate")
    if art_date is not None:
        y = _text(art_date.find("Year"))
        m = _text(art_date.find("Month"))
        if y:
            year, month = y, (m if m else "01")
    if not year:
        pub_date = art.find("Journal/JournalIssue/PubDate")
        if pub_date is not None:
            y = _text(pub_date.find("Year"))
            m = _text(pub_date.find("Month"))
            if not y:
                # MedlineDate fallback e.g. "2020 Jan-Feb"
                md = _text(pub_date.find("MedlineDate"))
                if md and len(md) >= 4 and md[:4].isdigit():
                    y = md[:4]
                    m = md[5:8] if len(md) >= 8 else ""
            if y:
                year = y
                month = m if m else "01"
    if not year:
        year, month = "0000", "01"
    if month and not month.isdigit():
        month = MONTHS.get(month[:3], "01")
    month = (month or "01").zfill(2)

    # DOI
    doi = ""
    for eloc in art.findall("ELocationID"):
        if eloc.get("EIdType") == "doi":
            doi = _text(eloc)
            break
    if not doi:
        for aid in article_el.findall("PubmedData/ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = _text(aid)
                break

    # Authors + per-author affiliation strings (for the Ilsan co-occurrence filter)
    authors = []
    author_affils = []  # list of (author_display_name, [affiliation strings])
    author_list = art.find("AuthorList")
    if author_list is not None:
        for a in author_list.findall("Author"):
            last = _text(a.find("LastName"))
            fore = _text(a.find("ForeName")) or _text(a.find("Initials"))
            if not last and not fore:
                collective = _text(a.find("CollectiveName"))
                if collective:
                    authors.append(collective)
                continue
            initials = "".join(p[0] for p in fore.split()) if fore else ""
            display = f"{last} {initials}".strip() if last else fore
            authors.append(display)
            affils = [ _text(aff.find("Affiliation")) for aff in a.findall("AffiliationInfo") ]
            affils = [x for x in affils if x]
            author_affils.append((display, affils))

    return {
        "pmid": pmid,
        "title": title,
        "journal": journal,
        "date": f"{year}-{month}",
        "year": year,
        "month": month,
        "doi": doi,
        "authors": ", ".join(authors),
        "_author_affils": author_affils,  # stripped before writing to disk
    }


def efetch_articles(pmids):
    """Fetch full metadata for a list of pmids, batching to stay well under NCBI limits."""
    records = []
    batch_size = 150
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        print(f"  fetching metadata {i+1}-{i+len(batch)} of {len(pmids)}...")
        raw = http_get(f"{EUTILS}/efetch.fcgi", {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
        })
        root = ET.fromstring(raw)
        for article_el in root.findall("PubmedArticle"):
            try:
                records.append(parse_article(article_el))
            except Exception as e:
                print(f"  WARNING: failed to parse one article: {e}", file=sys.stderr)
        time.sleep(0.34 if not API_KEY else 0.11)
    return records


def ilsan_same_author_filter(records):
    """Keep only records where a SINGLE author's affiliation contains both
    'ilsan hospital' and 'internal medicine' (case-insensitive)."""
    kept = []
    for r in records:
        ok = False
        for _name, affils in r.get("_author_affils", []):
            for aff in affils:
                al = aff.lower()
                if "ilsan hospital" in al and "internal medicine" in al:
                    ok = True
                    break
            if ok:
                break
        if ok:
            kept.append(r)
    return kept


def clean_record(r):
    r = dict(r)
    r.pop("_author_affils", None)
    return r


def load_existing(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return []


def update_dataset(key, query, existing_path, apply_ilsan_filter=False):
    print(f"[{key}] searching PubMed...")
    pmids = esearch_all_pmids(query)
    print(f"[{key}] {len(pmids)} candidate PMIDs")

    existing = load_existing(existing_path)
    existing_pmids = {r["pmid"] for r in existing}
    new_pmids = [p for p in pmids if p not in existing_pmids]
    print(f"[{key}] {len(new_pmids)} new PMIDs to fetch")

    if not new_pmids:
        return existing, 0

    fetched = efetch_articles(new_pmids)
    fetched = [r for r in fetched if r["year"] != "0000" and int(r["year"]) >= 2020]

    if apply_ilsan_filter:
        before = len(fetched)
        fetched = ilsan_same_author_filter(fetched)
        print(f"[{key}] same-author affiliation filter: {before} -> {len(fetched)}")

    added = len(fetched)
    merged = [clean_record(r) for r in fetched] + existing
    # dedupe defensively, keep first occurrence, sort desc by date
    seen = set()
    deduped = []
    for r in merged:
        if r["pmid"] in seen:
            continue
        seen.add(r["pmid"])
        deduped.append(r)
    deduped.sort(key=lambda r: (r["year"], r["month"], r["pmid"]), reverse=True)
    return deduped, added


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    yonsei_path = os.path.join(DATA_DIR, "yonsei.json")
    ilsan_path = os.path.join(DATA_DIR, "ilsan.json")
    meta_path = os.path.join(DATA_DIR, "meta.json")

    yonsei_records, yonsei_added = update_dataset("yonsei", QUERIES["yonsei"], yonsei_path, apply_ilsan_filter=False)
    ilsan_records, ilsan_added = update_dataset("ilsan", QUERIES["ilsan"], ilsan_path, apply_ilsan_filter=True)

    with open(yonsei_path, "w") as f:
        json.dump(yonsei_records, f, ensure_ascii=False, indent=2)
    with open(ilsan_path, "w") as f:
        json.dump(ilsan_records, f, ensure_ascii=False, indent=2)

    kst = timezone(timedelta(hours=9))
    now_kst = datetime.now(kst).strftime("%Y-%m-%d %H:%M KST")
    with open(meta_path, "w") as f:
        json.dump({"last_updated": now_kst}, f, ensure_ascii=False, indent=2)

    total_added = yonsei_added + ilsan_added
    print(f"\nDone. Yonsei +{yonsei_added} (total {len(yonsei_records)}), "
          f"Ilsan +{ilsan_added} (total {len(ilsan_records)}).")

    # Emit a simple flag file for the workflow to decide whether to commit.
    with open(os.path.join(DATA_DIR, ".last_run_added"), "w") as f:
        f.write(str(total_added))


if __name__ == "__main__":
    main()
