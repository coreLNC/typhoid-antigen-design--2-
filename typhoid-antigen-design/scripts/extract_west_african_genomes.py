"""List West African genomes from NCBI Datasets genome-summary JSON-lines files.

usage: python3 -I scripts/extract_west_african_genomes.py OUT.tsv SUMMARY1.jsonl [SUMMARY2.jsonl ...]

The .jsonl inputs come from `datasets summary genome taxon <taxon> --as-json-lines`
(one file each for S. Typhi, S. Paratyphi A and complete Salmonella genomes). They are
large (15-60 MB) and are not stored in this repository; OUT.tsv is the small derived table.
A genome is kept when the country part of the BioSample geo_loc_name (text before ':')
is a West African country. Matching is on the whole country token, so Niger != Nigeria.
"""
import csv, json, sys
WEST_AFRICA = {"Benin", "Burkina Faso", "Cabo Verde", "Cape Verde", "Cote d'Ivoire", "Ivory Coast",
               "Gambia", "The Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Mali",
               "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo"}
FIELDS = ["accession", "organism", "serovar", "strain", "country", "geo_loc_name", "collection_date",
          "assembly_level", "release_date", "biosample", "source_file"]

def row(d, src):
    bs = d.get("assembly_info", {}).get("biosample", {})
    geo = bs.get("geo_loc_name", "") or ""
    return {"accession": d.get("accession", ""), "organism": d.get("organism", {}).get("organism_name", ""),
            "serovar": bs.get("serovar", "") or "",
            "strain": (d.get("organism", {}).get("infraspecific_names", {}) or {}).get("strain", ""),
            "country": geo.split(":")[0].strip(), "geo_loc_name": geo,
            "collection_date": bs.get("collection_date", ""),
            "assembly_level": d.get("assembly_info", {}).get("assembly_level", ""),
            "release_date": d.get("assembly_info", {}).get("release_date", ""),
            "biosample": bs.get("accession", ""), "source_file": src}

def main(out, *files):
    seen, rows = set(), []
    for f in files:
        for line in open(f):
            r = row(json.loads(line), f.split("/")[-1])
            if r["country"] in WEST_AFRICA and (r["accession"], r["source_file"]) not in seen:
                seen.add((r["accession"], r["source_file"])); rows.append(r)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(len(rows), "West African records written to", out)

if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    main(*sys.argv[1:])
