import yaml

# Supplementary tables: compiled from pipeline outputs into supplementary_tables/
with open("config/supplementary_tables.yaml") as f:
  supplementary_tables = yaml.safe_load(f)["supplementary_tables"]

# Path template per source type ({comparison} is filled in from config/supplementary_tables.yaml)
supplementary_table_sources = {
  "rnaseq_degs": "output/rnaseq/analyses/{comparison}/degs/degs_{comparison}.csv",
  "rnaseq_fgsea": "output/rnaseq/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
  "olink_deps": "output/olink/analyses/{comparison}/deps/deps_{comparison}.csv",
  "melonardanaz2025_degs": "output/scrnaseq_melonardanaz2025/analyses/{comparison}/degs/degs_{comparison}.csv",
  "melonardanaz2025_fgsea": "output/scrnaseq_melonardanaz2025/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
  "toedter2011_degs": "output/microarray_toedter2011/analyses/{comparison}/degs/degs_{comparison}.csv",
  "toedter2011_fgsea": "output/microarray_toedter2011/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
}

for _table, _spec in supplementary_tables.items():
  if _spec["source"] not in supplementary_table_sources:
    raise ValueError(f"{_table}: unknown source '{_spec['source']}'")

def supplementary_table_input(wildcards):
  spec = supplementary_tables[wildcards.supplementary_table]
  return supplementary_table_sources[spec["source"]].format(comparison=spec["comparison"])

# Rules

rule supplementary_table:
  input:
    results_csv=supplementary_table_input,
  output:
    table_xlsx="supplementary_tables/{supplementary_table}.xlsx",
  params:
    source=lambda wc: supplementary_tables[wc.supplementary_table]["source"],
    comparison=lambda wc: supplementary_tables[wc.supplementary_table]["comparison"],
  log:
    "supplementary_tables/logs/{supplementary_table}.log",
  conda:
    "../envs/python_tables.yaml"
  message:
    "--- Supplementary table {wildcards.supplementary_table} ---"
  threads: 1
  script:
    "../scripts/supplementary_tables/compile_table.py"

rule supplementary_table_legends:
  input:
    "config/supplementary_tables.yaml",
  output:
    legends_txt="supplementary_tables/supplementary_table_legends.txt",
  message:
    "--- Supplementary table legends ---"
  run:
    with open(output.legends_txt, "w") as out:
      for table, spec in supplementary_tables.items():
        out.write(" ".join(spec["legend"].split()) + "\n\n")
