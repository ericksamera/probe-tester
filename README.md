# probe-tester

`probe-tester` is a modular Python toolkit for evaluating primers and probes against large collections of genomes.
It supports automated download of NCBI genomes, in silico PCR (`IPCRESS`) matching, summary/statistics reporting, and modern, user-friendly CLI progress indicators.

---

## **Features**

- **Automated download** of GenBank genomes by species, genus, or custom set (NCBI Datasets CLI required)
- **Flexible in silico PCR (IPCRESS)**: test any forward/reverse/probe set on downloaded genomes
- **Multiprocessing support**: scale up to thousands of genomes
- **Rich progress bars** (if available) or clean fallback output
- **Summary and specificity panelization**: separates “target” and “non-target” organisms for true molecular diagnostics benchmarking
- **CSV/Markdown export, table summaries, and metadata in results**
- **Graceful error handling, dependency checks, and quality-of-life CLI features**

---

## **Quick Start**

### **Requirements**

- `python >= 3.8`
- [`ncbi-datasets-cli`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- `Go >= 1.22` to build the pinned `ipcr` submodule (not needed when using an external `ipcr` via `--ipcr-bin`)
- Python packages: (optionally) `rich`,

**Easiest, install with conda**

```shell
mamba env create -f environment.yaml
```

### **ipcr submodule**

`probe-tester` uses [`KPU-AGC/ipcr`](https://github.com/KPU-AGC/ipcr) as a
Git submodule at `third_party/ipcr`. Clone with submodules enabled:

```shell
git clone --recurse-submodules <probe-tester-url>
```

For an existing checkout:

```shell
git submodule update --init --recursive
```

On the first `assay` run, `probe-tester` builds the pinned submodule revision
to `bin/ipcr`. If the submodule revision changes, the binary is rebuilt
automatically. `--ipcr-bin /path/to/ipcr` remains available as an explicit
override, and a system `ipcr` on `PATH` is used when the submodule is absent.

---

### **1. Download Genomes**

Download a specific taxonomic level of genomes:

Ex. Download 5 genomes for the species "Mycoplasmopsis bovis:

```sh
python main.py download --taxon "Mycoplasmopsis bovis" --mode species --max-genomes 5
```

Or knowing the target species, download genomes from the upper taxonomic level.

Ex. Download 5 genomes for all members of the genus that "Mycoplasmopsis bovis" belongs to.

```sh
python main.py download --taxon "Mycoplasmopsis bovis" --mode parent --max-genomes 5
```

Add `--dry-run` to preview downloads, or `--force` to skip confirmation for large jobs.

---

### **Marker Reference Panels (BOLD)**

Download public marker sequences for all identified species within a taxonomic query.
`COI` and `CO1` are normalized to BOLD's `COI-5P` marker code.

```sh
python main.py download-marker \
  --taxon "Ixodidae" \
  --marker COI \
  --min-length 500
```

This writes one multi-record FASTA per species below `markers/COI-5P/`. Test the
panel with the normal assay command by selecting `--markers-dir`:

```sh
python main.py assay \
  --markers-dir markers/COI-5P \
  --forward ... \
  --reverse ... \
  --probe ... \
  --threads 8
```

Marker mode currently requires the `ipcr` engine because its JSONL output
preserves FASTA record IDs. Each barcode record is counted independently in
the summary, including records with no amplification.

### **2. Run the Assay**

Test your primers and probe against all genomes:

```sh
python main.py assay \
  --forward TCTAATTTTTTCATCATCGCTAATGC \
  --reverse TCAGGCCTTTGCTACAATGAAC \
  --probe AACTGCATCATATCACATACT \
  --run-name Parker-2017 \
  --threads 4
```

- Uses multiprocessing (`--threads`) if desired for speed.
- `IPCRESS` is called under the hood.
- Results are written to a JSON file (named by run or timestamp).

---

### **3. Summarize Results**

Get a sensitivity/specificity table, with targets vs. non-targets, in pretty, CSV, or markdown format:

```sh
python main.py summarize --input results-Parker-2017.json \
  --target "Mycoplasmopsis-bovis" \
  --format text
```

Or for all _Mycoplasmopsis-_ species:

```sh
python main.py summarize --input results-Parker-2017.json \
  --target "Mycoplasmopsis-*"
```

Export to CSV/markdown for figures or publication:

```sh
python main.py summarize --input results.json --format csv --target "Mycoplasmopsis-bovis"
```

---

## **Advanced Usage**

- **List available species/genomes:**

  ```sh
  python main.py list --taxon "Mycoplasmopsis" --mode parent
  ```

---

## **Best Practices**

- Use `--dry-run` to preview large downloads before running them.
- Use `--run-name` for descriptive result files.
- Use `--threads N` on modern CPUs for much faster assay runs.
- Regularly update the NCBI Datasets CLI for the latest assembly info.

---

## **Example Full Workflow**

```sh
# Download up to 10 genomes for all organisms belonging to Mycoplasmopsis bovis's genus.
python main.py download --taxon "Mycoplasmopsis bovis" --max-genomes 10

# Run probe set against all genomes, using 8 CPU threads
python main.py assay --forward ... --reverse ... --probe ... --threads 8 --run-name my-panel

# Summarize, showing only Mycoplasmopsis-bovis as the target
python main.py summarize --input results-my-panel.json --target "Mycoplasmopsis-bovis" --format csv
```

---

**Happy analyzing!**
