***Design guides***

Web app for designing base editor guide RNAs.

***Setup*** 

```sh
conda env create -f environment.yml
conda activate grna_dash
git clone https://github.com/mhegde/base-editor-design-tool.git
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.g
python app.py
```