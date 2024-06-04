import pandas as pd
import subprocess
import time
import os
import zipfile

def event_table(df):
    event_rows =[]
    for i,r in df.iterrows():
        if pd.isna(r["Amino acid edits"]):
            continue
        amino_acid_edits = r["Amino acid edits"].split(";");
        mutation_categories = r["Mutation category"].split(";");
        nucleotide_edits = r["Nucleotide edits"].split(";"); 
        for j in range(len(amino_acid_edits)):
            edit = amino_acid_edits[j]
            if edit == "" or "utr" in edit or "UTR" in edit or "Exon" in edit:
                continue
            event_row = r.copy()
            event_row["edit"] = edit
            event_row["edit_location"] = int(edit[3:(len(edit)-3)])
            event_row["original_aa"] = edit[0] + edit[1] + edit[2]
            event_row["new_aa"] = edit[-3] + edit[-2] + edit[-1]
            event_row["mutation_category"] = mutation_categories[j]
            event_row["nucleotide_edits"] = nucleotide_edits[j]
            event_rows.append(event_row)
    return pd.DataFrame(event_rows)

def write_input_file(transcript_id, gene_symbol):
    with open("input.txt", "w") as myfile:
        myfile.write("Transcript ID" + "\t" +  "Gene Symbol" + "\n")
        myfile.write(transcript_id + "\t" +  gene_symbol)

def submit_jobs(pam_types, 
                edit_types,
                intron_buffer, 
                filter_gc_content, 
                edit_window,
                transcript_id,
                gene_symbol,
                sg_len=20
                ):
    output_names = {}
    ps = []
    for pam in pam_types:
        for edit in edit_types:
            write_input_file(transcript_id, gene_symbol),
            output_name = "_".join([pam, gene_symbol,transcript_id, edit[0], edit[-1]])
            output_names[output_name] = pam
            command = "python base_editing_guide_designs.py --input-file input.txt --input-type tid -" \
            "-variant-file variant_summary.txt --pam {} --edit-window {} --sg-len {} --edit {} " \
            "--intron-buffer {} --filter-gc {} --output-name " \
            "{}".format(pam, edit_window, sg_len, edit, intron_buffer, filter_gc_content, output_name)
            p = subprocess.check_call(command)
            time.sleep(3)
    exit_codes = [p.wait() for p in ps]
    return output_names


def process_output(output_names):
    output_dfs = []
    results_dir = ""
    dirs = os.listdir("./")
    with zipfile.ZipFile('grna_results.zip', 'w') as myzip:
        for output in output_names.keys():
            for dir in dirs:
                if output in dir:
                    results_dir = dir
                    break
            file_name = results_dir + "/" + "sgrna_designs_" + output + ".txt"
            myzip.write(f)
            df = pd.read_csv(file_name, delimiter="\t")
            df["pam_pattern"] = output_names[output]
            output_dfs.append(df)
    return pd.concat(output_dfs)