#!/usr/bin/env python3

"""
    Generates the HTML CaTCH report.
"""

import sys
import argparse
import pandas as pd
import base64
import pdb


LOGO_PATH = "/data/templates/img/CaTCH_logo.png"


def main():
    parser = argparse.ArgumentParser(description = "Creates the basic CaTCH report")

    parser.add_argument("--name", type = str, required = True, help = "experiment name")
    parser.add_argument("--tsne", type = str, required = True, help = "path to the file containing the tSNE data")
    parser.add_argument("--metadata", type = str, required = True, help = "path to the file containing the metadata")
    parser.add_argument("--template", type = str, required = True, help = "path to the template file")
    parser.add_argument("--out", type = str, required = True, help = "path to the output report file")


    try:
        opts = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    tsne_data = pd.read_csv(opts.tsne, sep = " ", header = 0)
    meta_data = pd.read_csv(opts.metadata, sep = " ", header = 0)
    
    template_content = ""
    with (open(opts.template, "rt")) as hInput:
        for line in hInput:
            template_content += line

    # CaTCH logo
    binary_content = open(LOGO_PATH, 'rb').read() 
    base64_utf8_str = base64.b64encode(binary_content).decode('utf-8')
    dataurl = f'data:image/png;base64,{base64_utf8_str}'
    template_content = template_content.replace("%LOGO_DATA%", dataurl)

    # Experiment name
    template_content = template_content.replace("%EXP_NAME%", opts.name)

    # tSNE with samples
    vals_x = [f"{x:.2f}" for x in tsne_data['tSNE1'].to_list()]
    vals_y = [f"{y:.2f}" for y in tsne_data['tSNE2'].to_list()]
    vals_samples = [f"'{s}'" for s in tsne_data['Sample'].to_list()]
    vals_clusters = [str(x) for x in meta_data['Cluster'].to_list()]
    vals_stages = [f"'{s}'" for s in meta_data['CellStage'].to_list()]
    vals_features = [str(x) for x in meta_data['detected'].to_list()]
    vals_reads = [str(x) for x in meta_data['sum'].to_list()]
    vals_qccategories = [f"'{s}'" for s in meta_data['Category'].to_list()]
    vals_catchbc_counts = [x.split(";")[0] for x in meta_data['CaTCH.BC.counts'].fillna("0").to_list()]
    vals_catch_status = [f"'{s}'" for s in meta_data['CaTCH.Multiplet'].fillna("No data").replace([True, False], ["Multiplet", "Singlet"]).to_list()]
    template_content = template_content.replace("%TSNE_X%", f"[{','.join(vals_x)}]")\
                                       .replace("%TSNE_Y%", f"[{','.join(vals_y)}]")\
                                       .replace("%TSNE_SAMPLES%", f"[{','.join(vals_samples)}]")\
                                       .replace("%TOTAL_CELLS%", f"{len(vals_x):,}")\
                                       .replace("%TSNE_CLUSTERS%", f"[{','.join(vals_clusters)}]")\
                                       .replace("%TSNE_CELLSTAGES%", f"[{','.join(vals_stages)}]")\
                                       .replace("%BOXPLOTRAW_FEATURES%", f"[{','.join(vals_features)}]")\
                                       .replace("%BOXPLOTRAW_READS%", f"[{','.join(vals_reads)}]")\
                                       .replace("%BOXPLOTRAW_QCCATEGORIES%", f"[{','.join(vals_qccategories)}]")\
                                       .replace("%BOXPLOTRAW_CATCHBCCOUNTS%", f"[{','.join(vals_catchbc_counts)}]")\
                                       .replace("%BOXPLOTRAW_CATCHSTATUS%", f"[{','.join(vals_catch_status)}]")
                                        

    # Add a separate table row for each sample
    html = ""
    tmp = tsne_data["Sample"].value_counts().sort_index()
    for sample in tmp.index.to_list():
        n = tmp[sample]
        html += f'<tr><td><a href="#" src="{sample}.html" class="popup-link">{sample}</a></td><td>{n:,}</td></tr>'
    template_content = template_content.replace("%SAMPLES_TABLE%", html)

    with (open(opts.out, "wt")) as hOutput:
        print(template_content, file = hOutput)



if __name__ == "__main__":
    main()