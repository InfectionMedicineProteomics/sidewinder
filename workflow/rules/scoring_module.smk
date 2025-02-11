#!/usr/bin/env python3
"""Based on the TX-MS web-app implementation Cheetah-MS.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


rule modeling:
    # TODO:
    #   - Deal with no-hit-scenarios
    #   - Save euclidean distances to file
    #       - Include aa-numbering and xyz-coord anchors
    #   - Set container image from config.
    input:
        PDB = rules.pdb_handling.output.PDB,
        top_xls = rules.ms2.output.top_xls,
        dock_file = rules.megadock.output
    output:
        top_model = os.path.join(OUTDIR_BASE,
                                 txms_files,
                                 "{sample}",
                                 "{binder}", "{target}", "best_model.pdb" ),
        top_xls = os.path.join(OUTDIR_BASE,
                               txms_files,
                               "{sample}",
                               "{binder}", "{target}", "best_model_xls.txt"),
        scores = os.path.join(OUTDIR_BASE,
                              txms_files,
                              "{sample}",
                              "{binder}", "{target}", "lr_models/scores.csv")
    params:
        cut_off = config['distance_cut_off'],
        n_top_f = config['num_of_top_filters'],
        n_docks = config['model_samples'],
        modeling_script = f'{WD}/scripts/sidewinder-model.py',
        outdir = os.path.join(OUTDIR_BASE,
                              txms_files, "{sample}", "{binder}", "{target}")
    conda:
        f'{WD}/envs/biopython_env.yml'
    threads:
        round(workflow.cores * 0.10)
    shell:
        "python3 {params.modeling_script} "
        "--complex_pdb {input.PDB} "
        "--output_dir {params.outdir} "
        "--top_xls_file {input.top_xls} "
        "--dock_file {input.dock_file} "
        "--n_models {params.n_docks} "
        "--n_filters {params.n_top_f} "
        "--cut_off {params.cut_off}"
