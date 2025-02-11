#!/usr/bin/env python3
"""Sidewinder Structure module.

Does this and that...

Input:
    - PDB file(s) or FASTA for the interacting proteins
Output:
    - PDB docking models
    - Theoretical cross-links
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


rule pdb_handling:
    input:
        # PDB_a = ancient(f'{config["targets_dir"]}/{{target}}.pdb'),
        PDB_a = lambda wildcards: str(config["targets_dir"])  + '/' + wildcards.target + '.pdb',
        # PDB_b = ancient(f'{config["binders_dir"]}/{{binder}}.pdb'),
        PDB_b = lambda wildcards: str(config["binders_dir"])  + '/' + wildcards.binder + '.pdb'
    output:
        PDB = os.path.join(OUTDIR_BASE, docking_files, "{binder}", "{target}", "complex_AB.pdb"),  ##   f'{{outdir}}/{{binder}}/{{target}}/complex_AB.pdb',
        seq_a = os.path.join(OUTDIR_BASE, docking_files, "{binder}", "{target}", "chain_A.txt"), # f'{{outdir}}/{{binder}}/{{target}}/chain_A.txt',
        seq_b = os.path.join(OUTDIR_BASE, docking_files, "{binder}", "{target}", "chain_B.txt"),  #f'{{outdir}}/{{binder}}/{{target}}/chain_B.txt',
        pdb_a = os.path.join(OUTDIR_BASE, docking_files, "{binder}", "{target}", "chain_A.pdb"),  #f'{{outdir}}/{{binder}}/{{target}}/chain_A.pdb',
        pdb_b = os.path.join(OUTDIR_BASE, docking_files, "{binder}", "{target}", "chain_B.pdb"), ##f'{{outdir}}/{{binder}}/{{target}}/chain_B.pdb'
    params:
        pdb_handling = f'{WD}/scripts/utils/pdb_handling.py',
        outdir = os.path.join(OUTDIR_BASE,
                              docking_files, "{binder}", "{target}")
    conda:
        f'{WD}/envs/biopython_env.yml'
    shell:
        "python3 {params.pdb_handling} "
        "--pdb1 {input.PDB_a} "
        "--pdb2 {input.PDB_b} "
        "--output_dir {params.outdir}"

if BLOCK_FV:
    # TODO:
    #   - Save process to SQLite DB and make output temporary.
    rule block_pdb:
        input:
            mc_pdb = rules.pdb_handling.input.PDB_a,
            sc_pdb = rules.pdb_handling.output.pdb_a
        output:
            pdb_a = os.path.join(OUTDIR_BASE,
                                 docking_files,
                                 "{binder}", "{target}", "chain_A_blocked.pdb")
        params:
            fv_blocking = f'{WD}/scripts/mdock-block_pdb.py',
            #outdir = f'{{outdir}}/{{binder}}/{{target}}'
            outdir = os.path.join(OUTDIR_BASE,
                                  docking_files, "{binder}", "{target}")
        conda:
            f'{WD}/envs/biopython_env.yml'
        shell:
            "python3 {params.fv_blocking} "
            "--multi_chain_pdb {input.mc_pdb} "
            "--single_chain_pdb {input.sc_pdb} "
            "--output_dir {params.outdir}"

rule seq2xl:
    input:
        PDB_a = rules.pdb_handling.input.PDB_a,
        PDB_b = rules.pdb_handling.input.PDB_b
    output:
        os.path.join(OUTDIR_BASE,
                     docking_files, "{binder}", "{target}", "all_xls.txt")
    params:
        seq2xl = f'{WD}/scripts/utils/seq2xl_v1.5.py',
        #outdir = f'{{outdir}}/{{binder}}/{{target}}'
        outdir = os.path.join(OUTDIR_BASE,
                              docking_files,"{binder}", "{target}")
    conda:
        f'{WD}/envs/biopython_env.yml'
    shell:
        "python3 {params.seq2xl} "
        "--pdb_file1 {input.PDB_a} "
        "--pdb_file2 {input.PDB_b} "
        "--output_dir {params.outdir}"

rule megadock:
    input:
        pdb_a = (rules.block_pdb.output.pdb_a
                 if BLOCK_FV else rules.pdb_handling.output.pdb_a),
        pdb_b = rules.pdb_handling.output.pdb_b
    output:
        os.path.join(OUTDIR_BASE,
                     docking_files, "{binder}_+_{target}_megadock.out")
    params:
        # megadock = '/opt/MEGADOCK/megadock-gpu',  # GPU
        megadock = '/usr/local/megadock/megadock-4.1.1/megadock',  # CPU
        predictions = config['docking_samples'],
        pdb_a = os.path.join(OUTDIR_BASE,
                             docking_files,
                             "{binder}", "{target}", "chain_A_blocked.pdb")
                if BLOCK_FV else
                os.path.join(OUTDIR_BASE,
                             docking_files,
                             "{binder}", "{target}", "chain_A.pdb"),
        pdb_b = os.path.join(OUTDIR_BASE,
                             docking_files,
                             "{binder}", "{target}", "chain_B.pdb"),
    singularity:
        # 'workflow/envs/megadock_4.1.4-gpu.sif'  # GPU
        f'{WD}/envs/megadock.sif'  # CPU
    threads:
        max(1, int(workflow.cores/1))
    resources:
        nvidia_gpu=4
    log:
        logfile = os.path.join(OUTDIR_BASE,
                               docking_files,
                               "{binder}", "{target}", 'megadock.log')
    shell:
        "{params.megadock} "
        "-R {params.pdb_a} "
        "-L {params.pdb_b} "
        "-o {output} "
        "-N {params.predictions} "
        "&> {log.logfile}"
