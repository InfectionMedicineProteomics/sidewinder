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


SCHEME = config['cdr_annotation_scheme']


rule pdb_handling:
    input:
        PDB_a = lambda wc: str(config["antibody_dir"]) + f'/{wc.antibody}.pdb',
        PDB_b = lambda wc: str(config["antigen_dir"]) + f'/{wc.antigen}.pdb'
    output:
        pdb_a = os.path.join(OUTDIR_BASE,
                             DOCKING_FILES,
                             "{antigen}_+_{antibody}", "{antibody}.pdb"),
        pdb_b = os.path.join(OUTDIR_BASE,
                             DOCKING_FILES,
                             "{antigen}_+_{antibody}", "{antigen}.pdb")
    params:
        pdb_handling = f'{WD}/scripts/utils/pdb_handling.py',
        outdir = os.path.join(OUTDIR_BASE,
                              DOCKING_FILES, "{antigen}_+_{antibody}"),
        seq_outdir = os.path.join(OUTDIR_BASE,
                                  SUPPORT_FILES)
    conda:
        f'{WD}/envs/biopython_env.yml'
    shell:
        "python3 {params.pdb_handling} "
        "--pdb1 {input.PDB_a} "
        "--pdb2 {input.PDB_b} "
        "--output_dir {params.outdir} "
        "--seq_output_dir {params.seq_outdir}"

if BLOCK_FV:
    # TODO:
    #   - Save process to SQLite DB and make output temporary.
    rule block_pdb:
        """Block Fab variable region (FV) in antibody PDB file."""
        input:
            mc_pdb = rules.pdb_handling.input.PDB_a,
            sc_pdb = rules.pdb_handling.output.pdb_a
        output:
            pdb_a = os.path.join(OUTDIR_BASE,
                                 DOCKING_FILES,
                                 "{antigen}_+_{antibody}",
                                 "{antibody}_blocked.pdb"),
            json = os.path.join(OUTDIR_BASE,
                                  DOCKING_FILES,
                                  "{antigen}_+_{antibody}",
                                  f'{{antibody}}_fv_{SCHEME}.json')
        params:
            fv_blocking = f'{WD}/scripts/mdock-block_pdb.py',
            annotation_scheme = SCHEME,
            outdir = os.path.join(OUTDIR_BASE,
                                  DOCKING_FILES, "{antigen}_+_{antibody}")
        conda:
            f'{WD}/envs/biopython_env.yml'
        shell:
            "python3 {params.fv_blocking} "
            "--multi_chain_pdb {input.mc_pdb} "
            "--single_chain_pdb {input.sc_pdb} "
            "--annotation_scheme {params.annotation_scheme} "
            "--output_dir {params.outdir}"

rule seq2xl:
    # TODO: Should make this run only once per unique antibody-antigen pair.
    input:
        PDB_a = lambda wc: os.path.join(OUTDIR_BASE,
                                        DOCKING_FILES,
                                        f'{wc.antigen}_+_{AB_ID_TO_REP[wc.antibody_id]}',
                                        f'{AB_ID_TO_REP[wc.antibody_id]}.pdb'),
        PDB_b = lambda wc: os.path.join(OUTDIR_BASE,
                                        DOCKING_FILES,
                                        f'{wc.antigen}_+_{AB_ID_TO_REP[wc.antibody_id]}',
                                        f'{wc.antigen}.pdb')
    output:
        os.path.join(OUTDIR_BASE,
                     SUPPORT_FILES,
                     "{antigen}_+_{antibody_id}.xls")
    params:
        seq2xl = f'{WD}/scripts/utils/seq2xl_v1.5.py',
    conda:
        f'{WD}/envs/biopython_env.yml'
    shell:
        "python3 {params.seq2xl} "
        "--pdb_file1 {input.PDB_a} "
        "--pdb_file2 {input.PDB_b} "
        "--output_file {output}"

rule megadock_docking:
    input:
        pdb_a = (rules.block_pdb.output.pdb_a
                 if BLOCK_FV else rules.pdb_handling.output.pdb_a),
        pdb_b = rules.pdb_handling.output.pdb_b
    output:
        os.path.join(OUTDIR_BASE,
                     DOCKING_FILES,
                     "{antigen}_+_{antibody}",
                     "{antigen}_+_{antibody}_megadock.out")
    params:
        megadock = '/opt/MEGADOCK/megadock-gpu',  # GPU
        # megadock = '/opt/MEGADOCK/megadock',  # CPU
        predictions = config['docking_samples'],
        outdir = os.path.join(OUTDIR_BASE,
                              DOCKING_FILES, "{antigen}_+_{antibody}")
    singularity:
        f'{WD}/envs/megadock_4.1.4-gpu.sif'  # GPU
        # f'{WD}/envs/megadock_4.1.4-cpu.sif'  # CPU
    threads:
        max(1, int(workflow.cores/2))
    resources:
        nvidia_gpu=1
    log:
        logfile = os.path.join(OUTDIR_BASE,
                               DOCKING_FILES,
                               '{antigen}_+_{antibody}',
                               '{antigen}_+_{antibody}_megadock.log')
    shell:
        "{params.megadock} "
        "-R {input.pdb_a} "
        "-L {input.pdb_b} "
        "-o {output} "
        "-N {params.predictions} "
        "-D "
        "-e 1.5 "
        "-d 1.5 "
        "&> {log.logfile}"

rule megadock_ensemble_generation:
    input:
        rules.megadock_docking.output
    output:
        os.path.join(OUTDIR_BASE,
                     DOCKING_FILES,
                     "{antigen}_+_{antibody}",
                     "{antigen}_+_{antibody}_megadock.done")
    params:
        decoygen = f'{WD}/../bin/decoygen',
        receptor = rules.pdb_handling.output.pdb_a,  # Antibody.
        ligand = rules.pdb_handling.output.pdb_b,  # Antigen.
        decoys = config['model_samples'],
        outdir = os.path.join(OUTDIR_BASE,
                              DOCKING_FILES, "{antigen}_+_{antibody}"),
        ensemble_prefix = "ensemble/{antigen}_+_{antibody}",
    shell:
        "mkdir -p {params.outdir}/ensemble "
        "&& "
        "for i in {{1..{params.decoys}}}; do"
        "  {params.decoygen}"
        "  {params.outdir}/tmp.pdb"  # Decoy filename.
        "  {params.ligand}"  # Ligand=antigen.
        "  {input}"
        "  $i"  # Model number.
        "  ; head -n-1 {params.receptor} "  # Receptor=antibody.
        "    > {params.outdir}/{params.ensemble_prefix}_$i.pdb"
        "  ; cat {params.outdir}/tmp.pdb"
        "    >> {params.outdir}/{params.ensemble_prefix}_$i.pdb"
        "  ; rm {params.outdir}/tmp.pdb; done"
        "&& touch {output}"
