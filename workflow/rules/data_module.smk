#!/usr/bin/env python3
"""Based on the TX-MS web-app implementation Cheetah-MS.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


rule msconvert:
    # TODO:
    #   - Sort proper thermoRawFileParser inclusion
    input:
        ms_file = os.path.join(config['ms_files'],
                           f'{{sample}}{config["ms_ext"]}')
    output:
        mzml = os.path.join(OUTDIR_BASE,
                            DATA_FILES, 'mzML_files', '{sample}.mzML')
    params:
        outdir = os.path.join(OUTDIR_BASE, DATA_FILES, 'mzML_files'),
        thermoRawFileParser = config['thermoRawFileParser']
    conda:
        f'{WD}/envs/thermoRawFileParser_env.yml'
    shell:
        "[[ $(sed 's/^.*\.//' <(echo {input.ms_file})) == mzML ]]"
        " &&"
        " ln -s {input.ms_file} {output.mzml}"
        " ||"
        " {params.thermoRawFileParser}"
        " -i={input.ms_file}"
        " -o={params.outdir}"
        " -f=2"

# Need to add this rule to make filtering more efficient, and only filter once
# per antibody-antigen-sample combination.
# rule filter_spectra:
#     input:
#         mzml = rules.msconvert.output.mzml,
#         xls = rules.seq2xl.output
#     output:
#         filt_mzml = os.path.join(OUTDIR_BASE,
#                                  DATA_FILES,
#                                  'filtered_mzML',
#                                  '{sample}_filtered.mzML')
#     params:
#         filter_script = f'{WD}/scripts/filter_mzml.py',
#         outdir = os.path.join(OUTDIR_BASE,
#                               DATA_FILES,
#                               'filtered_mzML')
#     conda:
#         f'{WD}/envs/pyteomics_env.yml'
#     shell:
#         "python3 {params.filter_script} "
#         "--mzml_file {input.mzml} "
#         "--output_file {output.filt_mzml}"

rule spectra_annotation:
    # TODO: Speed up!
    input:
        mzml = rules.msconvert.output.mzml,
        xls = rules.seq2xl.output
    output:
        sql = os.path.join(OUTDIR_BASE,
                           SCORE_FILES,
                           "{sample}",
                           "{antigen}", "{antibody}", "spectra_annotation.sql"),
        img_dir = directory(os.path.join(OUTDIR_BASE,
                                         SCORE_FILES,
                                         "{sample}",
                                         "{antigen}",
                                         "{antibody}", "top_spectra")),
        top_xls = os.path.join(OUTDIR_BASE,
                               SCORE_FILES,
                               "{sample}",
                               "{antigen}", "{antibody}", "top_xls.txt"),
        # Move to previous filtering rule to reduce redundancy:
        mgf_filt = os.path.join(OUTDIR_BASE,
                                DATA_FILES,
                                "{sample}_{antibody}_{antigen}_filtered.mzML")
    params:
        ms_script = f'{WD}/scripts/sidewinder-ms_v2.py',
        x_linker = get_linker_num,
        mass_delta_cutoff = 0.01,
        outdir = os.path.join(OUTDIR_BASE,
                              SCORE_FILES,
                              "{sample}", "{antigen}", "{antibody}")
    conda:
        f'{WD}/envs/pyteomics_env.yml'
    shell:
        "python3 {params.ms_script} "
        "--mzml_file {input.mzml} "
        # "--mgf_file {input.mgf} "
        "--xl_file {input.xls} "
        "--x_linker {params.x_linker} "
        "--mass_delta_cutoff {params.mass_delta_cutoff} "
        "--output_dir {params.outdir}"
