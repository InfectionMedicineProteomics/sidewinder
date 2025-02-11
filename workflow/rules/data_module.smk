#!/usr/bin/env python3
"""Based on the TX-MS web-app implementation Cheetah-MS.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


rule msconvert:
    # TODO:
    #   - Sort proper thermoRawFileParser inclusion
    input:
        ms2 = os.path.join(config["ms2_files"],
                           f"{{sample}}{config['ms2_ext']}")
    output:
        mzml = os.path.join(OUTDIR_BASE,
                           sample_files, "ms2_convert", "{sample}.mzML")
    params:
        outdir = os.path.join(OUTDIR_BASE, sample_files, "ms2_convert"),
        thermoRawFileParser = config['thermoRawFileParser']
    conda:
        f'{WD}/envs/mono_env.yml'
    shell:
        "[[ $(sed 's/^.*\.//' <(echo {input.ms2})) == mzML ]]"
        " &&"
        " ln -s {input.ms2} {output.mzml}"
        " ||"
        " mono {params.thermoRawFileParser}"
        " -i={input.ms2}"
        " -o={params.outdir}"
        " -f=2"

rule ms2:
    # TODO: Speed up!
    input:
        mzml = rules.msconvert.output.mzml,
        xls = rules.seq2xl.output
    output:
        sql = os.path.join(OUTDIR_BASE,
                           txms_files,
                           "{sample}",
                           "{binder}", "{target}", "ms2_results.sql"),
        img_dir = directory(os.path.join(OUTDIR_BASE,
                                         txms_files,
                                         "{sample}",
                                         "{binder}",
                                         "{target}", "top_spectra")),
        top_xls = os.path.join(OUTDIR_BASE,
                               txms_files,
                               "{sample}",
                               "{binder}", "{target}", "top_xls.txt"),
        mgf_filt = os.path.join(OUTDIR_BASE,
                                txms_files,
                                "{sample}",
                                "{binder}",
                                "{target}", "{sample}_filtered.mzML")
    params:
        ms2_script = f'{WD}/scripts/sidewinder-ms_v2.py',
        x_linker = get_linker_num,
        mass_delta_cutoff = 0.01,
        outdir = os.path.join(OUTDIR_BASE,
                              txms_files, "{sample}", "{binder}", "{target}")
    conda:
        f'{WD}/envs/pyteomics_env.yml'
    shell:
        "python3 {params.ms2_script} "
        "--mzml_file {input.mzml} "
        # "--mgf_file {input.mgf} "
        "--xl_file {input.xls} "
        "--x_linker {params.x_linker} "
        "--mass_delta_cutoff {params.mass_delta_cutoff} "
        "--output_dir {params.outdir}"
