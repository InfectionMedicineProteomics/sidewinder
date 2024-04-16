#!/usr/bin/env python3


from pathlib import Path
from typing import List

import matplotlib.pyplot as plt


def fig_maker(main_spectra: List[float],
              main_intensity: List[float],
              covered_frags: List[str],
              covered_mz: List[float],
              peptide1: str,
              peptide2: str,
              xl: str,
              num_mgf: str,
              delta: float,
              output_dir: Path):
    """Generates a plot visualizing identified fragments for a cross-linked peptide.

    This function takes information about the MS/MS spectrum, identified
    fragment ions, cross-linked peptides, and other analysis parameters, and
    creates a plot to visualize the matched fragments. The plot includes:

    - The original MS/MS spectrum (in yellow).
    - Identified fragment ions for the first peptide (in red).
    - Identified fragment ions for the second peptide (in blue).
    - Identified fragment ions potentially from both peptides (in green).
    - Fragment ion labels are displayed next to their corresponding m/z values.

    The plot is saved as a PNG image file named after the MGF spectrum number
    and the cross-link sequence.

    Args:
        main_spectra (List[float]): List of m/z values from the MS/MS spectrum.
        main_intensity (List[float]): List of intensity values from the MS/MS spectrum.
        covered_frags (List[str]): List of identified fragment ion sequences.
        covered_mz (List[float]): List of m/z values for the identified fragment ions.
        peptide1 (str): Amino acid sequence of the first peptide in the cross-link.
        peptide2 (str): Amino acid sequence of the second peptide in the cross-link.
        xl (str): The cross-linked peptide sequence in Kojak format.
        num_mgf (str): The MGF spectrum number (used for naming the output file).
        delta (float): The mass tolerance used for matching fragment ions.
        output_dir (Path): Path to the directory where the plot will be saved.
    """
    plt.switch_backend('agg')

    plot_mz_list = []

    plot_intensity_list = []

    plot_p1_mz = []

    plot_p1_int = []

    plot_p1_frag = []

    plot_p2_mz = []

    plot_p2_int = []

    plot_p2_frag = []

    plot_p1p2_mz = []

    plot_p1p2_int = []

    plot_p1p2_frag = []

    for y_num, ymz in enumerate(covered_mz):

        for mz_num, anymz in enumerate(main_spectra):

            min_MZ = abs(anymz - ymz)

            if min_MZ <= delta:

                if main_intensity[mz_num] != 0.0:

                    plot_intensity_list.append(main_intensity[mz_num])

                    plot_mz_list.append(anymz)

                    if covered_frags[y_num].replace("+","") in peptide1:

                        plot_p1_mz.append(anymz)

                        plot_p1_int.append(main_intensity[mz_num])

                        plot_p1_frag.append(covered_frags[y_num])

                    elif covered_frags[y_num].replace("+","") in peptide2:

                        plot_p2_mz.append(anymz)

                        plot_p2_int.append(main_intensity[mz_num])

                        plot_p2_frag.append(covered_frags[y_num])

                    else:

                        plot_p1p2_mz.append(anymz)

                        plot_p1p2_int.append(main_intensity[mz_num])

                        plot_p1p2_frag.append(covered_frags[y_num])

                    break

    spec_fig = plt.figure()

    ax = spec_fig.add_subplot(111)

    plt.title(xl)

    plt.xlim(0, max(covered_mz)+100)

    plt.ylim(0, max(main_intensity)+100)

    plt.xlabel('m/z')

    plt.ylabel('intensity')

    #plt.Axes.autoscale_view(tight=None, scalex=True, scaley=True)

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    plt.vlines(main_spectra,0,main_intensity,colors='y',linestyles='solid')

    plt.vlines(plot_p1_mz,0,plot_p1_int,color='r',linestyles='solid')

    plt.vlines(plot_p2_mz,0,plot_p2_int,color='b',linestyles='solid')

    plt.vlines(plot_p1p2_mz,0,plot_p1p2_int,color='g',linestyles='solid')

    for num1, frags1 in enumerate(plot_p1_frag):

        plt.text(plot_p1_mz[num1],
                 plot_p1_int[num1]+int(min(main_intensity)),
                 frags1,
                 rotation='vertical',
                 fontsize=4,
                 horizontalalignment='center',
                 verticalalignment='bottom', color='red')

    for num2, frags2 in enumerate(plot_p2_frag):

        plt.text(plot_p2_mz[num2],
                 plot_p2_int[num2]+int(min(main_intensity)),
                 frags2,
                 rotation='vertical',
                 fontsize=4,
                 horizontalalignment='center',
                 verticalalignment='bottom', color='blue')

    for num12, frags12 in enumerate(plot_p1p2_frag):

        plt.text(plot_p1p2_mz[num12],
                 plot_p1p2_int[num12]+int(min(main_intensity)),
                 frags12,
                 rotation='vertical',
                 fontsize=4,
                 horizontalalignment='center',
                 verticalalignment='bottom', color='green')

    # print ("************RESULT SECTION**************")

    # print ("covered_frags_list:\n", covered_frags)

    # print ("covered_mz_list:\n", covered_mz)

    # print ("hypo_intensity_list:\n", plot_intensity_list)

    # print ("Storing figure: ", str(num_mgf)+xl+".png")

    spec_image_file = output_dir / f'{num_mgf}{xl}.png'

    plt.savefig(spec_image_file, format='png', dpi=600)

    plt.close()
