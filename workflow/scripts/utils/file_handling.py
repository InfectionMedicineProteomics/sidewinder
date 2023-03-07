from pathlib import Path
from typing import List


def case_ins_ext_glob(directory: Path,
                      extensions: List[str] = ['.*']) -> List[Path]:
    """Case-insensitive file globing with specified extension.

    Takes a path object and a list of extensions (specify like
    ['.tsv', '.csv']) as input.

    Parameters
    ----------
    directory : pathlib.Path
    extensions : list of strings
        Default is to check for all extensions.

    Returns
    -------
    List of pathlib.Path objects
    """

    dir_files = [*directory.glob('*.*')]

    files = []

    for ext in extensions:

        files.extend(

            [x for x in dir_files if x.suffix.lower() == ext.lower()]

        )

    return files
