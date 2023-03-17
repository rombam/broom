import numpy as np

from copy import deepcopy


class Printable:
    """Helper class to add readable string and dictionary exports to any class.

    Note
    ----
    By default, includes all of the child class' attributes in the output string/dict.

    """

    def __str__(self):
        """Return a readable string containing attribute information.

        Returns
        -------
        class_str : str
            Readable string form of the class object.

        """
        clsdict = self.__dict__
        clsname = self.__class__.__name__
        return '\n'.join([f'{clsname}(']+[f'{k} = {v}' for k, v in clsdict.items()]
                         + [')'])

    def to_dict(self, ignore=[]):
        """Return the class' attributes as a dictionary.

        Returns
        -------
        clsdict : dict
            Dictionary with class attributes as keys and their contents as values.
        ignore : List[str]
            List of attribute names to ignore when converting to a dictionary.

        """
        clsdict = deepcopy(self.__dict__)
        for k, v in clsdict.items():
            if k not in ignore:
                if isinstance(v, Printable):
                    clsdict[k] = v.to_dict()
                elif isinstance(v, list) and isinstance(v[0], Printable):
                    clsdict[k] = [val.to_dict() for val in v]
        return clsdict

    def to_fort(self, ignore=[], indent=0):
        """Return a Fortran-formatted class string.

        Parameters
        ----------
        ignore : str, List[str], optional
            Name or list of attribute names to be ignored for printing.
        indent : int, optional
            Number of spaces to indent each string line by.

        Returns
        -------
        fort_str : str
            Fortran-formatted list of strings.

        TODO: clean this up. ignore and indent are not needed in to_fort, just in
        _fort_strs.
        """
        if isinstance(ignore, str):
            ignore = [ignore]

        fort_strs = self._fort_strs(ignore=ignore, indent=indent)

        connector = f'\n{indent*" "}'
        return indent*" "+connector.join(fort_strs)

    def _fort_strs(self, ignore=[], indent=0):
        """Return a list of Fortran-formatted strings for all the class attributes.

        Parameters
        ----------
        ignore : str, List[str], optional
            Name or list of attribute names to be ignored for printing.
        indent : int, optional
            Number of spaces to indent each string line by.

        Returns
        -------
        fort_strs : List[str]
            List of Fortran-formatted strings.

        """
        indent = f'{indent*" "}'
        fort_strs = []
        clsdict = self.__dict__
        for k, v in clsdict.items():
            if k not in ignore:
                if v is not None:
                    if isinstance(v, Printable):
                        fort_strs.append('')
                        fort_strs.extend(v._fort_strs())
                    elif isinstance(v, list):
                        if isinstance(v[0], Printable):
                            for val in v:
                                fort_strs.extend(val._fort_strs())
                        elif isinstance(v[0], str):
                            fort_strs.extend([f'{indent}{k} = {fortranslate(val)}'
                                              for val in v])
                    else:
                        fort_strs.append(f'{indent}{k} = {fortranslate(v)}')
        return fort_strs


def vec_fortstr(values):
    """Output a Fortran-formatted string representation of a 1D array.

    Parameters
    ----------
    values : List[numbers], np.ndarray
        (N, ) 1D number array.
    Returns
    -------
    string : str
        Fortran-formatted string representation of the array.

    """
    return f'(/ {", ".join([str(num) for num in values])} /)'


def bool_fortstr(value):
    """Output a Fortran-formatted string representation of a Boolean value.

    Parameters
    ----------
    value : bool
        Boolean value to translate.

    Returns
    -------
    string : str
        Fortran-formatted boolean value.

    """
    return 'T' if value else 'F'


def fortranslate(value):
    """Print a value in the Fortran format equivalent.

    Parameters
    ----------
    value : Any
        Variable to be converted to Fortran string format.

    Returns
    -------
    fort_str : str
        Fortran-formatted string.

    """
    if isinstance(value, (list, np.ndarray)):
        return vec_fortstr(value)
    elif isinstance(value, bool):
        return bool_fortstr(value)
    else:
        return str(value)
