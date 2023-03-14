import numpy as np
import pandas as pd

from itertools import groupby
from copy import deepcopy


class PolarTable():
    """An airfoil polar table. Contains aerodynamic coefficient information at several
    Reynolds and Mach numbers."""
    # TODO: think about the implementation with aerotools Polar

    def __init__(self, name='', desc=''):
        """Constructor method.

        Parameters
        ----------
        name : str, optional
            Name of the airfoil table.
        desc : str, optional
            Description of the airfoil table.

        """
        self.name = name
        self.desc = desc

        # Initialize data table
        self._cols = ['Re', 'M', 'alpha', 'Cl', 'Cd', 'Cm']
        self._table = pd.DataFrame(columns=self._cols)

    @property
    def name(self):
        """Return the name of the table.

        Returns
        -------
        name : str
            Name of the table.

        """
        return self._name

    @name.setter
    def name(self, value):
        """Set the name of the table.

        Parameters
        ----------
        value : str
            Name of the table.

        Raises
        ------
        TypeError
            If input value is not a string.

        """
        if not isinstance(value, str):
            raise TypeError(f'value is of type {type(value)}, should be str')

        self._name = value

    @property
    def desc(self):
        """Return the description of the table.

        Returns
        -------
        desc : str
            Description of the table.

        """
        return self._desc

    @desc.setter
    def desc(self, value):
        """Set the description of the table.

        Parameters
        ----------
        value : str
            Description of the table.

        Raises
        ------
        TypeError
            If input value is not a string.

        """
        if not isinstance(value, str):
            raise TypeError(f'value is of type {type(value)}, should be str')

        self._desc = value

    @property
    def table(self):
        """Return the table as a DataFrame.

        Returns
        -------
        table : pd.DataFrame
            DataFrame containing all the data of the C81 table.

        """
        return self._table

    @table.setter
    def table(self, df):
        """Setter method for the table attribute.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing lift coefficient Cl, drag coefficient Cd, 
            pitching moment coefficient Cm at different angles of attack alpha at
            one or several Mach and Reynolds numbers.
            Required columns: ['Re', 'M', 'alpha', 'Cl', 'Cd', 'Cm']
            Data types: [float, float, float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        cols = self._cols
        if not set(df.columns) == set(cols):
            raise KeyError(f'Input columns are {df.columns}, should be {cols}')

        self._table = df

    def append_point(self, re, m, alpha, cl, cd, cm):
        """Append lift, drag and moment coefficient curves for a single Re and M
        condition to the table.

        Parameters
        ----------
        re : float
            Reynolds number.
        m : float
            Mach number.
        alpha : array-like(float)
            Angles of attack. Units: degrees.
        cl : array-like(float)
            Lift coefficients.
        cd : array-like(float)
            Drag coefficients.
        cm : array-like(float)
            Moment coefficients.

        Raises
        ------
        ValueError
            If inputs do not have consistent dimensions.

        """

        def all_equal(iterable):
            """Return True if all the elements are equal to each other.

            Returns
            -------
            all_equal : bool
                Whether all of the elements of the iterable are equal.

            """
            g = groupby(iterable)
            return next(g, True) and not next(g, False)

        if not all_equal([len(alpha), len(cl), len(cd), len(cm)]):
            raise ValueError(f'Inconsistent input dimensions: len(alpha) = {len(alpha)}' /
                             f', len(cl) = {len(cl)}, len(cd) = {len(cd)}' /
                             f', len(cm) = {len(cm)}')
        n_points = len(alpha)
        re_data = np.array(n_points*[re])
        m_data = np.array(n_points*[m])

        # Assembling the input dictionary
        input_data = {'Re': re_data,
                      'M': m_data,
                      'alpha': alpha,
                      'Cl': cl,
                      'Cd': cd,
                      'Cm': cm}
        table = pd.DataFrame.from_dict(input_data)

        self.table = pd.concat([self.table, table], ignore_index=True)

    def append_table(self, table):
        """Append lift, drag and moment coefficient curves in a DataFrame, for an
        arbitrary amount of Reynolds and Mach numbers.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing lift coefficient Cl, drag coefficient Cd, 
            pitching moment coefficient Cm at different angles of attack alpha at
            one or several Mach and Reynolds numbers.
            Required columns: ['Re', 'M', 'alpha', 'Cl', 'Cd', 'Cm']
            Data types: [float, float, float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        self.table = pd.concat([self.table, table], ignore_index=True)

    def c81_text(self, re, comments=('', '')):
        """Return a single-Reynolds number table in C81 string format.

        Parameters
        ----------
        re : str, Path
            Path to the output file.
        comments : array-like(str)
            Two element iterable containing two string comments to identify the polar.

        Returns
        -------
        table : list(str)
            Table string lines in a list.

        """
        data = deepcopy(self._table.query(f"Re == {re}"))
        name = self.name
        table = []
        ml = md = mm = len(data['M'].unique())
        al = ad = am = len(data['alpha'].unique())
        m = data['M'].unique()

        table.append(f'{comments[0]}')
        table.append(f'{comments[1]}')
        table.append(f'{name:<33}   {ml:02}{al:02}{md:02}{ad:02}{mm:02}{am:02}')
        for c in ['Cl', 'Cd', 'Cm']:
            table.append(f'{"":>10}'+"".join([f'{mn:>10.4f}' for mn in m]))
            for a in data["alpha"].unique():
                coeffs = data.query(f"alpha == {a}")[c].to_list()
                table.append(f'{a:>10.2f}' + "".join([f'{ci:>10.4f}' for ci in coeffs]))

        return table

    def to_dust(self, filename, re):
        """Save the current PolarTable object as a DUST C81-formatted text file.

        Parameters
        ----------
        filename : str, Path
            Path to the output file.
        re : array-like(float/int)
            List or single Reynolds number to include in the C81 table. Should include
            at least two Reynolds numbers, in ascending order.

        """
        dust_table = []

        # Header: only the number of Reynolds number entries matters
        dust_table.extend([f'{len(re):<7n} {0:<7n} {0:<7n}',
                           f'{0:<7.2f} {1:<7.2f}',
                           f'{0.1:<7.4f} {0.1:<7.4f}'])

        # Body
        for reynolds in re:
            table = self.c81_text(reynolds, comments=('COMMENT #1', f'{reynolds} 0.0'))
            dust_table.extend(table)

        with open(filename, 'w+') as f:
            f.writelines([f'{line}\n' for line in dust_table])

        print(f'C81-DUST airfoil table output to {filename}.')
