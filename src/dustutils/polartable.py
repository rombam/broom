import numpy as np
import pandas as pd

from itertools import groupby
from functools import cached_property
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
        self._cols = ['re', 'm', 'alpha', 'cl', 'cd', 'cm']
        self._table = pd.DataFrame(columns=self._cols)

        self._cl = pd.DataFrame(columns=['re', 'm', 'alpha', 'cl'])
        self._cd = pd.DataFrame(columns=['re', 'm', 'alpha', 'cd'])
        self._cm = pd.DataFrame(columns=['re', 'm', 'alpha', 'cm'])

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
        """Update and return the table as a DataFrame.

        Returns
        -------
        table : pd.DataFrame
            DataFrame containing all the data of the C81 table.

        """
        temp = pd.concat([self.cl, self.cd, self.cm]).reset_index(drop=True)
        self._table = temp.groupby(['re', 'm', 'alpha']).mean().reset_index()
        return self._table

    @table.setter
    def table(self, df):
        """Setter method for the table attribute.

        Note
        ----
        Setting the table will also overwrite the cl, cd, and cm attributes.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing lift coefficient Cl, drag coefficient Cd, 
            pitching moment coefficient Cm at different angles of attack alpha at
            one or several Mach and Reynolds numbers.
            Required columns: ['re', 'm', 'alpha', 'cl', 'cd', 'cm']
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
        self.cl = df[['re', 'm', 'alpha', 'cl']].dropna()
        self.cd = df[['re', 'm', 'alpha', 'cd']].dropna()
        self.cm = df[['re', 'm', 'alpha', 'cm']].dropna()

    @property
    def cl(self):
        """Return the lift coefficient table.

        Returns
        -------
        cl : pd.DataFrame
            DataFrame containing the lift coefficient Cl at different angles of
            attack alpha at one or several Mach and Reynolds numbers.
            Columns: ['re', 'm', 'alpha', 'cl']
            Data types: [float, float, float, float]

        """
        return self._cl

    @cl.setter
    def cl(self, df):
        """Setter method for the cl attribute.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing lift coefficient Cl at different angles of
            attack alpha at one or several Mach and Reynolds numbers.
            Required columns: ['re', 'm', 'alpha', 'cl']
            Data types: [float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        cols = ['re', 'm', 'alpha', 'cl']
        if not set(df.columns) == set(cols):
            raise KeyError(f'Input columns are {df.columns}, should be {cols}')

        self._cl = df

    @property
    def cd(self):
        """Return the drag coefficient table.

        Returns
        -------
        cd : pd.DataFrame
            DataFrame containing the drag coefficient Cd at different angles of
            attack alpha at one or several Mach and Reynolds numbers.
            Columns: ['re', 'm', 'alpha', 'cd']
            Data types: [float, float, float, float]

        """
        return self._cd

    @cd.setter
    def cd(self, df):
        """Setter method for the cd attribute.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing drag coefficient Cd at different angles of
            attack alpha at one or several Mach and Reynolds numbers.
            Required columns: ['re', 'm', 'alpha', 'cd']
            Data types: [float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        cols = ['re', 'm', 'alpha', 'cd']
        if not set(df.columns) == set(cols):
            raise KeyError(f'Input columns are {df.columns}, should be {cols}')

        self._cd = df

    @property
    def cm(self):
        """Return the pitching moment coefficient table.

        Returns
        -------
        cm : pd.DataFrame
            DataFrame containing the pitching moment coefficient Cm at different
            angles of attack alpha at one or several Mach and Reynolds numbers.
            Columns: ['re', 'm', 'alpha', 'cm']
            Data types: [float, float, float, float]

        """
        return self._cm

    @cm.setter
    def cm(self, df):
        """Setter method for the cm attribute.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing pitching moment coefficient Cm at different
            angles of attack alpha at one or several Mach and Reynolds numbers.
            Required columns: ['re', 'm', 'alpha', 'cm']
            Data types: [float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        cols = ['re', 'm', 'alpha', 'cm']
        if not set(df.columns) == set(cols):
            raise KeyError(f'Input columns are {df.columns}, should be {cols}')

        self._cm = df

    def append_point(self, re, m, alpha, cl=None, cd=None, cm=None):
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
            If no input coefficient is given.
            If inputs do not have consistent dimensions.

        """
        vnames = ['cl', 'cd', 'cm']
        if all(coeff is None for coeff in [cl, cd, cm]):
            raise ValueError('At least one coefficient must be provided.')

        def all_equal(iterable):
            """Return True if all the elements are equal to each other.

            Returns
            -------
            all_equal : bool
                Whether all of the elements of the iterable are equal.

            """
            g = groupby(iterable)
            return next(g, True) and not next(g, False)

        n_points = len(alpha)
        re_data = np.array(n_points*[re])
        m_data = np.array(n_points*[m])

        for idx, coeff in enumerate([cl, cd, cm]):
            if coeff is not None:
                cname = vnames[idx]
                if not all_equal([len(alpha), len(coeff)]):
                    raise ValueError('Inconsistent input dimensions: ' /
                                     f'len(alpha) = {len(alpha)}' /
                                     f', len({cname}) = {len(coeff)}')

                input_data = {'re': re_data,
                              'm': m_data,
                              'alpha': alpha,
                              cname: coeff}
                table = pd.DataFrame.from_dict(input_data).dropna()
                if cname == 'cl':
                    self.cl = pd.concat([self.cl, table]).drop_duplicates(ignore_index=True)
                elif cname == 'cd':
                    self.cd = pd.concat([self.cd, table]).drop_duplicates(ignore_index=True)
                elif cname == 'cm':
                    self.cm = pd.concat([self.cm, table]).drop_duplicates(ignore_index=True)

    def append_table(self, table):
        """Append lift, drag and moment coefficient curves in a DataFrame, for an
        arbitrary amount of Reynolds and Mach numbers.

        Parameters
        ----------
        df : pd.DataFrame
            Pandas DataFrame containing lift coefficient Cl, drag coefficient Cd, 
            pitching moment coefficient Cm at different angles of attack alpha at
            one or several Mach and Reynolds numbers.
            Required columns: ['re', 'm', 'alpha', 'cl', 'cd', 'cm']
            Data types: [float, float, float, float, float, float]

        Raises
        ------
        KeyError
            If the input DataFrame does not have the correct columns.

        """
        self.table = pd.concat([self.table, table]).drop_duplicates(ignore_index=True)

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

        cl = deepcopy(self.cl.query(f"re == {re}"))
        cd = deepcopy(self.cd.query(f"re == {re}"))
        cm = deepcopy(self.cm.query(f"re == {re}"))
        name = self.name
        ml = len(cl['m'].unique())
        md = len(cd['m'].unique())
        mm = len(cm['m'].unique())
        al = len(cl['alpha'].unique())
        ad = len(cd['alpha'].unique())
        am = len(cm['alpha'].unique())
        m = [ml, md, mm]
        coeffs = [cl, cd, cm]
        table = []
        table.append(f'{comments[0]}')
        table.append(f'{comments[1]}')
        table.append(f'{name:<33}   {ml:02}{al:02}{md:02}{ad:02}{mm:02}{am:02}')
        for idx, c in enumerate(['cl', 'cd', 'cm']):
            machs = coeffs[idx]["m"].unique()
            table.append(f'{"":>10}'+"".join([f'{mn:>10.4f}' for mn in machs]))
            for a in coeffs[idx]["alpha"].unique():
                coefdata = coeffs[idx].query(f"alpha == {a}")[c].to_list()
                table.append(f'{a:>10.2f}' + "".join([f'{ci:>10.4f}' for ci in coefdata]))

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
