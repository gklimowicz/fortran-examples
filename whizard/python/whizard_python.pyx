"""
WHIZARD Python Interface.

This module provides a WHIZARD Python interface using the C API of WHIZARD.

Warning
-------
Do not instantiate several Whizard instances.

See also
--------
Manual : API section
api.nw

Note
----
We implement the different WHIZARD handles, i.e. handling WHIZARD or event samples, as class variables.
Therefore, different instances of WhizardSample or Whizard share a common handle object.

Furthermore, the Python interface makes heavy use of UTF-8 encoded strings and goes all the way to ensure that these strings are converted to binary strings (including to append a null character), i.e. ASCII-based, before handing them to the C API.

Example
-------
>>> import whizard
>>> wz = whizard.Whizard()
>>> wz.options('--model', 'sm')
"""
from cpython cimport array
cimport cwhizard

import array

cdef class WhizardSample:
    """
    WHIZARD event stream sampler.

    Provides access to WHIZARD's internal stream of event samples.
    Access to the event stream is granted either indirectly as Python iterator
    or directly using the WHIZARD API.

    Attributes
    ----------
    _c_sample : cwhizard.sample_handle_t
        WHIZARD Sample Handler.
    _begin : int
        Start of event stream iterator.
    _end : int
        End of event stream iterator.
    _it : int
        Current status of event stream iterator.

    Methods
    -------
    __iter__
        Returns itself as an iterator object.
    __next__
        Implement iterator and return current iterator index.
    open
        Open event stream.
    next_event
        Iterate and access next event.
    close
        Close event stream and destroy object.
    get_event_index
        Returns the current event index with regard to WHIZARD.
    get_process_index
        Returns the current process index with regard to WHIZARD.
    get_process_id
        Returns the process id as defined by the process command.
    get_fac_scale
        Returns the factorization scale.
    get_alpha_s
        Returns alpha_s.
    get_weight
        Returns the event weight, among others, includes flux factor and phasespace weight.
    get_sqme
        Returns the squared matrix element of the event.

    Examples
    --------
    The sample class provides itself as an iterator:
    >>> sample = wz.new_sample("process_name_1")
    >>> for it in sample:
    >>>		idx = sample.get_event_index()
    >>>     # [...]

    The direct approach using the WHIZARD API looks like:
    >>> sample = wz.new_sample("process_name_1")
    >>> it_begin, it_end = sample.open()
    >>> for it in range(it_begin, it_end + 1):
    >>>     idx = sample.get_event_index()
    >>>
    >>> sample.close() # alternative: del(sample)
    """

    cdef cwhizard.sample_handle_t _c_sample
    cdef readonly int _begin, _end, _it

    def __cinit__(self, Whizard wz, str name):
        cdef array.array c_name = array.array('b', name.encode('UTF-8') + b'\0')
        cwhizard.whizard_new_sample(
            &wz._c_whizard,
            c_name.data.as_chars,
            &self._c_sample)
        self.__check_handle()

    def __dealloc__(self):
        if self._c_sample is not NULL:
            cwhizard.whizard_sample_close(&self._c_sample)

    def __check_handle(self):
        """Check whether WHIZARD has been correctly allocated else raise a MemoryError.
        """
        if self._c_sample is NULL:
            raise MemoryError()

    def __iter__(self):
        it_begin, it_end = self.open()
        self._begin = it_begin
        self._end = it_end
        self._it = self._begin
        return self

    def __next__(self):
        """
        Next item.

        Yields
        ------
        it : int
            Current iteration index in the event stream.

        Raises
        ------
        StopIteration
        """
        if self._it <= self._end:
            self.next_event()
            it = self._it
            self._it = it + 1
            return it
        else:
            self.close()
            raise StopIteration

    cpdef (int, int) open(self):
        """
        Open event stream.

        Returns
        -------
        it_begin : int
            Start index in event stream.
        it_end : int
            End index in event stream.
        """
        cdef int it_begin, it_end
        self.__check_handle()
        cwhizard.whizard_sample_open(&self._c_sample, &it_begin, &it_end)
        return it_begin, it_end

    cpdef next_event(self):
        """Go to the next sample in the event stream.
        """
        self.__check_handle()
        cwhizard.whizard_sample_next_event(&self._c_sample)

    cpdef close(self):
        """
        Close event sampler.

        Delete self.
        """
        del(self)

    cpdef int get_event_index(self):
        """
        Get the event index.

        Returns
        -------
        idx : int
        """
        cdef int idx
        self.__check_handle()
        cwhizard.whizard_sample_get_event_index(&self._c_sample, &idx)
        return idx

    cpdef int get_process_index(self):
        """
        Get process index.

        Returns
        -------
        i_proc : int
        """
        cdef int i_proc
        self.__check_handle()
        cwhizard.whizard_sample_get_process_index(&self._c_sample, &i_proc)
        return i_proc

    cpdef str get_process_id(self):
        """
        Get the process id string.

        Returns
        -------
        proc_id : str
        """
        self.__check_handle()
        cdef int strlen = cwhizard.whizard_sample_get_process_id_len(&self._c_sample) + 1
        cdef array.array c_proc_id = array.array('b', [])
        array.resize(c_proc_id, strlen) # Caveat: inconsistency in the strlen handling!
        cwhizard.whizard_sample_get_process_id(
            &self._c_sample,
            c_proc_id.data.as_chars,
            strlen)
        return c_proc_id.tobytes().decode('UTF-8')[:-1]

    cpdef double get_fac_scale(self):
        """
        Get the factorization scale.

        Returns
        -------
        f_scale : double
            In GeV.
        """
        cdef double f_scale
        self.__check_handle()
        cwhizard.whizard_sample_get_fac_scale(&self._c_sample, &f_scale)
        return f_scale

    cpdef double get_alpha_s(self):
        """
        Get alpha_s.

        Returns
        -------
        alpha_s : double
        """
        cdef double alpha_s
        self.__check_handle()
        cwhizard.whizard_sample_get_alpha_s(&self._c_sample, &alpha_s)
        return alpha_s

    cpdef double get_weight(self):
        """
        Get the event weight.

        The weight is a product of squared matrix element, flux factor, phase-space volume
        and integration weight.

        Returns
        -------
        weight : double
        """
        cdef double weight
        self.__check_handle()
        cwhizard.whizard_sample_get_weight(&self._c_sample, &weight)
        return weight

    cpdef double get_sqme(self):
        """
        Get the squared matrix element.

        Returns
        -------
        sqme : double
        """
        cdef double sqme
        self.__check_handle()
        cwhizard.whizard_sample_get_sqme(&self._c_sample, &sqme)
        return sqme

cdef class Whizard:
    """Whizard Python Class

    Attributes
    ----------
    _c_whizard : cwhizard.whizard_t
        WHIZARD handle.

    Methods
    -------
    option(key: str, value: str)
        Provide commandline-like option, see =whizard --help= for more details on available options.
        Must-be called *before* init.
    init
        Initialize a WHIZARD execution.
    set_double(var: str, value: double)
        Assign a Sindarin variable of type double.
    set_int(var: str, value: int)
        Assign a Sindarin variable of type int.
    set_bool(var: str, value: int)
        Assign a Sindarin variable of type bool.
    set_string(var: str, value: str)
        Assign a Sindarin variable of type string.
    get_double(var: str) : double
        Return the value of a Sindarin variable of type double.
    get_int(var: str) : int
        Return the value of a Sindarin variable of type int.
    get_bool(var: str) : bool
        Return the value of a Sindarin variable of type bool.
    get_string(var: str) : str
        Return the value of a Sindarin variable of type string.
    flv_string(pdg: int) : str
        Convert a PDG code to the associated particle name in the WHIZARD model.
    flv_array_string(pdg: list(int)) : list(str)
        Convert a list of PDG codes to their associated particle names in the WHIZARD model, respectively.
    command(cmd: str)
        Execute a Sindarin command.
    get_integration_result(proc_id: str) : (double, double)
        Retrieve the estimate of the integral and the estimated Monte Carlo error for process with proc_id.
    new_sample(name: str) : WhizardSample
        Allocate a WHIZARD sample handler accessing (mulitple) event stream(s).

    Notes
    -----
    After instantiating a Whizard object, _commandline_ options can be set, however, only *before* a call to init.

    Examples
    --------
    >>> import whizard
    >>>
    >>> wz = whizard.Whizard()
    >>> wz.option("logfile", "python_example.log") # set commandline options before initialization
    >>> wz.init()
    >>>
    >>> wz.set_int("seed", 1234)
    >>> wz.command("sqrts = 1000")
    >>>
    >>> sample = wz.new_sample("process_name")

    In order to define a process, run the process command:
    >>> wz.command("process eett = e1, E1 => t, T")

    In all subsequent calls, we refer to the above process as _eett_.
    """
    cdef cwhizard.whizard_t _c_whizard

    def __cinit__(self):
        # Only access cdef fields in self
        cwhizard.whizard_create(&self._c_whizard)
        if self._c_whizard is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_whizard is not NULL:
            cwhizard.whizard_final(&self._c_whizard)

    cpdef option(self, str key, str value):
        """Set commandline options.

        See also
        --------
        For more details: whizard --help.
        """
        # https://docs.cython.org/en/latest/src/tutorial/strings.html#encoding-text-to-bytes
        # https://docs.python.org/3/library/array.html
        # https://cython.readthedocs.io/en/latest/src/tutorial/array.html
        # The string handling function string_c2f expects a '\0'-terminated character string.
        # Therefore, we need to append a null character (as binary).
        cdef array.array c_key = array.array('b', key.encode('UTF-8') + b'\0')
        cdef array.array c_value = array.array('b', value.encode('UTF-8') + b'\0')
        cwhizard.whizard_option(
            &self._c_whizard,
            c_key.data.as_chars,
            c_value.data.as_chars)

    cpdef init(self):
        """Initialize a WHIZARD run.
        """
        cwhizard.whizard_init(&self._c_whizard)

    cpdef set_double(self, str var, double value):
        """
        Set a Sindarin variable of type double.

        Parameters
        ---------
        var : str
            Sindarin variable name.
        value : double
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cwhizard.whizard_set_double(
            &self._c_whizard,
            c_var.data.as_chars,
            value)

    cpdef set_int(self, str var, int value):
        """
        Set a Sindarin variable of type int.

        Parameters
        ----------
        var : str
            Sindarin variable name.
        value : int
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cwhizard.whizard_set_int(
            &self._c_whizard,
            c_var.data.as_chars,
            value)

    cpdef set_bool(self, str var, bint value):
        """
        Set a Sindarin variable of type bool.

        Parameters
        ----------
        var : str
            Sindarin variable name.
        value : bool
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cwhizard.whizard_set_bool(
            &self._c_whizard,
            c_var.data.as_chars,
            value)

    cpdef set_string(self, str var, str value):
        """
        Set a Sindarin variable of type string.

        Parameters
        ----------
        var : str
            Sindarin variable name.
        valiue : str
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cdef array.array c_value = array.array('b', value.encode('UTF-8') + b'\0')
        cwhizard.whizard_set_char(
            &self._c_whizard,
            c_var.data.as_chars,
            c_value.data.as_chars)

    cpdef double get_double(self, str var) except? 0:
        """
        Get the value of a Sindarin variable of type double.

        Parameters
        ----------
        var : str
            Sindarin variable name.

        Returns
        -------
        double
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cdef double value
        cdef bint err = cwhizard.whizard_get_double(
            &self._c_whizard,
            c_var.data.as_chars,
            &value)
        if err:
            raise IndexError("Cannot retrieve double variable from WHIZARD.")
        return value

    cpdef int get_int(self, str var) except? 0:
        """
        Get the value of a Sindarin variable of type int.

        Parameters
        ----------
        var : str
            Sindarin variable name.

        Returns
        -------
        int
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cdef int value
        cdef bint err = cwhizard.whizard_get_int(
            &self._c_whizard,
            c_var.data.as_chars,
            &value)
        if err:
            raise IndexError("Cannot retrieve int variable from WHIZARD.")
        return value

    cpdef bint get_bool(self, str var) except? False:
        """
        Get the value of a Sindarin variable of type bool.

        Parameters
        ----------
        var : str
            Sindarin variable name.

        Returns
        -------
        bool
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cdef int value = 0
        cdef bint err = cwhizard.whizard_get_bool(
            &self._c_whizard,
            c_var.data.as_chars,
            &value)
        if err:
            raise IndexError("Cannot retrieve boolean variable from WHIZARD.")
        return value

    cpdef str get_string(self, str var):
        """
        Get the value of a Sindarin variable of type str.

        Parameters
        ----------
        var : str
            Sindarin variable name.

        Returns
        -------
        str
        """
        cdef array.array c_var = array.array('b', var.encode('UTF-8') + b'\0')
        cdef int strlen = cwhizard.whizard_get_char_len(
            &self._c_whizard,
            c_var.data.as_chars)
        cdef array.array c_value = array.array('b', [])
        array.resize(c_value, strlen)
        cdef bint err = cwhizard.whizard_get_char(
            &self._c_whizard,
            c_var.data.as_chars,
            c_value.data.as_chars,
            strlen)
        if err:
            raise IndexError("Cannot retrieve string variable from WHIZARD.")
        # Remove trailing null character from C.
        return c_value.tobytes().decode('UTF-8')[:-1]

    cpdef str flv_string(self, int pdg):
        """
        Retrieve the string representation of a PDG code.

        Returns
        -------
        flv_string : str
        """
        cdef int strlen = cwhizard.whizard_flv_string_len(&self._c_whizard, pdg)
        cdef array.array c_str = array.array('b', [])
        array.resize(c_str, strlen)
        cdef bint err = cwhizard.whizard_flv_string(
            &self._c_whizard,
            pdg,
            c_str.data.as_chars,
            strlen)
        if err:
            raise IndexError("Cannot retrieve string variable from WHIZARD.")
        # Remove trailing null character from C.
        return c_str.tobytes().decode('UTF-8')[:-1]

    cpdef str flv_array_string(self, list flavors):
        """
        Retrieve the string representation for a list of PDG codes.

        Returns
        -------
        flv_array_string : list(str)
        """
        cdef array.array fa = array.array('i', flavors)
        cdef array.array c_str = array.array('b', [])
        cdef int strlen = cwhizard.whizard_flv_array_string_len(
            &self._c_whizard,
            fa.data.as_ints,
            len(fa))
        array.resize(c_str, strlen)
        cdef bint err = cwhizard.whizard_flv_array_string(
            &self._c_whizard,
            &fa.data.as_ints[0],
            len(fa),
            c_str.data.as_chars, strlen)
        # Remove trailing null character from C.
        return c_str.tobytes().decode('UTF-8')[:-1]

    cpdef command(self, str cmd):
        """
        Execute a Sindarin command.

        Examples
        --------
        >>> wz.command("integrate(proc_id)")
        """
        cdef array.array c_cmd = array.array('b', cmd.encode('UTF-8') + b'\0')
        cwhizard.whizard_command(
            &self._c_whizard,
            c_cmd.data.as_chars)

    cpdef (double, double) get_integration_result(self, str proc_id):
        """
        Get the integration results for process with proc_id.

        Returns
        -------
        integral : double
            Integrated cross section or width in pb/fb.
        error : double
            Error estimate on the adaptive Monte Carlo integration.
        """
        cdef array.array c_proc_id = array.array('b', proc_id.encode('UTF-8') + b'\0')
        cdef double integral, error
        cdef bint err = cwhizard.whizard_get_integration_result(
            &self._c_whizard,
            c_proc_id.data.as_chars,
            &integral,
            &error)
        return integral, error

    cpdef WhizardSample new_sample(self, name):
        """
        Allocate a WhizardSample handle which accesses WHIZARD's event stream.

        Name may be a list of proc_id or a string where each proc_id is separated by a comma.

        Returns
        -------
        WhizardSample
        """
        if isinstance(name, list):
            name = ', '.join(name)
        return WhizardSample(self, name)
