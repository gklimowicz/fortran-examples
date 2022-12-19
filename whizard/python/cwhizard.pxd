cdef extern from "whizard.h":
    ctypedef void* whizard_t
    ctypedef void* sample_handle_t
    void whizard_create (whizard_t wh)
    void whizard_option (whizard_t wh, const char* key, const char* value)
    void whizard_init (whizard_t wh)
    void whizard_final (whizard_t wh)

    void whizard_set_double (whizard_t wh, const char* var, const double value)
    void whizard_set_int (whizard_t wh, const char* var, const int value)
    void whizard_set_bool (whizard_t wh, const char* var, const bint value)
    void whizard_set_char (whizard_t wh, const char* var, const char* value)
    int whizard_get_double (whizard_t wh, const char* var, double* value)
    int whizard_get_int (whizard_t wh, const char* var, int* value)
    int whizard_get_bool (whizard_t wh, const char* var, int* value)
    int whizard_get_char (whizard_t wh, const char* var, char* value, const int strlen)
    int whizard_get_char_len (whizard_t wh, const char* var)

    int whizard_flv_string (whizard_t wh, const int pdg, char* fstr, const int strlen)
    int whizard_flv_string_len (whizard_t wh, const int pdg)
    int whizard_flv_array_string (whizard_t wh, const int* pdg, const int nf, char* fstr, const int strlen)
    int whizard_flv_array_string_len (whizard_t wh, const int* pdg, const int nf)

    void whizard_command (whizard_t wh, const char* cmd)
    int whizard_get_integration_result (whizard_t wh, const char* proc_id, double* integral, double* error)

    void whizard_new_sample (whizard_t wh, const char* name, void* sample)
    void whizard_sample_open (sample_handle_t sample, int* it_begin, int* it_end)
    void whizard_sample_next_event (sample_handle_t sample)
    void whizard_sample_close (sample_handle_t sample)

    void whizard_sample_get_event_index (sample_handle_t sample, int* idx)
    void whizard_sample_get_process_index (sample_handle_t sample, int* i_proc)
    void whizard_sample_get_process_id (sample_handle_t sample, char* proc_id, const int strlen)
    int whizard_sample_get_process_id_len (sample_handle_t sample)
    void whizard_sample_get_fac_scale (sample_handle_t sample, double* f_scale)
    void whizard_sample_get_alpha_s (sample_handle_t sample, double* alpha_s)
    void whizard_sample_get_weight (sample_handle_t sample, double* weight)
    void whizard_sample_get_sqme (sample_handle_t sample, double* sqme)
