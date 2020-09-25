#ifndef INCLUDE_LICOM_ERROR
#define INCLUDE_LICOM_ERROR

//extern void LICOM_Error_mod_mp_LICOM_ErrorSet(int *, char []);
extern void licom_error_mod_mp_exit_licom_(int *, char []);

#define sigExit 0
#define sigAbort -1
#define LICOM_stdout 6

#endif // !INCLUDE_LICOM_ERROR
