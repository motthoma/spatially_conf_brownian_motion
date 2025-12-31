/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

typedef struct {
    char *destdir_name; 
    char *fullpath; 
    char *fname_simparams;
    char *fname_intparams;
    char *fname_confparams;

} T_DestPaths;

extern 
T_DestPaths DestPaths;
/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void CODEHAND_makedirectory(const char *confprfx, const char *intprfx); 

void CODEHAND_copy_file_to_dest(char *filename);

void CODEHAND_delerrorfiles(void);

void CODEHAND_copycode(void);

void CODEHAND_copy_comp_gen_header(void);
