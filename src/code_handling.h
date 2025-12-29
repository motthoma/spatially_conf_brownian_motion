

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
} T_DestPaths;

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void CODEHAND_makedirectory(const char *confprfx, const char *intprfx); 

void CODEHAND_copy_main(void);

void CODEHAND_delerrorfiles(void);

void CODEHAND_copycode(void);

void CODEHAND_copy_comp_gen_header(void);
