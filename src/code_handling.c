#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "sim_config.h"
#include "comp_gen_header.h"
#include "code_handling.h"

T_DestPaths DestPaths;

/* ------------------------------------------------------------ */
/* Helper: robust file copy                                     */
/* ------------------------------------------------------------ */
static int copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    if (!in) {
        fprintf(stderr, "copy_file: failed to open '%s': %s\n",
                src, strerror(errno));
        return -1;
    }

    FILE *out = fopen(dst, "wb");
    if (!out) {
        fprintf(stderr, "copy_file: failed to create '%s': %s\n",
                dst, strerror(errno));
        fclose(in);
        return -1;
    }

    char buf[8192];
    size_t nread;
    while ((nread = fread(buf, 1, sizeof buf, in)) > 0) {
        size_t nwritten = fwrite(buf, 1, nread, out);
        if (nwritten != nread) {
            fprintf(stderr, "copy_file: write error: %s\n", strerror(errno));
            fclose(in);
            fclose(out);
            return -1;
            }
    }
 
    fclose(in);
    fclose(out);
    return 0;
}

/* ------------------------------------------------------------ */
/* Create run-subdirectory name                                 */
/* ------------------------------------------------------------ */
void CODEHAND_makedirectory(const char *confprfx, const char *intprfx)
{
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);

    /* Allocate run directory name */
    DestPaths.destdir_name = malloc(256);
    if (!DestPaths.destdir_name) {
        fprintf(stderr, "malloc failed\n");
    }

#ifdef LJMIN
    snprintf(DestPaths.destdir_name, 256,
             "%s_%s_%d_%02d_%02d_%02dh_%02dmin_%02dsec_R_%.3lf_LJMIN_%.3lf_EPS_%.2lf_F_%.2lf_setn_%d",
             confprfx, intprfx,
             tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
             tm->tm_hour, tm->tm_min, tm->tm_sec,
             R_CONF, LJMIN, EPS_L, SimParams.F, SimParams.parts_per_set);

#else
    snprintf(DestPaths.destdir_name, 256,
             "%s_%s_%d_%02d_%02d_%02dh_%02dmin_%02dsec_R_%.3lf_F_%.2lf_setn_%d",
             confprfx, intprfx,
             tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
             tm->tm_hour, tm->tm_min, tm->tm_sec,
             R_CONF, SimParams.F, SimParams.parts_per_set);
#endif

    /* Ensure "runs" directory exists */
    mkdir(RUNS_DIR, 0755);

    DestPaths.fullpath = malloc(512);
    snprintf(DestPaths.fullpath, 512, "%s/%s", RUNS_DIR, DestPaths.destdir_name);

    if (mkdir(DestPaths.fullpath, 0700) != 0) {
        fprintf(stderr, "mkdir fullpath '%s' failed: %s\n",
                DestPaths.fullpath, strerror(errno));
        free(DestPaths.destdir_name);
    }

    /*if (chdir(DestPaths.fullpath) != 0) {
        fprintf(stderr, "chdir '%s' failed: %s\n",
                DestPaths.fullpath, strerror(errno));
        free(DestPaths.destdir_name);
    }*/
}

/* ------------------------------------------------------------ */
/* Copy needed files into run directory                         */
/* ------------------------------------------------------------ */
void CODEHAND_copy_file_to_dest(char *filename)
{  
    char *dest_path = malloc(800);
    char *src_path = malloc(800);
    snprintf(dest_path,
             800,
             "%s/%s",
             DestPaths.fullpath, filename);
    
    snprintf(src_path,
             800,
             "./%s",
             filename);

    copy_file(src_path, dest_path);
}

void CODEHAND_copycode(void)
{
    CODEHAND_copy_file_to_dest("code_handling.c");
    CODEHAND_copy_file_to_dest("code_handling.h");
}
