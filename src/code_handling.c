#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "comp_gen_header.h"
#include "code_handling.h"

#define RUNS_DIR "../runs"

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
char* CODEHAND_makedirectory(const char *confprfx, const char *intprfx)
{
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);

    /* Allocate run directory name */
    char *dirname = malloc(256);
    if (!dirname) {
        fprintf(stderr, "malloc failed\n");
        return NULL;
    }

#ifdef LJMIN
    snprintf(dirname, 256,
             "%s_%s_%d_%02d_%02d_%02dh_%02dmin_%02dsec_R_%.3lf_LJMIN_%.3lf_EPS_%.2lf_F_%.2lf_setn_%d",
             confprfx, intprfx,
             tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
             tm->tm_hour, tm->tm_min, tm->tm_sec,
             R_CONF, LJMIN, EPS_L, SimParams.F, SimParams.parts_per_set);
#else
    snprintf(dirname, 256,
             "%s_%s_%d_%02d_%02d_%02dh_%02dmin_%02dsec_R_%.3lf_F_%.2lf_setn_%d",
             confprfx, intprfx,
             tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
             tm->tm_hour, tm->tm_min, tm->tm_sec,
             R_CONF, SimParams.F, SimParams.parts_per_set);
#endif

    /* Ensure "runs" directory exists */
    mkdir(RUNS_DIR, 0755);

    /* Build full path: runs/<dirname> */
    char fullpath[512];
    snprintf(fullpath, sizeof fullpath, "%s/%s", RUNS_DIR, dirname);

    if (mkdir(fullpath, 0700) != 0) {
        fprintf(stderr, "mkdir '%s' failed: %s\n",
                fullpath, strerror(errno));
        free(dirname);
        return NULL;
    }

    if (chdir(fullpath) != 0) {
        fprintf(stderr, "chdir '%s' failed: %s\n",
                fullpath, strerror(errno));
        free(dirname);
        return NULL;
    }

    return dirname;
}

/* ------------------------------------------------------------ */
/* Copy needed files into run directory                         */
/* ------------------------------------------------------------ */
void CODEHAND_copy_main(void)
{
    copy_file("../main_brownconf.c", "main_brownconf.c");
    copy_file("../makefile",          "makefile");
}

void CODEHAND_copycode(void)
{
    copy_file("../code_handling.c", "code_handling.c");
    copy_file("../code_handling.h", "code_handling.h");
}

void CODEHAND_copy_comp_gen_header(void)
{
    copy_file("../comp_gen_header.h", "comp_gen_header.h");
}
