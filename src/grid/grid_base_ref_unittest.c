/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "grid_base_ref_replay.h"

static int run_test(const char task_file[]) {
    char filename[1024] = "";

    const char* cp2k_root_dir = getenv("CP2K_ROOT_DIR");
    if (cp2k_root_dir != NULL) {
        assert(strlen(cp2k_root_dir) < 512);
        assert(strcpy(filename, cp2k_root_dir) != NULL);
        if (filename[strlen(filename) - 1] != '/') {
            assert(strcat(filename, "/") != NULL);
        }
    }

    assert(strcat(filename, "src/grid/sample_tasks/") != NULL);
    assert(strcat(filename, task_file) != NULL);

    const double max_diff = grid_collocate_replay(filename, 1);
    if (max_diff > 1e-16) {
        printf("Max diff too high, test failed.\n");
        return 1;
    } else {
        printf("Max diff looks good, test passed.\n\n");
        return 0;
    }
}

int main(){
    int errors = 0;
    errors += run_test("collocate_ortho_density.task");
    errors += run_test("collocate_ortho_tau.task");
    errors += run_test("collocate_general_density.task");
    errors += run_test("collocate_general_tau.task");
    return errors;
}

//EOF
