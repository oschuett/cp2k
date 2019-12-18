/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "grid_base_ref_replay.h"

int main(){
    char filename[1024] = "";

    const char* cp2k_root_dir = getenv("CP2K_ROOT_DIR");
    if (cp2k_root_dir != NULL) {
        assert(strlen(cp2k_root_dir) < 512);
        assert(strcpy(filename, cp2k_root_dir) != NULL);
        if (filename[strlen(filename) - 1] != '/') {
            assert(strcat(filename, "/") != NULL);
        }
    }

    const char rel_path[] = "src/grid/sample_tasks/collocate_ortho_density.task";
    assert(strcat(filename, rel_path) != NULL);

    const double max_diff = grid_collocate_replay(filename, 1);
    if (max_diff > 1e-16) {
        printf("Max diff too high, test failed.\n");
        return 1;
    } else {
        printf("Max diff looks good, test passed.\n");
        return 0;
    }
}

//EOF
