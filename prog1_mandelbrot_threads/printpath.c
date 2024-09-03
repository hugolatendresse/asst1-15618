#include <stdio.h>
#include <unistd.h>
#include <limits.h>

int main() {
    char cwd[PATH_MAX];

    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        printf("new change, prog1_mandelbrot_threads dir has this working directory! %s\n", cwd);
    } else {
        perror("getcwd() error");
    }

    return 0;
}
