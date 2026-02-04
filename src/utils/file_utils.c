#include "file_utils.h"
#include <stdio.h>
#include <string.h>

int save_to_file(const char *filename, const void *data, size_t size) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) { perror("File open error"); return 0; }
    size_t written = fwrite(data, 1, size, fp);
    fclose(fp);
    return written == size;
}

int load_from_file(const char *filename, void *data, size_t size) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0; // 文件不存在
    size_t read = fread(data, 1, size, fp);
    fclose(fp);
    return read == size;
}