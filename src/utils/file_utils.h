#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <stddef.h>


int save_to_file(const char *filename, const void *data, size_t size);
int load_from_file(const char *filename, void *data, size_t size);


#endif /* FILE_UTILS_H */
