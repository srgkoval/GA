#ifndef SUPPORT_H
#define SUPPORT_H

#include <string>

void wait_and_exit();
void wait_and_exit(std::string error_message);

void f_read_table(double *x, double *y, int *n, const std::string &filename);
void f_write_table(double *x, double *y, int n, const std::string &filename);

#endif