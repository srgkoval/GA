#include "support.h"
#include <iostream>
#include <fstream>
#include <conio.h>

void wait()
{
	std::cout << "\npress Enter to continue...";
	std::cin.get();
}

void wait_and_exit()
{
	wait();
	exit(EXIT_FAILURE);
}

void wait_and_exit(std::string error_message)
{
	std::cout << error_message << "\n";
	wait_and_exit();
}


double round(double x)
{
    return x < 0.0 ? ceil(x - 0.5) : floor(x + 0.5);
}

// file import/export =============================================================================

void f_read_table(double *x, double *y, int *n, const std::string &filename)
{
	std::ifstream file_input(filename.c_str());
		
	if (!file_input) 
	{
		std::cout << "f_read_table: error reading " << filename.c_str() << "\n";
		wait_and_exit();
	}
		
	int input_length = 0;
	for( ; !file_input.eof(); input_length++)
	{
		file_input >> x[input_length] >> y[input_length];
		//std::cout << x[input_length] << "\t" << y[input_length] << "\n";
		//std::cin.get();
	}


	*n = input_length;

	file_input.close();
}


void f_write_table(double *x, double *y, int n, const std::string &filename)
{
	std::ofstream file_output(filename.c_str());
		
	if (!file_output) 
	{
		std::cout << "f_write_table: error reading " << filename.c_str() << "\n";
		wait_and_exit();
	}
		
	for(int i = 0; i < n; i++)
	{
		file_output << x[i] << "\t" << y[i];
		if (i != n - 1)
			file_output << "\n";
	}

	file_output.close();
}