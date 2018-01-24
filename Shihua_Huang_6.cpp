// Processing tick-data using the Fast Fourier Transform (FFT)
// 
#define _CRT_SECURE_NO_WARNINGS
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmat.h"
#include "newmatap.h"  // Need this for the FFT algorithm

using namespace std;

// function that is used for sorting a 2-col array based on the value of the entries in the 2nd column.
// taken from http://stackoverflow.com/questions/3041897/sorting-a-2-dimensional-array-on-multiple-columns
bool compareTwoRows2(double* rowA, double* rowB) {
	return ((rowA[1]>rowB[1]) || ((rowA[1] == rowB[1]) && (rowA[0]>rowB[0])));
}

class Filtering_Instance
{
	// Private
	int no_of_terms, no_of_data_points;

	// Private member function that computes the mean of a data array
	double compute_mean(ColumnVector data)
	{
		// write the code to compute the mean of "data"
		double sum;
		int n;
		n = data.Nrows();
		sum = 0;
		for (int i = 1; i <= n; i++)
		{
			sum = sum + data(i);
		}
		double mean = sum / n;
		return mean;
	}

	// Private member function that computes the magnitude of (an array of) complex #s
	void compute_magnitude(ColumnVector &magnitude, ColumnVector real_part, ColumnVector imag_part)
	{
		// write the code to compute sqrt(real_part(i)^2 + imag_part(i)^2)
		for (int i = 1; i <= real_part.Nrows(); i++)
		{
			magnitude(i) = sqrt(real_part(i)*real_part(i) + abs(imag_part(i)*imag_part(i)));
			//double A = 2 / magnitude.Nrows()*real_part(i + 1);
			//double B = 2 / magnitude.Nrows()*imag_part(i + 1);
			//magnitude(i) = sqrt(A * A + B * B);
		}
	}

	// Private member function that reads the data from the input file 
	// and stores it in "data" 
	void get_data(char* file_name, ColumnVector &data)
	{
		// write code that reads the ticker-data from the input file
		// and stores it in "data"
		ifstream input_filename(file_name);
		if (input_filename.is_open())
		{
			for (int i = 1; i <= data.Nrows(); i++) 
			{
				double value_just_read;
				input_filename >> value_just_read;
				data(i) = value_just_read;
			}
		}
		else
		{
			cout << "Input file missing" << endl;
			exit(0);
		}
	}

	// private member function that writes the data file into a file
	void write_data(char* file_name, ColumnVector &filtered_data)
	{
		
		// write code that writes "data" to file_name.
		ofstream out_filename(file_name);
		if (out_filename.is_open())
		{
			for (int i = 1; i <= filtered_data.Nrows(); i++)
			{
				double value_output = filtered_data(i);
				out_filename << value_output<< endl;
			}
			out_filename.close();
		}
	}

	// private member function that filters data using the FFT 
	// The filtered data is computed using the top "no_of_terms"-many 
	// magnitude-components of the orginal data
	void filter_the_data(ColumnVector &data, ColumnVector &filtered_data, int no_of_terms)
	{
		ColumnVector fft_real_part(data.Nrows()), fft_imag_part(data.Nrows());
		ColumnVector mean_adjusted_data(data.Nrows()), magnitude(data.Nrows());

		double mean = compute_mean(data);
		for (int i = 1; i <= data.Nrows(); i++)
			mean_adjusted_data(i) = data(i) - mean;

		RealFFT(mean_adjusted_data, fft_real_part, fft_imag_part);
		compute_magnitude(magnitude, fft_real_part, fft_imag_part);

		// creating a two dimensional array: first col is the index; second col is the 
		// magnitude.  The plan is to have this 2-D array sorted using the 2nd col (ie.
		// sorted based on magnitude). Then we pick the top "no_of_terms"-many of these
		// components to reconstitute/reconstruct the signal back. 

		double** two_dimensional_array = new double*[fft_imag_part.Nrows()];
		for (int i = 0; i < fft_imag_part.Nrows(); i++)
			two_dimensional_array[i] = new double[2];

		for (int i = 0; i < fft_imag_part.Nrows(); i++)
		{
			two_dimensional_array[i][0] = i;
			two_dimensional_array[i][1] = magnitude(i + 1);
		}
		std::sort(two_dimensional_array, two_dimensional_array + fft_imag_part.Nrows(), &compareTwoRows2);

		// if do_we_pick_this(i) == 1, then we keep that component for reconstruction
		// of the filtered signal.  The rest of the array-names should be self-explanatory
		ColumnVector do_we_pick_this(fft_imag_part.Nrows());
		ColumnVector filtered_fft_real_part(fft_imag_part.Nrows());
		ColumnVector filtered_fft_imag_part(fft_imag_part.Nrows());

		// write the code for picking the top "no_of_terms" many	magnitudes
		// put the real-part in "filtered_fft_real_part" and imaginary-part in 
		// "filtered_fft_imag_part" -- and reconstruct the filtered signal as 
		// shown below. 
		for (int i =1; i <= fft_imag_part.Nrows(); i++)
			do_we_pick_this(i) = 0;//initialize the value of all do_we_pick_this for n/2+1 points into 0
		for (int i = 0; i < no_of_terms; i++)
			do_we_pick_this(two_dimensional_array[i][0] + 1) = 1;//pick those top no_of_terms which index is also sorted along with magnitude;
		                                                         //therefore, we use two_dimensional_array[i][0] + 1 to record the top no_of_terms
		for (int i = 1; i <= fft_imag_part.Nrows(); i++)//for n/2+1 points
		{
			if (do_we_pick_this(i) == 1)
			{
				filtered_fft_real_part(i) = fft_real_part(i);//keep the original value for the new filtered fft real part
				filtered_fft_imag_part(i) = fft_imag_part(i);//keep the original value for the new filtered fft imaginary part
			}
			else
			{
				filtered_fft_real_part(i) = 0;//zero out remaininig entries of the new filtered fft real part
				filtered_fft_imag_part(i) = 0;//zero out remaininig entries of the new filtered fft imaginary part
			}
		}
		// reconstructed signal using just the "no_of_terms"-many top-magnitude 
		// components.
		RealFFTI(filtered_fft_real_part, filtered_fft_imag_part, filtered_data);

		// write code to add the mean-back to the reconstructed, filtered-signal
		for (int i = 1; i <= filtered_data.Nrows(); i++)
			filtered_data(i) = filtered_data(i) + mean;
	}

public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void run_the_filter(int argc, char * const argv[])
	{
		sscanf(argv[1], "%d", &no_of_terms);
		sscanf(argv[2], "%d", &no_of_data_points);

		std::cout << "Input File Name: " << argv[3] << std::endl;
		std::cout << "Number of data points in the input file = " << no_of_data_points << endl;
		std::cout << "Number of dominant terms in the FFT = " << no_of_terms << endl;
		std::cout << "Output File Name: " << argv[4] << std::endl;

		ColumnVector data(no_of_data_points), filtered_data(no_of_data_points);

		// get ticker data
		get_data(argv[3], data);

		// filter the ticker data
		filter_the_data(data, filtered_data, no_of_terms);

		// write the filtered data
		write_data(argv[4], filtered_data);
	}
};


int main (int argc, char* argv[])
{
	Filtering_Instance x;
	x.run_the_filter(argc, argv);
}