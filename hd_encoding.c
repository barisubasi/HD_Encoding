#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

/*
Author: Peter Sutor
Date : March 31st, 2022

Description :
	Library for encoding data as 1d arrays in a way that they can be decoded or
	looked up afterwards.This library is designed to handle the broad data types
	of arrays of numbers, text, and images.

	Usage :
	To use, import this file as a library to your project.Note that some global
	variables in this library need to be initialized in order for this library to
	work, as described below in the GLOBAL DECLARATIONS section.
	*/


	//---------------------------------------------------------------------------- -
	//INCLUDE STATEMENTS---------------------------------------------------------- -
	//---------------------------------------------------------------------------- -
#include <stdio.h>
#include <math.h>			//Lines between 23-30 in 1d_encoding.py file
#include <time.h>
//#include "testfile.c"

//---------------------------------------------------------------------------- -
//GLOBAL DECLARATIONS---------------------------------------------------------- -
//---------------------------------------------------------------------------- -
//Lines between 42-101 in 1d_encoding.py file

enum { size = 10000 };

const int native_lenght = 8192;	/*The native length of an array to use.Any
									operation that creates, manipulates, or
									performs operations on arrays will
									assume this length for arrays.Default
									value is 2 ^ 13, or 8, 192 length arrays,
									the closest power to 10,000, the typical
									length of array in Hyperdimensional Computing
									literature.*/

int** intensities;			/* A 2d array of intensity values for RGB
							images, with 256 intensity values, each
							represented as a 1d 1d array of native_length
							length.Each row in this variable is a 1d
							representing the intensity corresponding to
							the row value(row 0 = intensity 0, row 1 =intensity
							1, ... row 255 = intensity 255).
							MUST BE INITIALIZED BY :
							populate_intensities */


int row_permutation[size];     /*A 1d array of integers of native_length length,
								that represented a random permutation of a
								arrays's elements. This particular
								permutation represents row movements in an
								image, symbolically.
								MUST BE INITIALIZED BY :
								populate_row_col_permutations */


int row_permutations[size][size]; /*A 2d array of 1d arrays, where each row is
									a aray of length native_length.The row
									index corresponds to the power of the
									permutation, which is being applied to the
									array row_permutation, symbolically
									representing movements to the right in a
									sequence of data.Row 0 means no permutation
									(0 power), row 1 is a single application of
									row_permutation, ... row i is i permutations
									with the pattern row_permutation(i'th power).
									MUST BE INITIALIZED BY :
									populate_row_col_permutations
									where the rows argument signifies the number
									of rows to create in row_permutations.*/

int col_permutation[size];/*A array of integers of native_length length,
							that represented a random permutation of a
							array's elements. This particular
							permutation represents column movements in an
							image, symbolically.
							MUST BE INITIALIZED BY:
							   populate_row_col_permutations*/

int** col_permutations;      /*A 2d  array of 1d arrays, where each row is
							   a array of length native_length. The row
							   index corresponds to the power of the
							   permutation, which is being applied to the
							   array col_permutation, symbolically
							   representing movements to the right in a
							   sequence of data. Row 0 means no permutation
							   (0 power), row 1 is a single application of
							   row_permutation, ... row i is i permutations
							   with the pattern row_permutation (i'th power).
							   MUST BE INITIALIZED BY:
								  populate_row_col_permutations
							   where the cols argument signifies the number
							   of rows to create in col_permutations.  */

int** locations_basis;/*		A 2d  array of 1d arrays, where each row is
							   a array of length native_length. The row
							   index corresponds to a particular location of
							   sequential data, relative to the first thing
							   in the sequence. This is used to bind multiple
							   sequential elements together into a single
							   array, where the corresponding row in
							   the locations basis is that element in the
							   sequence. The computation looks as follows:
								  X = a_1 XOR l_1 + a_2 XOR l_2 + ...
							   where a_i is the element i in the sequence of
							   data, l_i is the corresponding identifying row
							   of the locations basis, and + is the consensus
							   summation (can use consensus_sum function to
							   compute this).
							   MUST BE INITIALIZED BY:
								  populate_locations
							   where the num argument is the number of row of
							   the locations basis to create.*/

							   //---------------------------------------------------------------------------- -
							   //END GLOBAL DECLARATIONS---------------------------------------------------- -
							   //---------------------------------------------------------------------------- -




							   //---------------------------------------------------------------------------- -
							   //FUNCTION DEFINITIONS--------------------------------------------------------
							   //---------------------------------------------------------------------------- -



							   /*
							   -----------------------------------------------------------------------------
							   Function: new_array
							   Creates a random array of native_length length, or of the length
							   provided by argument length.

							   Parameters:
							   length        The length of the random array to create. Defaults to
												native_length. An optional parameter.

							   Output:
							   Returns a  array of native_length or specified length, containing
							   random 0's or 1's, emulating a binary array.

							   */
							   //Lines between 131-145

int* new_array()
{
	int r[size]; //creates a 1d array 
	//Simply use random function to generate a array of given
	//length, by choosing randomly between 0 and 1.
	for (int i = 0; i < native_lenght; ++i)
	{
		r[i] = rand() % 2;
	}
	return r;
}

/*
-----------------------------------------------------------------------------
Function: new_array_set
Creates an array of random arrays in the form of a 2D  array,
where each row of the array is a array of length native_length, or what
is specified by the user in the length argument.

Parameters:
number        The number of array to create, i.e. the number of rows in the
				 output 2D array.
length        The length of the random array to create. Defaults to
				 native_length. An optional parameter.

Output:
augmented     A list of modified images corresponding to image_list, with
				 intensities divided by 2.

*/
//Lines between 150-167 in 
int** new_array_set(int number)
{
	int** r = malloc(number * sizeof(*r));
	// Uses random function to generate a 2Darray of number
	//rows and length cols, by choosing randomly between 0 and 1.
	for (int i = 0; i < number; ++i)
	{
		for (int j = 0; j < native_lenght; j++)
		{
			r[i][j] = rand() % 6;
		}
	}
}


/*
-----------------------------------------------------------------------------
Function: set_native_length
Alter the native_length global variable, which is used as the default length
of array in libary functions.

Parameters:
length        The new length to set the native_length to.

Output:
No output.
*/
//Lines between 172-184 in 1d_encoding.py file

void set_native_length(int length)
{
	const int native_lenght = length; //Set native_length to new length.

}

/*
-----------------------------------------------------------------------------
Function: swap
Randomizes array elements

Parameters:
a        The first element to swap by its adress
b        The second element to swap by its adress

Output:
No output.

*/
//Not in 1d.encoding.py file
void swap(int* a, int* b)
{
	int temp = *a;
	*a = *b;
	*b = temp;

}
/*
-----------------------------------------------------------------------------
Function: populate_row_col_permutations
Populates the row and col permutations. Namely, the global variables:
   row_permutation
   row_permutations
   col_permutation
   col_permutations
in order to pre-compute permutation arrays for encoding operations, to speed
up the run time. User specifies the max number of rows and the max number of
columns they expect in their data as arguments, and the number of rows in
row_permutations and col_permutations are dictated by this. These effectively
represent the max power of those permutations.

Parameters:
rows          The number of rows to have in row_permutations.
cols          The number of columns to have in col_permutations.

Output:
No output.
*/
//Lines between 188-248 in 1d_encoding.py file
void populate_row_col_permutations(int rows, int cols)
{
	//-------ROW PART---------

	// Set the row_permutations pattern, by generating a random permutation of
	//the array of integers[0, 1, 2, ..., native_length - 1].
	for (int i = 0; i < native_lenght; ++i)
		row_permutation[i] = i + 1;


	//Generate a array of zeros, then set the first row to be a
	//non - permutation, and the second to be a permutation of power 1, or just
	//row_permutation.
	for (int i = 0; i < rows; ++i)
		if (i == 0)
			for (int j = 0; j < native_lenght; ++j)
				row_permutations[i][j] = j + 1;
		else if (i == 1)
			for (int k = native_lenght - 1; k > 0; k--)
			{
				int j = rand() % (k + 1);
				swap(&row_permutation[k], &row_permutation[j]);
				row_permutations[i][k] = row_permutation[k];
			}
		else
			for (int j = 0; j < native_lenght; ++j)
				row_permutations[i][j] = 0;

	//Go through the rest of the rows of row_permutation, i.e., the powers of a
	//permutation, and set each to be the previous row permuted once again by
	//the pattern row_permutation.
	if (rows > 2)
		for (int i = 2; i < rows; i++)
			for (int j = 0; j < native_lenght; j++)
				row_permutations[i][j] = row_permutations[i - 1][row_permutation[j]];


	//---------COL PART-----------

	// Set the col_permutations pattern, by generating a random permutation of
	//the array of integers[0, 1, 2, ..., native_length - 1].

	for (int i = 0; i < native_lenght; ++i)
		col_permutation[i] = i + 1;

	// Generate a array of zeros, then set the first row to be a
	//non - permutation, and the second to be a permutation of power 1, or just
	//col_permutation.
	for (int i = 0; i < cols; ++i)
		if (i == 0)
			for (int j = 0; j < native_lenght; ++j)
				col_permutations[i][j] = j + 1;
		else if (i == 1)
			for (int k = native_lenght - 1; k > 0; k--)
			{
				int j = rand() % (k + 1);
				swap(&col_permutation[k], &col_permutation[j]);
				//Use swap function to shuffle
				col_permutations[i][k] = col_permutation[k];
			}
		else
			for (int j = 0; j < native_lenght; ++j)
				col_permutations[i][j] = 0;

	//Go through the rest of the rows of col_permutation, i.e., the powers of a
	//permutation, and set each to be the previous row permuted once again by
	//the pattern col_permutation.

	if (cols > 2)
		for (int i = 2; i < cols; i++)
			for (int j = 0; j < native_lenght; j++)
				col_permutations[i][j] = col_permutations[i - 1][col_permutation[j]];

}

/*
-----------------------------------------------------------------------------
Function: populate_intensities
Populate the intensities for RGB image data. This is a 2D array of dimensions
256 x native_length, of 0's and 1's, where each row is a array of
length native_length, representing the intensity of the same value as it's
row index. Thus, row 0 = intensity 0, row 1 = intensity 1, ..., row 255 =
intensity 255, each represented symbolically as a random array.

Parameters:
No input parameters.

Output:
No output.

*/
//Lines between 251-270 in 1d_encoding.py file

void populate_intensities()
{
	intensities = new_array_set(256);
	//No need to mark it as global.It's already a global value
	// Call new_array_set to generate the intensities and assign the result to
	//the global variable intensities.
}

/*
-----------------------------------------------------------------------------
Function: populate_locations
Populate the global locations_basis, which is a basis of arrays that
represents each location in a binding operation, so that sequential data can
be decoded later.

Parameters:
num           The number of locations to make in the locations basis, where
				 this value defines the number of rows in locations_basis.

Output:
No output.

*/
//Lines between 273-291 in 1d_encoding.py file
void populate_locations(int num)
{
	//No need to mark it as global.It's already a global value
	locations_basis = new_array_set(num);
	//Call new_array_set to generate the locations basis and assign the result
	//to the global variable locations_basis.
}


/*
-----------------------------------------------------------------------------
Function: XOR
Performs the bitwise XOR operation on arrays a and b, and returns the
result.

Parameters:
a             The first array to use in a XOR b, as a array.
b             The second array to use in a XOR b, as a array.
size_t        The size of the first two parameters.

Output:
Returns the result of a XOR b, the bitwise XOR of a and b, as a numpy array.
*/
//Lines between 293-310 in 1d_encoding.py file
int* XOR(int a[], int b[], int size_t)
{
	int new_array[size];
	//Creates a new array

	//if array a's i'th element matches with the
	//array b's i'th element, new_array array's i'th element
	//is 0. Else 1;
	for (int i = 0; i < size_t; i++)
	{
		if (a[i] == b[i])
		{
			new_array[i] = 0;
		}
		else
		{
			new_array[i] = 1;
		}
	}
	return new_array;
}

/*
# -----------------------------------------------------------------------------
# Function: hamming
# Measures the Hamming Distance between two arrays a and b, which is the sum
# of their XOR, divided by the native_length, in order to normalize the result
# as a ratio of 1's to 0's. The smaller the result, the closer the arrays are
# to each other.
#
# Parameters:
# a             The first array to measure Hamming Distance between it
#                  and another array.
# b             The second array to measure Hamming Distance between it
#                  and another array.
#
# Output:
# Returns the normalized Hamming Distance between a and b, as a value in the
# range [0, 1].
*/
//Lines between 313-333 in 1d_encoding.py file

float hamming(int a[], int b[])
{
	int n = sizeof b / sizeof * b; //n=array b's size
	int* k = XOR(a, b, n);
	//k value holds the value of XOR function for array a and b

	int total = 0;
	for (int i = 0; i < n; i++)
	{
		total += k[i];  //sum of the XOR values
	}
	return total / native_lenght; //divides total to native_lenght and returns.
}


/*
-----------------------------------------------------------------------------
Function: permute
Permutes array a's elements according to a given permutations pattern,
where both are arrays.

Parameters:
a             The array to permute, as a array.
pattern       The pattern to permute a with. This is represented as a numpy
				 array of length equal to a, which is a permutation of the
				 array [0, 1, 2, ..., native_length-1], signifying how to
				 reorder the elements of a.

Output:
Returns the numpy array a, with elements reordered according to pattern.

*/
//Lines between 336-353 in 1d_encoding.py file
int* permute(int a[], int pattern[])
{
	//Creates 1d array.
	int r[size];

	//applies the pattern in pattern array.
	for (int i = 0; i < sizeof(a); i++)
	{
		r[i] = a[pattern[i]];

	}
	return r;
}

/*
-----------------------------------------------------------------------------
Function: permute_power
Permutes array a repeatedly using a permutation pattern, power-many
times, which is the pattern specified by a permutation power.

Parameters:
a             The array to permute, as a array.
pattern       The pattern to permute a with. This is represented as a numpy
				 array of length equal to a, which is a permutation of the
				 array [0, 1, 2, ..., native_length-1], signifying how to
				 reorder the elements of a.
power         The number of times to permute a by the pattern, or the power
				 of the permutation.

Output:
Returns the numpy array a, with elements reordered according to pattern to
the power'th power, done by repeatedly applying the permutation specified by
pattern.
*/
//Lines between 356-382 in 1d_encoding.py file
int* permute_power(int a[], int pattern[], int power)
{
	//# For each power, we need to permute again.
	for (int i = 0; i < power; i++)
	{
		a = permute(a, pattern);
	}
	//Return our permuted array.
	return a;
}
/*
-----------------------------------------------------------------------------
Function: consensus_sum
Computes the consensus sum of a set of arrays. The consensus summation
is defined as a array, where each output bit is the bit value that has
more representation across the same component in the set of arrays. So
if the arrays have more 1's in that component, that is the results, or
vice versa if there's more 0's. If there is a tie, the tie breaker is a
random selection between 0 and 1. This is achieved by adding a random
array to the set.

Parameters:
array      The 2D array that specifies the set of arrays to
				 compute the consensus sum of. The rows of this array are
				 arrays of length native_length, which are to be added.

Output:
Returns a single array as a array, which is the consensus sum of
the row arrays of array.
*/
//Lines between 385-426 in 1d_encoding.py file
int* consensus_sum(int hdarray[][size], int array_size)
{
	int rows = array_size;
	//If the number of rows is even, a tie-breaker is needed, so we add another
	//row to array, which is a random array.
	if (rows % 2 == 0)
		rows += 1;

	/*
	 Get the logical array where 1's are TRUE and 0's are FALSE. Then, sum
	this array along its rows to get a array of counts of 1's in each
	component of the set of arrays to sum. Divide this by rows to get the
	ratio of 1's to 0's in the terms of the sum, in each component.
	*/
	int logical[size][size];
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < array_size; j++)
		{
			if (hdarray[i][j] > 0)
			{
				logical[i][j] = 1;
			}
			else
			{
				logical[i][j] = 0;
			}
		}
	}

	//returns the sum of cols divided by rows
	//if it's bigger than 0.5 it rounds to 1 
	//if its smaller than 0.5 it rounds to 0
	int l[size];
	float total;
	for (int i = 0; i < array_size; i++)
	{
		total = 0;
		for (int j = 0; j < rows; j++)
		{
			total += logical[j][i];
		}

		if (total / rows > 0.5)
			l[i] = 1;
		else
			l[i] = 0;
	}

	return l;

}

/*
-----------------------------------------------------------------------------
Function: encode_image_cpu
Encodes an image I as a array by computing the following on each row of
pixels, r:
   R^r I(r, 0) XOR R^r C^1 I(r, 1) XOR ...
where R is the pattern row_permutation, R^r is the pattern in R permuted r
times, C is the col_permutation, and applying 2 permutation patterns after
another permutes them accordingly. The final permutation is applied to the
hyperdimensional representation of pixel I(r, c), which, if this pixel has
intensity x, is the x index in the global variable intensities. This gets
permuted by the permutation pattern, which represents the row and column
permutation. The XOR of these represents the row of pixels as a array.
Finally, each row is XOR'd together to get the final image encoding as a
array.THIS FUNCTION USES THE CPU EXPLICITLY

MUST INITIALIZE BY FIRST CALLING:
   populate_intensities
   populate_row_col_permutations
   populate_locations

Parameters:
image         A 2D array of intensities, for one channel, of an image. You
				 must call encode_image_cpu 3 times for each RGB channel.

Output:
Returns a array encoding of the input image in one channel.
*/
//Lines between 429-534 in 1d_encoding.py file
int* encode_image_cpu(int image[][size])
{
	//Dimensions of the image.
	int shape[2] = { native_lenght,2 };

	// The output of encoding; contains the buckets of each row of pixels in the
	//image, and the original shape of the image, which we can pre - populate.

	int buckets[size][size];
	int encoding[size];

	for (int i = 0; i < shape[0]; i++)
		for (int j = 0; j < native_lenght; i++)
			buckets[i][j] = 0;

	//Initialize the encoding of the image as a array, as all 0's for now.
	for (int i = 0; i < native_lenght; i++)
		encoding[i] = 0;

	//Counts of the number of bindings, bundles, and permutations done.
	int bindings = 0, bundles = 0, permutations = 0;

	//For each row of pixels, encode the row as a array, then XOR that 
	//with the encoding array.
	for (int row = 0; row < shape[0]; row++)
	{
		//Initialize the bucket that will contain this row's pixels as a set
		//of arrays for each pixel.
		int bucket[size][size];

		for (int i = 0; i < shape[1]; i++)
			for (int j = 0; j < native_lenght; i++)
				bucket[i][j] = 0;

		//For each column of the row, i.e., the pixel on that row, encode it and add it
		for (int col = 0; col < shape[1]; col++)
		{
			//The current pixel is encoded as the array of the intensity
			//of the pixel, permuted by the row permutation of this row,
			//permuted by the column permutation of this column in the image.

			int tempintensities[size];
			int temprow_permutations[size];
			int tempcol_permutations[size];

			for (int i = 0; i < sizeof(bucket); i++)
			{
				for (int j = 0; j < sizeof(intensities); j++)
				{
					int x = image[row, col];
					tempintensities[j] = intensities[x][j];
				}
				for (int j = 0; j < sizeof(row_permutation); j++)
				{
					temprow_permutations[j] = tempintensities[row_permutations[row][j]];
				}
				for (int j = 0; j < sizeof(col_permutation); j++)
				{
					tempcol_permutations[j] = temprow_permutations[col_permutations[col][j]];

				}
			}

			permutations += 2; // Count the permutations done here.

			//Add this pixel to the encoding by XORing it with encoding.
			//Remember: the encoding is each pixel, permuted to match the
			//row / column location, XOR'd together.


			int temp[size];
			for (int i = 0; i < native_lenght; i++)
			{
				temp[i] = bucket[col][i];
			}
			int tmp;
			int n = sizeof encoding / sizeof * encoding;
			tmp = XOR(encoding, temp, n);
			free(encoding);
			int encoding = tmp;

			bindings += 1; //Count the binding done here.

			/*
			The bucket must be XOR'd with the correct location from the
			locations basis, to indicate it's sequential location. The
			final results will eventually be
			   p_1 XOR l_1 + p_2 XOR l_2 + ...
			where p_i is the pixel in the row sequence, and l_i is the
			corresponding location array from the locations basis.
			*/

			int tempbucket[size];
			int templocal_basis[size];
			for (int i = 0; i < sizeof(bucket); i++)
			{
				tempbucket[i] = bucket[col][i];
			}
			for (int i = 0; i < sizeof(locations_basis); i++)
			{
				templocal_basis[i] = locations_basis[col][i];
			}
			for (int i = 0; i < sizeof(bucket); i++)
			{
				bucket[col][i] = XOR(tempbucket, templocal_basis, sizeof(bucket))[i];
			}
			bindings += 1;
		}


		/*
		Collapse the bucket into a single array compression by
		performing consensus sum on the rows of the bucket. This is the
		hyper array representing:
		   p_1 XOR l_1 + p_2 XOR l_2 + ...
		where p_i is the pixel in the row sequence, and l_i is the
		corresponding location array from the locations basis.
		*/
		for (int i = 0; i < sizeof(buckets); i++)
			buckets[row][i] = consensus_sum(bucket, sizeof(bucket))[i];

		bundles += shape[1];  // Count this bundle done here

	}
	int to_return[size][size];
	//Set the encoding in our returned object.
	//buckets_index = 0;
	//shape_index = 1;
	//encoding_index = 2;
	for (int i = 0; i < sizeof(buckets); i++)
		to_return[0][i] = buckets[i];

	for (int i = 0; i < sizeof(shape); i++)
		to_return[1][i] = shape[i];

	for (int i = 0; i < sizeof(encoding); i++)
		to_return[2][i] = encoding[i];

	// Print how many of each operation we had to do.
	printf("Bindings: %d\n", bindings);
	printf("Bundles: %d\n", bundles);
	printf("Permutations: %d\n", permutations);

	//Return our encoding and it's buckets, along with any other info.
	return to_return;
}

/*
-----------------------------------------------------------------------------
Function: decode_image_cpu
Decodes an image encoded by encode_image_cpu. This is done by probing for
each intensity value in each row bucket, to find which pixel goes where. The
result of this will build the full image in one channel, if each pixel is
correctly found. Note, there will be some error in finding the correct pixel,
which is dictated by randomness, the odds increase with the number of rows in
the location basis.THIS FUNCTION USES THE CPU EXPLICITLY

MUST INITIALIZE BY FIRST CALLING:
   populate_intensities
   populate_row_col_permutations
   populate_locations
or by populating the global variables with what you used to encode.

Parameters:
buckets       The buckets for each row generated during encoding.
encoding      The array that represents the encoding of the image.
shape         The dimensions of the image, as a tuple of rows and columns.

Output:
Returns a 2D array with values between 0 and 255, signifying the intensities
of each pixel, in one channel.
*/
//Lines between 537-638 in 1d_encoding.py file

int** decode_image_cpu(int buckets[][size], int encoding[], int shape[])
{
	//What will contain our image, eventually, built pixel by pixel.
	int* image[size][size];
	for (int i = 0; i < shape[0]; i)
		for (int j = 0; j < shape[1]; j)
			image[i][j] = 0;

	//Counts of each operation we had to do to encode the image.
	int bindings = 0, permutations = 0, hammings = 0;

	//For each row of pixels, take bucket representing that rowand decode it.
	for (int row = 0; row < shape[0]; row++)
	{
		int bucket[size];

		for (int i = 0; i < sizeof(buckets) / sizeof buckets[0]; i++)
		{
			bucket[i] = buckets[row][i];
		}

		//For each column, predict the pixel that must have been at this row/
		//col location in the original image, from the encoding.
		for (int col = 0; col < shape[0]; col++)
		{
			//Our final permutation for this location, using the row and col
			//permutations that are cached.
			int permuted[size];
			int temprow[size];
			int tempcol[size];

			for (int i = 0; i < sizeof(row_permutations); i++)
				temprow[i] = row_permutations[row, i];
			for (int i = 0; i < sizeof(col_permutations); i++)
				tempcol[i] = col_permutations[row, i];

			for (int i = 0; i < sizeof temprow; i++)
				permuted[i] = temprow[tempcol[i]];

			permutations += 1; //# Count this permutation.

			//Need to keep track of the closest distance seen so far (closest
			//to 0, the best matching intensity, and the corresponding
			//intensity array.

			float best_distance = 1.0;
			int best_match = -1;
			int best_encoding[size];
			for (int i = 0; i < native_lenght; i++)
				best_encoding[i] = 0;

			int tempencode = best_encoding;

			//For each possible intensity, probe the bucket to see how good of
			//a match it is in terms of Hamming Distance.If we find a better
			//one, update our best matches.
			for (int intensity = 0; intensity < sizeof(intensities); intensity++)
			{
				//Our current hypethesis (the currenty intensity array).
				int tempinten[size];
				int hypothesis[size];

				for (int i = 0; i < sizeof(intensities); i++)
					tempinten[i] = intensities[intensity][i];

				for (int i = 0; i < sizeof tempinten; i++)
					hypothesis[i] = tempinten[permuted[i]];

				permutations += 1;//# Count this permutation.

				/*
				Compute the Hamming Distance between our location basis for
				this spot in the sequence of row pixels, and the XOR of our
				hypothesis and bucket. The XOR of the hypothesis and bucket
				collapses the superposed row pixel arrays from encoding into
				either pure random noise, or a slightly meaningful deviation
				for the matching location. By comparing the distance between
				the location and this result, our best match becomes our
				prediction for the pixel intensity.
				*/


				int templocal[size];
				for (int i = 0; i < sizeof locations_basis; i++)
					templocal[i] = locations_basis[col][i];

				float distance = hamming(XOR(hypothesis, bucket, sizeof hypothesis), templocal);
				hammings += 1; //Count this Hamming distance.
				bindings += 2; //Count the bindings we had do to here.

				//If this is a better match (smaller distance), update our best
				//match information.
				if (distance < best_distance)
				{
					best_distance = distance;
					best_match = intensity;
					tempencode = hypothesis;
				}
			}

			//Now that we have our prediction for the pixel intensity at this
			//location, set it in the image.
			*image[row, col] = best_match;

			//# XOR our pixel prediction array with our encoding to
			//"remove" it from the encoding.
			encoding = XOR(encoding, best_encoding, (sizeof encoding / sizeof * encoding));
			bindings += 1;   //Count the binding we had to do here.
		}
	}
	//Print how many operations we had to do.
	printf("Bindings: %d\n", bindings);
	printf("Permutations: %d\n", permutations);
	printf("Hamming Distances: \n", hammings);
	int val = 0;
	for (int i = 0; i < sizeof encoding / sizeof encoding[0]; i++)
	{
		if (encoding[i] > 0)
			val = 1;
		break;
	}
	if (val == 0)
		printf("Pixels have been lost...");

	//Return our predicted image.
	return image;
}
//--------------------------TEST FILE------------------------------
#include <time.h>
//#include <stb_image.h>


#define FALSE (1==0)
#define TRUE  (1==1)


char data[] = "images";
char device[] = "cpu";
char path[] = "/C:User";


int run_cpu_images(int file_list)
{
	/*
	index            element
   --------		   --------------
	 0				  images
	 1				 encodings
	 2               location_basis
	 3				 row_permutation
	 4				 row_permutations
	 5				 col_permutation
	 6				 col_permutations
	 7				 intensities
	*/

	int result_map[8];

	result_map[0] = file_list;
	result_map[1];
	result_map[2] = locations_basis;
	result_map[3] = row_permutation;
	result_map[4] = row_permutations;
	result_map[5] = col_permutation;
	result_map[6] = col_permutations;
	result_map[7] = intensities;


	for (int i = 0; i < sizeof(file_list); i++)
	{
		int temp_r[size][3];
		int temp_g[size][3];
		int temp_b[size][3];
		int templist = file_list;
		clock_t t = clock();

		int tempr_channel[size][size];

		int tempg_channel[size][size];

		int tempb_channel[size][size];

		tempr_channel[size][size] = encode_image_cpu(temp_r);

		tempg_channel[size][size] = encode_image_cpu(temp_g);

		tempb_channel[size][size] = encode_image_cpu(temp_b);

		t = clock() - t;
		double time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("Elapsed Time Encoding: %f\n", time_taken);

		//result_map[1] = {1,templist[i] , tempr_channel,tempg_channel,tempb_channel };

		t = clock();

		int r_decoded = decode_image_cpu(tempr_channel[0], tempr_channel[1], tempr_channel[2]);

		int g_decoded = decode_image_cpu(tempg_channel[0], tempg_channel[1], tempg_channel[2]);

		int b_decoded = decode_image_cpu(tempb_channel[0], tempb_channel[1], tempb_channel[2]);
	
		t = clock() - t;
		time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("Elapsed Time Encoding: %f\n", time_taken);
	
	}
	return result_map;


}


int run_gpu_images(file_list)
{
	return TRUE;
}

int run_cpu_text(file_list)
{
	return TRUE;
}

int run_gpu_text(file_list)
{
	return TRUE;
}

int run_cpu_numbers(file_list)
{
	return TRUE;
}
int run_gpu_numbers(file_list)
{
	return TRUE;
}
int* find_max_dimensions(file_list)
{
	int rows = 0;
	int cols = 0;
	int temp[2];
	for (int i = 0; i < sizeof(file_list); i++)
	{
		if (rows < sizeof(temp[0]))
		{
			rows = sizeof(temp[0]);
		}

		if (cols < sizeof(temp[1]))
		{
			cols = sizeof(temp[1]);
		}
	}
	int res[2] = { rows,cols };

	return res;
}
	
//Checks for the path if it is exist.
int canCreateFile(char* path)
{
	FILE* file = fopen(path, "w");
	if (file)
	{
		fclose(file);
		return 1;
	}
	return 0;
}


int main()
{ 
	printf("%s\n", path);
	printf("%s\n", data);
	printf("%s\n", device);


	if (canCreateFile(path) == 1)
	{
		//151 159

	}
	char files = "skip";

	if (data == "images")
	{
		int max_rows = find_max_dimensions(files)[0];
		int max_cols = find_max_dimensions(files)[1];

		populate_intensities();
		populate_row_col_permutations(max_rows, max_cols);
		populate_locations(max_cols);

		if (device == "gpu")
		{
			int encoded= run_gpu_images(files);
		}
		else
		{
			int encoded = run_cpu_images(files);
		}
	}
	else if (data == "text")
	{
		if (device == "gpu")
		{
			int encoded = run_gpu_images(files);
		}
		else
		{
			int encoded= run_cpu_images(files);
		}

	}
	else if (data == "numbers")
	{
		if (device == "gpu")
		{
			int encoded= run_gpu_images(files);
		}
		else
		{
			int encoded= run_cpu_images(files);
		}
	}
	else
	{
		int max_rows = find_max_dimensions(files)[0];
		int max_cols = find_max_dimensions(files)[1];

		populate_intensities();
		populate_row_col_permutations(max_rows, max_cols);
		populate_locations(max_cols);

		if (device == "gpu")
		{
			int encoded = run_gpu_images(files);
		}
		else
		{
			int encoded= run_cpu_images(files);
		}

	}

	return 0;
}
