#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Author: Peter Sutor
# Date: March 31st, 2022
#
# Description:
# Library for encoding data as HD vectors in a way that they can be decoded or
# looked up afterwards. This library is designed to handle the broad data types
# of arrays of numbers, text, and images.
#
# Usage:
# To use, import this file as a library to your project. Note that some global
# variables in this library need to be initialized in order for this library to
# work, as described below in the GLOBAL DECLARATIONS section.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# IMPORT STATEMENTS -----------------------------------------------------------
# -----------------------------------------------------------------------------

# Full imports.
# import os
# import argparse
import numpy as np
# import cupy as cp

# Imported files from packages.
from numpy import random

# -----------------------------------------------------------------------------
# END IMPORT STATEMENTS -------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# GLOBAL DECLARATIONS ---------------------------------------------------------
# -----------------------------------------------------------------------------

# Global variables.
native_length = 2 ** 13     # The native length of a hypervector to use. Any
                            #    operation that creates, manipulates, or
                            #    performs operations on hypervectors will
                            #    assume this length for hypervectors. Default
                            #    value is 2^13, or 8,192 length hypervectors,
                            #    the closest power to 10,000, the typical
                            #    length of vector in Hyperdimensional Computing
                            #    literature.
intensities = None          # A 2d numpy array of intensity values for RGB
                            #    images, with 256 intensity values, each
                            #    represented as a HD vector of native_length
                            #    length. Each row in this variable is a HD
                            #    representing the intensity corresponding to
                            #    the row value (row 0 = intensity 0, row 1 =
                            #    intensity 1, ... row 255 = intensity 255).
                            #    MUST BE INITIALIZED BY:
                            #       populate_intensities
row_permutation = None      # A vector of integers of native_length length,
                            #    that represented a random permutation of a
                            #    hypervector's elements. This particular
                            #    permutation represents row movements in an
                            #    image, symbolically.
                            #    MUST BE INITIALIZED BY:
                            #       populate_row_col_permutations
row_permutations = None     # A 2d numpy array of HD vectors, where each row is
                            #    a hypervector of length native_length. The row
                            #    index corresponds to the power of the
                            #    permutation, which is being applied to the
                            #    hypervector row_permutation, symbolically
                            #    representing movements to the right in a
                            #    sequence of data. Row 0 means no permutation
                            #    (0 power), row 1 is a single application of
                            #    row_permutation, ... row i is i permutations
                            #    with the pattern row_permutation (i'th power).
                            #    MUST BE INITIALIZED BY:
                            #       populate_row_col_permutations
                            #    where the rows argument signifies the number
                            #    of rows to create in row_permutations.
col_permutation = None      # A vector of integers of native_length length,
                            #    that represented a random permutation of a
                            #    hypervector's elements. This particular
                            #    permutation represents column movements in an
                            #    image, symbolically.
                            #    MUST BE INITIALIZED BY:
                            #       populate_row_col_permutations
col_permutations = None     # A 2d numpy array of HD vectors, where each row is
                            #    a hypervector of length native_length. The row
                            #    index corresponds to the power of the
                            #    permutation, which is being applied to the
                            #    hypervector col_permutation, symbolically
                            #    representing movements to the right in a
                            #    sequence of data. Row 0 means no permutation
                            #    (0 power), row 1 is a single application of
                            #    row_permutation, ... row i is i permutations
                            #    with the pattern row_permutation (i'th power).
                            #    MUST BE INITIALIZED BY:
                            #       populate_row_col_permutations
                            #    where the cols argument signifies the number
                            #    of rows to create in col_permutations.
locations_basis = None      # A 2d numpy array of HD vectors, where each row is
                            #    a hypervector of length native_length. The row
                            #    index corresponds to a particular location of
                            #    sequential data, relative to the first thing
                            #    in the sequence. This is used to bind multiple
                            #    sequential elements together into a single
                            #    hypervector, where the corresponding row in
                            #    the locations basis is that element in the
                            #    sequence. The computation looks as follows:
                            #       X = a_1 XOR l_1 + a_2 XOR l_2 + ...
                            #    where a_i is the element i in the sequence of
                            #    data, l_i is the corresponding identifying row
                            #    of the locations basis, and + is the consensus
                            #    summation (can use consensus_sum function to
                            #    compute this).
                            #    MUST BE INITIALIZED BY:
                            #       populate_locations
                            #    where the num argument is the number of row of
                            #    the locations basis to create.

# -----------------------------------------------------------------------------
# END GLOBAL DECLARATIONS -----------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS --------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function: new_vector
# Creates a random hypervector of native_length length, or of the length
# provided by argument length.
#
# Parameters:
# length        The length of the random hypervector to create. Defaults to
#                  native_length. An optional parameter.
#
# Output:
# Returns a numpy vector of native_length or specified length, containing
# random 0's or 1's, emulating a binary hypervector.
def new_vector(length=native_length):
    # Simply use numpy's random function to generate a numpy array of given
    # length, by choosing randomly between 0 and 1.
    return np.random.randint(2, size=length)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: new_vector_set
# Creates an array of random hypervectors in the form of a 2D numpy array,
# where each row of the array is a hypervector of length native_length, or what
# is specified by the user in the length argument.
#
# Parameters:
# number        The number of vector to create, i.e. the number of rows in the
#                  output 2D numpy array.
# length        The length of the random hypervector to create. Defaults to
#                  native_length. An optional parameter.
#
# Output:
# augmented     A list of modified images corresponding to image_list, with
#                  intensities divided by 2.
def new_vector_set(number, length=native_length):
    # Simple use numpy's random function to generate a 2D matrix of number
    # rows and length cols, by choosing randomly between 0 and 1.
    return np.random.randint(2, size=(number, length))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: set_native_length
# Alter the native_length global variable, which is used as the default length
# of hypervector in libary functions.
#
# Parameters:
# length        The new length to set the native_length to.
#
# Output:
# No output.
def set_native_length(length):
    global native_length        # Mark as global.

    native_length = length      # Set native_length to new length.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: populate_row_col_permutations
# Populates the row and col permutations. Namely, the global variables:
#    row_permutation
#    row_permutations
#    col_permutation
#    col_permutations
# in order to pre-compute permutation vectors for encoding operations, to speed
# up the run time. User specifies the max number of rows and the max number of
# columns they expect in their data as arguments, and the number of rows in
# row_permutations and col_permutations are dictated by this. These effectively
# represent the max power of those permutations.
#
# Parameters:
# rows          The number of rows to have in row_permutations.
# cols          The number of columns to have in col_permutations.
#
# Output:
# No output.
def populate_row_col_permutations(rows, cols):
    # Mark these as global variables to alter.
    global row_permutation, row_permutations, col_permutation, col_permutations

    # Set the row_permutations pattern, by generating a random permutation of
    # the vector of integers [0, 1, 2, ..., native_length-1].
    row_permutation = np.arange(native_length)

    # Generate a matrix of zeros, then set the first row to be a
    # non-permutation, and the second to be a permutation of power 1, or just
    # row_permutation.
    row_permutations = np.zeros((rows, native_length), dtype='int32')
    row_permutations[0, :] = np.arange(native_length)
    random.shuffle(row_permutation)
    row_permutations[1, :] = row_permutation

    # Go through the rest of the rows of row_permutation, i.e., the powers of a
    # permutation, and set each to be the previous row permuted once again by
    # the pattern row_permutation.
    if rows > 2:
        for i in range(2, rows):
            row_permutations[i, :] = (row_permutations[i - 1, :])[row_permutation]

    # Set the col_permutations pattern, by generating a random permutation of
    # the vector of integers [0, 1, 2, ..., native_length-1].
    col_permutation = np.arange(native_length)

    # Generate a matrix of zeros, then set the first row to be a
    # non-permutation, and the second to be a permutation of power 1, or just
    # col_permutation.
    col_permutations = np.zeros((cols, native_length), dtype='int32')
    col_permutations[0, :] = np.arange(native_length)
    random.shuffle(col_permutation)
    col_permutations[1, :] = col_permutation

    # Go through the rest of the rows of col_permutation, i.e., the powers of a
    # permutation, and set each to be the previous row permuted once again by
    # the pattern col_permutation.
    if cols > 2:
        for i in range(2, cols):
            col_permutations[i, :] = (col_permutations[i - 1, :])[col_permutation]
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: populate_intensities
# Populate the intensities for RGB image data. This is a 2D array of dimensions
# 256 x native_length, of 0's and 1's, where each row is a hypervector of
# length native_length, representing the intensity of the same value as it's
# row index. Thus, row 0 = intensity 0, row 1 = intensity 1, ..., row 255 =
# intensity 255, each represented symbolically as a random hypervector.
#
# Parameters:
# No input parameters.
#
# Output:
# No output.
def populate_intensities():
    global intensities                  # Mark as a global variable.

    # Call new_vector_set to generate the intensities and assign the result to
    # the global variable intensities.
    intensities = new_vector_set(256)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: populate_locations
# Populate the global locations_basis, which is a basis of hypervectors that
# represents each location in a binding operation, so that sequential data can
# be decoded later.
#
# Parameters:
# num           The number of locations to make in the locations basis, where
#                  this value defines the number of rows in locations_basis.
#
# Output:
# No output.
def populate_locations(num):
    global locations_basis                  # Mark as global variable.

    # Call new_vector_set to generate the locations basis and assign the result
    # to the global variable locations_basis.
    locations_basis = new_vector_set(num)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: xor
# Performs the bitwise XOR operation on hypervectors a and b, and returns the
# result.
#
# Parameters:
# a             The first hypervector to use in a XOR b, as a numpy array.
# b             The second hypervector to use in a XOR b, as a numpy array.
#
# Output:
# Returns the result of a XOR b, the bitwise XOR of a and b, as a numpy vector.
def xor(a, b):
    # Compute the logical vector a != b, which is TRUE where the elements of a
    # and b do not match, and FALSE where they do. Setting this to an int32
    # type transforms the numpy vector into a vector of 1's and 0's.
    return (a != b).astype('int32')
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: hamming
# Measures the Hamming Distance between two vectors a and b, which is the sum
# of their XOR, divided by the native_length, in order to normalize the result
# as a ratio of 1's to 0's. The smaller the result, the closer the vectors are
# to each other.
#
# Parameters:
# a             The first hypervector to measure Hamming Distance between it
#                  and another hypervector.
# b             The second hypervector to measure Hamming Distance between it
#                  and another hypervector.
#
# Output:
# Returns the normalized Hamming Distance between a and b, as a value in the
# range [0, 1].
def hamming(a, b):
    # The hamming distance is defined as the sum of bitwise XOR. To normalize,
    # we divide the result by native_length.
    return xor(a, b).sum() / native_length
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: permute
# Permutes hypervector a's elements according to a given permutations pattern,
# where both are numpy arrays.
#
# Parameters:
# a             The hypervector to permute, as a numpy array.
# pattern       The pattern to permute a with. This is represented as a numpy
#                  array of length equal to a, which is a permutation of the
#                  vector [0, 1, 2, ..., native_length-1], signifying how to
#                  reorder the elements of a.
#
# Output:
# Returns the numpy vector a, with elements reordered according to pattern.
def permute(a, pattern):
    # Easily performed by simply numpy indexing with pattern.
    return a[pattern]
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: permute_power
# Permutes hypervector a repeatedly using a permutation pattern, power-many
# times, which is the pattern specified by a permutation power.
#
# Parameters:
# a             The hypervector to permute, as a numpy array.
# pattern       The pattern to permute a with. This is represented as a numpy
#                  array of length equal to a, which is a permutation of the
#                  vector [0, 1, 2, ..., native_length-1], signifying how to
#                  reorder the elements of a.
# power         The number of times to permute a by the pattern, or the power
#                  of the permutation.
#
# Output:
# Returns the numpy vector a, with elements reordered according to pattern to
# the power'th power, done by repeatedly applying the permutation specified by
# pattern.
def permute_power(a, pattern, power):
    # For each power, we need to permute again.
    for i in range(power):
        # Permute a by the pattern, repeatedly.
        a = permute(a, pattern)

    # Return our permuted vector.
    return a
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: consensus_sum
# Computes the consensus sum of a set of hypervectors. The consensus summation
# is defined as a hypervector, where each output bit is the bit value that has
# more representation across the same component in the set of hypervectors. So
# if the hypervectors have more 1's in that component, that is the results, or
# vice versa if there's more 0's. If there is a tie, the tie breaker is a
# random selection between 0 and 1. This is achieved by adding a random
# hypervector to the set.
#
# Parameters:
# hdmatrix      The 2D numpy matrix that specifies the set of hypervectors to
#                  compute the consensus sum of. The rows of this matrix are
#                  hypervectors of length native_length, which are to be added.
#
# Output:
# Returns a single hypervector as a numpy array, which is the consensus sum of
# the row vectors of hdmatrix.
def consensus_sum(hdmatrix):
    # Get the number of rows in the matrix.
    rows = hdmatrix.shape[0]

    # If the number of rows is even, a tie-breaker is needed, so we add another
    # row to hdmatrix, which is a random hypervector.
    if rows % 2 == 0:
        # Vertically stack the matrix with another hypervector in the last row.
        # This new vector is randomly created. Increase the row count.
        hdmatrix = np.vstack([hdmatrix, new_vector()])
        rows += 1

    # Get the logical vector where 1's are TRUE and 0's are FALSE. Then, sum
    # this matrix along its rows to get a vector of counts of 1's in each
    # component of the set of vectors to sum. Divide this by rows to get the
    # ratio of 1's to 0's in the terms of the sum, in each component.
    logical = (hdmatrix > 0)
    logical = logical.sum(axis=0) / rows

    # The consensus sum is simply a hypervector that is TRUE where the
    # components had more 1's than 0's, and FALSE otherwise. We reinterpret
    # this as a vector of 1's and 0's by changing the type to int32.
    return (logical > 0.5).astype('int32')
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: encode_image_cpu
# Encodes an image I as a hypervector by computing the following on each row of
# pixels, r:
#    R^r I(r, 0) XOR R^r C^1 I(r, 1) XOR ...
# where R is the pattern row_permutation, R^r is the pattern in R permuted r
# times, C is the col_permutation, and applying 2 permutation patterns after
# another permutes them accordingly. The final permutation is applied to the
# hyperdimensional representation of pixel I(r, c), which, if this pixel has
# intensity x, is the x index in the global variable intensities. This gets
# permuted by the permutation pattern, which represents the row and column
# permutation. The XOR of these represents the row of pixels as a hypervector.
# Finally, each row is XOR'd together to get the final image encoding as a
# hypervector.THIS FUNCTION USES THE CPU EXPLICITLY
#
# MUST INITIALIZE BY FIRST CALLING:
#    populate_intensities
#    populate_row_col_permutations
#    populate_locations
#
# Parameters:
# image         A 2D matrix of intensities, for one channel, of an image. You
#                  must call encode_image_cpu 3 times for each RGB channel.
#
# Output:
# Returns a hypervector encoding of the input image in one channel.
def encode_image_cpu(image):
    # Dimensions of the image.
    shape = image.shape

    # The output of encoding; contains the buckets of each row of pixels in the
    # image, and the original shape of the image, which we can pre-populate.
    to_return = {
        'buckets': np.zeros((shape[0], native_length)),
        'shape': shape
                 }

    # Initialize the encoding of the image as a hypervector, as all 0's for
    # now.
    encoding = np.zeros(native_length)

    # Counts of the number of bindings, bundles, and permutations done.
    bindings = 0
    bundles = 0
    permutations = 0

    # For each row of pixels, encode the row as a hypervector, then XOR that
    # with the encoding vector.
    for row in range(shape[0]):
        # The row of pixels.
        r = image[row, :]

        # Initialize the bucket that will contain this row's pixels as a set
        # of hypervectors for each pixel.
        bucket = np.zeros((shape[1], native_length))

        # For each column of the row, i.e., the pixel on that row, encode it
        # and add it
        for col in range(shape[1]):
            # The current pixel is encoded as the hypervector of the intensity
            # of the pixel, permuted by the row permutation of this row,
            # permuted by the column permutation of this column in the image.
            bucket[col, :] = (( intensities[image[row, col], :])[row_permutations[row, :]])[col_permutations[col, :]]

            permutations += 2           # Count the permutations done here.

            # Add this pixel to the encoding by XORing it with encoding.
            # Remember: the encoding is each pixel, permuted to match the
            # row/column location, XOR'd together.
            encoding = xor(encoding, bucket[col, :])
            bindings += 1               # Count the binding done here.

            # The bucket must be XOR'd with the correct location from the
            # locations basis, to indicate it's sequential location. The
            # final results will eventually be
            #    p_1 XOR l_1 + p_2 XOR l_2 + ...
            # where p_i is the pixel in the row sequence, and l_i is the
            # corresponding location hypervector from the locations basis.
            bucket[col, :] = xor(bucket[col, :], locations_basis[col, :])
            bindings += 1               # Count the binding done here.

        # Collapse the bucket into a single hypervector compression by
        # performing consensus sum on the rows of the bucket. This is the
        # hyper vector representing:
        #    p_1 XOR l_1 + p_2 XOR l_2 + ...
        # where p_i is the pixel in the row sequence, and l_i is the
        # corresponding location hypervector from the locations basis.
        to_return['buckets'][row, :] = consensus_sum(bucket)
        bundles += shape[1]             # Count this bundle done here.

    # Set the encoding in our returned object.
    to_return['encoding'] = encoding

    # Print how many of each operation we had to do.
    print("Bindings: ", bindings)
    print("Bundles: ", bundles)
    print("Permutations: ", permutations)

    # Return our encoding and it's buckets, along with any other info.
    return to_return
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Function: decode_image_cpu
# Decodes an image encoded by encode_image_cpu. This is done by probing for
# each intensity value in each row bucket, to find which pixel goes where. The
# result of this will build the full image in one channel, if each pixel is
# correctly found. Note, there will be some error in finding the correct pixel,
# which is dictated by randomness, the odds increase with the number of rows in
# the location basis.THIS FUNCTION USES THE CPU EXPLICITLY
#
# MUST INITIALIZE BY FIRST CALLING:
#    populate_intensities
#    populate_row_col_permutations
#    populate_locations
# or by populating the global variables with what you used to encode.
#
# Parameters:
# buckets       The buckets for each row generated during encoding.
# encoding      The hypervector that represents the encoding of the image.
# shape         The dimensions of the image, as a tuple of rows and columns.
#
# Output:
# Returns a 2D matrix with values between 0 and 255, signifying the intensities
# of each pixel, in one channel.
def decode_image_cpu(buckets, encoding, shape):
    # What will contain our image, eventually, built pixel by pixel.
    image = np.zeros(shape)

    # Counts of each operation we had to do to encode the image.
    bindings = 0
    permutations = 0
    hammings = 0

    # For each row of pixels, take bucket representing that row and decode it.
    for row in range(shape[0]):
        bucket = buckets[row, :]        # Obtain the corresponding bucket.

        # For each column, predict the pixel that must have been at this row/
        # col location in the original image, from the encoding.
        for col in range(shape[1]):
            # Our final permutation for this location, using the row and col
            # permutations that are cached.
            permuted = (row_permutations[row, :])[col_permutations[col, :]]
            permutations += 1           # Count this permutation.

            # Need to keep track of the closest distance seen so far (closest
            # to 0, the best matching intensity, and the corresponding
            # intensity hypervector.
            best_distance = 1.0
            best_match = -1
            best_encoding = np.zeros(native_length)

            # For each possible intensity, probe the bucket to see how good of
            # a match it is in terms of Hamming Distance. If we find a better
            # one, update our best matches.
            for intensity in range(intensities.shape[0]):
                # Our current hypethesis (the currenty intensity hypervector).
                hypothesis = (intensities[intensity, :])[permuted]
                permutations += 1       # Count this permutation.

                # Compute the Hamming Distance between our location basis for
                # this spot in the sequence of row pixels, and the XOR of our
                # hypothesis and bucket. The XOR of the hypothesis and bucket
                # collapses the superposed row pixel vectors from encoding into
                # either pure random noise, or a slightly meaningful deviation
                # for the matching location. By comparing the distance between
                # the location and this result, our best match becomes our
                # prediction for the pixel intensity.
                distance = hamming(
                    xor(hypothesis, bucket),
                    locations_basis[col, :]
                )
                hammings += 1           # Count this Hamming distance.
                bindings += 2           # Count the bindings we had do to here.

                # If this is a better match (smaller distance), update our best
                # match information.
                if distance < best_distance:
                    best_distance = distance
                    best_match = intensity
                    best_encoding = hypothesis

            # Now that we have our prediction for the pixel intensity at this
            # location, set it in the image.
            image[row, col] = best_match

            # XOR our pixel prediction hypervector with our encoding to
            # "remove" it from the encoding.
            encoding = xor(encoding, best_encoding)
            bindings += 1               # Count the binding we had to do here.

    # Print how many operations we had to do.
    print("Bindings: ", bindings)
    print("Permutations: ", permutations)
    print("Hamming Distances: ", hammings)

    # If our encoding has not been decoded during this process to a vector of
    # all zeros, we must have lost at least one pixel.
    if (encoding > 0).sum() > 0:
        print("Pixels have been lost...")

    # Return our predicted image.
    return image
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# END FUNCTION DEFINITIONS ----------------------------------------------------
# -----------------------------------------------------------------------------
