import os
import argparse
import numpy as np
import hd_encoding

from PIL import Image
from timeit import default_timer as timer


data = "images"
device = "cpu"

PARSER = argparse.ArgumentParser()

PARSER.add_argument("path", type=str)
PARSER.add_argument("--data", nargs='?', type=str, const="images")
PARSER.add_argument("--device", "-d", nargs='?', type=str, const="cpu")


def run_cpu_images(file_list):
    result_map = {
        'images': file_list,
        'encodings': {},
        'locations_basis': hd_encoding.locations_basis,
        'row_permutation': hd_encoding.row_permutation,
        'row_permutations': hd_encoding.row_permutations,
        'col_permutation': hd_encoding.col_permutation,
        'col_permutations': hd_encoding.col_permutations,
        'intensities': hd_encoding.intensities
    }

    for image in file_list:
        im = Image.open(image)

        I = np.array(im)

        r = I[:, :, 0]
        g = I[:, :, 1]
        b = I[:, :, 2]

        timed = timer()
        r_channel = hd_encoding.encode_image_cpu(r)
        g_channel = hd_encoding.encode_image_cpu(g)
        b_channel = hd_encoding.encode_image_cpu(b)
        print("Elapsed Time Encoding: ", timer() - timed)

        result_map['encodings'][image] = {
            'r': r_channel,
            'g': g_channel,
            'b': b_channel
        }

        timed = timer()
        r_decoded = hd_encoding.decode_image_cpu(
            r_channel['buckets'],
            r_channel['encoding'],
            r_channel['shape']
        )
        g_decoded = hd_encoding.decode_image_cpu(
            g_channel['buckets'],
            g_channel['encoding'],
            g_channel['shape']
        )
        b_decoded = hd_encoding.decode_image_cpu(
            b_channel['buckets'],
            b_channel['encoding'],
            b_channel['shape']
        )
        print("Elapsed Time Decoding: ", timer() - timed)

        print("Missed pixels in r channel", ((r - r_decoded.astype("int32")) > 0).sum())
        print("Missed pixels in g channel", ((g - g_decoded.astype("int32")) > 0).sum())
        print("Missed pixels in b channel", ((b - b_decoded.astype("int32")) > 0).sum())

        im.close()

    return result_map


def run_gpu_images(file_list):
    return True


def run_cpu_text(file_list):
    return True


def run_gpu_text(file_list):
    return True


def run_cpu_numbers(file_list):
    return True


def run_gpu_numbers(file_list):
    return True


def find_max_dimensions(file_list):
    rows = 0
    cols = 0

    for image in file_list:
        im = Image.open(image)

        I = np.array(im)

        shape = I.shape

        if rows < shape[0]:
            rows = shape[0]

        if cols < shape[1]:
            cols = shape[1]

        im.close()

    return [rows, cols]


if __name__ == "__main__":
    args = PARSER.parse_args()
    path = args.path

    if args.data is not None:
        if args.data == "images":
            data = "images"
        elif args.data == "text":
            data = "text"
        elif args.data == "numbers":
            data = "numbers"
        else:
            data = "images"

    if args.device is not None:
        if args.device == "cpu":
            device = "cpu"
        elif args.device == "gpu":
            device = "gpu"
        else:
            device = "cpu"


    print("path: ", path)
    print("data: ", data)
    print("device: ", device)
    print("\n")


    if os.path.isdir(path):
        files = os.listdir(path)
        files = [os.path.join(path, x) for x in files]
    elif os.path.isfile(path):
        files = [path]
    else:
        raise Exception("Path provided is an invalid folder or file: ", path)

    print(files)

    if data == "images":
        [max_rows, max_cols] = find_max_dimensions(files)

        hd_encoding.populate_intensities()
        hd_encoding.populate_row_col_permutations(max_rows, max_cols)
        hd_encoding.populate_locations(max_cols)

        if device == "gpu":
            encoded = run_gpu_images(files)
        else:
            encoded = run_cpu_images(files)

    elif data == "text":
        if device == "gpu":
            encoded = run_gpu_text(files)
        else:
            encoded = run_cpu_text(files)

    elif data == "numbers":
        if device == "gpu":
            encoded = run_gpu_numbers(files)
        else:
            encoded = run_cpu_numbers(files)
    else:
        [max_rows, max_cols] = find_max_dimensions(files)

        hd_encoding.populate_intensities()
        hd_encoding.populate_row_col_permutations(max_rows, max_cols)
        hd_encoding.populate_locations(max_cols)

        if device == "gpu":
            encoded = run_gpu_images(files)
        else:
            encoded = run_cpu_images(files)
