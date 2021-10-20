#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from runners import run

def merge_all_sub_files(sub_files, output, outdir, chunk_size, ref_dict, pjo) :
    """
    While there are more than `chunk_size` files to merge,
    create chunks of `chunk_size` files to merge
    """

    chunks_to_remove = []
    n_chunk = 0 # outside of loop because if more than 1 iteration of while
    # there would be conflicting filenames appearing
    while len(sub_files) >= chunk_size :
        new_sub_files = [] # used for resetting sub_files list at each iteration of while loop

        sub_chunks = chunks_of_n(sub_files, chunk_size) # each chunk contains chunk_size filenames to merge
        for chunk in sub_chunks :
            n_chunk += 1
            chunk_list_file = os.path.join(outdir, "chunk_{}.list".format(n_chunk))
            merged_chunk = os.path.join(outdir, "chunk_{}.merged.g.vcf".format(n_chunk))
            chunks_to_remove.append(chunk_list_file)
            chunks_to_remove.append(merged_chunk)
            new_sub_files.append(merged_chunk)

            f = open(chunk_list_file, "w")
            for file in chunk :
                f.write(file + "\n")
            f.close()

            # Making and running merge command
            merge_vcfs(chunk_list_file, ref_dict, merged_chunk, pjo)

        sub_files = new_sub_files # reset sub_files with the merged_chunks now

    # Create a list with remaining sub files (chunks if there were > chunk_size files to merge)
    list_file = os.path.join(outdir, "merge.list")
    f = open(list_file, "w")
    for file in sub_files :
        f.write(file + "\n")
    f.close()

    # Making and running merge command
    merge_vcfs(list_file, ref_dict, output, pjo)
    return chunks_to_remove

def merge_vcfs(input, reference_dict, output, java_options) :
    """Create a picard MergeVcfs command and run it"""
    dc_merge = {"j":java_options, "l":input, "d":reference_dict, "o":output}
    cmd = "picard MergeVcfs {j} I={l} D={d} O={o}"
    cmd = cmd.format(**dc_merge)
    print(cmd + "\n")
    run(cmd)

def chunks(l, n):
    """Yield n number of sequential chunks from l."""
    d, r = divmod(len(l), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        yield l[si:si+(d+1 if i < r else d)]

def chunks_of_n(l, n):
    """Yield chunks from l of size n."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
