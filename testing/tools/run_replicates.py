#!/usr/bin/env python3 

import sys
import os
import shutil
import subprocess
import re
import multiprocessing
import hashlib
import glob

def hash_text(text) -> int:
    # Create hash. To make this resistant to line-end changes we will convert all whitespace to text,
    # and collapse repeated spaces. Eventually we get a single line
    text=re.sub("RNGSeed\s+[-]?[0-9]+",f"RNGSeed 0", text)
    text_prev = None
    while text != text_prev:
        text_prev = text
        text = text.replace("  "," ")
        text = text.replace("\t"," ")
        text = text.replace("\n"," ")
        text = text.replace("\r"," ")
    # https://stackoverflow.com/a/16008760
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[0:6]

def calc_dmpci_file_hash(directory, dmpci_name=None) -> int:
    """
    calc_dmpci_file_hash(dmpci_file_path)       : calculate hash of given file
    calc_dmpci_file_hash(directory, dmpci_name) : calculate hash of {directory}/dmpci.{dmpci_name}
    """
    if dmpci_name==None:
        src_path=directory
    else:
        src_path=f"{directory}/dmpci.{dmpci_name}"
    with open(src_path, "r") as src:
        return hash_text(src.read())

def make_dir_name(n):
    if n<=10:
        return f"0000-0009/{n}"
    elif n<=100:
        return f"0010-0099/{n}"
    elif n<1000:
        return f"0100-0999/{n}"
    else:
        return f"1000-9999+/{n}"
    
def do_task(params):
    (i,vars)=params
    working_root_dir=vars["working_root_dir"]
    dmpci_name=vars["dmpci_name"]
    dmpci_contents=vars["dmpci_contents"]
    dmpci_hash=vars["dmpci_hash"]
    osprey_path=vars["osprey_path"]

    replicate_dir=f"{working_root_dir}/{make_dir_name(i)}"
    os.makedirs(replicate_dir, exist_ok=True)
    replicate_dmpci_file=f"{replicate_dir}/dmpci.{dmpci_name}-{dmpci_hash}"

    if os.path.exists(f"{replicate_dir}/finished.{dmpci_name}-{dmpci_hash}"):
        return

    os.makedirs(replicate_dir,exist_ok=True)
    replica_dmpci_contents=re.sub("RNGSeed\s+[-]?[0-9]+",f"RNGSeed {i}", dmpci_contents)

    with open(replicate_dmpci_file,"w") as dst:
        dst.write(replica_dmpci_contents)
    #print(replica_dmpci_contents)

    subprocess.run([osprey_path,f"{dmpci_name}-{dmpci_hash}"], cwd=replicate_dir, capture_output=True)

    if i != 0:
        for x in glob.glob(os.path.join(replicate_dir,f"{dmpci_name}-{dmpci_hash}.*.vtk")):
            os.remove(x)
        for x in glob.glob(os.path.join(replicate_dir,f"{dmpci_name}-{dmpci_hash}.*.dat")):
            os.remove(x)

    with open(f"{replicate_dir}/finished.{dmpci_name}-{dmpci_hash}", "w") as dst:
        dst.write("ok")

if __name__=="__main__":
    osprey_path=sys.argv[1]
    osprey_path=os.path.abspath(osprey_path)
    assert os.path.isfile(osprey_path)

    dmpci_file=sys.argv[2]
    assert os.path.isfile(dmpci_file)
    m=re.match(".*dmpci[.]([-_a-zA-Z0-9]+)$", dmpci_file)
    assert m
    dmpci_name=m.group(1)
    sys.stderr.write(f"Source = {dmpci_file}\n")

    with open(dmpci_file, "r") as src:
        dmpci_contents=src.read()
    # Create normalised version on hash 0
    dmpci_hash=hash_text(dmpci_contents)

    if len(sys.argv) < 4:
        replicates=1000
    else:
        replicates=int(sys.argv[3])

    if len(sys.argv) < 5:
        working_root_dir=f"./working"
        os.makedirs(working_root_dir,exist_ok=True)
    else:
        working_root_dir=sys.argv[4]


    assert os.path.isdir(working_root_dir)


    sys.stderr.write(f"dmpci-file=X{dmpci_file}X, dmpci-name=X{dmpci_name}X")

    vars = {}
    vars["working_root_dir"]=working_root_dir
    vars["dmpci_name"]=dmpci_name
    vars["dmpci_contents"]=dmpci_contents
    vars["dmpci_hash"]=dmpci_hash
    vars["osprey_path"]=osprey_path

    with multiprocessing.Pool() as pool:
        done=0
        for _ in pool.imap_unordered(do_task,  [(i,vars) for i in range(replicates)]):
            done+=1
            if (done%10)==0:
                sys.stderr.write(f"Done {done} of {replicates}\n")
