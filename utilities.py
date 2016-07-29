#!/usr/bin/python
import subprocess, datetime, hashlib, os

from constants import *


def run_cmd(cmd, returnout=True, shell=False, pwd='./'):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # TODO remove when done testing
    print cmd

    # has the command so we can match the logs together
    m = hashlib.md5()
    m.update(str(cmd))

    # log the command
    with open(pwd + 'log/luigi.cmd.log', 'a+') as fout:
        fout.write("{time} {hash}# {actn}\n".format(time=str(datetime.datetime.now()),
                                                    hash=m.hexdigest(),
                                                    actn=" ".join(cmd)))

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=shell,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

    if returnout:
        return stdout
    else:
        # TODO send output to a log file
        # with open('', 'w') as fout:
        pass


def unzip_file(gzip):
    """
    Unzip a gzipped file, using multi-threading when available

    :param gzip: The path to the gzip file
    :return: The uncompressed data stream
    """
    try:
        # use unpigz for multithreaded unzipping (if installed)
        return run_cmd(["unpigz",
                        "-c",                 # output to stdout
                        "-p", MAX_CPU_CORES,  # use the maximum cores
                        gzip])                # gzip file

    except OSError as e:
        # if it's not installed
        if e.errno == os.errno.ENOENT:

            # use single-threaded gunzip
            return run_cmd(["gunzip",
                            "-c",             # output to stdout
                            gzip])            # gzip file
        else:
            # escalate the exception
            raise e


def curl_download(url, filename):
    """
    Downloads a remote url to a local file path using cURL
    """

    # download the file
    run_cmd(["curl",
             "-s",                  # download silently
             "--output", filename,  # output path
             url])                  # from this url