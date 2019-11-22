"""Script to run stacking scripts on the DESY cluster.

Through use of argparse, a given configuration for the code can be selected.
This can be given from the command line, in the form:

python RunCluster.py -c Desired_Configuration_Name -n Number_Of_Tasks -s

Each available configuration must be listed in "config.ini", and controls
options for fitting, such as which catalogue is to be used, and which seasons
of data should be included. If -x is included, then a new job is submitted
to the cluster. Having submitted the job to the cluster it will be run in
parallel Number_of_Tasks times. The shell script SubmitOne.sh is called for
each task, which in turn calls RunLocal.py with the given configuration setting.

The Wait function will periodically query the cluster
to check on the status of the job, and will output the job status occasionally.

Once all sub-tasks are completed, the script will proceed to call
MergeFiles.run() for the given configuration, combining results.

"""
import subprocess
import time
import os
import os.path
import argparse
from flarestack.shared import log_dir, fs_dir
from flarestack.cluster.make_desy_cluster_script import make_desy_submit_file, submit_file

username = os.path.basename(os.environ['HOME'])

cmd = 'qstat -u ' + username


def wait_for_cluster():
    """Runs the command cmd, which queries the status of the job on the
    cluster, and reads the output. While the output is not an empty
    string (indicating job completion), the cluster is re-queried
    every 30 seconds. Occasionally outputs the number of remaining sub-tasks
    on cluster, and outputs full table result every ~ 8 minutes. On
    completion of job, terminates function process and allows the script to
    continue.
    """
    time.sleep(10)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    tmp = process.stdout.read().decode()
    i = 31
    j = 6
    while tmp != '':
        if i > 3:

            n_total = len(str(tmp).split('\n')) - 3

            running_process = subprocess.Popen(
                cmd + " -s r", stdout=subprocess.PIPE, shell=True)
            running_tmp = running_process.stdout.read().decode()

            if running_tmp != '':
                n_running = len(running_tmp.split('\n')) - 3
            else:
                n_running = 0

            print(time.asctime(time.localtime()), n_total, "entries in queue.")
            print("Of these,", n_running, "are running tasks, and")
            print(n_total-n_running, "are jobs still waiting to be executed.")
            print(time.asctime(time.localtime()), "Waiting for Cluster")
            i = 0
            j += 1

        time.sleep(30)
        i += 1
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        tmp = process.stdout.read().decode()


def submit_to_cluster(path, n_cpu=2, n_jobs=10):

    for file in os.listdir(log_dir):
        os.remove(log_dir + file)

    # Submits job to the cluster

    submit_cmd = "qsub "

    if n_cpu > 1:
        submit_cmd += " -pe multicore {0} -R y ".format(n_cpu)

    ram_per_core = "{0:.1f}G".format(6./float(n_cpu) + 2.)
    print("Ram per core:", ram_per_core)

    submit_cmd += "-t 1-{0}:1 {1} {2} {3}".format(
        n_jobs, submit_file, path, n_cpu
    )

    make_desy_submit_file(ram_per_core)

    print(time.asctime(time.localtime()), submit_cmd, "\n")
    os.system(submit_cmd)


if not os.path.isfile(submit_file):
    make_desy_submit_file()
