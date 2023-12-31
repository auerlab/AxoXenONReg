#!/bin/sh -e

##########################################################################
#   Run before all other scripts on supported platforms.
#   Must be run by a systems manager.
##########################################################################

cluster_run=cluster-run
srun="srun --ntasks=1 --mem=1g"
node_spec=compute

case $(uname) in
FreeBSD)
    # Install ports on all compute nodes
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y rna-seq' $node_spec"
    ;;

*)
    # Check for pkgsrc installed via auto-pkgsrc-setup
    if which sbatch; then
	cat << EOM

$0: You appear to be using a non-FreeBSD cluster.

You can use pkgsrc and install the biology/rna-seq package on all
compute nodes or in a shared location that compute nodes can access.

EOM
    fi
    exit 1
    ;;

esac
