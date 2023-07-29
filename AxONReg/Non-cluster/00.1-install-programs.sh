#!/bin/sh -e

##########################################################################
#   Run before all other scripts on supported platforms.
#   Must be run by a systems manager.
##########################################################################

case $(uname) in
FreeBSD)
    # Install ports on all compute nodes
    printf "Root "
    su -m root -c 'pkg install -y rna-seq'
    ;;

*)
    if which auto-pkgsrc-dir; then
	cd $(auto-pkgsrc-dir)/biology/rna-seq
	bmake deinstall clean clean-depends install
    else
	cat << EOM

$0: No pkgsrc installation found.

If you have a pkgsrc tree installed, install sysutils/auto-admin so
that $0 can use auto-pkgsrc-prefix to find it.

Otherwise, install a pkgsrc tree using auto-pkgsrc-setup.  See the
documentation and script at:

    http://netbsd.org/~bacon/

EOM
	exit 1
    fi
    ;;

esac
