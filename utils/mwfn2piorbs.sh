#!/bin/env bash
# general dependencies:
#    bash       (to run this script)
#    util-linux (for getopt)

VERSION="220425"

###### Modify according to your system
multiwfn_dir="~/software/multiwfn/Multiwfn_3.8_dev_bin_Linux"
multiwfn_exe="Multiwfn"
nprocdef=1
def_wfn_exten="wfn"
date=$(date +%g%m%d%H%M%S)
#######################################

PROGNAME="$(basename "$0")"
MIN_REQUIRED_ARGS="1"
# make sure that all command outputs are in English
# so we can parse them correctly
export LC_ALL=C

multiwfn_dir="$(eval readlink -f "${multiwfn_dir}")"
multiwfn="${multiwfn_dir}/${multiwfn_exe}"
[ ! -f "${multiwfn}" -o ! -x "${multiwfn}" ] && echo -e "\n\t**** ${multiwfn} NOT FOUND OR NOT EXECUTABLE !!! ****\n" && exit 1

##################### FUNCTIONS ##########################
usage() {
    cat << AAL

NAME
	$PROGNAME - Extracts from Multiwfn output the pi spin orbitals and generate a .sg file.

SYNOPSIS
	$PROGNAME [OPTIONS]...

DESCRIPTION
	This script extracts from the Multiwfn output the identification of the pi spin orbitals
	and creates a .sg file including the numbering of these orbitals.

	File specification is mandatory.

	-h, --help                      Show this help and exit.

  	-v, --version                   Print version number and exit.

  	-p, --nthreads                  Set number of procs in Multiwfn (default: ${nprocdef})

  	-w, --wfn                       .wfn file necessary to obtain the pi spin orbitals.

AUTHOR
	Written by Nicolás Otero Martínez

AAL
}

remove_ext() {                   ### Remove extension from a file
    echo "$1" | awk -F . 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
}

extract_ext() {
    echo "$1" | awk -F. ' { print $NF }'
}

is_integer() {
    return $(test "$@" -eq "$@" > /dev/null 2>&1);
} 

is_integer_and_gt_x() {
    local input="$1"
    local gtinteger="$2"
    if is_integer "${input}";
        then
            [[ ${input} -gt ${gtinteger} ]] && echo "1" || echo "0"
        else
             echo "0"
    fi
}
proc_output() {
	awk 	'
			/The number of orbitals:/	{	nmo = substr($5,1,length($5)-1)	}
			/Expected pi orbitals/		{
								getline
								a   = 1
							}
			/Total number of pi orbitals/	{	a   = 0 	}
			a				{
								pi[n++]=$1
							}
			END 				{
								nmo++
								nmo--   # It is strange ... Without this trick, nmo is treated uncorrectly
								n = 0
								for(i=1;i<=nmo;i++) 		{ mo[i] = i }
								print nmo,length(pi),nmo-length(pi),"0"
								for(i=0;i<length(pi);i++) 	{ delete mo[pi[i]] }
								for(i=0;i<length(pi)-1;i++) 	{ printf "%s ",pi[i] }
								printf "%s\n",pi[length(pi)-1]
								for(i in mo) 			{ sigma[n++] = mo[i] }
								if (length(sigma) != nmo-length(pi)) print "Error detecting pi orbitals"
								for(i in sigma) 		{ printf "%s ",sigma[i] }
								print "\n"
							}
		'
}
################## MAIN SCRIPT #########################
ARGS=( "$@" )
GETOPT_ARGS=$(getopt -o hvw:p: -l "help","version","wfn:","nthreads:" -n "${PROGNAME}" -- "$@")
[[ $? -ne 0 ]] && exit 1
eval set -- "${GETOPT_ARGS}"
while :
	do
		case "$1" in
			-h|--help)
				usage
				exit 0
				;;
			-v|--version)
				echo -e "\n\n\t\t${PROGNAME} -- version ${VERSION}\n"
				shift
				;;
 			-w|--wfn)
				shift
				file="$1"
				shift
				;;
 			-p|--nthreads)
				shift
				nproc="$1"
                 		[ "$(is_integer_and_gt_x "${nproc}" 0)" != 1 ] && \\
					echo -e "\n\t**** WRONG # THREADS ****\n" && exit 2

				shift
				;;
			--)
				shift
				break
				;;
		esac
done

if [[ ${#ARGS[@]} -lt ${MIN_REQUIRED_ARGS} ]]
	then
		echo -e "\n\n\t\t** INSUFFICIENT NUMBER OF ARGUMENTS (${#ARGS[@]} < ${MIN_REQUIRED_ARGS}) **"
		usage >&2
		exit 0
fi
 
if [ ! -n "${file}" ]
	then
        	echo -e "\n\n\t\t** FILE NOT SPECIFIED. Exiting ... **\n"
        	usage >&2
        	exit 0
fi

if [ -f "${file}" ]
	then
        	wfn_exten="$(extract_ext "${file}")"
		filename_wo_exten="$(remove_ext "${file}")"
	else
		wfn_exten="${def_wfn_exten}"
		filename_wo_exten="${file}"
		[ ! -f "${filename_wo_exten}"."${wfn_exten}" ] && echo -e "\n\t**** File ${filename_wo_exten}.${wfn_exten} not found ****\n" && exit 1
fi

output_file="${filename_wo_exten}.sg"
[ -f "${output_file}" ] && echo " ** Renaming old ${output_file} to ${output_file}.${date}" && mv "${output_file}" "${output_file}.${date}"

nproc="${nproc:-${nprocdef}}"

echo -e " ** Reading ${file} and printing in ${output_file}:\n----------------------"
# 1000 10 #threads 100 22 0 0 0 q   <<< To browse Multiwfn menus
"${multiwfn}" "${file}" << EOF | proc_output | tee "${output_file}"
1000
10
${nproc}
100
22
0
0
0
q
EOF
echo -e "----------------------\n >> Done! <<"

# vim: et sts=4 sw=4
