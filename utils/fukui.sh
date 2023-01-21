#!/bin/bash
#FUKUI_DIRECTORY="/home/users/nico/software/bin/brabo/"
FUKUI_DIRECTORY="/home/nicux/Documentos/research/code/fortran/fukui/build4"
FUKUI2B_PROGRAM="${FUKUI_DIRECTORY}/fukui.x"
remove_ext() {                   ### Remove extension from a file
	echo "$1" | awk -F . 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
}
atom_number() {
	local var="$1"
#	echo ${var} | awk -F'[(_,.)]' '{gsub(/[[:alpha:]]/, "",$(NF-1)); printf "%i\n",$(NF-1)}'
	echo "${var}" | awk -F_ '{gsub(/[[:alpha:]]/, "",$(NF)); printf "%i\n",$(NF)}'
}
if_not_exist-exit() {
	local file="$1"
	local  ext="$2"
	local text="$3"
	if [ ! -f "${file}" ] 
		then
			if [ ! -f "${file}.${ext}" ] 
				then
					echo " ${file}.${ext} not found"
					exit 1
				else
					eval ${text}="\"${file}\""
			fi
		else
			eval ${text}="\"${file}\""
	fi
}
remove_atom_and_ext() {
	local file="$1"
	echo "${file}" | awk -F _ 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
}
ARGUMS=( $* )
case "${#ARGUMS[@]}" in
	1 ) 
		if_not_exist-exit "${ARGUMS[0]}" 			stk STK_FILE
		if_not_exist-exit "$(remove_atom_and_ext ${ARGUMS[0]})" wfn WFN_FILE
		ATOM_NMB="$(atom_number ${ARGUMS[0]})"					;;
	2 ) 
		if_not_exist-exit "${ARGUMS[1]}" stk STK_FILE
		if_not_exist-exit "${ARGUMS[0]}" wfn WFN_FILE
		ATOM_NMB="$(atom_number ${ARGUMS[1]})"					;;
	* )
		echo " ** Wrong number of arguments ** "
		exit 1									;;
esac
if [ ! -f "${WFN_FILE}.mxal" ]
	then
		awk '/GAUSSIAN/{printf "%s\n%s\n1000000\n%s\n",$7,$2,$5;exit}' "${WFN_FILE}.wfn" > "${WFN_FILE}.mxal"
fi
echo " >> ${FUKUI2B_PROGRAM} : Calculation using wf:${WFN_FILE} and grd:${STK_FILE} for atom:${ATOM_NMB}"
"${FUKUI2B_PROGRAM}" "${WFN_FILE}" "${STK_FILE}" "${ATOM_NMB}"
