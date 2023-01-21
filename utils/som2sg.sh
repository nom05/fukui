#!/bin/bash
### .som -> .sg
inpfile=${1:?Input file needed}

remove_ext() {                   ### Remove extension from a file
	echo "$1" | awk -F . 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
	#echo ${1%.*}
}

outfile="$(remove_ext "${inpfile}").sg"

awk_check() {
	local inpfile=$1
	awk 	'BEGIN	{
        			pidim=0
        			sigdim=0
        			cerdim=0
			}
		/PI/	{
        			pidim++
        			a[pidim]=$2
			}
		/SIGMA/	{
        			sigdim++
        			b[sigdim]=$2
			}
		/CERO/	{
        			cerdim++
        			c[cerdim]=$2
			}
		END	{
        			print "TOTMO=" NR/6
        			print "PIMO="  pidim
        			print "SIGMO=" sigdim
        			print "CERMO=" cerdim
        			print "PIMAT=("
        			for(i=1 ; i<=pidim ; i++)	{
                							print a[i]
        							}
        			print ")"
        			print "SIGMAT=("
        			for(i=1 ; i<=sigdim ; i++)	{
                							print b[i]
        							}
        			print ")"
        			print "CERMAT=("
        			for(i=1 ; i<=cerdim ; i++)	{
                							print c[i]
        							}
        			print ")"
}
	' "${inpfile}"
}
eval $(awk_check $inpfile)
echo "$TOTMO $PIMO $SIGMO $CERMO" | tee "${outfile}"
echo ${PIMAT[@]} | tee -a "${outfile}"
echo ${SIGMAT[@]} | tee -a "${outfile}"
echo ${CERMAT[@]} | tee -a "${outfile}"
