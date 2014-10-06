#!/bin/bash -e

BOOST_DIR=/Developer/usr/lib/boost_1_36_0/
YASMIC_DIR=.

source ccfiles.sh
OFILES=`echo ${CCFILES} | sed -e 's/\.cc/\.o/g'`

CFLAGS="-O2 -arch x86_64 -DMATLAB_BGL_LARGE_ARRAYS -fPIC -c -I${BOOST_DIR} -I${YASMIC_DIR}"
#CFLAGS="-g -W -DMATLAB_BGL_LARGE_ARRAYS -fPIC -c -I${BOOST_DIR} -I${YASMIC_DIR}"

function echocmd {
	echo $@
	$@
}

for file in ${CCFILES}; do
	echocmd g++ $CFLAGS $file
done

echocmd ar rc libmbgl-macosx-intel-64-large.a ${OFILES} 
	
echocmd rm ${OFILES}
