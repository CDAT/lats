#export CFLAGS="-Wall -g -m64 -pipe -O2  -fPIC"
#export CXXLAGS="${CFLAGS}"
#export CPPFLAGS="-I${PREFIX}/include"
#export LDFLAGS="-L${PREFIX}/lib"

if [ `uname` == Linux ]; then
    echo "Linux  "${PREFIX}
else
    echo "Mac  "${PREFIX}
fi
sed "s#@cdat_EXTERNALS@#${PREFIX}#g;" Makefile.gfortran.in > Makefile.gfortran
make  -f Makefile.gfortran
make -f Makefile.gfortran install
