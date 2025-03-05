#!/bin/tcsh

# Setup AMS environment
setenv AMSSRC $AMSWD

setenv FLAG ""
if ( "$1" == "debug" ) then
    setenv FLAG $1
else if ( "$1" == "clean" ) then
    make clean && exit
endif

setenv USEPRHEION 1
unsetenv USEEVENTORDER
unsetenv USEHEL1
unsetenv USENEWL1L9G
unsetenv USEMCTKRAW
unsetenv USEONEEV
unsetenv USEADDTKHIT
unsetenv USENOLINEARCOR
unsetenv USECALIB
setenv USEHEINNER 1
setenv USESAVENEG 1
setenv USEPRL1 1

# Compile AMSVX static library first
echo
echo "build AMSVX..."
echo
set current_dir=`pwd`
cd $AMSWD/install/
if ( "$FLAG" == "debug" ) then
    make debug_static PGTRACK=1 -j8
else
    make static PGTRACK=1 -j8
endif
cd $current_dir

# Handle command line arguments
echo
echo "build Pass8..."
echo
setenv USEPASS7 1

make $FLAG -j4