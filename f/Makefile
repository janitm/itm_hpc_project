# Makefile for After post-processing programs
#
#  Modified by Ben Murphy - 3/11/08
#
#################################################################
#
#  Portland Group Compiler Flags (pgf77)
#
#    -tp p6		specify the target processor 
#			  (p6=pentium pro/pentium II)
#    -pc 64		precision option (64=double precision)
#    -Kieee		Perform floating-point operations 
#			  in strict conformance with the IEEE 
#			  754 standard
#    -Mdalign		Align doubles in structures on 8-byte
#			  boundaries
#    -Mextend		Allow 132 column source lines
#    -Mnoframe		Do not set up a true stack frame pointer 
#			  for fxns
#    -byteswapio	Swap from big-endian to little-endian
#			  or vice-versa on input/output of
#			  unformatted FORTRAN data
#    -W1,		Pass each comma-delimeted option to 
#			  the linker
#       -Bstatic	Static linking
#    -Minfo		Emit useful info to stderr
#
#    -c			skip the link step; compile and assemble 
#			  only
#    -o			output file
#
#  Debugging Options
#
#    -g			generate symbolic debug info.
#			sets optimizations level to zero.
#
################################################################




FC = pgf90


#FLGS = -g -tp p6 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -Wl,-Bstatic -Minfo
#note: it is not sensible to use -wl,-Bstatic, since this links everything statically and
#and that means that the libraries need to be in the correct order

#second note: the -tp p6 flag is really important, but the default is actually pretty sensible

FLGS = -g -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio  -Minfo

runme: runme.f hadvppm.f
	$(FC) $(FLGS) -c runme.f hadvppm.f
	$(FC) $(FLGS) runme.o hadvppm.o -o runme


clean: runme
	rm *.o runme
