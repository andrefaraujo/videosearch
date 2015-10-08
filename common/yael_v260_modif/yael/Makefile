.SUFFIXES: .c .o .swg 

include ../makefile.inc


all: libyael.a libyael.$(SHAREDEXT) _yael.so


ifeq "$(USEARPACK)" "yes"
  EXTRAYAELLDFLAG=$(ARPACKLDFLAGS)
  EXTRAMATRIXFLAG=-DHAVE_ARPACK
  EXTRAEIGSFLAG=-DHAVE_ARPACK
endif

ifeq "$(USEOPENMP)" "yes"
  OPENMPFLAG = -fopenmp
  EXTRAMATRIXFLAG+=-fopenmp
  EXTRAYAELLDFLAG+=-fopenmp
endif


#############################
# Various  

LIBOBJ=matrix.o vector.o nn.o kmeans.o hkm.o gmm.o eigs.o \
	spectral_clustering.o sorting.o binheap.o machinedeps.o vlad.o kmlsh.o

libyael.a: $(LIBOBJ)
	ar r libyael.a $^

libyael.$(SHAREDEXT): $(LIBOBJ)
	gcc $(LDFLAGS) $(YAELSHAREDFLAGS) -o libyael.$(SHAREDEXT) $^ $(LAPACKLDFLAGS) $(THREADLDFLAGS) $(EXTRAYAELLDFLAG)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@ $(FLAGS) $(EXTRACFLAGS)

nn.o: EXTRACFLAGS=$(LAPACKCFLAGS) $(THREADCFLAGS) $(OPENMPFLAG)

matrix.o: EXTRACFLAGS=$(LAPACKCFLAGS) $(THREADCFLAGS) $(EXTRAMATRIXFLAG)

kmeans.o: EXTRACFLAGS=$(THREADCFLAGS)

hkm.o: EXTRACFLAGS=$(THREADCFLAGS)

gmm.o: EXTRACFLAGS=$(THREADCFLAGS) $(LAPACKCFLAGS)

kmlsh.o: EXTRACFLAGS=$(THREADCFLAGS) $(LAPACKCFLAGS)

eigs.o: EXTRACFLAGS=$(LAPACKCFLAGS) $(LAPACKCFLAGS) $(EXTRAEIGSFLAG)

yael_wrap.o: EXTRACFLAGS=$(PYTHONCFLAGS) $(NUMPYCFLAGS)

mmap_and_lock: mmap_and_lock.c
	$(CC) $(CFLAGS) -o $@ $<



#############################
# Python wrapper 


yael_wrap.c : yael.swg common.swg nn.h kmeans.h gmm.h eigs.h hkm.h \
	sorting.h matrix.h vector.h machinedeps.h vlad.h kmlsh.h
	$(SWIG) -I.. $(NUMPYSWIGFLAGS) $< 

_yael.so: yael_wrap.o libyael.$(SHAREDEXT)
	$(CC) $(LDFLAGS) -o $@ $(WRAPLDFLAGS) $^ $(PYTHONLDFLAGS) $(LAPACKLDFLAGS) $(THREADLDFLAGS) $(YAELLDFLAGS)

#############################
# Dependencies  

# for i in *.c ; do cpp -I.. -MM $i; done

binheap.o: binheap.c binheap.h sorting.h
eigs.o: eigs.c eigs.h vector.h sorting.h machinedeps.h 
gmm.o: gmm.c vector.h matrix.h kmeans.h nn.h gmm.h sorting.h
kmeans.o: kmeans.c vector.h kmeans.h nn.h sorting.h
hkm.o: hkm.c hkm.h kmeans.h vector.h kmeans.h nn.h
machinedeps.o: machinedeps.c machinedeps.h
matrix.o: matrix.c matrix.h vector.h sorting.h machinedeps.h eigs.h
nn.o: nn.c machinedeps.h vector.h nn.h binheap.h sorting.h
sorting.o: sorting.c sorting.h machinedeps.h binheap.h vector.h
spectral_clustering.o: spectral_clustering.c spectral_clustering.h nn.h eigs.h vector.h \
  matrix.h kmeans.h spectral_clustering.h
vector.o: vector.c vector.h
vlad.o: vlad.c vlad.h nn.h vector.h sorting.h
kmlsh.o: kmlsh.c ../yael/vector.h ../yael/nn.h ../yael/hkm.h \
 ../yael/sorting.h ../yael/kmeans.h ../yael/machinedeps.h ../yael/kmlsh.h


clean:
	rm -f libyael.a libyael.$(SHAREDEXT)* *.o yael.py *.pyc yael_wrap.c _yael.so 
