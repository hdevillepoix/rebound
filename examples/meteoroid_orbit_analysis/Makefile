##----------------------
CFLAGS = -O3 -g -hfp3 -homp -h stdc
CPPFLAGS = -DLIBREBOUND -DOPENMP -D_GNU_SOURCE


all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) $(CFLAGS) -I../../src/ $(CPPFLAGS) -L/group/pawsey0152/hdevillepoix/rebound/src problem.c -lrebound -lm -lrt -o rebound.x
	@echo ""
	@echo "REBOUND compiled successfully."

librebound:
	@echo "Compiling static library librebound.a ..."
	$(MAKE) -C ../../src/

clean:
	@echo "Cleaning up static library librebound.so ..."
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound.x


