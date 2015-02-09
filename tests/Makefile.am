noinst_HEADERS = radar.h polar.h odim.h compo.h

check_PROGRAMS = decbufr build-composite build-compo-cumul build-local build-cumul \
	build-wrwp build-wrwp-OMM build-polar build-antilles build-odim-polar bufr2hdf5 hdf2bufr \
	build-odim-compo comp2hdf5 hdf2comp pdh2comp otan2comp

decbufr_SOURCES = decbufr.c
build_composite_SOURCES = build-composite.c
build_compo_cumul_SOURCES = build-compo-cumul.c
build_local_SOURCES = build-local.c
build_cumul_SOURCES = build-cumul.c
build_wrwp_SOURCES = build-wrwp.c
build_wrwp_OMM_SOURCES = build-wrwp-OMM.c
build_polar_SOURCES = build-polar.c
build_antilles_SOURCES = build-antilles.c
build_odim_polar_SOURCES = build-odim-polar.c odim.c
bufr2hdf5_SOURCES = bufr2hdf5.c odim.c
hdf2bufr_SOURCES = hdf2bufr.c odim.c
build_odim_compo_SOURCES = build-odim-compo.c odim.c
comp2hdf5_SOURCES = comp2hdf5.c odim.c
hdf2comp_SOURCES = hdf2comp.c odim.c
pdh2comp_SOURCES = pdh2comp.c odim.c
otan2comp_SOURCES = otan2comp.c odim.c

LDADD = ../lib/libOperaBufr.a

REF = -ref

REFERENCES = composite$(REF).bfr cumul$(REF).bfr local$(REF).bfr wrwp-opera$(REF).bfr wrwp-OMM$(REF).bfr \
	polar$(REF).bfr antilles-CMO$(REF).bfr compo-cumul$(REF).bfr odim-polar$(REF).bfr odim-compo$(REF).bfr

EXTRA_DIST = $(REFERENCES) \
	french-200010301300.bfr french-cumul-200605180605.bfr cumul-abbeville-200512271800.bfr \
	polar-20040619.bfr abbeville-200512271800.bfr wrwp-wmo.bfr antilles.bfr \
	200906150321_t1 200906150323_t2 200906150325_t3 RAD_NL60_VOL_OP_200705211200.h5 \
	Euroc_200911291715.bfr T_PAGZ43_C_EFKL_20101217061505.h5_v2.1

check: 
	./build-compo-cumul -d ../tables french-cumul-200605180605.bfr compo-cumul.bfr 2>log 
	(cmp compo-cumul.bfr compo-cumul$(REF).bfr && echo "PASS compo-cumul") || echo "FAIL compo-cumul"
	./build-composite -d ../tables french-200010301300.bfr composite.bfr 2>>log 
	(cmp composite.bfr composite$(REF).bfr && echo "PASS composite") || echo "FAIL composite"
	./build-local -d ../tables abbeville-200512271800.bfr local.bfr 2>>log 
	(cmp local.bfr local$(REF).bfr && echo "PASS local") || echo "FAIL local"
	./build-cumul -d ../tables cumul-abbeville-200512271800.bfr cumul.bfr 2>>log 
	(cmp cumul.bfr cumul$(REF).bfr && echo "PASS cumul") || echo "FAIL cumul"
	./build-wrwp -d ../tables wrwp-wmo.bfr wrwp-opera.bfr 2>>log 
	(cmp wrwp-opera.bfr wrwp-opera$(REF).bfr && echo "PASS wrwp-opera") || echo "FAIL wrwp-opera"
	./build-wrwp-OMM -d ../tables wrwp-wmo.bfr wrwp-OMM.bfr 2>>log 
	(cmp wrwp-OMM.bfr wrwp-OMM$(REF).bfr && echo "PASS wrwp-OMM") || echo "FAIL wrwp-OMM"
	./build-polar -d ../tables polar-20040619.bfr polar.bfr 2>>log 
	(cmp polar.bfr polar$(REF).bfr && echo "PASS polar") || echo "FAIL polar"
	./build-antilles -d ../tables antilles.bfr antilles-CMO.bfr 2>>log
	(cmp antilles-CMO.bfr antilles-CMO$(REF).bfr && echo "PASS antilles-CMO") || echo "FAIL antilles-CMO"

odim: build-odim-polar
	./build-odim-polar -d ../tables 200906150321_t1 200906150323_t2 200906150325_t3 odim-polar.bfr 2>log 
	(cmp odim-polar.bfr odim-polar$(REF).bfr && echo "PASS odim-polar") || echo "FAIL odim-polar"
bufr2hdf: bufr2hdf5
	./bufr2hdf5 -d ../tables odim-polar.bfr odim-polar.hdf5 2>>log 
hdf2bfr: hdf2bufr
	./hdf2bufr -d ../tables odim-polar.hdf5 odim-polar.bfr2 2>>log 

odim2: bufr2hdf5 hdf2bufr
	./hdf2bufr -d ../tables RAD_NL60_VOL_OP_200705211200.h5 odim-suede.bfr 2>>log 
	./bufr2hdf5 -d ../tables odim-suede.bfr odim-suede.h5 2>>log 
	./hdf2bufr -d ../tables odim-suede.h5 odim-suede.bfr2 2>>log 
	(cmp odim-suede.bfr odim-suede.bfr2 && echo "PASS odim-suede") || echo "FAIL odim-suede"

compo: build-odim-compo
	./build-odim-compo -d ../tables Euroc_200911291715.bfr odim-compo.bfr 2>>log 
	(cmp odim-compo.bfr odim-compo$(REF).bfr && echo "PASS odim-compo") || echo "FAIL odim-compo"
comp2hdf: comp2hdf5
	./comp2hdf5 -d ../tables odim-compo.bfr odim-compo.hdf5 2>>log 
hdf2com: hdf2comp
	./hdf2comp -d ../tables odim-compo.hdf5 odim-compo.bfr2 2>>log 

# test for ODIM v2.1
hdf21: hdf2bufr bufr2hdf5
	./hdf2bufr -d ../tables ./T_PAGZ43_C_EFKL_20101217061505.h5_v2.1 odim-v2.1.bfr 2>>log 
	./bufr2hdf5 -d ../tables/ odim-v2.1.bfr odim-v2.1.h5 2>>log 
	./hdf2bufr -d ../tables odim-v2.1.h5 odim-v2.1.bfr2 2>>log 
	(cmp odim-v2.1.bfr odim-v2.1.bfr2 && echo "PASS odim-v2.1") || echo "FAIL odim-v2.1"

idempotent: odim compo bufr2hdf hdf2bfr comp2hdf hdf2com odim2 hdf21
	(cmp odim-polar.bfr odim-polar.bfr2 && echo "PASS idem odim-polar") || echo "FAIL idem odim-polar"
	(cmp odim-compo.bfr odim-compo.bfr2 && echo "PASS idem odim-compo") || echo "FAIL idem odim-compo"
