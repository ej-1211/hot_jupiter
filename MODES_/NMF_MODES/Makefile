include ./Mkinclude

all:
	cd share && $(MAKE)
	pwd
	cd hough && $(MAKE)
	pwd
	cd VSF && $(MAKE)
	pwd
	cd normal && $(MAKE)
	pwd
	cd main && $(MAKE)
	pwd
	cd gauss && $(MAKE)

clean:
	cd share && $(MAKE) clean
	cd hough && $(MAKE) clean
	cd VSF && $(MAKE) clean
	cd normal && $(MAKE) clean
	cd main && $(MAKE) clean
	cd gauss && $(MAKE) clean
	cd bin && rm *
