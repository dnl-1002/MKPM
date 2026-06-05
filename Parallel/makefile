FLAG = -O3 -std=c++14 -msse4.2 -Ofast -fopenmp -DOPENMP

#FLAG = -g --std=c++14 -O2 
ifdef PROMPT
	FLAG := $(FLAG) -DPROMPT
endif

ifdef OPT
	FLAG := $(FLAG) -Ofast
endif
#-fcilkplus

all: 
	g++ $(FLAG) kPlexEnum.cpp -o PlexEnum
	g++ toBin.cpp -o toBin
clean:
	rm PlexEnum_par
