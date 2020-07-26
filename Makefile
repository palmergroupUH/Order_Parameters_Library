SHELL = /bin/sh

all: 

	$(MAKE) --directory=general  
	$(MAKE) --directory=IO 
	
	cd IO; $(MAKE) clean  
