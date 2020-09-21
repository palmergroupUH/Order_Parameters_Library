SHELL = /bin/sh

all: 

	$(MAKE) --directory=general  
	$(MAKE) --directory=IO 
	$(MAKE) --directory=algorithm
	$(MAKE) --directory=mathlib 
	$(MAKE) --directory=order_parameter 
	$(MAKE) --directory=pair_correlation
	$(MAKE) --directory=pub1
 	
	cd IO; $(MAKE) clean  
	cd algorithm; $(MAKE) clean  
	cd mathlib; $(MAKE) clean  
	cd order_parameter; $(MAKE) clean  
	cd pair_correlation; $(MAKE) clean  
	cd pub1; $(MAKE) clean  


