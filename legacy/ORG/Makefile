lego:  mk_config_param.o mk_config_calc.o mk_config.o
	gcc  -O  -o lego -lm mk_config.o mk_config_param.o mk_config_calc.o 

mk_config_param.o:  mk_config_param.c 
	gcc  -O  -static -Wall -c  mk_config_param.c 

mk_config_calc.o:  mk_config_calc.c  mk_config_globals.h 
	gcc  -O  -static -Wall -c mk_config_calc.c 

mk_config.o: mk_config.c  mk_config_globals.h
	gcc  -O  -static -Wall -c mk_config.c
clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O