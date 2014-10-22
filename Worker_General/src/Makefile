COMPILER= g++ -g -I../eigen/ -I../include/

all: imcomp_analysis corrtest karytest difficulty_models difficulty_model_test realdatatests kary3_realdatatest

datasets: temporal_data.o rte_data.o imcomp.o imana_data.o mooc_data.o wsd_data.o wordsim_data.o

difficulty_models: worker.o difficulty_models.cpp datasets pool.o test.o
	$(COMPILER) -o difficulty_models difficulty_models.cpp temporal_data.o rte_data.o imcomp_data.o imana_data.o mooc_data.o wsd_data.o wordsim_data.o worker.o pool.o test.o

difficulty_model_test: worker.o difficulty_model_test.cpp datasets pool.o test.o
	$(COMPILER) -o difficulty_model_test difficulty_model_test.cpp temporal_data.o rte_data.o imcomp_data.o imana_data.o mooc_data.o wsd_data.o wordsim_data.o worker.o pool.o test.o

realdatatests: worker.o realdatatests.cpp datasets pool.o test.o
	$(COMPILER) -o  realdatatests realdatatests.cpp temporal_data.o rte_data.o imcomp_data.o imana_data.o mooc_data.o wsd_data.o wordsim_data.o worker.o pool.o test.o
	
kary3_realdatatest: worker.o kary3_realdatatest.cpp datasets pool.o test.o
	$(COMPILER) -o kary3_realdatatest kary3_realdatatest.cpp temporal_data.o rte_data.o imcomp_data.o imana_data.o mooc_data.o wsd_data.o wordsim_data.o worker.o pool.o test.o	

corrtest: worker.o corrtest.cpp pool.o test.o
	$(COMPILER) -o corrtest corrtest.cpp worker.o pool.o test.o

karytest: worker.o karytest.cpp pool.o test.o
	$(COMPILER) -o karytest karytest.cpp worker.o pool.o test.o

imcomp_analysis: worker.o imcomp_analysis.cpp imcomp_data.o pool.o test.o
	$(COMPILER) -o imcomp_analysis imcomp_analysis.cpp imcomp_data.o worker.o pool.o test.o

worker.o: worker.cpp ../include/worker.h
	$(COMPILER) -c worker.cpp

pool.o: pool.cpp ../include/pool.h
	$(COMPILER) -c pool.cpp

test.o: test.cpp ../include/test.h
	$(COMPILER) -c test.cpp

imcomp_data.o: imcomp_data.cpp ../include/imcomp_data.h
	$(COMPILER) -c imcomp_data.cpp
	
temporal_data.o: temporal_data.cpp ../include/temporal_data.h
	$(COMPILER) -c temporal_data.cpp

rte_data.o: rte_data.cpp ../include/rte_data.h
	$(COMPILER) -c rte_data.cpp

imana_data.o: imana_data.cpp ../include/imana_data.h
	$(COMPILER) -c imana_data.cpp

mooc_data.o: mooc_data.cpp ../include/mooc_data.h
	$(COMPILER) -c mooc_data.cpp

wsd_data.o: wsd_data.cpp ../include/wsd_data.h
	$(COMPILER) -c wsd_data.cpp
	
wordsim_data.o: wordsim_data.cpp ../include/wordsim_data.h
	$(COMPILER) -c wordsim_data.cpp

try: try.cpp
	$(COMPILER) -o try try.cpp

clean: 
	rm imcomp_analysis corrtest karytest difficulty_models *.o 
	
	
