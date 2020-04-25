task_simple_1cty = one_country/src
task_simple_2cty = two_countries/src

all :
	$(MAKE) -C $(task_simple_1cty)
	$(MAKE) -C $(task_simple_2cty)

.PHONY: cleanall
cleanall:
	-rm -r \
	one_country/input \
	two_countries/input \
	\
	one_country/output \
	two_countries/output \
