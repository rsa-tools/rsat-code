usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' makefile

all: dirs listing_1998 listing program_description new_scripts

dirs:
	mkdir -p material

program_description:
	get-program-description -all -v 1	\
		| add-linenb			\
		-o material/program_descriptions.tab 
	text-to-html -font variable -chunk 200 -border '1 cellpadding=3'	\
		-i material/program_descriptions.tab			\
		-o program_descriptions.html 

REF_DATE=1998_0305
listing_1998:
	ls -lR ../rsat_${REF_DATE}/perl-scripts/ > material/listing_perl_scripts_${REF_DATE}.txt

DATE=`date +%Y_%m%d`
listing:
	ls -lR ../perl-scripts/ > material/listing_perl_scripts_${DATE}.txt

new_scripts:
	ls -1 ../rsat_${REF_DATE}/perl-scripts/ > material/${REF_DATE}
	ls -1  ../perl-scripts/  > material/${DATE}
	diff material/${REF_DATE} material/${DATE} | grep "^>" | perl -pe "s/^> //" > material/listing_new_scripts.txt

archive:
	cd .. ; zip -ry archives/doc_${DATE}.zip doc