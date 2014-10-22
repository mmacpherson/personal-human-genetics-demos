all: nancy_g.23andMe.txt marcia_g.23andMe.txt joseph_g.23andMe.txt

nancy_g.23andMe.txt:
	curl https://opensnp.org/data/662.23andme.340?1352935156 -o nancy_g.23andMe.txt

marcia_g.23andMe.txt:
	curl https://opensnp.org/data/747.23andme.353?1353906096 -o marcia_g.23andMe.txt

joseph_g.23andMe.txt:
	curl https://opensnp.org/data/663.23andme.305?1350666749 -o joseph_g.23andMe.txt

clean:
	rm -f nancy_g.23andMe.txt marcia_g.23andMe.txt joseph_g.23andMe.txt
