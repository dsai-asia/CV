
SRC=02-Projective.tex

TARGDISP=$(SRC:.tex=.pdf)

all: $(TARGDISP)

$(TARGDISP): $(SRC)
	pdflatex $(SRC)

zip:
	zip 02-Projective *.tex *.pdf *.jpg Makefile

clean:
	rm -f *~ *log *nav *aux $(TARGPRINT) $(TARGDISP) *.out *.snm *.toc \
       	      *.vrb

