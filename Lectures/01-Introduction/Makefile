
SRC=01-Introduction.tex

TARGDISP=$(SRC:.tex=.pdf)

all: $(TARGDISP)

$(TARGDISP): $(SRC)
	pdflatex $(SRC)

zip:
	zip 01-Introduction *.tex *.pdf *.jpg Makefile

clean:
	rm -f *~ *log *nav *aux $(TARGPRINT) $(TARGDISP) *.out *.snm *.toc \
       	      *.vrb

