
SRC=06-NView.tex

TARGDISP=$(SRC:.tex=.pdf)

all: $(TARGDISP)

$(TARGDISP): $(SRC)
	pdflatex $(SRC)

clean:
	rm -f *~ *log *nav *aux $(TARGDISP) *.out *.snm *.toc \
              *.vrb

