
SRC=04-Cameras.tex

TARGDISP=$(SRC:.tex=.pdf)

all: $(TARGDISP)

$(TARGDISP): $(SRC)
	pdflatex $(SRC)

zip:
	zip P_Demo *.m

clean:
	rm -f *~ *log *nav *aux $(TARGDISP) *.out *.snm *.toc \
              *.vrb

