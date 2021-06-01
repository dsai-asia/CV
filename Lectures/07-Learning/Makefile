
SRC=07-Learning.tex

TARGDISP=$(SRC:.tex=.pdf)

all: $(TARGDISP) #prog1

$(TARGDISP): $(SRC) *.tex
	pdflatex $(SRC)

zip:
	zip -r PCA-Demo pca_test.m pca_train.m run.m Jigme.png

clean:
	rm -f *~ *log *nav *aux $(TARGDISP) *.out *.snm *.toc \
              *.vrb

