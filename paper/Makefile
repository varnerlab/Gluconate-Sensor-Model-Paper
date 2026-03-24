TEX = pdflatex
BIB = bibtex
PAPER = Paper

.PHONY: all clean distclean

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex References.bib
	$(TEX) $(PAPER)
	$(BIB) $(PAPER)
	$(TEX) $(PAPER)
	$(TEX) $(PAPER)

clean:
	rm -f $(PAPER).aux $(PAPER).bbl $(PAPER).blg $(PAPER).log \
	      $(PAPER).out $(PAPER).toc $(PAPER).lof $(PAPER).lot \
	      $(PAPER).fls $(PAPER).fdb_latexmk $(PAPER).synctex.gz \
	      $(PAPER).run.xml $(PAPER)-blx.bib

distclean: clean
	rm -f $(PAPER).pdf
