sources := 		RcppEigen-Introduction.Rnw
vignettes := 		RcppEigen-Introduction.pdf

all:			${vignettes}

${vignettes}: 		${sources}
			Rscript -e 'utils::Sweave("RcppEigen-Introduction.Rnw")'
			pdflatex RcppEigen-Introduction.tex
			bibtex RcppEigen-Introduction.aux
			pdflatex RcppEigen-Introduction.tex
			pdflatex RcppEigen-Introduction.tex
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			cp $@ ..

clean:
			@rm -rf *.aux *.log *.out *.tex *.pdf *.bbl *.blg pinp.cls auto/
