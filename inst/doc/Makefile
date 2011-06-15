
all: 	RcppEigen-unitTests.pdf

clean:	
	rm -f RcppEigen*.log RcppEigen*.aux RcppEigen*.out RcppEigen*.tex

pdfclean:
	rm -f RcppEigen-unitTests.pdf

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

RcppEigen-unitTests.pdf:
	rm -fr unitTests-results/*
	$(RSCRIPT) unitTests/RcppEigen-unitTests.R
	$(RPROG) CMD Sweave RcppEigen-unitTests.Rnw
	$(RSCRIPT) -e "tools::texi2dvi( 'RcppEigen-unitTests.tex', pdf = TRUE, clean = TRUE )"
	rm -fr RcppEigen-unitTests.tex

