# Makefile to generate fig_noise.pdf
# Reinhard Heckel, 2013


all: fig_noise.pdf

fig_noise.pdf:	CE_huge_noise.dat 
	pdflatex -halt-on-error -interaction=batchmode fig_noise.tex 

CE_huge_noise.dat:
	matlab -nojvm < run_huge_noise.m

clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux

