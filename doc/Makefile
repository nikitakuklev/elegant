
all:: elegant.ps elegant.pdf html

elegant.ps:
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc
	latex elegant
	latex elegant
	dvips elegant -o elegant.ps

elegant.pdf: elegant.ps
	ps2pdf13 elegant.ps elegant.pdf

html:
	latex2html -local_icons elegant

clean::
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc elegant.ps elegant.pdf
	rm -rf elegant