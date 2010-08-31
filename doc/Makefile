
all:: elegant.ps elegant.pdf

elegant.ps:
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc dictionary.tex
	../O.linux-x86_64/elegant print_dict.ele
	latex elegant
	latex elegant
	dvips elegant -o elegant.ps

elegant.pdf: elegant.ps
	ps2pdf13 elegant.ps elegant.pdf

clean::
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc elegant.ps elegant.pdf dictionary.tex
	rm -rf elegant
