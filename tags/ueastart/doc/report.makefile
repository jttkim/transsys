# makefile for typesetting transsys_report
# (yes, this means I haven't mastered document handling with
# automake -- yet...)

REPORTEXAMPLES = tr_ex00.tra tr_ex01.tra tr_ex02.tra tr_ex03.tra tr_ex04.tra tr_ex05.tra tr_ex06.tra \
trl_ex00.trl trl_ex01.trl trl_ex02.trl trl_ex03.trl trl_ex04.trl trl_ex05.trl trl_ex06.trl trl_ex07.trl trl_ex08.trl trl_ex09.trl trl_ex10.trl
REPORTFULLEXAMPLES = tr_ex02full.tra  tr_ex04full.tra  tr_ex05full.tra tr_ex06full.tra \
trl_ex04full.trl trl_ex05full.trl trl_ex06full.trl trl_ex09full.trl trl_ex10full.trl
REPORTGNUPLOTEPS = tr_ex02.eps tr_ex04.eps tr_ex05.eps tr_ex06.eps
REPORTTRANSSYSEPS = tr_ex06_netgraph.eps
REPORTLTRANSSYSEPS = trl_ex05_glgraph.eps trl_ex06_glgraph.eps trl_ex09_glgraph.eps

transsys_report.pdf : transsys_report.ps
	ps2pdf13 transsys_report.ps transsys_report.pdf

transsys_report.ps : transsys_report.dvi
	dvips -t a4 -P amz -P cmz -o transsys_report.ps transsys_report

transsys_report.dvi : transsys_report.tex transsys_report.gpc $(REPORTEXAMPLES) $(REPORTFULLEXAMPLES) $(REPORTTRANSSYSEPS) $(REPORTLTRANSSYSEPS)
	gnuplot transsys_report.gpc
	latex transsys_report
	latex transsys_report
	latex transsys_report
	latex transsys_report

tr_ex06_netgraph.eps : tr_ex06full.tra
	transps -e tr_ex06full.tra tr_ex06_netgraph.eps

trl_ex05_glgraph.ppm : trl_ex05full.trl
	@echo 'Please press "i", then "<Esc>" in ltransgl window'
	ltransgl -p -5,5.5,11 -r 90,0,0 -d 4 -i trl_ex05_glgraph.ppm trl_ex05full.trl 

trl_ex06_glgraph.ppm : trl_ex06full.trl
	@echo 'Please press "i", then "<Esc>" in ltransgl window'
	ltransgl -p 0,8,15 -r 90,0,0 -d 4 -i trl_ex06_glgraph.ppm trl_ex06full.trl 

trl_ex09_glgraph.ppm : trl_ex09full.trl
	@echo 'Please press "i", then "<Esc>" in ltransgl window'
	ltransgl -p 0,2.6,5.2 -r 90,0,0 -d 500 -i trl_ex09_glgraph.ppm trl_ex09full.trl 

%.eps : %.ppm
	pnmtops -noturn $*.ppm > $*.eps
