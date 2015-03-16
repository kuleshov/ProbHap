SIH:=path/to/SingleIndividualHaplotyper/SIH.jar
INPUTDIR:=../data/input/
VALIDATIONDIR:=../data/validation/
MIXSIH:=path/to/MixSIH
PROBHAP:=../probhap.py
PROBHAP_POSTPROCESS:=../probhap-postprocess.py
EVALUATOR:=../libprism/eval/evaluator.py
COMBINER:=../libprism/eval/combine_stats.py

CHR:=chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22

run:
	make $(CHR:=.refhap.stats)
	make $(CHR:=.probhap.corrected.stats)
	# make refhap.combined-stats probhap.corrected.combined-stats dgs.combined-stats fasthare.combined-stats mixsih.combined-stats

other:
	make $(CHR:=.dgs.out) $(CHR:=.fasthare.out) $(CHR:=.mixsih.out)

$(CHR:=.refhap.out): %.refhap.out: $(INPUTDIR)/%.matrix.SORTED
	time java -cp $(SIH) mpg.molgen.sih.main.SIH $< $@

$(CHR:=.probhap.out): %.probhap.out: $(INPUTDIR)/%.matrix.SORTED
	time python $(PROBHAP) \
				--reads $< \
				--parsed-reads $*.probhap.reads \
				--phase $@ \
				--assignments $*.probhap.assignments

$(CHR:=.probhap.corrected.out): %.probhap.corrected.out: %.probhap.out
	python $(PROBHAP_POSTPROCESS) \
		--filtered-reads $*.probhap.reads \
		--assignments $*.probhap.assignments \
		--blocks $< \
		--corrected-blocks $@.tmp
	mv $@.tmp $@

$(CHR:=.mixsih.out): %.mixsih.out: $(INPUTDIR)/%.matrix.SORTED
	time $(MIXSIH) $< /dev/null $@

$(CHR:=.dgs.out): %.dgs.out: $(INPUTDIR)/%.matrix.SORTED
	time java -cp $(SIH) mpg.molgen.sih.main.SIH -a DGS $< $@

$(CHR:=.fasthare.out): %.fasthare.out: $(INPUTDIR)/%.matrix.SORTED
	time java -cp $(SIH) mpg.molgen.sih.main.SIH -a FastHare $< $@

%.stats: %.out
	python $(EVALUATOR) \
			--master $(VALIDATIONDIR)/$(firstword $(subst ., ,$*)).valid.master \
			--chr $(firstword $(subst ., ,$*)) \
			--alg-phase $< > $@.tmp
	mv $@.tmp $@

%.combined-stats: $(CHR:=.%.stats)
	python $(COMBINER) \
		--stat-files $^ \
		> $@.tmp
	mv $@.tmp $@

# %.stats: %.out
# 	python $(EVALUATOR) \
# 			--master ~/probhap/probhap-local/data/chr20.master \
# 			--alg-phase $< > $@.tmp
# 	mv $@.tmp $@
