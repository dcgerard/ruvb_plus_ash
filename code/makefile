all : ex1 ex2

.PHONY : ex1
ex1 : exchangeable_mouthwash_sims.R
	mkdir -p ../output/alpha_1_sims_out
	Rscript exchangeable_mouthwash_sims.R

.PHONY : ex2
ex2 : exchangeable_mouthwash_sims_asifalpha0.R
	mkdir -p ../output/alpha_1_sims_out
	Rscript exchangeable_mouthwash_sims_asifalpha0.R

.PHONY : ex3
ex3 : exchangeable_mouthwash_sims_pen.R
	mkdir -p ../output/alpha_1_sims_out
	Rscript exchangeable_mouthwash_sims_pen.R
